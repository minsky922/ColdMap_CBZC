//  Copyright (c) 2011-present, Facebook, Inc.  All rights reserved.
//  This source code is licensed under both the GPLv2 (found in the
//  COPYING file in the root directory) and Apache 2.0 License
//  (found in the LICENSE.Apache file in the root directory).
//
// Copyright (c) 2011 The LevelDB Authors. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file. See the AUTHORS file for names of contributors.

#include "db/builder.h"

#include <algorithm>
#include <deque>
#include <vector>

#include "db/blob/blob_file_builder.h"
#include "db/compaction/compaction_iterator.h"
#include "db/event_helpers.h"
#include "db/internal_stats.h"
#include "db/merge_helper.h"
#include "db/output_validator.h"
#include "db/range_del_aggregator.h"
#include "db/table_cache.h"
#include "db/version_edit.h"
#include "file/file_util.h"
#include "file/filename.h"
#include "file/read_write_util.h"
#include "file/writable_file_writer.h"
#include "monitoring/iostats_context_imp.h"
#include "monitoring/thread_status_util.h"
#include "options/options_helper.h"
#include "rocksdb/db.h"
#include "rocksdb/env.h"
#include "rocksdb/iterator.h"
#include "rocksdb/options.h"
#include "rocksdb/table.h"
#include "table/block_based/block_based_table_builder.h"
#include "table/format.h"
#include "table/internal_iterator.h"
#include "table/unique_id_impl.h"
#include "test_util/sync_point.h"
#include "util/stop_watch.h"

namespace ROCKSDB_NAMESPACE {

class TableFactory;

TableBuilder* NewTableBuilder(const TableBuilderOptions& tboptions,
                              WritableFileWriter* file) {
  assert((tboptions.column_family_id ==
          TablePropertiesCollectorFactory::Context::kUnknownColumnFamily) ==
         tboptions.column_family_name.empty());
  return tboptions.ioptions.table_factory->NewTableBuilder(tboptions, file);
}

/*  memtable에 저장된 데이터를 디스크에 SST 파일 형식으로 저장하는 과정에서
 * 사용 - Flush 작업의 일환*/
Status BuildTable(
    const std::string& dbname, VersionSet* versions,
    const ImmutableDBOptions& db_options, const TableBuilderOptions& tboptions,
    const FileOptions& file_options, TableCache* table_cache,
    InternalIterator* iter,
    std::vector<std::unique_ptr<FragmentedRangeTombstoneIterator>>
        range_del_iters,
    FileMetaData* meta, std::vector<BlobFileAddition>* blob_file_additions,
    std::vector<SequenceNumber> snapshots,
    SequenceNumber earliest_write_conflict_snapshot,
    SequenceNumber job_snapshot, SnapshotChecker* snapshot_checker,
    bool paranoid_file_checks, InternalStats* internal_stats,
    IOStatus* io_status, const std::shared_ptr<IOTracer>& io_tracer,
    BlobFileCreationReason blob_creation_reason, EventLogger* event_logger,
    int job_id, const Env::IOPriority io_priority,
    TableProperties* table_properties, Env::WriteLifeTimeHint write_hint,
    const std::string* full_history_ts_low,
    BlobFileCompletionCallback* blob_callback, uint64_t* num_input_entries,
    uint64_t* memtable_payload_bytes, uint64_t* memtable_garbage_bytes) {
  assert((tboptions.column_family_id ==
          TablePropertiesCollectorFactory::Context::kUnknownColumnFamily) ==
         tboptions.column_family_name.empty());
  auto& mutable_cf_options = tboptions.moptions;
  auto& ioptions = tboptions.ioptions;
  // Reports the IOStats for flush for every following bytes.
  const size_t kReportFlushIOStatsEvery = 1048576;
  OutputValidator output_validator(
      tboptions.internal_comparator,
      /*enable_order_check=*/
      mutable_cf_options.check_flush_compaction_key_order,
      /*enable_hash=*/paranoid_file_checks);
  Status s;
  meta->fd.file_size = 0;
  iter->SeekToFirst();
  std::unique_ptr<CompactionRangeDelAggregator> range_del_agg(
      new CompactionRangeDelAggregator(&tboptions.internal_comparator,
                                       snapshots));
  uint64_t num_unfragmented_tombstones = 0;
  uint64_t total_tombstone_payload_bytes = 0;
  for (auto& range_del_iter : range_del_iters) {
    num_unfragmented_tombstones +=
        range_del_iter->num_unfragmented_tombstones();
    total_tombstone_payload_bytes +=
        range_del_iter->total_tombstone_payload_bytes();
    range_del_agg->AddTombstones(std::move(range_del_iter));
  }
  // sst 파일이름 생성(컬럼패밀리경로와 파일id 기반으로)
  std::string fname = TableFileName(ioptions.cf_paths, meta->fd.GetNumber(),
                                    meta->fd.GetPathId());
  std::vector<std::string> blob_file_paths;
  std::string file_checksum = kUnknownFileChecksum;
  std::string file_checksum_func_name = kUnknownFileChecksumFuncName;
#ifndef ROCKSDB_LITE
  EventHelpers::NotifyTableFileCreationStarted(ioptions.listeners, dbname,
                                               tboptions.column_family_name,
                                               fname, job_id, tboptions.reason);
#endif  // !ROCKSDB_LITE
  Env* env = db_options.env;
  assert(env);
  FileSystem* fs = db_options.fs.get();
  assert(fs);

  TableProperties tp;
  bool table_file_created = false;
  if (iter->Valid() || !range_del_agg->IsEmpty()) {
    std::unique_ptr<CompactionFilter> compaction_filter;
    if (ioptions.compaction_filter_factory != nullptr &&
        ioptions.compaction_filter_factory->ShouldFilterTableFileCreation(
            tboptions.reason)) {
      CompactionFilter::Context context;
      context.is_full_compaction = false;
      context.is_manual_compaction = false;
      context.column_family_id = tboptions.column_family_id;
      context.reason = tboptions.reason;
      compaction_filter =
          ioptions.compaction_filter_factory->CreateCompactionFilter(context);
      if (compaction_filter != nullptr &&
          !compaction_filter->IgnoreSnapshots()) {
        s.PermitUncheckedError();
        return Status::NotSupported(
            "CompactionFilter::IgnoreSnapshots() = false is not supported "
            "anymore.");
      }
    }
    // SST 파일을 작성하는 역할
    TableBuilder* builder;
    // 파일에 데이터를 실제로 쓰는 역할
    std::unique_ptr<WritableFileWriter> file_writer;
    {
      std::unique_ptr<FSWritableFile> file;
#ifndef NDEBUG
      bool use_direct_writes = file_options.use_direct_writes;
      TEST_SYNC_POINT_CALLBACK("BuildTable:create_file", &use_direct_writes);
#endif  // !NDEBUG
      IOStatus io_s = NewWritableFile(fs, fname, &file, file_options);
      assert(s.ok());
      s = io_s;
      if (io_status->ok()) {
        *io_status = io_s;
      }
      if (!s.ok()) {
        EventHelpers::LogAndNotifyTableFileCreationFinished(
            event_logger, ioptions.listeners, dbname,
            tboptions.column_family_name, fname, job_id, meta->fd,
            kInvalidBlobFileNumber, tp, tboptions.reason, s, file_checksum,
            file_checksum_func_name);
        return s;
      }
      // 파일이 성공적으로 생성되었음을 기록
      table_file_created = true;
      FileTypeSet tmp_set = ioptions.checksum_handoff_file_types;
      // 파일에 I/O 우선순위와 수명 힌트를 설정
      file->SetIOPriority(io_priority);
      //
      std::string sst(".sst");
      if (std::equal(sst.rbegin(), sst.rend(), fname.rbegin())) {
        file->level_ = tboptions.level_at_creation;
        // builder->SetFileNumber(meta->fd.GetNumber());
        file->fno_ = meta->fd.GetNumber();  // 파일 번호 설정
        // file->input_fno_.clear();  // 입력 파일 번호 목록 초기화
      }
      //
      file->SetWriteLifeTimeHint(write_hint);
      file_writer.reset(new WritableFileWriter(
          std::move(file), fname, file_options, ioptions.clock, io_tracer,
          ioptions.stats, ioptions.listeners,
          ioptions.file_checksum_gen_factory.get(),
          tmp_set.Contains(FileType::kTableFile), false));
      // TableBuilder를 통해 SST 파일을 생성
      builder = NewTableBuilder(tboptions, file_writer.get());
    }
    // MergeHelper는 컴팩션 과정에서 데이터 병합을 관리
    MergeHelper merge(
        env, tboptions.internal_comparator.user_comparator(),
        ioptions.merge_operator.get(), compaction_filter.get(), ioptions.logger,
        true /* internal key corruption is not ok */,
        snapshots.empty() ? 0 : snapshots.back(), snapshot_checker);
    // Blob 파일 빌더 초기화 (Blob 파일이 활성화된 경우에만)
    std::unique_ptr<BlobFileBuilder> blob_file_builder(
        (mutable_cf_options.enable_blob_files &&
         tboptions.level_at_creation >=
             mutable_cf_options.blob_file_starting_level &&
         blob_file_additions)
            ? new BlobFileBuilder(
                  versions, fs, &ioptions, &mutable_cf_options, &file_options,
                  job_id, tboptions.column_family_id,
                  tboptions.column_family_name, io_priority, write_hint,
                  io_tracer, blob_callback, blob_creation_reason,
                  &blob_file_paths, blob_file_additions)
            : nullptr);

    // 컴팩션 작업 중 병합, 삭제 및 스냅샷 처리
    const std::atomic<bool> kManualCompactionCanceledFalse{false};
    CompactionIterator c_iter(
        iter, tboptions.internal_comparator.user_comparator(), &merge,
        kMaxSequenceNumber, &snapshots, earliest_write_conflict_snapshot,
        job_snapshot, snapshot_checker, env,
        ShouldReportDetailedTime(env, ioptions.stats),
        true /* internal key corruption is not ok */, range_del_agg.get(),
        blob_file_builder.get(), ioptions.allow_data_in_errors,
        ioptions.enforce_single_del_contracts,
        /*manual_compaction_canceled=*/kManualCompactionCanceledFalse,
        /*compaction=*/nullptr, compaction_filter.get(),
        /*shutting_down=*/nullptr, db_options.info_log, full_history_ts_low);
    // 이터레이터를 첫 번째 키로 이동하고, 데이터를 순차적으로 처리
    c_iter.SeekToFirst();
    for (; c_iter.Valid(); c_iter.Next()) {
      const Slice& key = c_iter.key();                // 현재 키 가져오기
      const Slice& value = c_iter.value();            // 현재 값 가져오기
      const ParsedInternalKey& ikey = c_iter.ikey();  // 내부 키 구문 분석
      // Generate a rolling 64-bit hash of the key and values
      // Note :
      // Here "key" integrates 'sequence_number'+'kType'+'user key'.
      // 키와 값에 대한 해시 생성 및 테이블 빌더에 추가
      s = output_validator.Add(key, value);
      if (!s.ok()) {
        break;
      }
      // 키-값 쌍을 테이블에 추가
      builder->Add(key, value);
      /* // TableBuilder는 추상 클래스.
      // NewTableBuilder() 함수는 BlockBasedTableBuilder 인스턴스를 반환.
      // builder->Add() 호출 시 BlockBasedTableBuilder::Add()가 실행.*/

      // 메타데이터 경계 업데이트
      s = meta->UpdateBoundaries(key, value, ikey.sequence, ikey.type);
      if (!s.ok()) {
        break;
      }

      // TODO(noetzli): Update stats after flush, too.
      // 쓰기 중 통계 업데이트
      if (io_priority == Env::IO_HIGH &&
          IOSTATS(bytes_written) >= kReportFlushIOStatsEvery) {
        ThreadStatusUtil::SetThreadOperationProperty(
            ThreadStatus::FLUSH_BYTES_WRITTEN, IOSTATS(bytes_written));
      }
    }
    if (!s.ok()) {
      c_iter.status().PermitUncheckedError();
    } else if (!c_iter.status().ok()) {
      s = c_iter.status();
    }
    // 삭제 범위를 처리하는 이터레이터를 통해 삭제 범위 추가
    if (s.ok()) {
      auto range_del_it = range_del_agg->NewIterator();
      for (range_del_it->SeekToFirst(); range_del_it->Valid();
           range_del_it->Next()) {
        auto tombstone = range_del_it->Tombstone();
        auto kv = tombstone.Serialize();
        builder->Add(kv.first.Encode(), kv.second);
        meta->UpdateBoundariesForRange(kv.first, tombstone.SerializeEndKey(),
                                       tombstone.seq_,
                                       tboptions.internal_comparator);
      }
    }
    // 테이블 빌더가 비어있지 않으면 테이블을 완성
    TEST_SYNC_POINT("BuildTable:BeforeFinishBuildTable");
    const bool empty = builder->IsEmpty();
    if (num_input_entries != nullptr) {
      *num_input_entries =
          c_iter.num_input_entry_scanned() + num_unfragmented_tombstones;
    }
    if (!s.ok() || empty) {
      builder->Abandon();
    } else {
      s = builder->Finish();
    }
    if (io_status->ok()) {
      *io_status = builder->io_status();
    }
    // 테이블 빌드가 성공적으로 완료되었으면 메타데이터 업데이트
    if (s.ok() && !empty) {
      uint64_t file_size = builder->FileSize();
      meta->fd.file_size = file_size;
      meta->marked_for_compaction = builder->NeedCompact();
      assert(meta->fd.GetFileSize() > 0);
      tp = builder
               ->GetTableProperties();  // refresh now that builder is finished
      // 메모리 테이블에서 실제 데이터 크기와 쓰레기 데이터 크기를 계산
      if (memtable_payload_bytes != nullptr &&
          memtable_garbage_bytes != nullptr) {
        const CompactionIterationStats& ci_stats = c_iter.iter_stats();
        uint64_t total_payload_bytes = ci_stats.total_input_raw_key_bytes +
                                       ci_stats.total_input_raw_value_bytes +
                                       total_tombstone_payload_bytes;
        uint64_t total_payload_bytes_written =
            (tp.raw_key_size + tp.raw_value_size);
        // Prevent underflow, which may still happen at this point
        // since we only support inserts, deletes, and deleteRanges.
        if (total_payload_bytes_written <= total_payload_bytes) {
          *memtable_payload_bytes = total_payload_bytes;
          *memtable_garbage_bytes =
              total_payload_bytes - total_payload_bytes_written;
        } else {
          *memtable_payload_bytes = 0;
          *memtable_garbage_bytes = 0;
        }
      }
      if (table_properties) {
        *table_properties = tp;
      }
    }
    delete builder;

    // Finish and check for file errors
    TEST_SYNC_POINT("BuildTable:BeforeSyncTable");
    if (s.ok() && !empty) {
      // 파일을 동기화 (쓰기 완료를 보장하기 위한 단계)
      StopWatch sw(ioptions.clock, ioptions.stats, TABLE_SYNC_MICROS);
      *io_status = file_writer->Sync(ioptions.use_fsync);
    }
    TEST_SYNC_POINT("BuildTable:BeforeCloseTableFile");
    if (s.ok() && io_status->ok() && !empty) {
      // 파일을 닫음 (파일을 닫기 전에 동기화가 완료되어야 함)
      *io_status = file_writer->Close();
    }
    if (s.ok() && io_status->ok() && !empty) {
      // 체크섬 정보를 메타데이터에 추가 (데이터 무결성을 보장)
      // Add the checksum information to file metadata.
      meta->file_checksum = file_writer->GetFileChecksum();
      meta->file_checksum_func_name = file_writer->GetFileChecksumFuncName();
      file_checksum = meta->file_checksum;
      file_checksum_func_name = meta->file_checksum_func_name;
      // Set unique_id only if db_id and db_session_id exist
      // 파일에 고유 ID 설정
      if (!tboptions.db_id.empty() && !tboptions.db_session_id.empty()) {
        if (!GetSstInternalUniqueId(tboptions.db_id, tboptions.db_session_id,
                                    meta->fd.GetNumber(), &(meta->unique_id))
                 .ok()) {
          // if failed to get unique id, just set it Null
          meta->unique_id = kNullUniqueId64x2;
        }
      }
    }

    if (s.ok()) {
      s = *io_status;
    }

    if (blob_file_builder) {
      // Blob 파일 빌더가 존재하면 마무리 (Finish) 작업 진행
      if (s.ok()) {
        s = blob_file_builder->Finish();
      } else {
        blob_file_builder->Abandon(s);
      }
      blob_file_builder.reset();
    }

    // TODO Also check the IO status when create the Iterator.

    TEST_SYNC_POINT("BuildTable:BeforeOutputValidation");
    if (s.ok() && !empty) {
      // 파일이 유효한지 검사 (테이블 검증)
      // Verify that the table is usable
      // We set for_compaction to false and don't OptimizeForCompactionTableRead
      // here because this is a special case after we finish the table building.
      // No matter whether use_direct_io_for_flush_and_compaction is true,
      // the goal is to cache it here for further user reads.
      ReadOptions read_options;
      std::unique_ptr<InternalIterator> it(table_cache->NewIterator(
          read_options, file_options, tboptions.internal_comparator, *meta,
          nullptr /* range_del_agg */, mutable_cf_options.prefix_extractor,
          nullptr,
          (internal_stats == nullptr) ? nullptr
                                      : internal_stats->GetFileReadHist(0),
          TableReaderCaller::kFlush, /*arena=*/nullptr,
          /*skip_filter=*/false, tboptions.level_at_creation,
          MaxFileSizeForL0MetaPin(mutable_cf_options),
          /*smallest_compaction_key=*/nullptr,
          /*largest_compaction_key*/ nullptr,
          /*allow_unprepared_value*/ false));
      s = it->status();
      if (s.ok() && paranoid_file_checks) {
        OutputValidator file_validator(tboptions.internal_comparator,
                                       /*enable_order_check=*/true,
                                       /*enable_hash=*/true);
        for (it->SeekToFirst(); it->Valid(); it->Next()) {
          // Generate a rolling 64-bit hash of the key and values
          file_validator.Add(it->key(), it->value()).PermitUncheckedError();
        }
        s = it->status();
        if (s.ok() && !output_validator.CompareValidator(file_validator)) {
          s = Status::Corruption("Paranoid checksums do not match");
        }
      }
    }
  }

  // Check for input iterator errors
  if (!iter->status().ok()) {
    s = iter->status();
  }

  if (!s.ok() || meta->fd.GetFileSize() == 0) {
    TEST_SYNC_POINT("BuildTable:BeforeDeleteFile");

    constexpr IODebugContext* dbg = nullptr;

    if (table_file_created) {
      Status ignored = fs->DeleteFile(fname, IOOptions(), dbg);
      ignored.PermitUncheckedError();
    }

    assert(blob_file_additions || blob_file_paths.empty());

    if (blob_file_additions) {
      for (const std::string& blob_file_path : blob_file_paths) {
        Status ignored = DeleteDBFile(&db_options, blob_file_path, dbname,
                                      /*force_bg=*/false, /*force_fg=*/false);
        ignored.PermitUncheckedError();
        TEST_SYNC_POINT("BuildTable::AfterDeleteFile");
      }
    }
  }

  Status status_for_listener = s;
  if (meta->fd.GetFileSize() == 0) {
    fname = "(nil)";
    if (s.ok()) {
      status_for_listener = Status::Aborted("Empty SST file not kept");
    }
  }
  // Output to event logger and fire events.
  EventHelpers::LogAndNotifyTableFileCreationFinished(
      event_logger, ioptions.listeners, dbname, tboptions.column_family_name,
      fname, job_id, meta->fd, meta->oldest_blob_file_number, tp,
      tboptions.reason, status_for_listener, file_checksum,
      file_checksum_func_name);

  return s;
}

}  // namespace ROCKSDB_NAMESPACE
