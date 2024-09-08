// Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved.
// Copyright (c) 2019-present, Western Digital Corporation
//  This source code is licensed under both the GPLv2 (found in the
//  COPYING file in the root directory) and Apache 2.0 License
//  (found in the LICENSE.Apache file in the root directory).
/* ZenFS의 주요 기능을 구현하는 파일로,
 파일 시스템의 마운트, 언마운트, 읽기, 쓰기 등의 기능을 포함합니다. 이 파일은
 ZenFS의 핵심적인 로직을 담고 있습니다.*/
#if !defined(ROCKSDB_LITE) && defined(OS_LINUX)

#include "fs_zenfs.h"
///
#include <chrono>
#include <ctime>
#include <iostream>
///
#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <mntent.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <set>
#include <sstream>
#include <utility>
#include <vector>

#ifdef ZENFS_EXPORT_PROMETHEUS
#include "metrics_prometheus.h"
#endif
#include "rocksdb/utilities/object_registry.h"
#include "snapshot.h"
#include "util/coding.h"
#include "util/crc32c.h"

#define DEFAULT_ZENV_LOG_PATH "/tmp/"

namespace ROCKSDB_NAMESPACE {
/* Slice 객체에서 슈퍼블록(superblock)을 디코딩하여 데이터를 읽어오는 작업을
 * 수행합니다. 주어진 입력이 예상된 크기와 일치하는지 확인하고, 다양한 필드를
 * 추출하여 슈퍼블록 구조체에 저장합니다. 또한, 슈퍼블록의 무결성을 검증하여
 * 올바른 포맷인지 확인합니다.*/
Status Superblock::DecodeFrom(Slice* input) {
  if (input->size() != ENCODED_SIZE) {
    return Status::Corruption("ZenFS Superblock",
                              "Error: Superblock size missmatch");
  }

  GetFixed32(input, &magic_);
  memcpy(&uuid_, input->data(), sizeof(uuid_));
  input->remove_prefix(sizeof(uuid_));
  GetFixed32(input, &sequence_);
  GetFixed32(input, &superblock_version_);
  GetFixed32(input, &flags_);
  GetFixed32(input, &block_size_);
  GetFixed32(input, &zone_size_);
  GetFixed32(input, &nr_zones_);
  GetFixed32(input, &finish_treshold_);
  memcpy(&aux_fs_path_, input->data(), sizeof(aux_fs_path_));
  input->remove_prefix(sizeof(aux_fs_path_));
  memcpy(&zenfs_version_, input->data(), sizeof(zenfs_version_));
  input->remove_prefix(sizeof(zenfs_version_));
  memcpy(&reserved_, input->data(), sizeof(reserved_));
  input->remove_prefix(sizeof(reserved_));
  assert(input->size() == 0);

  if (magic_ != MAGIC)
    return Status::Corruption("ZenFS Superblock", "Error: Magic missmatch");
  if (superblock_version_ != CURRENT_SUPERBLOCK_VERSION) {
    return Status::Corruption(
        "ZenFS Superblock",
        "Error: Incompatible ZenFS on-disk format version, "
        "please migrate data or switch to previously used ZenFS version. "
        "See the ZenFS README for instructions.");
  }

  return Status::OK();
}

/*슈퍼블록(superblock) 의 데이터를 주어진 std::string 객체에
    인코딩하는 작업을 수행합니다 .이는 슈퍼블록의 모든 필드를
    시리얼라이즈(serialized) 하여 std::string 객체에 저장합니다.이 작업은
    슈퍼블록의 현재 상태를 저장하거나 전송하기 위해 필요합니다*/
void Superblock::EncodeTo(std::string* output) {
  sequence_++; /* Ensure that this superblock representation is unique */
  output->clear();
  PutFixed32(output, magic_);
  output->append(uuid_, sizeof(uuid_));
  PutFixed32(output, sequence_);
  PutFixed32(output, superblock_version_);
  PutFixed32(output, flags_);
  PutFixed32(output, block_size_);
  PutFixed32(output, zone_size_);
  PutFixed32(output, nr_zones_);
  PutFixed32(output, finish_treshold_);
  output->append(aux_fs_path_, sizeof(aux_fs_path_));
  output->append(zenfs_version_, sizeof(zenfs_version_));
  output->append(reserved_, sizeof(reserved_));
  assert(output->length() == ENCODED_SIZE);
}

void Superblock::GetReport(std::string* reportString) {
  reportString->append("Magic:\t\t\t\t");
  PutFixed32(reportString, magic_);
  reportString->append("\nUUID:\t\t\t\t");
  reportString->append(uuid_);
  reportString->append("\nSequence Number:\t\t");
  reportString->append(std::to_string(sequence_));
  reportString->append("\nSuperblock Version:\t\t");
  reportString->append(std::to_string(superblock_version_));
  reportString->append("\nFlags [Decimal]:\t\t");
  reportString->append(std::to_string(flags_));
  reportString->append("\nBlock Size [Bytes]:\t\t");
  reportString->append(std::to_string(block_size_));
  reportString->append("\nZone Size [Blocks]:\t\t");
  reportString->append(std::to_string(zone_size_));
  reportString->append("\nNumber of Zones:\t\t");
  reportString->append(std::to_string(nr_zones_));
  reportString->append("\nFinish Threshold [%]:\t\t");
  reportString->append(std::to_string(finish_treshold_));
  reportString->append("\nGarbage Collection Enabled:\t");
  reportString->append(std::to_string(!!(flags_ & FLAGS_ENABLE_GC)));
  reportString->append("\nAuxiliary FS Path:\t\t");
  reportString->append(aux_fs_path_);
  reportString->append("\nZenFS Version:\t\t\t");
  std::string zenfs_version = zenfs_version_;
  if (zenfs_version.length() == 0) {
    zenfs_version = "Not Available";
  }
  reportString->append(zenfs_version);
}

/* ZonedBlockDevice 객체와 슈퍼블록이 호환되는지 확인하는 작업을 수행*/
Status Superblock::CompatibleWith(ZonedBlockDevice* zbd) {
  if (block_size_ != zbd->GetBlockSize())
    return Status::Corruption("ZenFS Superblock",
                              "Error: block size missmatch");
  if (zone_size_ != (zbd->GetZoneSize() / block_size_))
    return Status::Corruption("ZenFS Superblock", "Error: zone size missmatch");
  if (nr_zones_ > zbd->GetNrZones())
    return Status::Corruption("ZenFS Superblock",
                              "Error: nr of zones missmatch");

  return Status::OK();
}

/*주어진 Slice 데이터를 메타 로그에 기록합니다. 물리적 저장 크기를 계산하고,
 * 메모리를 할당하여 데이터를 준비한 후, CRC32C를 사용하여 오류 검사를
 * 수행합니다. 준비된 데이터를 존에 추가한 후 메모리를 해제하고 결과 상태를
 * 반환합니다. 이 과정을 통해 데이터 무결성을 보장하고 안정적으로 저장할 수
 * 있습니다*/
IOStatus ZenMetaLog::AddRecord(const Slice& slice) {
  uint32_t record_sz = slice.size();
  const char* data = slice.data();
  size_t phys_sz;
  uint32_t crc = 0;
  char* buffer;
  int ret;
  IOStatus s;

  phys_sz = record_sz + zMetaHeaderSize;

  if (phys_sz % bs_) phys_sz += bs_ - phys_sz % bs_;

  assert(data != nullptr);
  assert((phys_sz % bs_) == 0);

  ret = posix_memalign((void**)&buffer, sysconf(_SC_PAGESIZE), phys_sz);
  if (ret) return IOStatus::IOError("Failed to allocate memory");

  memset(buffer, 0, phys_sz);

  crc = crc32c::Extend(crc, (const char*)&record_sz, sizeof(uint32_t));
  crc = crc32c::Extend(crc, data, record_sz);
  crc = crc32c::Mask(crc);

  EncodeFixed32(buffer, crc);
  EncodeFixed32(buffer + sizeof(uint32_t), record_sz);
  memcpy(buffer + sizeof(uint32_t) * 2, data, record_sz);

  s = zone_->Append(buffer, phys_sz);

  free(buffer);
  return s;
}

/*Slice 객체를 사용하여 메타 로그에서 데이터를 읽어오는 작업을 수행합니다.
 * 함수는 특정 위치에서 데이터를 읽고, 읽기 위치를 업데이트합니다. 또한,
 * EOF(파일 끝)와 영역 경계를 체크하여 적절한 오류를 반환합니다*/
IOStatus ZenMetaLog::Read(Slice* slice) {
  char* data = (char*)slice->data();  // 읽을 데이터를 저장할 버퍼
  size_t read = 0;                    // 현재까지 읽은 데이터 크기
  size_t to_read = slice->size();     // 읽어야 할 전체 데이터 크기
  int ret;                            // 읽기 작업의 반환값

  // EOF 체크
  if (read_pos_ >= zone_->wp_) {
    slice->clear();
    return IOStatus::OK();
  }

  // 존 경계 체크
  if ((read_pos_ + to_read) > (zone_->start_ + zone_->max_capacity_)) {
    return IOStatus::IOError("Read across zone");
  }

  // 데이터 읽기
  while (read < to_read) {
    ret = zbd_->Read(data + read, read_pos_, to_read - read, false);

    if (ret == -1 && errno == EINTR) continue;  // EINTR 오류가 발생하면 재시도
    if (ret < 0)
      return IOStatus::IOError("Read failed");  // 읽기 실패 시 오류 반환

    read += ret;       // 읽은 데이터 크기 업데이트
    read_pos_ += ret;  // 읽기 위치 업데이트
  }

  return IOStatus::OK();
}

/*메타 로그에서 레코드를 읽어와 주어진 Slice 객체와 std::string 객체에 저장하는
 * 작업을 수행합니다. 주요 작업 단계는 레코드의 헤더를 읽고, 레코드 크기 및 CRC
 * 값을 확인한 다음, 데이터를 읽고 CRC 검사를 통해 유효성을 확인하는 것입니다*/
IOStatus ZenMetaLog::ReadRecord(Slice* record, std::string* scratch) {
  Slice header;
  uint32_t record_sz = 0;
  uint32_t record_crc = 0;
  uint32_t actual_crc;
  IOStatus s;

  scratch->clear();
  record->clear();

  scratch->append(zMetaHeaderSize, 0);
  header = Slice(scratch->c_str(), zMetaHeaderSize);

  s = Read(&header);
  if (!s.ok()) return s;

  // EOF?
  if (header.size() == 0) {
    record->clear();
    return IOStatus::OK();
  }

  GetFixed32(&header, &record_crc);
  GetFixed32(&header, &record_sz);

  scratch->clear();
  scratch->append(record_sz, 0);

  *record = Slice(scratch->c_str(), record_sz);
  s = Read(record);
  if (!s.ok()) return s;

  actual_crc = crc32c::Value((const char*)&record_sz, sizeof(uint32_t));
  actual_crc = crc32c::Extend(actual_crc, record->data(), record->size());

  if (actual_crc != crc32c::Unmask(record_crc)) {
    return IOStatus::IOError("Not a valid record");
  }

  /* Next record starts on a block boundary */
  if (read_pos_ % bs_) read_pos_ += bs_ - (read_pos_ % bs_);

  return IOStatus::OK();
}

/*이 생성자는 ZonedBlockDevice, 보조 파일 시스템(aux_fs), 로거(logger)를 받아
 * ZenFS 파일 시스템을 초기화하는 역할을 합니다*/
ZenFS::ZenFS(ZonedBlockDevice* zbd, std::shared_ptr<FileSystem> aux_fs,
             std::shared_ptr<Logger> logger)
    : FileSystemWrapper(aux_fs), zbd_(zbd), logger_(logger) {
  // ZenFS 초기화 시작
  Info(logger_, "ZenFS initializing");

  // ZenFS 초기화 매개변수 로깅
  Info(logger_, "ZenFS parameters: block device: %s, aux filesystem: %s",
       zbd_->GetFilename().c_str(), target()->Name());

  // 다음 파일 ID 초기화
  next_file_id_ = 1;

  // 메타데이터 작성기 초기화
  metadata_writer_.zenFS = this;
}

/*이 소멸자는 ZenFS 객체가 소멸될 때 호출되며, ZenFS 파일 시스템의 정리 작업을
 * 수행합니다. 주요 작업은 로깅, 가비지 컬렉션 스레드 종료, 메타 로그 정리, 파일
 * 맵 정리, ZonedBlockDevice 객체 삭제 등입니다*/
ZenFS::~ZenFS() {
  Status s;                              // 상태 변수 초기화
  Info(logger_, "ZenFS shutting down");  // 종료 메시지 로깅
  zbd_->LogZoneUsage();                  // 존 사용량 로깅
  LogFiles();                            // 파일 정보 로깅

  if (gc_worker_) {          // 가비지 컬렉션 스레드 종료
    run_gc_worker_ = false;  // 가비지 컬렉션 실행 중지
    gc_worker_->join();      // 가비지 컬렉션 스레드 종료 대기
  }

  if (bg_reset_worker_) {
    run_bg_reset_worker_ = false;
    // cv_.notify_all();
    bg_reset_worker_->join();
  }

  std::cout << "GetTotalBytesWritten :: " << zbd_->GetTotalBytesWritten()
            << "\n"
            << "UserByteWritten : " << zbd_->GetUserBytesWritten() << "\n";
  std::cout << "FAR STAT :: WA_zc (mb) : "
            << (zbd_->GetTotalBytesWritten() - zbd_->GetUserBytesWritten()) /
                   (1 << 20)
            << "\n";

  meta_log_.reset(nullptr);  // 메타 로그 정리
  uint64_t non_free = zbd_->GetUsedSpace() + zbd_->GetReclaimableSpace();
  uint64_t free = zbd_->GetFreeSpace();
  uint64_t free_percent = (100 * free) / (free + non_free);

  printf("@@~Zenfs Last Free percent freepercent %ld \n", free_percent);
  ClearFiles();  // 파일 맵 정리
  delete zbd_;   // ZonedBlockDevice 객체 삭제
}

void ZenFS::BackgroundStatTimeLapse() {
  while (run_bg_reset_worker_) {
    free_percent_ = zbd_->CalculateFreePercent();
    zbd_->AddTimeLapse(mount_time_);
    sleep(1);
    mount_time_.fetch_add(1);
  }
}

uint64_t ZenFS::EstimateFileAge(Env::WriteLifeTimeHint hint) {
  switch (hint) {
    case Env::WLTH_SHORT:
      return 1;
    case Env::WLTH_MEDIUM:
      return 2;
    case Env::WLTH_LONG:
      return 3;
    case Env::WLTH_EXTREME:
      return 4;
    default:
      return 1;  // 기본적으로 짧은 수명으로 간주
  }
}

void ZenFS::ZoneCleaning(bool forced) {
  uint64_t zone_size = zbd_->GetZoneSize();
  size_t should_be_copied = 0;
  auto start_chrono =
      std::chrono::high_resolution_clock::now();  // valid data copy time
  auto start_time = std::chrono::system_clock::now();
  auto start_time_t = std::chrono::system_clock::to_time_t(start_time);
  std::cout << "ZoneCleaning started at: " << std::ctime(&start_time_t)
            << std::endl;
  zc_lock_.lock();  // 락을 걸어 동시 접근 방지

  ZenFSSnapshot snapshot;
  ZenFSSnapshotOptions options;
  options.zone_ = 1;
  options.zone_file_ = 1;
  options.log_garbage_ = 1;

  GetZenFSSnapshot(snapshot, options);  // 스냅샷 가져오기
  size_t all_inval_zone_n = 0;
  std::vector<std::pair<uint64_t, uint64_t>> victim_candidate;
  std::set<uint64_t> migrate_zones_start;

  // 모든 파일의 WriteLifeTimeHint 가져오기 - LIZC
  auto lifetime_hints = GetWriteLifeTimeHints();

  // 스냅샷에서 모든 존을 순회하며 상태 확인
  for (const auto& zone : snapshot.zones_) {
    if (zone.capacity != 0) {
      continue;
    }
    if (zone.used_capacity == zone.max_capacity) {
      continue;
    }
    uint64_t garbage_percent_approx =
        100 - 100 * zone.used_capacity / zone.max_capacity;  // 가비지 비율
    if (zone.used_capacity > 0) {  // 유효 데이터(valid data)가 있는 경우
      // 현재 시간을 얻습니다. - CBZC1
      // auto current_time = std::chrono::system_clock::now();
      uint64_t total_age = 0;

      // IOOptions 객체를 생성합니다.
      rocksdb::IOOptions io_options;

      // zone_files_를 사용하여 해당 존에 속한 파일들에 접근합니다.
      for (const auto& zone_file : snapshot.zone_files_) {
        for (const auto& extent : zone_file.extents) {
          if (extent.zone_start ==
              zone.start) {  // 현재 존에 속하는 익스텐트만 처리
            // 익스텐트 정보 출력
            // std::cout << "  Zone File ID: " << zone_file.file_id
            //           << " | Filename: " << zone_file.filename
            //           << " | Extent Start: " << extent.start
            //           << " | Extent Length: " << extent.length << std::endl;
            // std::cout << "  extent.File: " << extent.filename
            //           << " | Zone Start: " << zone.start << std::endl;

            /* CBZC1*/
            // uint64_t file_mod_time = 0;

            // // 파일의 수정 시간을 가져옵니다.
            // IOStatus s = GetFileModificationTime(extent.filename, io_options,
            //                                      &file_mod_time, nullptr);
            // // 수정 시간을 제대로 가져왔다면 출력합니다.
            // if (s.ok()) {
            //   uint64_t file_age =
            //       std::chrono::duration_cast<std::chrono::seconds>(
            //           current_time.time_since_epoch())
            //           .count() -
            //       file_mod_time;
            //   // std::cout << "File_age: " << file_age << std::endl;
            //   total_age += file_age;
            // }

            /*CBZC2 - LIZC*/
            // 수명 힌트를 기반으로 나이 추정
            auto it = lifetime_hints.find(extent.filename);
            if (it != lifetime_hints.end()) {
              Env::WriteLifeTimeHint hint = it->second;
              // // WriteLifeTimeHint 값 출력
              // std::cout << "Filename: " << extent.filename
              //           << ", WriteLifeTimeHint: " << static_cast<int>(hint)
              //           << std::endl;

              uint64_t file_age = EstimateFileAge(hint);  // 추정된 나이 사용
              // // 계산된 file_age 값 출력
              // std::cout << "Estimated file age for " << extent.filename << ":
              // "
              //           << file_age << std::endl;
              total_age += file_age;
            }
          }
        }
      }
      std::cout << "Total_age: " << total_age << std::endl;
      uint64_t denominator = (100 - garbage_percent_approx) * 2;
      // std::cout << "  Denominator: " << denominator << std::endl;
      /* greedy */
      // victim_candidate.push_back(
      //     {garbage_percent_approx, zone.start});

      /* cost-benefit */
      if (denominator != 0) {
        // u = valid
        // garbage_percent_approx(free space+invalid)=1-valid= 1-u ex) 80 %
        // cost = 2u = (100-gpa)*2

        uint64_t cost_benefit_score =
            garbage_percent_approx * total_age / denominator;

        // std::cout << "  Calculated cost-benefit score: " <<
        // cost_benefit_score
        //           << std::endl;

        victim_candidate.push_back({cost_benefit_score, zone.start});
        // std::cout << "  Added to victim_candidate: Zone Start: " <<
        // zone.start
        //           << std::endl;
      }
    } else {  // 유효 데이터가 없는 경우
      all_inval_zone_n++;
      std::cout << "all_inal_zone..." << std::endl;
    }
  }

  /* cost-benefit */
  // benefit/cost = free space generated*age of data / cost
  // = (1-u) * age / (1+u)
  // cost = u(유효데이터(zone.used_capacity)(비율)를 읽는 비용)
  // + u(valid data copy)
  // = 2u
  // benefit = 1-u(free space(garbage_percent_approx(?)))
  // * age(creation file(생성시점)을 기준으로)
  // age = 세그먼트 내 가장 최근에 수정된 블록의 시간을 공간이 여유 상태를
  // 유지할 시간의 추정치로 사용했습니다(즉, 가장 최근에 수정된 블록의 나이)

  // 가비지 비율에 따라 후보 존 정렬(greedy)
  std::cout << "Sorting victim candidates..." << std::endl;
  sort(victim_candidate.rbegin(), victim_candidate.rend());

  // victim_candidate 출력
  // std::cout << "Victim candidates:" << std::endl;
  // for (const auto& candidate : victim_candidate) {
  //   std::cout << "cost-benefit score: " << candidate.first
  //             << ", Zone Start: " << candidate.second << std::endl;
  // }

  uint64_t threshold = 0;
  uint64_t reclaimed_zone_n = 1;

  // if (forced) {
  //   reclaimed_zone_n += 1;
  // }

  // 청소할 존 수 계산
  reclaimed_zone_n = reclaimed_zone_n > victim_candidate.size()
                         ? victim_candidate.size()
                         : reclaimed_zone_n;

  // 청소 대상 존 선택
  for (size_t i = 0;
       (i < reclaimed_zone_n && migrate_zones_start.size() < reclaimed_zone_n);
       i++) {
    if (victim_candidate[i].first > threshold) {
      should_be_copied +=
          (zone_size - (victim_candidate[i].first * zone_size / 100));
      std::cout << "cost-benefit score: " << victim_candidate[i].first
                << ", Zone Start: " << victim_candidate[i].second << std::endl;
      migrate_zones_start.emplace(victim_candidate[i].second);
    }
  }

  std::cout << "ZoneCleaning::reclaimed_zone_n: " << reclaimed_zone_n << "\n";

  std::vector<ZoneExtentSnapshot*> migrate_exts;
  for (auto& ext : snapshot.extents_) {
    if (migrate_zones_start.find(ext.zone_start) != migrate_zones_start.end()) {
      migrate_exts.push_back(&ext);
    }
  }

  if (migrate_zones_start.size() > 0) {
    IOStatus s;
    Info(logger_, "Garbage collecting %d extents \n", (int)migrate_exts.size());
    s = MigrateExtents(migrate_exts);  // 익스텐트 마이그레이션
    zc_triggerd_count_.fetch_add(1);
    if (!s.ok()) {
      Error(logger_, "Garbage collection failed");
    }
    // 종료 시간 기록
    auto end_time = std::chrono::system_clock::now();
    auto end_time_t = std::chrono::system_clock::to_time_t(end_time);
    std::cout << "ZoneCleaning ended at: " << std::ctime(&end_time_t)
              << std::endl;
    if (should_be_copied > 0) {
      auto elapsed = std::chrono::high_resolution_clock::now() - start_chrono;
      long long microseconds =
          std::chrono::duration_cast<std::chrono::microseconds>(elapsed)
              .count();
      zbd_->AddZCTimeLapse(start_time_t, end_time_t, microseconds,
                           migrate_zones_start.size(), should_be_copied,
                           forced);
    }

    zc_lock_.unlock();  // 락 해제
    // return migrate_zones_start.size() + all_inval_zone_n;
  }
}

/* 이 함수는 백그라운드에서 주기적으로 실행되며,
 여유 공간 비율이 설정된 임계값 이하로 떨어지면 가비지 컬렉션을 트리거합니다.
 스냅샷을 통해 현재 파일 시스템 상태를 확인하고,
 가비지가 많은 존을 찾아 데이터를 새로운 존으로 마이그레이션합니다.
 이를 통해 ZenFS의 여유 공간을 확보하고 성능을 유지합니다.*/
void ZenFS::GCWorker() {
  while (run_gc_worker_) {
    usleep(100 * 1000);
    /* 여유 공간 계산*/
    // printf("#######GCWorker start\n");
    // // 사용된 공간과 회수 가능한 공간을 더한 값을 non_free에 저장
    // uint64_t non_free = zbd_->GetUsedSpace() + zbd_->GetReclaimableSpace();
    // // 여유 공간을 free에 저장
    // uint64_t free = zbd_->GetFreeSpace();
    // // 총 공간 대비 여유 공간의 비율을 free_percent에 계산
    // uint64_t free_percent = (100 * free) / (free + non_free);
    // std::cout << "GCWorker : free_percent : " << free_percent << "\n";
    //////////////////
    free_percent_ = zbd_->CalculateFreePercent();
    std::cout << "GCWorker : free_percent_ : " << free_percent_ << "\n";
    if (free_percent_ < 20) {  // 20% 이하일 때만 ZoneCleaning 실행
      ZoneCleaning(true);
    }

    //   usleep(1000 * 1000 * 10);  // 10초 동안 대기
    //   /* 여유 공간 계산*/
    //   // 사용된 공간과 회수 가능한 공간을 더한 값을 non_free에 저장
    //   uint64_t non_free = zbd_->GetUsedSpace() +
    //   zbd_->GetReclaimableSpace();
    //   // 여유 공간을 free에 저장
    //   uint64_t free = zbd_->GetFreeSpace();
    //   // 총 공간 대비 여유 공간의 비율을 free_percent에 계산
    //   uint64_t free_percent = (100 * free) / (free + non_free);
    //   //////////////////
    //   ZenFSSnapshot snapshot;
    //   ZenFSSnapshotOptions options;
    //   // 여유 공간 비율이 GC_START_LEVEL(20)보다 크면 루프의 다음 반복으로
    //   // 넘어갑니다. 즉, 여유 공간이 충분한 경우 GC를 수행하지 않습니다
    //   if (free_percent > GC_START_LEVEL) continue;
    //   /* 스냅샷 옵션을 설정하여 존, 파일, 가비지 정보를 포함하도록 합니다.
    //   */ options.zone_ = 1; options.zone_file_ = 1; options.log_garbage_ =
    //   1;
    //   // GetZenFSSnapshot 함수를 호출하여 파일 시스템의 현재 상태를
    //   스냅샷으로
    //   // 가져옵니다
    //   GetZenFSSnapshot(snapshot, options);
    //   // 가비지 컬렉션 임계값을 계산하여 저장합니다
    //   uint64_t threshold = (100 - GC_SLOPE * (GC_START_LEVEL -
    //   free_percent));
    //   /* 마이그레이션할 존 선택 */
    //   // 스냅샷에서 각 존을 반복하면서, 가비지 비율이 임계값을 초과하는
    //   존을 선택
    //   // 선택된 존의 시작 위치를 migrate_zones_start 집합에 저장
    //   std::set<uint64_t> migrate_zones_start;  // 마이그레이션할 존의 시작
    //   위치를
    //                                            // 저장할 집합 선언
    //   for (const auto& zone : snapshot.zones_) {  // 스냅샷의 각 존을 반복
    //     if (zone.capacity == 0) {  // 현재 존이 사용되지 않는 경우
    //       uint64_t garbage_percent_approx =  // 가비지 비율 계산
    //           100 - 100 * zone.used_capacity / zone.max_capacity;
    //       if (garbage_percent_approx >
    //               threshold &&  // 가비지 비율이 임계값을 초과하고
    //           garbage_percent_approx < 100) {  // 100% 미만인 경우
    //         migrate_zones_start.emplace(
    //             zone.start);  // 존의 시작 위치를 집합에 추가
    //       }
    //     }
    //   }
    //   /* 마이그레이션할 익스텐트 선택 */
    //   // 스냅샷에서 각 익스텐트를 반복하면서, 선택된 존에 속하는 익스텐트를
    //   // migrate_exts 벡터에 저장
    //   std::vector<ZoneExtentSnapshot*> migrate_exts;
    //   for (auto& ext : snapshot.extents_) {
    //     if (migrate_zones_start.find(ext.zone_start) !=
    //         migrate_zones_start.end()) {
    //       migrate_exts.push_back(&ext);
    //     }
    //   }
    //   /* 익스텐트 마이그레이션 */
    //   // 마이그레이션할 익스텐트가 있는 경우, MigrateExtents 함수를
    //   호출하여
    //   // 익스텐트를 마이그레이션 마이그레이션이 실패하면 오류를 로깅
    //   if (migrate_exts.size() > 0) {
    //     IOStatus s;
    //     Info(logger_, "Garbage collecting %d extents \n",
    //          (int)migrate_exts.size());
    //     s = MigrateExtents(migrate_exts);
    //     if (!s.ok()) {
    //       Error(logger_, "Garbage collection failed");
    //     }
    //   }
    // }
  }
}

IOStatus ZenFS::Repair() {
  std::map<std::string, std::shared_ptr<ZoneFile>>::iterator it;
  for (it = files_.begin(); it != files_.end(); it++) {
    std::shared_ptr<ZoneFile> zFile = it->second;
    if (zFile->HasActiveExtent()) {
      IOStatus s = zFile->Recover();
      if (!s.ok()) return s;
    }
  }

  return IOStatus::OK();
}

std::string ZenFS::FormatPathLexically(fs::path filepath) {
  fs::path ret = fs::path("/") / filepath.lexically_normal();
  return ret.string();
}

void ZenFS::LogFiles() {
  std::map<std::string, std::shared_ptr<ZoneFile>>::iterator it;
  uint64_t total_size = 0;

  Info(logger_, "  Files:\n");
  for (it = files_.begin(); it != files_.end(); it++) {
    std::shared_ptr<ZoneFile> zFile = it->second;
    std::vector<ZoneExtent*> extents = zFile->GetExtents();

    Info(logger_, "    %-45s sz: %lu lh: %d sparse: %u", it->first.c_str(),
         zFile->GetFileSize(), zFile->GetWriteLifeTimeHint(),
         zFile->IsSparse());
    for (unsigned int i = 0; i < extents.size(); i++) {
      ZoneExtent* extent = extents[i];
      Info(logger_, "          Extent %u {start=0x%lx, zone=%u, len=%lu} ", i,
           extent->start_,
           (uint32_t)(extent->zone_->start_ / zbd_->GetZoneSize()),
           extent->length_);

      total_size += extent->length_;
    }
  }
  Info(logger_, "Sum of all files: %lu MB of data \n",
       total_size / (1024 * 1024));
}

void ZenFS::ClearFiles() {
  std::map<std::string, std::shared_ptr<ZoneFile>>::iterator it;
  std::lock_guard<std::mutex> file_lock(files_mtx_);
  for (it = files_.begin(); it != files_.end(); it++) it->second.reset();
  files_.clear();
}

/* Assumes that files_mutex_ is held */
IOStatus ZenFS::WriteSnapshotLocked(ZenMetaLog* meta_log) {
  IOStatus s;
  std::string snapshot;

  EncodeSnapshotTo(&snapshot);
  s = meta_log->AddRecord(snapshot);
  if (s.ok()) {
    for (auto it = files_.begin(); it != files_.end(); it++) {
      std::shared_ptr<ZoneFile> zoneFile = it->second;
      zoneFile->MetadataSynced();
    }
  }
  return s;
}

IOStatus ZenFS::WriteEndRecord(ZenMetaLog* meta_log) {
  std::string endRecord;

  PutFixed32(&endRecord, kEndRecord);
  return meta_log->AddRecord(endRecord);
}

/* Assumes the files_mtx_ is held */
IOStatus ZenFS::RollMetaZoneLocked() {
  std::unique_ptr<ZenMetaLog> new_meta_log, old_meta_log;
  Zone* new_meta_zone = nullptr;
  IOStatus s;

  ZenFSMetricsLatencyGuard guard(zbd_->GetMetrics(), ZENFS_ROLL_LATENCY,
                                 Env::Default());
  zbd_->GetMetrics()->ReportQPS(ZENFS_ROLL_QPS, 1);

  IOStatus status = zbd_->AllocateMetaZone(&new_meta_zone);
  if (!status.ok()) return status;

  if (!new_meta_zone) {
    assert(false);  // TMP
    Error(logger_, "Out of metadata zones, we should go to read only now.");
    return IOStatus::NoSpace("Out of metadata zones");
  }

  Info(logger_, "Rolling to metazone %d\n", (int)new_meta_zone->GetZoneNr());
  new_meta_log.reset(new ZenMetaLog(zbd_, new_meta_zone));

  old_meta_log.swap(meta_log_);
  meta_log_.swap(new_meta_log);

  /* Write an end record and finish the meta data zone if there is space left
   */
  if (old_meta_log->GetZone()->GetCapacityLeft())
    WriteEndRecord(old_meta_log.get());
  if (old_meta_log->GetZone()->GetCapacityLeft())
    old_meta_log->GetZone()->Finish();

  std::string super_string;
  superblock_->EncodeTo(&super_string);

  s = meta_log_->AddRecord(super_string);
  if (!s.ok()) {
    Error(logger_,
          "Could not write super block when rolling to a new meta zone");
    return IOStatus::IOError("Failed writing a new superblock");
  }

  s = WriteSnapshotLocked(meta_log_.get());

  /* We've rolled successfully, we can reset the old zone now */
  if (s.ok()) old_meta_log->GetZone()->Reset();

  return s;
}

IOStatus ZenFS::PersistSnapshot(ZenMetaLog* meta_writer) {
  IOStatus s;

  std::lock_guard<std::mutex> file_lock(files_mtx_);
  std::lock_guard<std::mutex> metadata_lock(metadata_sync_mtx_);

  s = WriteSnapshotLocked(meta_writer);
  if (s == IOStatus::NoSpace()) {
    Info(logger_, "Current meta zone full, rolling to next meta zone");
    s = RollMetaZoneLocked();
  }

  if (!s.ok()) {
    Error(logger_,
          "Failed persisting a snapshot, we should go to read only now!");
  }

  return s;
}

IOStatus ZenFS::PersistRecord(std::string record) {
  IOStatus s;

  std::lock_guard<std::mutex> lock(metadata_sync_mtx_);
  s = meta_log_->AddRecord(record);
  if (s == IOStatus::NoSpace()) {
    Info(logger_, "Current meta zone full, rolling to next meta zone");
    s = RollMetaZoneLocked();
    /* After a successfull roll, a complete snapshot has been persisted
     * - no need to write the record update */
  }

  return s;
}

IOStatus ZenFS::SyncFileExtents(ZoneFile* zoneFile,
                                std::vector<ZoneExtent*> new_extents) {
  IOStatus s;

  std::vector<ZoneExtent*> old_extents = zoneFile->GetExtents();
  zoneFile->ReplaceExtentList(new_extents);
  zoneFile->MetadataUnsynced();
  s = SyncFileMetadata(zoneFile, true);

  if (!s.ok()) {
    return s;
  }

  // Clear changed extents' zone stats
  for (size_t i = 0; i < new_extents.size(); ++i) {
    ZoneExtent* old_ext = old_extents[i];
    if (old_ext->start_ != new_extents[i]->start_) {
      old_ext->zone_->used_capacity_ -= old_ext->length_;
    }
    delete old_ext;
  }

  return IOStatus::OK();
}

/* Must hold files_mtx_ */
IOStatus ZenFS::SyncFileMetadataNoLock(ZoneFile* zoneFile, bool replace) {
  std::string fileRecord;
  std::string output;
  IOStatus s;
  ZenFSMetricsLatencyGuard guard(zbd_->GetMetrics(), ZENFS_META_SYNC_LATENCY,
                                 Env::Default());

  if (zoneFile->IsDeleted()) {
    Info(logger_, "File %s has been deleted, skip sync file metadata!",
         zoneFile->GetFilename().c_str());
    return IOStatus::OK();
  }

  if (replace) {
    PutFixed32(&output, kFileReplace);
  } else {
    zoneFile->SetFileModificationTime(time(0));
    PutFixed32(&output, kFileUpdate);
  }
  zoneFile->EncodeUpdateTo(&fileRecord);
  PutLengthPrefixedSlice(&output, Slice(fileRecord));

  s = PersistRecord(output);
  if (s.ok()) zoneFile->MetadataSynced();

  return s;
}

IOStatus ZenFS::SyncFileMetadata(ZoneFile* zoneFile, bool replace) {
  std::lock_guard<std::mutex> lock(files_mtx_);
  return SyncFileMetadataNoLock(zoneFile, replace);
}

/* Must hold files_mtx_ */
std::shared_ptr<ZoneFile> ZenFS::GetFileNoLock(std::string fname) {
  std::shared_ptr<ZoneFile> zoneFile(nullptr);
  fname = FormatPathLexically(fname);
  if (files_.find(fname) != files_.end()) {
    zoneFile = files_[fname];
  }
  return zoneFile;
}

std::shared_ptr<ZoneFile> ZenFS::GetFile(std::string fname) {
  std::shared_ptr<ZoneFile> zoneFile(nullptr);
  std::lock_guard<std::mutex> lock(files_mtx_);
  zoneFile = GetFileNoLock(fname);
  return zoneFile;
}

/* Must hold files_mtx_ */
IOStatus ZenFS::DeleteFileNoLock(std::string fname, const IOOptions& options,
                                 IODebugContext* dbg) {
  std::shared_ptr<ZoneFile> zoneFile(nullptr);
  IOStatus s;

  fname = FormatPathLexically(fname);
  zoneFile = GetFileNoLock(fname);
  if (zoneFile != nullptr) {
    std::string record;

    files_.erase(fname);
    s = zoneFile->RemoveLinkName(fname);
    if (!s.ok()) return s;
    EncodeFileDeletionTo(zoneFile, &record, fname);
    s = PersistRecord(record);
    if (!s.ok()) {
      /* Failed to persist the delete, return to a consistent state */
      files_.insert(std::make_pair(fname.c_str(), zoneFile));
      zoneFile->AddLinkName(fname);
    } else {
      if (zoneFile->GetNrLinks() > 0) return s;
      /* Mark up the file as deleted so it won't be migrated by GC */
      zoneFile->SetDeleted();
      zoneFile.reset();
    }
  } else {
    s = target()->DeleteFile(ToAuxPath(fname), options, dbg);
  }

  return s;
}

IOStatus ZenFS::NewSequentialFile(const std::string& filename,
                                  const FileOptions& file_opts,
                                  std::unique_ptr<FSSequentialFile>* result,
                                  IODebugContext* dbg) {
  std::string fname = FormatPathLexically(filename);
  std::shared_ptr<ZoneFile> zoneFile = GetFile(fname);

  Debug(logger_, "New sequential file: %s direct: %d\n", fname.c_str(),
        file_opts.use_direct_reads);

  if (zoneFile == nullptr) {
    return target()->NewSequentialFile(ToAuxPath(fname), file_opts, result,
                                       dbg);
  }

  result->reset(new ZonedSequentialFile(zoneFile, file_opts));
  return IOStatus::OK();
}

IOStatus ZenFS::NewRandomAccessFile(const std::string& filename,
                                    const FileOptions& file_opts,
                                    std::unique_ptr<FSRandomAccessFile>* result,
                                    IODebugContext* dbg) {
  std::string fname = FormatPathLexically(filename);
  std::shared_ptr<ZoneFile> zoneFile = GetFile(fname);

  Debug(logger_, "New random access file: %s direct: %d\n", fname.c_str(),
        file_opts.use_direct_reads);

  if (zoneFile == nullptr) {
    return target()->NewRandomAccessFile(ToAuxPath(fname), file_opts, result,
                                         dbg);
  }

  result->reset(new ZonedRandomAccessFile(files_[fname], file_opts));
  return IOStatus::OK();
}

inline bool ends_with(std::string const& value, std::string const& ending) {
  if (ending.size() > value.size()) return false;
  return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

IOStatus ZenFS::NewWritableFile(const std::string& filename,
                                const FileOptions& file_opts,
                                std::unique_ptr<FSWritableFile>* result,
                                IODebugContext* /*dbg*/) {
  std::string fname = FormatPathLexically(filename);
  Debug(logger_, "New writable file: %s direct: %d\n", fname.c_str(),
        file_opts.use_direct_writes);

  return OpenWritableFile(fname, file_opts, result, nullptr, false);
}

IOStatus ZenFS::ReuseWritableFile(const std::string& filename,
                                  const std::string& old_filename,
                                  const FileOptions& file_opts,
                                  std::unique_ptr<FSWritableFile>* result,
                                  IODebugContext* dbg) {
  IOStatus s;
  std::string fname = FormatPathLexically(filename);
  std::string old_fname = FormatPathLexically(old_filename);
  Debug(logger_, "Reuse writable file: %s old name: %s\n", fname.c_str(),
        old_fname.c_str());

  if (GetFile(old_fname) == nullptr)
    return IOStatus::NotFound("Old file does not exist");

  /*
   * Delete the old file as it cannot be written from start of file
   * and create a new file with fname
   */
  s = DeleteFile(old_fname, file_opts.io_options, dbg);
  if (!s.ok()) {
    Error(logger_, "Failed to delete file %s\n", old_fname.c_str());
    return s;
  }

  return OpenWritableFile(fname, file_opts, result, dbg, false);
}

IOStatus ZenFS::FileExists(const std::string& filename,
                           const IOOptions& options, IODebugContext* dbg) {
  std::string fname = FormatPathLexically(filename);
  Debug(logger_, "FileExists: %s \n", fname.c_str());

  if (GetFile(fname) == nullptr) {
    return target()->FileExists(ToAuxPath(fname), options, dbg);
  } else {
    return IOStatus::OK();
  }
}

/* If the file does not exist, create a new one,
 * else return the existing file
 */
IOStatus ZenFS::ReopenWritableFile(const std::string& filename,
                                   const FileOptions& file_opts,
                                   std::unique_ptr<FSWritableFile>* result,
                                   IODebugContext* dbg) {
  std::string fname = FormatPathLexically(filename);
  Debug(logger_, "Reopen writable file: %s \n", fname.c_str());

  return OpenWritableFile(fname, file_opts, result, dbg, true);
}

/* Must hold files_mtx_ */
void ZenFS::GetZenFSChildrenNoLock(const std::string& dir,
                                   bool include_grandchildren,
                                   std::vector<std::string>* result) {
  auto path_as_string_with_separator_at_end = [](fs::path const& path) {
    fs::path with_sep = path / fs::path("");
    return with_sep.lexically_normal().string();
  };

  auto string_starts_with = [](std::string const& string,
                               std::string const& needle) {
    return string.rfind(needle, 0) == 0;
  };

  std::string dir_with_terminating_seperator =
      path_as_string_with_separator_at_end(fs::path(dir));

  auto relative_child_path =
      [&dir_with_terminating_seperator](std::string const& full_path) {
        return full_path.substr(dir_with_terminating_seperator.length());
      };

  for (auto const& it : files_) {
    fs::path file_path(it.first);
    assert(file_path.has_filename());

    std::string file_dir =
        path_as_string_with_separator_at_end(file_path.parent_path());

    if (string_starts_with(file_dir, dir_with_terminating_seperator)) {
      if (include_grandchildren ||
          file_dir.length() == dir_with_terminating_seperator.length()) {
        result->push_back(relative_child_path(file_path.string()));
      }
    }
  }
}

/* Must hold files_mtx_ */
IOStatus ZenFS::GetChildrenNoLock(const std::string& dir_path,
                                  const IOOptions& options,
                                  std::vector<std::string>* result,
                                  IODebugContext* dbg) {
  std::vector<std::string> auxfiles;
  std::string dir = FormatPathLexically(dir_path);
  IOStatus s;

  Debug(logger_, "GetChildrenNoLock: %s \n", dir.c_str());

  s = target()->GetChildren(ToAuxPath(dir), options, &auxfiles, dbg);
  if (!s.ok()) {
    /* On ZenFS empty directories cannot be created, therefore we cannot
       distinguish between "Directory not found" and "Directory is empty"
       and always return empty lists with OK status in both cases. */
    if (s.IsNotFound()) {
      return IOStatus::OK();
    }
    return s;
  }

  for (const auto& f : auxfiles) {
    if (f != "." && f != "..") result->push_back(f);
  }

  GetZenFSChildrenNoLock(dir, false, result);

  return s;
}

IOStatus ZenFS::GetChildren(const std::string& dir, const IOOptions& options,
                            std::vector<std::string>* result,
                            IODebugContext* dbg) {
  std::lock_guard<std::mutex> lock(files_mtx_);
  return GetChildrenNoLock(dir, options, result, dbg);
}

/* Must hold files_mtx_ */
IOStatus ZenFS::DeleteDirRecursiveNoLock(const std::string& dir,
                                         const IOOptions& options,
                                         IODebugContext* dbg) {
  std::vector<std::string> children;
  std::string d = FormatPathLexically(dir);
  IOStatus s;

  Debug(logger_, "DeleteDirRecursiveNoLock: %s aux: %s\n", d.c_str(),
        ToAuxPath(d).c_str());

  s = GetChildrenNoLock(d, options, &children, dbg);
  if (!s.ok()) {
    return s;
  }

  for (const auto& child : children) {
    std::string file_to_delete = (fs::path(d) / fs::path(child)).string();
    bool is_dir;

    s = IsDirectoryNoLock(file_to_delete, options, &is_dir, dbg);
    if (!s.ok()) {
      return s;
    }

    if (is_dir) {
      s = DeleteDirRecursiveNoLock(file_to_delete, options, dbg);
    } else {
      s = DeleteFileNoLock(file_to_delete, options, dbg);
    }
    if (!s.ok()) {
      return s;
    }
  }

  return target()->DeleteDir(ToAuxPath(d), options, dbg);
}

IOStatus ZenFS::DeleteDirRecursive(const std::string& d,
                                   const IOOptions& options,
                                   IODebugContext* dbg) {
  IOStatus s;
  {
    std::lock_guard<std::mutex> lock(files_mtx_);
    s = DeleteDirRecursiveNoLock(d, options, dbg);
  }
  if (s.ok()) s = zbd_->ResetUnusedIOZones();
  return s;
}

IOStatus ZenFS::OpenWritableFile(const std::string& filename,
                                 const FileOptions& file_opts,
                                 std::unique_ptr<FSWritableFile>* result,
                                 IODebugContext* dbg, bool reopen) {
  IOStatus s;
  std::string fname = FormatPathLexically(filename);
  bool resetIOZones = false;
  {
    std::lock_guard<std::mutex> file_lock(files_mtx_);
    std::shared_ptr<ZoneFile> zoneFile = GetFileNoLock(fname);

    /* if reopen is true and the file exists, return it */
    if (reopen && zoneFile != nullptr) {
      zoneFile->AcquireWRLock();
      result->reset(
          new ZonedWritableFile(zbd_, !file_opts.use_direct_writes, zoneFile));
      return IOStatus::OK();
    }

    if (zoneFile != nullptr) {
      s = DeleteFileNoLock(fname, file_opts.io_options, dbg);
      if (!s.ok()) return s;
      resetIOZones = true;
    }

    zoneFile = std::make_shared<ZoneFile>(zbd_, next_file_id_++,
                                          &metadata_writer_, this);
    zoneFile->SetFileModificationTime(time(0));
    zoneFile->AddLinkName(fname);

    /* RocksDB does not set the right io type(!)*/
    if (ends_with(fname, ".log")) {
      zoneFile->SetIOType(IOType::kWAL);
      zoneFile->SetSparse(!file_opts.use_direct_writes);
    } else {
      zoneFile->SetIOType(IOType::kUnknown);
    }

    /* Persist the creation of the file */
    s = SyncFileMetadataNoLock(zoneFile);
    if (!s.ok()) {
      zoneFile.reset();
      return s;
    }

    zoneFile->AcquireWRLock();
    files_.insert(std::make_pair(fname.c_str(), zoneFile));
    result->reset(
        new ZonedWritableFile(zbd_, !file_opts.use_direct_writes, zoneFile));
  }

  if (resetIOZones) s = zbd_->ResetUnusedIOZones();

  return s;
}

IOStatus ZenFS::DeleteFile(const std::string& fname, const IOOptions& options,
                           IODebugContext* dbg) {
  IOStatus s;

  Debug(logger_, "DeleteFile: %s \n", fname.c_str());

  files_mtx_.lock();
  s = DeleteFileNoLock(fname, options, dbg);
  files_mtx_.unlock();
  if (s.ok()) s = zbd_->ResetUnusedIOZones();
  zbd_->LogZoneStats();

  return s;
}

IOStatus ZenFS::GetFileModificationTime(const std::string& filename,
                                        const IOOptions& options,
                                        uint64_t* mtime, IODebugContext* dbg) {
  std::shared_ptr<ZoneFile> zoneFile(nullptr);
  std::string f = FormatPathLexically(filename);
  IOStatus s;

  Debug(logger_, "GetFileModificationTime: %s \n", f.c_str());
  std::lock_guard<std::mutex> lock(files_mtx_);
  if (files_.find(f) != files_.end()) {
    zoneFile = files_[f];
    *mtime = (uint64_t)zoneFile->GetFileModificationTime();
  } else {
    s = target()->GetFileModificationTime(ToAuxPath(f), options, mtime, dbg);
  }
  return s;
}

IOStatus ZenFS::GetFileSize(const std::string& filename,
                            const IOOptions& options, uint64_t* size,
                            IODebugContext* dbg) {
  std::shared_ptr<ZoneFile> zoneFile(nullptr);
  std::string f = FormatPathLexically(filename);
  IOStatus s;

  Debug(logger_, "GetFileSize: %s \n", f.c_str());

  std::lock_guard<std::mutex> lock(files_mtx_);
  if (files_.find(f) != files_.end()) {
    zoneFile = files_[f];
    *size = zoneFile->GetFileSize();
  } else {
    s = target()->GetFileSize(ToAuxPath(f), options, size, dbg);
  }

  return s;
}

void ZenFS::SetResetScheme(uint32_t r, bool f, uint64_t T) {
  std::cout << "ZenFS::SetResetScheme: r = " << r << ", f = " << f
            << ", T = " << T << std::endl;
  zbd_->SetResetScheme(r, f, T);
  run_bg_reset_worker_ = true;
  if (gc_worker_ != nullptr) {
    if (bg_reset_worker_ == nullptr) {
      bg_reset_worker_.reset(
          new std::thread(&ZenFS::BackgroundStatTimeLapse, this));
    }
  }
}

/* Must hold files_mtx_ */
IOStatus ZenFS::RenameChildNoLock(std::string const& source_dir,
                                  std::string const& dest_dir,
                                  std::string const& child,
                                  const IOOptions& options,
                                  IODebugContext* dbg) {
  std::string source_child = (fs::path(source_dir) / fs::path(child)).string();
  std::string dest_child = (fs::path(dest_dir) / fs::path(child)).string();
  return RenameFileNoLock(source_child, dest_child, options, dbg);
}

/* Must hold files_mtx_ */
IOStatus ZenFS::RollbackAuxDirRenameNoLock(
    const std::string& source_path, const std::string& dest_path,
    const std::vector<std::string>& renamed_children, const IOOptions& options,
    IODebugContext* dbg) {
  IOStatus s;

  for (const auto& rollback_child : renamed_children) {
    s = RenameChildNoLock(dest_path, source_path, rollback_child, options, dbg);
    if (!s.ok()) {
      return IOStatus::Corruption(
          "RollbackAuxDirRenameNoLock: Failed to roll back directory rename");
    }
  }

  s = target()->RenameFile(ToAuxPath(dest_path), ToAuxPath(source_path),
                           options, dbg);
  if (!s.ok()) {
    return IOStatus::Corruption(
        "RollbackAuxDirRenameNoLock: Failed to roll back auxiliary path "
        "renaming");
  }

  return s;
}

/* Must hold files_mtx_ */
IOStatus ZenFS::RenameAuxPathNoLock(const std::string& source_path,
                                    const std::string& dest_path,
                                    const IOOptions& options,
                                    IODebugContext* dbg) {
  IOStatus s;
  std::vector<std::string> children;
  std::vector<std::string> renamed_children;

  s = target()->RenameFile(ToAuxPath(source_path), ToAuxPath(dest_path),
                           options, dbg);
  if (!s.ok()) {
    return s;
  }

  GetZenFSChildrenNoLock(source_path, true, &children);

  for (const auto& child : children) {
    s = RenameChildNoLock(source_path, dest_path, child, options, dbg);
    if (!s.ok()) {
      IOStatus failed_rename = s;
      s = RollbackAuxDirRenameNoLock(source_path, dest_path, renamed_children,
                                     options, dbg);
      if (!s.ok()) {
        return s;
      }
      return failed_rename;
    }
    renamed_children.push_back(child);
  }

  return s;
}

/* Must hold files_mtx_ */
IOStatus ZenFS::RenameFileNoLock(const std::string& src_path,
                                 const std::string& dst_path,
                                 const IOOptions& options,
                                 IODebugContext* dbg) {
  std::shared_ptr<ZoneFile> source_file(nullptr);
  std::shared_ptr<ZoneFile> existing_dest_file(nullptr);
  std::string source_path = FormatPathLexically(src_path);
  std::string dest_path = FormatPathLexically(dst_path);
  IOStatus s;

  Debug(logger_, "Rename file: %s to : %s\n", source_path.c_str(),
        dest_path.c_str());

  source_file = GetFileNoLock(source_path);
  if (source_file != nullptr) {
    existing_dest_file = GetFileNoLock(dest_path);
    if (existing_dest_file != nullptr) {
      s = DeleteFileNoLock(dest_path, options, dbg);
      if (!s.ok()) {
        return s;
      }
    }

    s = source_file->RenameLink(source_path, dest_path);
    if (!s.ok()) return s;
    files_.erase(source_path);

    files_.insert(std::make_pair(dest_path, source_file));

    s = SyncFileMetadataNoLock(source_file);
    if (!s.ok()) {
      /* Failed to persist the rename, roll back */
      files_.erase(dest_path);
      s = source_file->RenameLink(dest_path, source_path);
      if (!s.ok()) return s;
      files_.insert(std::make_pair(source_path, source_file));
    }
  } else {
    s = RenameAuxPathNoLock(source_path, dest_path, options, dbg);
  }

  return s;
}

IOStatus ZenFS::RenameFile(const std::string& source_path,
                           const std::string& dest_path,
                           const IOOptions& options, IODebugContext* dbg) {
  IOStatus s;
  {
    std::lock_guard<std::mutex> lock(files_mtx_);
    s = RenameFileNoLock(source_path, dest_path, options, dbg);
  }
  if (s.ok()) s = zbd_->ResetUnusedIOZones();
  return s;
}

IOStatus ZenFS::LinkFile(const std::string& file, const std::string& link,
                         const IOOptions& options, IODebugContext* dbg) {
  std::shared_ptr<ZoneFile> src_file(nullptr);
  std::string fname = FormatPathLexically(file);
  std::string lname = FormatPathLexically(link);
  IOStatus s;

  Debug(logger_, "LinkFile: %s to %s\n", fname.c_str(), lname.c_str());
  {
    std::lock_guard<std::mutex> lock(files_mtx_);

    if (GetFileNoLock(lname) != nullptr)
      return IOStatus::InvalidArgument("Failed to create link, target exists");

    src_file = GetFileNoLock(fname);
    if (src_file != nullptr) {
      src_file->AddLinkName(lname);
      files_.insert(std::make_pair(lname, src_file));
      s = SyncFileMetadataNoLock(src_file);
      if (!s.ok()) {
        s = src_file->RemoveLinkName(lname);
        if (!s.ok()) return s;
        files_.erase(lname);
      }
      return s;
    }
  }
  s = target()->LinkFile(ToAuxPath(fname), ToAuxPath(lname), options, dbg);
  return s;
}

IOStatus ZenFS::NumFileLinks(const std::string& file, const IOOptions& options,
                             uint64_t* nr_links, IODebugContext* dbg) {
  std::shared_ptr<ZoneFile> src_file(nullptr);
  std::string fname = FormatPathLexically(file);
  IOStatus s;

  Debug(logger_, "NumFileLinks: %s\n", fname.c_str());
  {
    std::lock_guard<std::mutex> lock(files_mtx_);

    src_file = GetFileNoLock(fname);
    if (src_file != nullptr) {
      *nr_links = (uint64_t)src_file->GetNrLinks();
      return IOStatus::OK();
    }
  }
  s = target()->NumFileLinks(ToAuxPath(fname), options, nr_links, dbg);
  return s;
}

IOStatus ZenFS::AreFilesSame(const std::string& file, const std::string& linkf,
                             const IOOptions& options, bool* res,
                             IODebugContext* dbg) {
  std::shared_ptr<ZoneFile> src_file(nullptr);
  std::shared_ptr<ZoneFile> dst_file(nullptr);
  std::string fname = FormatPathLexically(file);
  std::string link = FormatPathLexically(linkf);
  IOStatus s;

  Debug(logger_, "AreFilesSame: %s, %s\n", fname.c_str(), link.c_str());

  {
    std::lock_guard<std::mutex> lock(files_mtx_);
    src_file = GetFileNoLock(fname);
    dst_file = GetFileNoLock(link);
    if (src_file != nullptr && dst_file != nullptr) {
      if (src_file->GetID() == dst_file->GetID())
        *res = true;
      else
        *res = false;
      return IOStatus::OK();
    }
  }
  s = target()->AreFilesSame(fname, link, options, res, dbg);
  return s;
}

void ZenFS::EncodeSnapshotTo(std::string* output) {
  std::map<std::string, std::shared_ptr<ZoneFile>>::iterator it;
  std::string files_string;
  PutFixed32(output, kCompleteFilesSnapshot);
  for (it = files_.begin(); it != files_.end(); it++) {
    std::string file_string;
    std::shared_ptr<ZoneFile> zFile = it->second;

    zFile->EncodeSnapshotTo(&file_string);
    PutLengthPrefixedSlice(&files_string, Slice(file_string));
  }
  PutLengthPrefixedSlice(output, Slice(files_string));
}

void ZenFS::EncodeJson(std::ostream& json_stream) {
  bool first_element = true;
  json_stream << "[";
  for (const auto& file : files_) {
    if (first_element) {
      first_element = false;
    } else {
      json_stream << ",";
    }
    file.second->EncodeJson(json_stream);
  }
  json_stream << "]";
}

Status ZenFS::DecodeFileUpdateFrom(Slice* slice, bool replace) {
  std::shared_ptr<ZoneFile> update(
      new ZoneFile(zbd_, 0, &metadata_writer_, this));
  uint64_t id;
  Status s;

  s = update->DecodeFrom(slice);
  if (!s.ok()) return s;

  id = update->GetID();
  if (id >= next_file_id_) next_file_id_ = id + 1;

  /* Check if this is an update or an replace to an existing file */
  for (auto it = files_.begin(); it != files_.end(); it++) {
    std::shared_ptr<ZoneFile> zFile = it->second;
    if (id == zFile->GetID()) {
      for (const auto& name : zFile->GetLinkFiles()) {
        if (files_.find(name) != files_.end())
          files_.erase(name);
        else
          return Status::Corruption("DecodeFileUpdateFrom: missing link file");
      }

      s = zFile->MergeUpdate(update, replace);
      update.reset();

      if (!s.ok()) return s;

      for (const auto& name : zFile->GetLinkFiles())
        files_.insert(std::make_pair(name, zFile));

      return Status::OK();
    }
  }

  /* The update is a new file */
  assert(GetFile(update->GetFilename()) == nullptr);
  files_.insert(std::make_pair(update->GetFilename(), update));

  return Status::OK();
}

Status ZenFS::DecodeSnapshotFrom(Slice* input) {
  Slice slice;

  assert(files_.size() == 0);

  while (GetLengthPrefixedSlice(input, &slice)) {
    std::shared_ptr<ZoneFile> zoneFile(
        new ZoneFile(zbd_, 0, &metadata_writer_, this));
    Status s = zoneFile->DecodeFrom(&slice);
    if (!s.ok()) return s;

    if (zoneFile->GetID() >= next_file_id_)
      next_file_id_ = zoneFile->GetID() + 1;

    for (const auto& name : zoneFile->GetLinkFiles())
      files_.insert(std::make_pair(name, zoneFile));
  }

  return Status::OK();
}

void ZenFS::EncodeFileDeletionTo(std::shared_ptr<ZoneFile> zoneFile,
                                 std::string* output, std::string linkf) {
  std::string file_string;

  PutFixed64(&file_string, zoneFile->GetID());
  PutLengthPrefixedSlice(&file_string, Slice(linkf));

  PutFixed32(output, kFileDeletion);
  PutLengthPrefixedSlice(output, Slice(file_string));
}

Status ZenFS::DecodeFileDeletionFrom(Slice* input) {
  uint64_t fileID;
  std::string fileName;
  Slice slice;
  IOStatus s;

  if (!GetFixed64(input, &fileID))
    return Status::Corruption("Zone file deletion: file id missing");

  if (!GetLengthPrefixedSlice(input, &slice))
    return Status::Corruption("Zone file deletion: file name missing");

  fileName = slice.ToString();
  if (files_.find(fileName) == files_.end())
    return Status::Corruption("Zone file deletion: no such file");

  std::shared_ptr<ZoneFile> zoneFile = files_[fileName];
  if (zoneFile->GetID() != fileID)
    return Status::Corruption("Zone file deletion: file ID missmatch");

  files_.erase(fileName);
  s = zoneFile->RemoveLinkName(fileName);
  if (!s.ok())
    return Status::Corruption("Zone file deletion: file links missmatch");

  return Status::OK();
}

Status ZenFS::RecoverFrom(ZenMetaLog* log) {
  bool at_least_one_snapshot = false;
  std::string scratch;
  uint32_t tag = 0;
  Slice record;
  Slice data;
  Status s;
  bool done = false;

  while (!done) {
    IOStatus rs = log->ReadRecord(&record, &scratch);
    if (!rs.ok()) {
      Error(logger_, "Read recovery record failed with error: %s",
            rs.ToString().c_str());
      return Status::Corruption("ZenFS", "Metadata corruption");
    }

    if (!GetFixed32(&record, &tag)) break;

    if (tag == kEndRecord) break;

    if (!GetLengthPrefixedSlice(&record, &data)) {
      return Status::Corruption("ZenFS", "No recovery record data");
    }

    switch (tag) {
      case kCompleteFilesSnapshot:
        ClearFiles();
        s = DecodeSnapshotFrom(&data);
        if (!s.ok()) {
          Warn(logger_, "Could not decode complete snapshot: %s",
               s.ToString().c_str());
          return s;
        }
        at_least_one_snapshot = true;
        break;

      case kFileUpdate:
        s = DecodeFileUpdateFrom(&data);
        if (!s.ok()) {
          Warn(logger_, "Could not decode file snapshot: %s",
               s.ToString().c_str());
          return s;
        }
        break;

      case kFileReplace:
        s = DecodeFileUpdateFrom(&data, true);
        if (!s.ok()) {
          Warn(logger_, "Could not decode file snapshot: %s",
               s.ToString().c_str());
          return s;
        }
        break;

      case kFileDeletion:
        s = DecodeFileDeletionFrom(&data);
        if (!s.ok()) {
          Warn(logger_, "Could not decode file deletion: %s",
               s.ToString().c_str());
          return s;
        }
        break;

      default:
        Warn(logger_, "Unexpected metadata record tag: %u", tag);
        return Status::Corruption("ZenFS", "Unexpected tag");
    }
  }

  if (at_least_one_snapshot)
    return Status::OK();
  else
    return Status::NotFound("ZenFS", "No snapshot found");
}

/* 파일 시스템을 마운트하는 것은 운영 체제와 저장 장치 간의 인터페이스를
  설정하여 사용자가 저장 장치의 데이터를 접근하고 관리할 수 있게 하는
  과정입니다. 파일 시스템이 마운트되면, 운영 체제는 특정 디렉토리(마운트
  포인트)를 통해 파일 시스템의 내용을 볼 수 있게 됩니다. */
/* Mount the filesystem by recovering form the latest valid metadata zone */
Status ZenFS::Mount(bool readonly) {
  std::vector<Zone*> metazones = zbd_->GetMetaZones();
  std::vector<std::unique_ptr<Superblock>> valid_superblocks;
  std::vector<std::unique_ptr<ZenMetaLog>> valid_logs;
  std::vector<Zone*> valid_zones;
  std::vector<std::pair<uint32_t, uint32_t>> seq_map;

  Status s;

  /* We need a minimum of two non-offline meta data zones */
  if (metazones.size() < 2) {
    Error(logger_,
          "Need at least two non-offline meta zones to open for write");
    return Status::NotSupported();
  }
  /* 각 메타존에서 유효한 슈퍼블록을 찾아서 valid_superblocks, valid_logs,
   * valid_zones 벡터에 저장 */
  /* Find all valid superblocks */
  for (const auto z : metazones) {
    std::unique_ptr<ZenMetaLog> log;
    std::string scratch;
    Slice super_record;

    if (!z->Acquire()) {
      assert(false);
      return Status::Aborted("Could not aquire busy flag of zone" +
                             std::to_string(z->GetZoneNr()));
    }

    // log takes the ownership of z's busy flag.
    log.reset(new ZenMetaLog(zbd_, z));

    if (!log->ReadRecord(&super_record, &scratch).ok()) continue;

    if (super_record.size() == 0) continue;

    std::unique_ptr<Superblock> super_block;

    super_block.reset(new Superblock());
    s = super_block->DecodeFrom(&super_record);
    if (s.ok()) s = super_block->CompatibleWith(zbd_);
    if (!s.ok()) return s;

    Info(logger_, "Found OK superblock in zone %lu seq: %u\n", z->GetZoneNr(),
         super_block->GetSeq());

    seq_map.push_back(std::make_pair(super_block->GetSeq(), seq_map.size()));
    valid_superblocks.push_back(std::move(super_block));
    valid_logs.push_back(std::move(log));
    valid_zones.push_back(z);
  }

  if (!seq_map.size()) return Status::NotFound("No valid superblock found");

  /* Sort superblocks by descending sequence number */
  std::sort(seq_map.begin(), seq_map.end(),
            std::greater<std::pair<uint32_t, uint32_t>>());

  bool recovery_ok = false;
  unsigned int r = 0;

  /* Recover from the zone with the highest superblock sequence number.
     If that fails go to the previous as we might have crashed when rolling
     metadata zone.
  */
  for (const auto& sm : seq_map) {
    uint32_t i = sm.second;
    std::string scratch;
    std::unique_ptr<ZenMetaLog> log = std::move(valid_logs[i]);

    s = RecoverFrom(log.get());
    if (!s.ok()) {
      if (s.IsNotFound()) {
        Warn(logger_,
             "Did not find a valid snapshot, trying next meta zone. Error: %s",
             s.ToString().c_str());
        continue;
      }

      Error(logger_, "Metadata corruption. Error: %s", s.ToString().c_str());
      return s;
    }

    r = i;
    recovery_ok = true;
    meta_log_ = std::move(log);
    break;
  }

  if (!recovery_ok) {
    return Status::IOError("Failed to mount filesystem");
  }

  Info(logger_, "Recovered from zone: %d", (int)valid_zones[r]->GetZoneNr());
  superblock_ = std::move(valid_superblocks[r]);
  zbd_->SetFinishTreshold(superblock_->GetFinishTreshold());

  IOOptions foo;
  IODebugContext bar;
  s = target()->CreateDirIfMissing(superblock_->GetAuxFsPath(), foo, &bar);
  if (!s.ok()) {
    Error(logger_, "Failed to create aux filesystem directory.");
    return s;
  }

  /* Free up old metadata zones, to get ready to roll */
  for (const auto& sm : seq_map) {
    uint32_t i = sm.second;
    /* Don't reset the current metadata zone */
    if (i != r) {
      /* Metadata zones are not marked as having valid data, so they can be
       * reset */
      valid_logs[i].reset();
    }
  }

  if (!readonly) {
    s = Repair();
    if (!s.ok()) return s;
  }

  if (readonly) {
    Info(logger_, "Mounting READ ONLY");
  } else {
    std::lock_guard<std::mutex> lock(files_mtx_);
    s = RollMetaZoneLocked();
    if (!s.ok()) {
      Error(logger_, "Failed to roll metadata zone.");
      return s;
    }
  }

  Info(logger_, "Superblock sequence %d", (int)superblock_->GetSeq());
  Info(logger_, "Finish threshold %u", superblock_->GetFinishTreshold());
  Info(logger_, "Filesystem mount OK");
  /* 사용하지 않는 IO 존을 리셋하여 공간을 회수합니다. */
  if (!readonly) {
    Info(logger_, "Resetting unused IO Zones..");
    IOStatus status = zbd_->ResetUnusedIOZones();
    if (!status.ok()) return status;
    Info(logger_, "  Done");
    /* 가비지 컬렉션이 활성화된 경우, 가비지 컬렉션 작업자를 시작합니다. */
    if (superblock_->IsGCEnabled()) {
      Info(logger_, "Starting garbage collection worker");
      run_gc_worker_ = true;
      gc_worker_.reset(new std::thread(&ZenFS::GCWorker, this));
      // run_bg_reset_worker_ = true;
      // if (bg_reset_worker_ == nullptr) {
      //   printf("starting bg_reset_worker");
      //   bg_reset_worker_.reset(
      //       new std::thread(&ZenFS::BackgroundStatTimeLapse, this));
      // }
    }
  }

  LogFiles();

  return Status::OK();
}

/* ZenFS 파일 시스템을 초기화하고, 메타 데이터 존을 설정하며, 슈퍼블록을
 * 기록하여 파일 시스템을 사용할 준비를 합니다*/
Status ZenFS::MkFS(std::string aux_fs_p, uint32_t finish_threshold,
                   bool enable_gc) {
  std::vector<Zone*> metazones = zbd_->GetMetaZones();
  std::unique_ptr<ZenMetaLog> log;
  Zone* meta_zone = nullptr;
  std::string aux_fs_path = FormatPathLexically(aux_fs_p);
  IOStatus s;

  if (aux_fs_path.length() > 255) {
    return Status::InvalidArgument(
        "Aux filesystem path must be less than 256 bytes\n");
  }
  ClearFiles();
  IOStatus status = zbd_->ResetUnusedIOZones();
  if (!status.ok()) return status;

  for (const auto mz : metazones) {
    if (!mz->Acquire()) {
      assert(false);
      return Status::Aborted("Could not aquire busy flag of zone " +
                             std::to_string(mz->GetZoneNr()));
    }

    if (mz->Reset().ok()) {
      if (!meta_zone) meta_zone = mz;
    } else {
      Warn(logger_, "Failed to reset meta zone\n");
    }

    if (meta_zone != mz) {
      // for meta_zone == mz the ownership of mz's busy flag is passed to log.
      if (!mz->Release()) {
        assert(false);
        return Status::Aborted("Could not unset busy flag of zone " +
                               std::to_string(mz->GetZoneNr()));
      }
    }
  }

  if (!meta_zone) {
    return Status::IOError("Not available meta zones\n");
  }

  log.reset(new ZenMetaLog(zbd_, meta_zone));

  Superblock super(zbd_, aux_fs_path, finish_threshold, enable_gc);
  std::string super_string;
  super.EncodeTo(&super_string);

  s = log->AddRecord(super_string);
  if (!s.ok()) return static_cast<Status>(s);

  /* Write an empty snapshot to make the metadata zone valid */
  s = PersistSnapshot(log.get());
  if (!s.ok()) {
    Error(logger_, "Failed to persist snapshot: %s", s.ToString().c_str());
    return Status::IOError("Failed persist snapshot");
  }

  Info(logger_, "Empty filesystem created");
  return Status::OK();
}

/*ZenFS에서 관리하는 모든 파일의 이름과 해당 파일의 쓰기 수명 힌트를 맵 형태로
   반환합니다. 이 정보를 통해 파일 시스템 내의 각 파일이 어떤 수명 힌트를 가지고
   있는지 쉽게 확인할 수 있습니다.
*/
std::map<std::string, Env::WriteLifeTimeHint> ZenFS::GetWriteLifeTimeHints() {
  std::map<std::string, Env::WriteLifeTimeHint> hint_map;

  for (auto it = files_.begin(); it != files_.end(); it++) {
    std::shared_ptr<ZoneFile> zoneFile = it->second;
    std::string filename = it->first;
    hint_map.insert(std::make_pair(filename, zoneFile->GetWriteLifeTimeHint()));
  }

  return hint_map;
}

#if !defined(NDEBUG) || defined(WITH_TERARKDB)
static std::string GetLogFilename(std::string bdev) {
  std::ostringstream ss;
  time_t t = time(0);
  struct tm* log_start = std::localtime(&t);
  char buf[40];

  std::strftime(buf, sizeof(buf), "%Y-%m-%d_%H:%M:%S.log", log_start);
  ss << DEFAULT_ZENV_LOG_PATH << std::string("zenfs_") << bdev << "_" << buf;

  return ss.str();
}
#endif

Status NewZenFS(FileSystem** fs, const std::string& bdevname,
                std::shared_ptr<ZenFSMetrics> metrics) {
  return NewZenFS(fs, ZbdBackendType::kBlockDev, bdevname, metrics);
}

Status NewZenFS(FileSystem** fs, const ZbdBackendType backend_type,
                const std::string& backend_name,
                std::shared_ptr<ZenFSMetrics> metrics) {
  std::shared_ptr<Logger> logger;
  Status s;

  // TerarkDB needs to log important information in production while ZenFS
  // doesn't (currently).
  //
  // TODO(guokuankuan@bytedance.com) We need to figure out how to reuse
  // RocksDB's logger in the future.
#if !defined(NDEBUG) || defined(WITH_TERARKDB)
  s = Env::Default()->NewLogger(GetLogFilename(backend_name), &logger);
  if (!s.ok()) {
    fprintf(stderr, "ZenFS: Could not create logger");
  } else {
    logger->SetInfoLogLevel(DEBUG_LEVEL);
#ifdef WITH_TERARKDB
    logger->SetInfoLogLevel(INFO_LEVEL);
#endif
  }
#endif

  ZonedBlockDevice* zbd =
      new ZonedBlockDevice(backend_name, backend_type, logger, metrics);
  IOStatus zbd_status = zbd->Open(false, true);
  if (!zbd_status.ok()) {
    Error(logger, "mkfs: Failed to open zoned block device: %s",
          zbd_status.ToString().c_str());
    return Status::IOError(zbd_status.ToString());
  }

  ZenFS* zenFS = new ZenFS(zbd, FileSystem::Default(), logger);
  s = zenFS->Mount(false);
  if (!s.ok()) {
    delete zenFS;
    return s;
  }

  *fs = zenFS;
  return Status::OK();
}

Status AppendZenFileSystem(
    std::string path, ZbdBackendType backend,
    std::map<std::string, std::pair<std::string, ZbdBackendType>>& fs_map) {
  std::unique_ptr<ZonedBlockDevice> zbd{
      new ZonedBlockDevice(path, backend, nullptr)};
  IOStatus zbd_status = zbd->Open(true, false);

  if (zbd_status.ok()) {
    std::vector<Zone*> metazones = zbd->GetMetaZones();
    std::string scratch;
    Slice super_record;
    Status s;

    for (const auto z : metazones) {
      Superblock super_block;
      std::unique_ptr<ZenMetaLog> log;
      if (!z->Acquire()) {
        return Status::Aborted("Could not aquire busy flag of zone" +
                               std::to_string(z->GetZoneNr()));
      }
      log.reset(new ZenMetaLog(zbd.get(), z));

      if (!log->ReadRecord(&super_record, &scratch).ok()) continue;
      s = super_block.DecodeFrom(&super_record);
      if (s.ok()) {
        /* Map the uuid to the device-mapped (i.g dm-linear) block device to
           avoid trying to mount the whole block device in case of a split
           device */
        if (fs_map.find(super_block.GetUUID()) != fs_map.end() &&
            fs_map[super_block.GetUUID()].first.rfind("dm-", 0) == 0) {
          break;
        }
        fs_map[super_block.GetUUID()] = std::make_pair(path, backend);
        break;
      }
    }
  }
  return Status::OK();
}

Status ListZenFileSystems(
    std::map<std::string, std::pair<std::string, ZbdBackendType>>& out_list) {
  std::map<std::string, std::pair<std::string, ZbdBackendType>> zenFileSystems;

  auto closedirDeleter = [](DIR* d) {
    if (d != nullptr) closedir(d);
  };
  std::unique_ptr<DIR, decltype(closedirDeleter)> dir{
      opendir("/sys/class/block"), std::move(closedirDeleter)};
  struct dirent* entry;

  while (NULL != (entry = readdir(dir.get()))) {
    if (entry->d_type == DT_LNK) {
      Status status =
          AppendZenFileSystem(std::string(entry->d_name),
                              ZbdBackendType::kBlockDev, zenFileSystems);
      if (!status.ok()) return status;
    }
  }

  struct mntent* mnt = NULL;
  FILE* file = NULL;

  file = setmntent("/proc/mounts", "r");
  if (file != NULL) {
    while ((mnt = getmntent(file)) != NULL) {
      if (!strcmp(mnt->mnt_type, "zonefs")) {
        Status status = AppendZenFileSystem(
            std::string(mnt->mnt_dir), ZbdBackendType::kZoneFS, zenFileSystems);
        if (!status.ok()) return status;
      }
    }
  }

  out_list = std::move(zenFileSystems);
  return Status::OK();
}

void ZenFS::GetZenFSSnapshot(ZenFSSnapshot& snapshot,
                             const ZenFSSnapshotOptions& options) {
  if (options.zbd_) {
    snapshot.zbd_ = ZBDSnapshot(*zbd_);
  }
  if (options.zone_) {
    zbd_->GetZoneSnapshot(snapshot.zones_);
  }
  if (options.zone_file_) {
    std::lock_guard<std::mutex> file_lock(files_mtx_);
    for (const auto& file_it : files_) {
      ZoneFile& file = *(file_it.second);

      /* Skip files open for writing, as extents are being updated */
      if (!file.TryAcquireWRLock()) continue;

      // file -> extents mapping
      snapshot.zone_files_.emplace_back(file);
      // extent -> file mapping
      for (auto* ext : file.GetExtents()) {
        snapshot.extents_.emplace_back(*ext, file.GetFilename());
      }

      file.ReleaseWRLock();
    }
  }

  if (options.trigger_report_) {
    zbd_->GetMetrics()->ReportSnapshot(snapshot);
  }

  if (options.log_garbage_) {
    zbd_->LogGarbageInfo();
  }
}

/*  주어진 익스텐트 목록을 파일 이름별로 그룹화하고,
 .sst 파일에 대해 유효 데이터를 마이그레이션합니다.
  마이그레이션이 성공적으로 완료되면 사용되지 않는 IO 존을 재설정*/
IOStatus ZenFS::MigrateExtents(
    const std::vector<ZoneExtentSnapshot*>& extents) {
  IOStatus s;
  // Group extents by their filename
  // file_extents는 파일 이름을 키로 하고, 해당 파일의 익스텐트 목록을 값으로
  // 가지는 맵 fname이 .sst로 끝나는 경우에만 file_extents 맵에 추가
  std::map<std::string, std::vector<ZoneExtentSnapshot*>> file_extents;
  for (auto* ext : extents) {
    std::string fname = ext->filename;
    // We only migrate SST file extents
    // if (ends_with(fname, ".sst")) {
    file_extents[fname].emplace_back(ext);
    // }
  }
  // 파일 익스텐트 마이그레이션 및 존 재설정
  for (const auto& it : file_extents) {
    s = MigrateFileExtents(it.first, it.second);
    if (!s.ok()) break;              // 마이 그레이션 실패하면 break
    s = zbd_->ResetUnusedIOZones();  // 사용되지 않는 IO 존을 재설정
    if (!s.ok()) break;  // 재설정이 실패하면 루프를 중단합니다
  }
  return s;
}
/* 주어진 파일의 익스텐트를 새로운 존으로 마이그레이션하는 작업을 수행 */
// 함수는 여러 단계를 통해 유효 데이터를 새로운 존으로 이동시키고, 파일
// 시스템의 무결성을 유지
IOStatus ZenFS::MigrateFileExtents(
    const std::string& fname,
    const std::vector<ZoneExtentSnapshot*>& migrate_exts) {
  IOStatus s = IOStatus::OK();
  uint64_t copied = 0;
  // 파일 이름과 익스텐트 개수를 로깅
  Info(logger_, "MigrateFileExtents, fname: %s, extent count: %lu",
       fname.data(), migrate_exts.size());

  // The file may be deleted by other threads, better double check.
  // 파일을 가져옵니다
  auto zfile = GetFile(fname);
  if (zfile == nullptr) {
    return IOStatus::OK();
  }

  // Don't migrate open for write files and prevent write reopens while we
  // migrate
  if (!zfile->TryAcquireWRLock()) {
    return IOStatus::OK();
  }
  // 새로운 익스텐트 리스트 생성
  std::vector<ZoneExtent*> new_extent_list;
  std::vector<ZoneExtent*> extents = zfile->GetExtents();
  for (const auto* ext : extents) {
    new_extent_list.push_back(
        new ZoneExtent(ext->start_, ext->length_, ext->zone_));
  }
  // 익스텐트 마이그레이션
  // Modify the new extent list
  for (ZoneExtent* ext : new_extent_list) {
    // Check if current extent need to be migrated
    // 유효 데이터(즉, 마이그레이션할 필요가 있는 익스텐트)를 확인
    auto it = std::find_if(migrate_exts.begin(), migrate_exts.end(),
                           [&](const ZoneExtentSnapshot* ext_snapshot) {
                             return ext_snapshot->start == ext->start_ &&
                                    ext_snapshot->length == ext->length_;
                           });

    if (it == migrate_exts.end()) {
      Info(logger_, "Migrate extent not found, ext_start: %lu", ext->start_);
      continue;
    }

    Zone* target_zone = nullptr;

    // Allocate a new migration zone.
    s = zbd_->TakeMigrateZone(&target_zone, zfile->GetWriteLifeTimeHint(),
                              ext->length_);
    if (!s.ok()) {
      continue;
    }

    if (target_zone == nullptr) {
      zbd_->ReleaseMigrateZone(target_zone);
      Info(logger_, "Migrate Zone Acquire Failed, Ignore Task.");
      continue;
    }

    uint64_t target_start = target_zone->wp_;
    copied += ext->length_;
    if (zfile->IsSparse()) {
      // For buffered write, ZenFS use inlined metadata for extents and each
      // extent has a SPARSE_HEADER_SIZE.
      target_start = target_zone->wp_ + ZoneFile::SPARSE_HEADER_SIZE;
      zfile->MigrateData(ext->start_ - ZoneFile::SPARSE_HEADER_SIZE,
                         ext->length_ + ZoneFile::SPARSE_HEADER_SIZE,
                         target_zone);
      // zbd_->AddGCBytesWritten(ext->length_ + ZoneFile::SPARSE_HEADER_SIZE);
      copied += ZoneFile::SPARSE_HEADER_SIZE;
    } else {
      zfile->MigrateData(ext->start_, ext->length_, target_zone);
      // zbd_->AddGCBytesWritten(ext->length_);
    }

    // If the file doesn't exist, skip
    if (GetFileNoLock(fname) == nullptr) {
      Info(logger_, "Migrate file not exist anymore.");
      zbd_->ReleaseMigrateZone(target_zone);
      break;
    }

    ext->start_ = target_start;
    ext->zone_ = target_zone;
    ext->zone_->used_capacity_ += ext->length_;

    zbd_->ReleaseMigrateZone(target_zone);
  }
  zbd_->AddGCBytesWritten(copied);
  SyncFileExtents(zfile.get(), new_extent_list);
  zfile->ReleaseWRLock();

  Info(logger_, "MigrateFileExtents Finished, fname: %s, extent count: %lu",
       fname.data(), migrate_exts.size());
  return IOStatus::OK();
}

extern "C" FactoryFunc<FileSystem> zenfs_filesystem_reg;

FactoryFunc<FileSystem> zenfs_filesystem_reg =
#if (ROCKSDB_MAJOR < 6) || (ROCKSDB_MAJOR <= 6 && ROCKSDB_MINOR < 28)
    ObjectLibrary::Default()->Register<FileSystem>(
        "zenfs://.*", [](const std::string& uri, std::unique_ptr<FileSystem>* f,
                         std::string* errmsg) {
#else
    ObjectLibrary::Default()->AddFactory<FileSystem>(
        ObjectLibrary::PatternEntry("zenfs", false)
            .AddSeparator("://"), /* "zenfs://.+" */
        [](const std::string& uri, std::unique_ptr<FileSystem>* f,
           std::string* errmsg) {
#endif
          std::string devID = uri;
          FileSystem* fs = nullptr;
          Status s;

          devID.replace(0, strlen("zenfs://"), "");
          if (devID.rfind("dev:") == 0) {
            devID.replace(0, strlen("dev:"), "");
#ifdef ZENFS_EXPORT_PROMETHEUS
            s = NewZenFS(&fs, ZbdBackendType::kBlockDev, devID,
                         std::make_shared<ZenFSPrometheusMetrics>());
#else
            s = NewZenFS(&fs, ZbdBackendType::kBlockDev, devID);
#endif
            if (!s.ok()) {
              *errmsg = s.ToString();
            }
          } else if (devID.rfind("uuid:") == 0) {
            std::map<std::string, std::pair<std::string, ZbdBackendType>>
                zenFileSystems;
            s = ListZenFileSystems(zenFileSystems);
            if (!s.ok()) {
              *errmsg = s.ToString();
            } else {
              devID.replace(0, strlen("uuid:"), "");

              if (zenFileSystems.find(devID) == zenFileSystems.end()) {
                *errmsg = "UUID not found";
              } else {

#ifdef ZENFS_EXPORT_PROMETHEUS
                s = NewZenFS(&fs, zenFileSystems[devID].second,
                             zenFileSystems[devID].first,
                             std::make_shared<ZenFSPrometheusMetrics>());
#else
                s = NewZenFS(&fs, zenFileSystems[devID].second,
                             zenFileSystems[devID].first);
#endif
                if (!s.ok()) {
                  *errmsg = s.ToString();
                }
              }
            }
          } else if (devID.rfind("zonefs:") == 0) {
            devID.replace(0, strlen("zonefs:"), "");
            s = NewZenFS(&fs, ZbdBackendType::kZoneFS, devID);
            if (!s.ok()) {
              *errmsg = s.ToString();
            }
          } else {
            *errmsg = "Malformed URI";
          }
          f->reset(fs);
          return f->get();
        });
};  // namespace ROCKSDB_NAMESPACE

#else

#include "rocksdb/env.h"

namespace ROCKSDB_NAMESPACE {
Status NewZenFS(FileSystem** /*fs*/, const ZbdBackendType /*backend_type*/,
                const std::string& /*backend_name*/,
                ZenFSMetrics* /*metrics*/) {
  return Status::NotSupported("Not built with ZenFS support\n");
}
std::map<std::string, std::string> ListZenFileSystems() {
  std::map<std::string, std::pair<std::string, ZbdBackendType>> zenFileSystems;
  return zenFileSystems;
}
}  // namespace ROCKSDB_NAMESPACE

#endif  // !defined(ROCKSDB_LITE) && defined(OS_LINUX)
