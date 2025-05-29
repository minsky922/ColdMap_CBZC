// Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved.
// Copyright (c) 2019-present, Western Digital Corporation
//  This source code is licensed under both the GPLv2 (found in the
//  COPYING file in the root directory) and Apache 2.0 License
//  (found in the LICENSE.Apache file in the root directory).
/*ZenFS의 주요 데이터 구조와 인터페이스를 정의하는 헤더 파일입니다. ZenFS의
 * 클래스와 함수 선언을 포함하고 있습니다*/
#pragma once

#if __cplusplus < 201703L
#include "filesystem_utility.h"
namespace fs = filesystem_utility;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif

#include <map>
#include <memory>
#include <thread>

#include "io_zenfs.h"
#include "metrics.h"
#include "rocksdb/db.h"
#include "rocksdb/env.h"
#include "rocksdb/file_system.h"
#include "rocksdb/status.h"
#include "snapshot.h"
#include "version.h"
#include "zbd_zenfs.h"

namespace ROCKSDB_NAMESPACE {

#if !defined(ROCKSDB_LITE) && defined(OS_LINUX)

class ZoneSnapshot;
class ZoneFileSnapshot;
class ZenFSSnapshot;
class ZenFSSnapshotOptions;

/* 파일 시스템의 메타데이터 및 설정 정보를 포함*/
class Superblock {
  uint32_t magic_ = 0;
  char uuid_[37] = {0};
  uint32_t sequence_ = 0;
  uint32_t superblock_version_ = 0;
  uint32_t flags_ = 0;
  uint32_t block_size_ = 0; /* in bytes */
  uint32_t zone_size_ = 0;  /* in blocks */
  uint32_t nr_zones_ = 0;   /* 총 존의 수 */
  char aux_fs_path_[256] = {0};
  uint32_t finish_treshold_ = 0;
  char zenfs_version_[64]{0};
  char reserved_[123] = {0};

 public:
  const uint32_t MAGIC = 0x5a454e46; /* ZENF */
  const uint32_t ENCODED_SIZE = 512;
  const uint32_t CURRENT_SUPERBLOCK_VERSION = 2;
  const uint32_t DEFAULT_FLAGS = 0;
  const uint32_t FLAGS_ENABLE_GC =
      1 << 0; /* 가비지 컬렉션(GC)을 활성화하는 플래그 */

  Superblock() {}

  /* Create a superblock for a filesystem covering the entire zoned block device
   */
  Superblock(ZonedBlockDevice* zbd, std::string aux_fs_path = "",
             uint32_t finish_threshold = 0, bool enable_gc = false) {
    std::string uuid = Env::Default()->GenerateUniqueId();
    int uuid_len =
        std::min(uuid.length(),
                 sizeof(uuid_) - 1); /* make sure uuid is nullterminated */
    memcpy((void*)uuid_, uuid.c_str(), uuid_len);
    magic_ = MAGIC;
    superblock_version_ = CURRENT_SUPERBLOCK_VERSION;
    flags_ = DEFAULT_FLAGS;
    if (enable_gc) flags_ |= FLAGS_ENABLE_GC;

    finish_treshold_ = finish_threshold;

    block_size_ = zbd->GetBlockSize();
    zone_size_ = zbd->GetZoneSize() / block_size_;
    nr_zones_ = zbd->GetNrZones();

    strncpy(aux_fs_path_, aux_fs_path.c_str(), sizeof(aux_fs_path_) - 1);

    std::string zenfs_version = ZENFS_VERSION;
    strncpy(zenfs_version_, zenfs_version.c_str(), sizeof(zenfs_version_) - 1);
  }

  Status DecodeFrom(Slice* input);
  void EncodeTo(std::string* output);
  Status CompatibleWith(ZonedBlockDevice* zbd);

  void GetReport(std::string* reportString);

  uint32_t GetSeq() { return sequence_; }
  std::string GetAuxFsPath() { return std::string(aux_fs_path_); }
  uint32_t GetFinishTreshold() { return finish_treshold_; }
  std::string GetUUID() { return std::string(uuid_); }
  bool IsGCEnabled() {
    return flags_ & FLAGS_ENABLE_GC;
  }; /* 가비지 컬렉션이 활성화되었는지 여부를 확인합니다.*/
};

class ZenMetaLog {
  uint64_t read_pos_;
  Zone* zone_;
  ZonedBlockDevice* zbd_;
  size_t bs_;

  /* Every meta log record is prefixed with a CRC(32 bits) and record length (32
   * bits) */
  const size_t zMetaHeaderSize = sizeof(uint32_t) * 2;

 public:
  ZenMetaLog(ZonedBlockDevice* zbd, Zone* zone) {
    // assert(zone->IsBusy());
    zbd_ = zbd;
    zone_ = zone;
    bs_ = zbd_->GetBlockSize();
    read_pos_ = zone->start_;
  }

  virtual ~ZenMetaLog() {
    // TODO: report async error status
    bool ok = zone_->Release();
    assert(ok);
    (void)ok;
  }

  IOStatus AddRecord(const Slice& slice);
  IOStatus ReadRecord(Slice* record, std::string* scratch);

  Zone* GetZone() { return zone_; };

 private:
  IOStatus Read(Slice* slice);
};

class ZenFS : public FileSystemWrapper {
  ZonedBlockDevice* zbd_;  // 존 기반 블록 디바이스
  std::map<std::string, std::shared_ptr<ZoneFile>> files_;  // 파일 맵
  std::mutex files_mtx_;                // 파일 맵 접근 제어 뮤텍스
  std::shared_ptr<Logger> logger_;      // 로거 객체
  std::atomic<uint64_t> next_file_id_;  // 다음 파일 ID

  Zone* cur_meta_zone_ = nullptr;           // 현재 메타 존
  std::unique_ptr<ZenMetaLog> meta_log_;    // 메타 로그 객체
  std::mutex metadata_sync_mtx_;            // 메타데이터 동기화 뮤텍스
  std::unique_ptr<Superblock> superblock_;  // 슈퍼블록 객체

  std::shared_ptr<Logger> GetLogger() { return logger_; }  // 로거 반환

  std::unique_ptr<std::thread> gc_worker_ = nullptr;  // 가비지 컬렉션 스레드
  //
  std::unique_ptr<std::thread> bg_reset_worker_ =
      nullptr;  // 매초마다 빈공간 계산하는 백그라운드 스레드
  //
  uint64_t free_percent_ = 100;
  bool run_gc_worker_ = false;  // 가비지 컬렉션 작업 실행 여부
  //
  bool run_bg_reset_worker_ = false;
  //
  std::atomic<int> zc_triggerd_count_{0};
  std::mutex zc_lock_;

  std::atomic<int> mount_time_{0};

  // bool coldest_type_set_  =false;
  // int coldest_type_ = -1;
  // std::mutex coldest_type_lock_;

  // std::atomic<int> file_operation_sequence_{0};

  uint64_t file_size_dist[5];

  DB* db_ptr_ = nullptr;

  struct ZenFSMetadataWriter : public MetadataWriter {
    ZenFS* zenFS;
    IOStatus Persist(ZoneFile* zoneFile) {
      Debug(zenFS->GetLogger(), "Syncing metadata for: %s",
            zoneFile->GetFilename().c_str());
      return zenFS->SyncFileMetadata(zoneFile);
    }
  };

  ZenFSMetadataWriter metadata_writer_;

  enum ZenFSTag : uint32_t {
    kCompleteFilesSnapshot = 1,
    kFileUpdate = 2,
    kFileDeletion = 3,
    kEndRecord = 4,
    kFileReplace = 5,
  };

  void LogFiles();
  void ClearFiles();
  std::string FormatPathLexically(fs::path filepath);
  IOStatus WriteSnapshotLocked(ZenMetaLog* meta_log);
  IOStatus WriteEndRecord(ZenMetaLog* meta_log);
  IOStatus RollMetaZoneLocked();
  IOStatus PersistSnapshot(ZenMetaLog* meta_writer);
  IOStatus PersistRecord(std::string record);
  IOStatus SyncFileExtents(ZoneFile* zoneFile,
                           std::vector<ZoneExtent*> new_extents);
  /* Must hold files_mtx_ */
  IOStatus SyncFileMetadataNoLock(ZoneFile* zoneFile, bool replace = false);
  /* Must hold files_mtx_ */
  IOStatus SyncFileMetadataNoLock(std::shared_ptr<ZoneFile> zoneFile,
                                  bool replace = false) {
    return SyncFileMetadataNoLock(zoneFile.get(), replace);
  }
  IOStatus SyncFileMetadata(ZoneFile* zoneFile, bool replace = false);
  IOStatus SyncFileMetadata(std::shared_ptr<ZoneFile> zoneFile,
                            bool replace = false) {
    return SyncFileMetadata(zoneFile.get(), replace);
  }

  void EncodeSnapshotTo(std::string* output);
  void EncodeFileDeletionTo(std::shared_ptr<ZoneFile> zoneFile,
                            std::string* output, std::string linkf);

  Status DecodeSnapshotFrom(Slice* input);
  Status DecodeFileUpdateFrom(Slice* slice, bool replace = false);
  Status DecodeFileDeletionFrom(Slice* slice);

  Status RecoverFrom(ZenMetaLog* log);

  std::string ToAuxPath(std::string path) {
    return superblock_->GetAuxFsPath() + path;
  }

  std::string ToZenFSPath(std::string aux_path) {
    std::string path = aux_path;
    path.erase(0, superblock_->GetAuxFsPath().length());
    return path;
  }

  /* Must hold files_mtx_ */
  std::shared_ptr<ZoneFile> GetFileNoLock(std::string fname);
  /* Must hold files_mtx_ */
  void GetZenFSChildrenNoLock(const std::string& dir,
                              bool include_grandchildren,
                              std::vector<std::string>* result);
  /* Must hold files_mtx_ */
  IOStatus GetChildrenNoLock(const std::string& dir, const IOOptions& options,
                             std::vector<std::string>* result,
                             IODebugContext* dbg);

  /* Must hold files_mtx_ */
  IOStatus RenameChildNoLock(std::string const& source_dir,
                             std::string const& dest_dir,
                             std::string const& child, const IOOptions& options,
                             IODebugContext* dbg);

  /* Must hold files_mtx_ */
  IOStatus RollbackAuxDirRenameNoLock(
      const std::string& source_path, const std::string& dest_path,
      const std::vector<std::string>& renamed_children,
      const IOOptions& options, IODebugContext* dbg);

  /* Must hold files_mtx_ */
  IOStatus RenameAuxPathNoLock(const std::string& source_path,
                               const std::string& dest_path,
                               const IOOptions& options, IODebugContext* dbg);

  /* Must hold files_mtx_ */
  IOStatus RenameFileNoLock(const std::string& f, const std::string& t,
                            const IOOptions& options, IODebugContext* dbg);

  std::shared_ptr<ZoneFile> GetFile(std::string fname);

  /* Must hold files_mtx_, On successful return,
   * caller must release files_mtx_ and call ResetUnusedIOZones() */
  IOStatus DeleteFileNoLock(std::string fname, const IOOptions& options,
                            IODebugContext* dbg);

  IOStatus Repair();

  /* Must hold files_mtx_ */
  IOStatus DeleteDirRecursiveNoLock(const std::string& d,
                                    const IOOptions& options,
                                    IODebugContext* dbg);

  /* Must hold files_mtx_ */
  IOStatus IsDirectoryNoLock(const std::string& path, const IOOptions& options,
                             bool* is_dir, IODebugContext* dbg) {
    if (GetFileNoLock(path) != nullptr) {
      *is_dir = false;
      return IOStatus::OK();
    }
    return target()->IsDirectory(ToAuxPath(path), options, is_dir, dbg);
  }

 protected:
  IOStatus OpenWritableFile(const std::string& fname,
                            const FileOptions& file_opts,
                            std::unique_ptr<FSWritableFile>* result,
                            IODebugContext* dbg, bool reopen);

 public:
  explicit ZenFS(ZonedBlockDevice* zbd, std::shared_ptr<FileSystem> aux_fs,
                 std::shared_ptr<Logger> logger);
  virtual ~ZenFS();

  Status Mount(bool readonly);
  Status MkFS(std::string aux_fs_path, uint32_t finish_threshold,
              bool enable_gc);
  std::map<std::string, Env::WriteLifeTimeHint> GetWriteLifeTimeHints();

  const char* Name() const override {
    return "ZenFS - The Zoned-enabled File System";
  }

  void EncodeJson(std::ostream& json_stream);

  void ReportSuperblock(std::string* report) { superblock_->GetReport(report); }

  virtual IOStatus NewSequentialFile(const std::string& fname,
                                     const FileOptions& file_opts,
                                     std::unique_ptr<FSSequentialFile>* result,
                                     IODebugContext* dbg) override;
  virtual IOStatus NewRandomAccessFile(
      const std::string& fname, const FileOptions& file_opts,
      std::unique_ptr<FSRandomAccessFile>* result,
      IODebugContext* dbg) override;
  virtual IOStatus NewWritableFile(const std::string& fname,
                                   const FileOptions& file_opts,
                                   std::unique_ptr<FSWritableFile>* result,
                                   IODebugContext* dbg) override;
  virtual IOStatus ReuseWritableFile(const std::string& fname,
                                     const std::string& old_fname,
                                     const FileOptions& file_opts,
                                     std::unique_ptr<FSWritableFile>* result,
                                     IODebugContext* dbg) override;
  virtual IOStatus ReopenWritableFile(const std::string& fname,
                                      const FileOptions& options,
                                      std::unique_ptr<FSWritableFile>* result,
                                      IODebugContext* dbg) override;
  virtual IOStatus FileExists(const std::string& fname,
                              const IOOptions& options,
                              IODebugContext* dbg) override;
  virtual IOStatus GetChildren(const std::string& dir, const IOOptions& options,
                               std::vector<std::string>* result,
                               IODebugContext* dbg) override;
  virtual IOStatus DeleteFile(const std::string& fname,
                              const IOOptions& options,
                              IODebugContext* dbg) override;
  virtual IOStatus LinkFile(const std::string& fname, const std::string& lname,
                            const IOOptions& options,
                            IODebugContext* dbg) override;
  virtual IOStatus NumFileLinks(const std::string& fname,
                                const IOOptions& options, uint64_t* nr_links,
                                IODebugContext* dbg) override;
  virtual IOStatus AreFilesSame(const std::string& fname,
                                const std::string& link,
                                const IOOptions& options, bool* res,
                                IODebugContext* dbg) override;

  IOStatus GetFileSize(const std::string& f, const IOOptions& options,
                       uint64_t* size, IODebugContext* dbg) override;
  IOStatus RenameFile(const std::string& f, const std::string& t,
                      const IOOptions& options, IODebugContext* dbg) override;

  IOStatus GetFreeSpace(const std::string& /*path*/,
                        const IOOptions& /*options*/, uint64_t* diskfree,
                        uint64_t* free_percent,
                        IODebugContext* /*dbg*/) override {
    if (diskfree != nullptr) {
      *diskfree = zbd_->GetFreeSpace();
    } else {  // if nullptr
      goto ret;
    }
    if (zbd_ && zbd_->GetZCRunning()) {
      while (zbd_->GetZCRunning()) {
        // std::cout << "GetFreeSpace->free_percent_: " << free_percent_
        //           << std::endl;
      }
    }
  ret:
    *free_percent = free_percent_;

    return IOStatus::OK();
  }

  IOStatus GetFileModificationTime(const std::string& fname,
                                   const IOOptions& options, uint64_t* mtime,
                                   IODebugContext* dbg) override;

  // The directory structure is stored in the aux file system

  IOStatus IsDirectory(const std::string& path, const IOOptions& options,
                       bool* is_dir, IODebugContext* dbg) override {
    std::lock_guard<std::mutex> lock(files_mtx_);
    return IsDirectoryNoLock(path, options, is_dir, dbg);
  }

  IOStatus NewDirectory(const std::string& name, const IOOptions& io_opts,
                        std::unique_ptr<FSDirectory>* result,
                        IODebugContext* dbg) override {
    Debug(logger_, "NewDirectory: %s to aux: %s\n", name.c_str(),
          ToAuxPath(name).c_str());
    return target()->NewDirectory(ToAuxPath(name), io_opts, result, dbg);
  }

  IOStatus CreateDir(const std::string& d, const IOOptions& options,
                     IODebugContext* dbg) override {
    Debug(logger_, "CreatDir: %s to aux: %s\n", d.c_str(),
          ToAuxPath(d).c_str());
    return target()->CreateDir(ToAuxPath(d), options, dbg);
  }

  IOStatus CreateDirIfMissing(const std::string& d, const IOOptions& options,
                              IODebugContext* dbg) override {
    Debug(logger_, "CreatDirIfMissing: %s to aux: %s\n", d.c_str(),
          ToAuxPath(d).c_str());
    return target()->CreateDirIfMissing(ToAuxPath(d), options, dbg);
  }

  IOStatus DeleteDir(const std::string& d, const IOOptions& options,
                     IODebugContext* dbg) override {
    std::vector<std::string> children;
    IOStatus s;

    Debug(logger_, "DeleteDir: %s aux: %s\n", d.c_str(), ToAuxPath(d).c_str());

    s = GetChildren(d, options, &children, dbg);
    if (children.size() != 0)
      return IOStatus::IOError("Directory has children");

    return target()->DeleteDir(ToAuxPath(d), options, dbg);
  }

  IOStatus DeleteDirRecursive(const std::string& d, const IOOptions& options,
                              IODebugContext* dbg);

  // We might want to override these in the future
  IOStatus GetAbsolutePath(const std::string& db_path, const IOOptions& options,
                           std::string* output_path,
                           IODebugContext* dbg) override {
    return target()->GetAbsolutePath(ToAuxPath(db_path), options, output_path,
                                     dbg);
  }

  IOStatus LockFile(const std::string& f, const IOOptions& options,
                    FileLock** l, IODebugContext* dbg) override {
    return target()->LockFile(ToAuxPath(f), options, l, dbg);
  }

  IOStatus UnlockFile(FileLock* l, const IOOptions& options,
                      IODebugContext* dbg) override {
    return target()->UnlockFile(l, options, dbg);
  }

  IOStatus GetTestDirectory(const IOOptions& options, std::string* path,
                            IODebugContext* dbg) override {
    *path = "rocksdbtest";
    Debug(logger_, "GetTestDirectory: %s aux: %s\n", path->c_str(),
          ToAuxPath(*path).c_str());
    return target()->CreateDirIfMissing(ToAuxPath(*path), options, dbg);
  }

  IOStatus NewLogger(const std::string& fname, const IOOptions& options,
                     std::shared_ptr<Logger>* result,
                     IODebugContext* dbg) override {
    return target()->NewLogger(ToAuxPath(fname), options, result, dbg);
  }

  // Not supported (at least not yet)
  IOStatus Truncate(const std::string& /*fname*/, size_t /*size*/,
                    const IOOptions& /*options*/,
                    IODebugContext* /*dbg*/) override {
    return IOStatus::NotSupported("Truncate is not implemented in ZenFS");
  }

  virtual IOStatus NewRandomRWFile(const std::string& /*fname*/,
                                   const FileOptions& /*options*/,
                                   std::unique_ptr<FSRandomRWFile>* /*result*/,
                                   IODebugContext* /*dbg*/) override {
    return IOStatus::NotSupported("RandomRWFile is not implemented in ZenFS");
  }

  virtual IOStatus NewMemoryMappedFileBuffer(
      const std::string& /*fname*/,
      std::unique_ptr<MemoryMappedFileBuffer>* /*result*/) override {
    return IOStatus::NotSupported(
        "MemoryMappedFileBuffer is not implemented in ZenFS");
  }
  //////////////////////////////////////
  void SetDBPtr(DB* ptr) override {
    db_ptr_ = ptr;
    zbd_->SetDBPtr(ptr);
  }
  void SetResetScheme(uint32_t r, uint32_t partial_reset_scheme, uint64_t T,
                      uint64_t zc, uint64_t until, uint64_t allocation_scheme,
                      uint64_t zc_scheme, double alpha_value,
                      double sigma_value, uint64_t finish_scheme,
                      uint64_t predict_cnt,
                      std::vector<uint64_t>& other_options) override;
  void GiveZenFStoLSMTreeHint(
      std::vector<uint64_t>& compaction_inputs_input_level_fno,
      std::vector<uint64_t>& compaction_inputs_output_level_fno,
      int output_level, bool trivial_move) override {


    int seq = zbd_->file_operation_sequence_.fetch_add(1);
    if(!trivial_move){
      int type=0;
      switch (output_level)
      {
      case 0:
        zbd_->latest_file_operation_sequence_[SeqL0L1andFlush] = seq;
        type=SeqL0L1andFlush;
        break;
      case 1:
        zbd_->latest_file_operation_sequence_[SeqL0L1andFlush] = seq;
        type=SeqL0L1andFlush;
        break;
      case 2:
        zbd_->latest_file_operation_sequence_[SeqL1L2] = seq;
        type=SeqL1L2;
        break;
      case 3:
        zbd_->latest_file_operation_sequence_[SeqL2L3] = seq;
        type=SeqL2L3;
        break;
      case 4:
        zbd_->latest_file_operation_sequence_[SeqL3L4] = seq;
        type=SeqL3L4;
        break;
      default:
        break;
      }

      {
        
        if(zbd_->coldest_type_set_== true){
          std::lock_guard<std::mutex> lg(zbd_->coldest_type_lock_);
          // todo
          bool ok = true;
          zbd_->check_coldest_[type]=true;
          for(int i =0;i<10;i++){
            if(zbd_->latest_file_operation_sequence_[i]==0){
              continue;
            }
            if(zbd_->check_coldest_[i]==true){
              continue;
            }
            ok=false;
            break;
          }


          if(ok==true){
            if(zbd_->coldest_type_!=type){
              zbd_->CBSC_mispredict_stats_[zbd_->coldest_type_].fetch_add(1);
            }
            zbd_->CBSC_total_predict_stats_[zbd_->coldest_type_].fetch_add(1);
            zbd_->coldest_type_set_=false;
          }
          else{
            if(zbd_->coldest_type_==type){
              zbd_->CBSC_mispredict_stats_[zbd_->coldest_type_].fetch_add(1);
              zbd_->CBSC_total_predict_stats_[zbd_->coldest_type_].fetch_add(1);
              zbd_->coldest_type_set_=false;
            }
          }

        }
      }
    }


    zbd_->GiveZenFStoLSMTreeHint(compaction_inputs_input_level_fno,
                                 compaction_inputs_output_level_fno,
                                 output_level, trivial_move);
  }
  double GetMaxInvalidateCompactionScore(std::vector<uint64_t>& file_candidates,
                                         uint64_t* candidate_size) override;
  bool IsZoneDevice() { return true; }
  // void ZoneCleaningWorker(bool run_once=false) override;
  void ZoneCleaning(bool forced);
  void ReCalculateLifetimes();
  void Adv_ReCalculateLifetimes();
  void NormalizeZoneLifetimes();
  double CalculateZoneLifetimeVariance();
  void CalculateHorizontalLifetimes(
      std::map<int, std::vector<std::pair<uint64_t, double>>>& level_file_map);
  struct FileInfo_ {
    uint64_t fno;
    double horizontal_lifetime;
    bool is_compacting;
    bool is_trivial;
    double sst_lifetime_value_;
  };
  struct ZoneLifetimeData {
    double total_lifetime;
    int file_count;
    // fno -> 해당 파일의 lifetime
    std::map<uint64_t, double> file_lifetimes;

    ZoneLifetimeData() : total_lifetime(0.0), file_count(0) {}
  };
  std::map<int, std::vector<FileInfo_>> level_file_map_;
  std::map<uint64_t, ZoneLifetimeData> zone_lifetime_map_;
  // std::unordered_map<int, std::vector<FileInfo_>> simulated_file_map;
  std::set<uint64_t> fno_already_propagated;
  std::set<uint64_t> fno_not_should_selected_as_pivot_again;

  void PredictCompaction(int step);

  uint64_t GetMaxLevelScoreLevel(std::array<uint64_t, 10>& tmp_lsm_tree,
                                 int initial_l0_files_n,
                                 std::unordered_set<int>& excluded_levels);

  uint64_t GetMaxHorizontalFno(int pivot_level);
  void PredictCompactionImpl(uint64_t& pivot_level,
                             std::array<uint64_t, 10>& tmp_lsm_tree,
                             uint64_t& pivot_fno,
                             std::vector<uint64_t>& unpivot_fno_list,
                             int initial_l0_files_n);
  void GetOverlappingFno(uint64_t pivot_fno, uint64_t pivot_level,
                         std::vector<uint64_t>& unpivot_fno_list);
  void Propagation(uint64_t pivot_fno, std::vector<uint64_t> unpivot_fno_list);

  void CalculateHorizontalLifetimes(
      std::map<int, std::vector<FileInfo_>>& level_file_map);
  int GetMountTime(void) override { return mount_time_.load(); }
  // bool IsZCRunning(void) { return run_gc_worker_; }
  void ZCLock(void) override { zc_lock_.lock(); }
  void ZCUnLock(void) override { zc_lock_.unlock(); }
  void BackgroundStatTimeLapse();
  uint64_t EstimateFileAge(Env::WriteLifeTimeHint hint);
  // uint64_t GetRecentModificationTime(ZenFSZone& zone);
  //////////////////////////////////
  void GetZenFSSnapshot(ZenFSSnapshot& snapshot,
                        const ZenFSSnapshotOptions& options);

  IOStatus MigrateExtents(const std::vector<ZoneExtentSnapshot*>& extents);

  IOStatus MigrateFileExtents(
      const std::string& fname,
      const std::vector<ZoneExtentSnapshot*>& migrate_exts);

  // std::map<uint64_t, std::tuple<double, int, std::vector<double>>>
  //     zone_lifetime_map_;

 private:
  // std::map<uint64_t, std::pair<double, int>> zone_lifetime_map_;

  const uint64_t GC_START_LEVEL =
      20;                      /* Enable GC when < 20% free space available */
  const uint64_t GC_SLOPE = 3; /* GC agressiveness */
  void GCWorker();
};
#endif  // !defined(ROCKSDB_LITE) && defined(OS_LINUX)

Status NewZenFS(
    FileSystem** fs, const std::string& bdevname,
    std::shared_ptr<ZenFSMetrics> metrics = std::make_shared<NoZenFSMetrics>());
Status NewZenFS(
    FileSystem** fs, const ZbdBackendType backend_type,
    const std::string& backend_name,
    std::shared_ptr<ZenFSMetrics> metrics = std::make_shared<NoZenFSMetrics>());
Status AppendZenFileSystem(
    std::string path, ZbdBackendType backend,
    std::map<std::string, std::pair<std::string, ZbdBackendType>>& fs_list);
Status ListZenFileSystems(
    std::map<std::string, std::pair<std::string, ZbdBackendType>>& out_list);

}  // namespace ROCKSDB_NAMESPACE
