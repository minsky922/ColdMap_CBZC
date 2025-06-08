// Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved.
// Copyright (c) 2019-present, Western Digital Corporation
//  This source code is licensed under both the GPLv2 (found in the
//  COPYING file in the root directory) and Apache 2.0 License
//  (found in the LICENSE.Apache file in the root directory).

#pragma once

#include <cstdint>
#if !defined(ROCKSDB_LITE) && defined(OS_LINUX)

#include <errno.h>
#include <libzbd/zbd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <iostream>
#include <mutex>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "metrics.h"
#include "rocksdb/comparator.h"
#include "rocksdb/db.h"
#include "rocksdb/env.h"
#include "rocksdb/file_system.h"
#include "rocksdb/io_status.h"
namespace ROCKSDB_NAMESPACE {

class ZonedBlockDevice;
class ZonedBlockDeviceBackend;
class ZoneSnapshot;
class ZenFSSnapshotOptions;
class Zone;
class ZoneFile;

// #define ZONE_CLEANING_KICKING_POINT (20)
#define ZEU_SIZE = 128 << 20;
#define READ_FD 0
#define READ_DIRECT_FD 1
#define WRITE_DIRECT_FD 2
#define KB (1024)

#define MB (1024 * KB)

#define ZENFS_SPARE_ZONES (0)

#define ZENFS_META_ZONES (3)
#define log2_DEVICE_IO_CAPACITY (6)  // 64GB

#define ZC_COMPACTION_IO_PRIORITY IOPRIO_PRIO_VALUE(IOPRIO_CLASS_RT, 4)

#define READ_DISK_COST 0
#define READ_PAGE_COST 1
#define WRITE_COST 2
#define FREE_SPACE_COST 3
#define SEQ_DIST_MAX 5000
// #define ZENFS_IO_ZONES (80)

#define ZONE_SIZE 1024

#define DEVICE_SIZE ((ZENFS_IO_ZONES) * (ZONE_SIZE))

// #define ZONE_SIZE_PER_DEVICE_SIZE (100 / (ZENFS_IO_ZONES))

#define WP_TO_RELATIVE_WP(wp, zone_sz, zidx) ((wp) - (zone_sz * zidx))

#define BYTES_TO_MB(bytes) (bytes >> 20)

#define PASS 0
#define SPINLOCK 1

#define LIZA 0
#define CAZA 1
#define CAZA_ADV 2

#define GREEDY 0
#define CBZC1 1
#define CBZC2 2
#define CBZC3 3
#define CBZC4 4
#define CBZC5 5
#define CBZC6 6
#define PARTIAL_RESET_KICKED_THRESHOLD 40
// | zone-reset | partial-reset |
#define RUNTIME_ZONE_RESET_DISABLED 0    // |      x     |       x       |
#define RUNTIME_ZONE_RESET_ONLY 1        // |      o     |       x       |
#define PARTIAL_RESET_WITH_ZONE_RESET 2  // |      o     |       o       |
#define PARTIAL_RESET_ONLY 3             // |      x     |       o       |
#define PARTIAL_RESET_AT_BACKGROUND 4    // |      ?     |       o       |
#define PARTIAL_RESET_BACKGROUND_T_WITH_ZONE_RESET \
  5                               // |      o     |       o       |
#define PROACTIVE_ZONECLEANING 6  // |      x     |       x       |

#define IS_BIG_SSTABLE(file_size) \
  (bool)((uint64_t)(file_size) > (uint64_t)(63 << 20))

#define FINISH_ENABLE 0
#define FINISH_DISABLE 1
#define FINISH_PROPOSAL 2
enum{
  SeqWAL=0,
  SeqL0L1andFlush=1,
  SeqL1L2=2,
  SeqL2L3=3,
  SeqL3L4=4,
  SeqL4L5=5,
};
enum WaitForOpenZoneClass {
  WAL = 0,
  ZC = 1,
  L0 = 2,
  L1 = 3,
  L2 = 4,
  L3 = 5,
  L4 = 6,
  L5 = 7,
  L6 = 8,
  L7 = 9
};
class ZoneList {
 private:
  void *data_;
  unsigned int zone_count_;

 public:
  ZoneList(void *data, unsigned int zone_count)
      : data_(data), zone_count_(zone_count) {};
  void *GetData() { return data_; };
  unsigned int ZoneCount() { return zone_count_; };
  ~ZoneList() { free(data_); };
};

// inline bool ends_with(std::string const &value, std::string const &ending) {
//   if (ending.size() > value.size()) return false;
//   return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
// }

struct start_end_to_be_sorted {
  long start_milliseconds;
  long end_milliseconds;
  bool zcied;
  static bool compare_start_end_to_be_sorted(const start_end_to_be_sorted &a,
                                             const start_end_to_be_sorted &b) {
    if (a.start_milliseconds == b.start_milliseconds) {
      return a.end_milliseconds < b.end_milliseconds;
    }
    return a.start_milliseconds < b.start_milliseconds;
  }
};

struct motivation_lifetime_diff {
  struct timespec file_created_time_;
  struct timespec file_deleted_time_;
  bool zcied_ = false;
};

struct motivation_zone_lifetime_diff {
  struct timespec zone_allocated_time_;
  struct timespec zone_resetted_time_;
  struct std::vector<motivation_lifetime_diff> motivation_lifetime_diffs;
};

class Zone {
  ZonedBlockDevice *zbd_;
  ZonedBlockDeviceBackend *zbd_be_;
  std::atomic_bool busy_;

 public:
  explicit Zone(ZonedBlockDevice *zbd, ZonedBlockDeviceBackend *zbd_be,
                std::unique_ptr<ZoneList> &zones, unsigned int idx);

  uint64_t start_;
  uint64_t capacity_; /* remaining capacity */
  uint64_t max_capacity_;
  uint64_t wp_;
  Env::WriteLifeTimeHint lifetime_;
  //
  uint64_t zidx_;  // not changed
  uint64_t zone_sz_;
  unsigned int log2_erase_unit_size_ = 0;
  uint64_t erase_unit_size_ = 0;
  uint64_t block_sz_;
  bool is_finished_ = false;
  uint64_t reset_count_ = 0;
  uint64_t finish_count_ = 0;
  enum State { EMPTY, OPEN, CLOSE, FINISH, RO, OFFLINE };
  State state_ = EMPTY;
  std::chrono::time_point<std::chrono::system_clock> recent_inval_time_;
  //
  std::atomic<uint64_t> used_capacity_;

  bool is_allocated_ = false;
  struct timespec allocated_time_;
  bool this_zone_motivation_check_ = false;

  struct std::vector<motivation_lifetime_diff> motivation_lifetime_diffs;
  std::mutex motivation_lifetime_diffs_lock_;

  IOStatus Reset();
  IOStatus Finish();
  IOStatus Close();

  IOStatus Append(char *data, uint32_t size);
  bool IsUsed();
  bool IsFull();
  bool IsEmpty();
  bool IsFinished();
  uint64_t GetZoneNr();
  uint64_t GetCapacityLeft();
  bool IsBusy() { return this->busy_.load(std::memory_order_relaxed); }
  bool Acquire() {
    bool expected = false;
    return this->busy_.compare_exchange_strong(expected, true,
                                               std::memory_order_acq_rel);
  }
  bool Release() {
    bool expected = true;
    return this->busy_.compare_exchange_strong(expected, false,
                                               std::memory_order_acq_rel);
  }

  void EncodeJson(std::ostream &json_stream);

  inline IOStatus CheckRelease();
};

class ZonedBlockDeviceBackend {
 public:
  uint32_t block_sz_ = 0;
  uint64_t zone_sz_ = 0;
  uint32_t nr_zones_ = 0;

 public:
  virtual IOStatus Open(bool readonly, bool exclusive,
                        unsigned int *max_active_zones,
                        unsigned int *max_open_zones) = 0;

  virtual std::unique_ptr<ZoneList> ListZones() = 0;
  virtual IOStatus Reset(uint64_t start, bool *offline,
                         uint64_t *max_capacity) = 0;
  virtual IOStatus Finish(uint64_t start) = 0;
  virtual IOStatus Close(uint64_t start) = 0;
  virtual int Read(char *buf, int size, uint64_t pos, bool direct) = 0;
  virtual int Write(char *data, uint32_t size, uint64_t pos) = 0;
  virtual int InvalidateCache(uint64_t pos, uint64_t size) = 0;
  virtual bool ZoneIsSwr(std::unique_ptr<ZoneList> &zones,
                         unsigned int idx) = 0;
  virtual bool ZoneIsOffline(std::unique_ptr<ZoneList> &zones,
                             unsigned int idx) = 0;
  virtual bool ZoneIsWritable(std::unique_ptr<ZoneList> &zones,
                              unsigned int idx) = 0;
  virtual bool ZoneIsActive(std::unique_ptr<ZoneList> &zones,
                            unsigned int idx) = 0;
  virtual bool ZoneIsOpen(std::unique_ptr<ZoneList> &zones,
                          unsigned int idx) = 0;
  virtual uint64_t ZoneStart(std::unique_ptr<ZoneList> &zones,
                             unsigned int idx) = 0;
  virtual uint64_t ZoneMaxCapacity(std::unique_ptr<ZoneList> &zones,
                                   unsigned int idx) = 0;
  virtual uint64_t ZoneWp(std::unique_ptr<ZoneList> &zones,
                          unsigned int idx) = 0;
  virtual std::string GetFilename() = 0;
  uint32_t GetBlockSize() { return block_sz_; };
  //  uint64_t GetBlockSize() { return 4096; };
  uint64_t GetZoneSize() { return zone_sz_; };
  uint32_t GetNrZones() { return nr_zones_; };
  virtual ~ZonedBlockDeviceBackend() {};
};

enum class ZbdBackendType {
  kBlockDev,
  kZoneFS,
};

class ZonedBlockDevice {
 private:
  FileSystemWrapper *zenfs_;
  std::unique_ptr<ZonedBlockDeviceBackend> zbd_be_;
  std::vector<Zone *> io_zones;
  std::vector<Zone *> meta_zones;
  time_t start_time_;
  std::shared_ptr<Logger> logger_;
  uint32_t finish_threshold_ = 0;
  //
  uint64_t zone_sz_;
  //
  /* FAR STATS*/
  std::atomic<uint64_t> bytes_written_{0};
  std::atomic<uint64_t> gc_bytes_written_{0};
  uint64_t zc_read_amp_ = 0;
  std::atomic<bool> force_zc_should_triggered_{false};
  uint64_t reset_threshold_ = 0;
  uint64_t reset_threshold_arr_[101];

  uint64_t finish_threshold_arr_[101];
  ///

  ///
  std::atomic<size_t> reset_count_{0};
  std::atomic<size_t> reset_count_before_full_{0};
  std::atomic<size_t> reset_size_before_full_{0};
  std::atomic<size_t> finish_count_{0};

  std::atomic<size_t> erase_size_{0};
  std::atomic<size_t> erase_size_zc_{0};
  std::atomic<size_t> erase_size_proactive_zc_{0};
  std::atomic<size_t> partial_erase_size_{0};
  std::atomic<size_t> partial_to_valid_erase_size_{0};

  std::atomic<size_t> partial_reset_count_{0};
  std::atomic<size_t> partial_reset_to_valid_count_{0};
  std::atomic<size_t> reset_count_zc_{0};
  std::atomic<size_t> reset_resigned_{0};
  std::atomic<size_t> reset_tried_{0};
  std::atomic<size_t> ratio_sum_{0};
  std::atomic<size_t> ratio_sum2_{0};

  size_t candidate_ratio_sum_ = 0;
  size_t candidate_valid_ratio_sum_ = 0;
  size_t no_candidate_valid_ratio_sum_ = 0;

  size_t candidate_ratio_sum_before_zc_ = 0;
  size_t candidate_valid_ratio_sum_before_zc_ = 0;
  size_t no_candidate_valid_ratio_sum_before_zc_ = 0;
  size_t before_zc_T_ = -1;
  bool before_zc_ = true;
  std::atomic<uint64_t> read_latency_sum_{0};
  std::atomic<uint64_t> read_n_{0};

  std::atomic<uint64_t> wasted_wp_{0};
  std::atomic<uint64_t> new_wasted_wp_{0};
  std::atomic<uint64_t> finished_wasted_wp_{0};
  std::atomic<uint64_t> finished_valid_data_{0};
  std::atomic<clock_t> runtime_reset_reset_latency_{0};
  std::atomic<clock_t> runtime_reset_latency_{0};


  std::atomic<uint64_t> predict_right_vertical_{0};
  std::atomic<uint64_t> predict_false_vertical_{0};
  std::atomic<uint64_t> predict_right_horizontal_{0};
  std::atomic<uint64_t> predict_false_horizontal_{0};

  std::atomic<uint64_t> device_free_space_;

  std::mutex compaction_refused_lock_;
  std::atomic<int> compaction_refused_by_zone_interface_{0};
  std::set<int> compaction_blocked_at_;
  std::vector<int> compaction_blocked_at_amount_;

  /* time_lapse */
  // int zc_io_block_ = 0;
  int zone_cleaning_io_block_ = 0;
  clock_t ZBD_mount_time_;
  bool zone_allocation_state_ = true;

  DB *db_ptr_;

  ZoneFile **sst_file_bitmap_;
  struct ZCStat {
    size_t zc_z;
    int s;
    int e;
    long long us;
    size_t copied;
    bool forced;
    double zlv;
    uint64_t invalid_data_size;
    uint64_t valid_data_size;
    uint64_t invalid_ratio;
    uint64_t page_cache_hit_size;
  };
  std::vector<ZCStat> zc_timelapse_;
  std::vector<uint64_t> zc_copied_timelapse_;

  std::mutex io_lock_;
  struct IOBlockStat {
    pid_t tid;
    int s;
    int e;
  };
  std::vector<IOBlockStat> io_block_timelapse_;
  int io_blocked_thread_n_ = 0;

  /* Protects zone_resuorces_  condition variable, used
   for notifying changes in open_io_zones_ */
  std::mutex zone_resources_mtx_;
  std::condition_variable zone_resources_;
  std::mutex zone_deferred_status_mutex_;
  IOStatus zone_deferred_status_;

  std::condition_variable migrate_resource_;
  std::mutex migrate_zone_mtx_;
  std::atomic<bool> migrating_{false};

  // unsigned int max_nr_active_io_zones_;
  std::atomic<int> max_nr_active_io_zones_;
  int max_nr_open_io_zones_;
  // unsigned int max_nr_open_io_zones_;
  //
  uint64_t cur_free_percent_ = 100;
  //

  std::shared_ptr<ZenFSMetrics> metrics_;

  void EncodeJsonZone(std::ostream &json_stream,
                      const std::vector<Zone *> zones);
  void CalculateResetThreshold(uint64_t free_percent);
  void CalculateFinishThreshold(uint64_t free_percent);
  uint32_t reset_scheme_;
  // bool reset_at_foreground_;
  uint64_t allocation_scheme_;
  uint64_t zc_scheme_;
  double alpha_value_;
  double sigma_value_;
  uint64_t finish_scheme_;
  uint64_t predict_cnt_;
  uint32_t partial_reset_scheme_;
  uint64_t input_aware_scheme_;
  uint64_t tuning_point_;
  uint64_t cbzc_enabled_;
  uint64_t cumulative_io_blocking_ = 0;  // ms
  uint64_t cumulative_1_ = 0;
  uint64_t cumulative_2_ = 0;
  uint64_t cumulative_3_ = 0;
  uint64_t cumulative_4_ = 0;
  uint64_t cumulative_5_ = 0;
  uint64_t cumulative_6_ = 0;
  uint64_t cumulative_7_ = 0;

  uint64_t calculate_lapse = 0;  // ms

  uint64_t default_extent_size_ = 256 << 20;
  enum {
    kEager = 0,
    kLazy = 1,
    kFAR = 2,
    kLazy_Log = 3,
    kLazy_Linear = 4,
    kCustom = 5,
    kLogLinear = 6,
    kNoRuntimeReset = 7,
    kNoRuntimeLinear = 8,
    kLazyExponential = 9
  };
  struct FARStat {
    uint64_t free_percent_;
    size_t reset_count_;
    size_t reset_count_zc_;
    size_t erase_size_ = 0;
    size_t erase_size_zc_ = 0;
    size_t erase_size_proactive_zc_ = 0;

    size_t RC_;
    int T_;
    uint64_t R_wp_;  // (%)
    uint64_t RT_;
    size_t candidate_ratio_;

    std::vector<int> num_files_levels_;

    std::vector<double> compaction_scores_;

    std::vector<uint64_t> levels_size_;
    uint64_t compaction_triggered_[10];
    double avg_same_zone_score_[10];
    double avg_inval_score_[10];
    double avg_invalid_ratio_;
    std::vector<uint64_t> invalid_percent_per_zone_;
    uint64_t cur_ops_;
    uint64_t cur_gc_written_;
    uint64_t valid_data_size_;
    uint64_t invalid_data_size_;
    uint64_t cumulative_io_blocking_;
    FARStat(uint64_t fr, size_t rc, uint64_t wwp, int T, uint64_t rt)
        : free_percent_(fr), RC_(rc), T_(T), RT_(rt) {
      if (RC_ == 0) {
        R_wp_ = 100;
      } else
        R_wp_ = (ZONE_SIZE * 100 - wwp * 100 / RC_) / ZONE_SIZE;
    }
  };
  //   void PrintStat(void) {
  //     printf("[%3d] | %3ld  | %3ld |  %3ld | [%3ld] |", T_, free_percent_,
  //     RC_,
  //            R_wp_, (RT_ >> 20));
  //   }
  std::vector<FARStat> far_stats_;
  std::mutex same_zone_score_mutex_;
  std::vector<double> same_zone_score_[10];
  std::vector<double> same_zone_score_for_timelapse_[10];

  std::vector<double> invalidate_score_[10];
  std::vector<double> invalidate_score_for_timelapse_[10];

  std::atomic<uint64_t> same_zone_score_atomic_[10];
  std::atomic<uint64_t> invalidate_score_atomic_[10];
  std::atomic<uint64_t> compaction_triggered_[10];

 public:
  uint64_t GetZCScheme() const { return zc_scheme_; }
  double GetAlphaValue() const { return alpha_value_; }
  double GetSigmaValue() const { return sigma_value_; }
  uint64_t GetDisableFinish() const { return finish_scheme_; }
  uint64_t GetPredictCnt() const { return predict_cnt_; }

  std::atomic<long> active_io_zones_;
  std::atomic<long> open_io_zones_;
  std::atomic<long> migration_io_zones_{0};

  std::atomic<uint64_t> propagation_count_{0};
  void AddPropagationCount(uint64_t val) {
    propagation_count_.fetch_add(val, std::memory_order_relaxed);
  }
  struct timespec motivation_mount_time_;
  std::vector<motivation_zone_lifetime_diff> motivation_zone_lifetime_diffs_;
  struct OverlappingStat {
    std::atomic<uint64_t> numerator;
    std::atomic<uint64_t> denominator;
  };
  OverlappingStat stats_[6];
  std::atomic<int> file_operation_sequence_{0};

  
  int latest_file_operation_sequence_[10];
  bool check_coldest_[10];
  
  std::atomic<int> CBSC_mispredict_stats_[10];
  std::atomic<int> CBSC_total_predict_stats_[10];

  bool coldest_type_set_  =false;
  int coldest_type_ = -1;
  std::mutex coldest_type_lock_;

  void SetCBSCColdestType(){
    int i;
    int min_seq=100000000;
    int ret=0;
    for(i = 0;i<10;i++){
      if(latest_file_operation_sequence_[i]==0){
        continue;
      }
      int tmp = latest_file_operation_sequence_[i];
      if(tmp < min_seq){
        min_seq= tmp;
        ret = i;
      } 
    }
    // return ret
    coldest_type_=ret;
  }

  // int GetCBSCColdestType(){
  //   return cur_coldest_type_;
  // }
  void PrintMisPredictStats(){
    int i;
    printf("PrintMisPredictStats\n");
    for(i = 0;i<10;i++){
      if(CBSC_total_predict_stats_[i]==0){
        continue;
      }
      int a = CBSC_mispredict_stats_[i].load();
      int b = CBSC_total_predict_stats_[i].load();
      printf("[%d] %d / %d = %d\n",i,a,b,
      a*10000/b);
    }
  }
  
  std::atomic<uint64_t> lsm_tree_[10];
  std::array<uint64_t, 10> GetCurrentLSMTree();
  uint64_t max_bytes_for_level_base_ = 256 << 20;
  std::atomic<uint64_t> total_deletion_after_copy_time_{0};
  std::atomic<uint64_t> total_deletion_after_copy_n_{0};
  std::atomic<uint64_t> total_deletion_after_copy_size_{0};
  std::atomic<uint64_t> actual_cost_benefit_score_{0};


  std::atomic<uint64_t> total_deletion_after_copy_seq_{0};
  std::atomic<uint64_t> total_deletion_after_copy_seq_distribution_[SEQ_DIST_MAX];
  std::atomic<uint64_t> cost_benefit_score_sum_sequence_mb_{0};

  struct LSM_Tree_State{
    bool set= false;
    uint64_t lsm_tree_[10];
    std::vector<uint64_t> ssts[10];
    std::unordered_map<uint64_t, uint64_t> oscore_map[10];
  };

  LSM_Tree_State prev_state_;
  int cur_max_level_ = 0;
  

  std::atomic<uint64_t> file_size_distribution_[1077];

  bool zc_until_set_ = false;
  uint64_t zc_ = 20;
  uint64_t until_ = 20;
  std::atomic<bool> zc_running_ = false;

  uint64_t GetLevelSizeLimit(int level) {
    uint64_t max_bytes_for_level = max_bytes_for_level_base_;
    for (int l = 1; l < level; l++) {
      max_bytes_for_level *= 10;
    }
    return max_bytes_for_level;
  }

  double PredictCompactionScore(int level) {
    if (cur_free_percent_ > 97) {
      return 0.0;
    }
    if (db_ptr_ == nullptr) {
      return 0.0;
    }
    // if (allocation_scheme_ == CAZA_ADV) {
    if (level == 0) {
      if (db_ptr_ == nullptr) {
        return 0.0;
      }
      return db_ptr_->ReCalculateCompactionScore(0);
    }

    if (level == 1) {
      return (double)((double)(lsm_tree_[level].load()) /
                      (double)(max_bytes_for_level_base_));
    }

    return (double)((double)(lsm_tree_[level].load()) /
                    (double)GetLevelSizeLimit(level));
    // }

    // return db_ptr_->ReCalculateCompactionScore(level);
  }

  double PredictCompactionScoreTmp(int level,
                                   std::array<uint64_t, 10> &tmp_lsm_tree,
                                   int initial_l0_files_n) {
    if (cur_free_percent_ > 97) {
      return 0.0;
    }
    if (db_ptr_ == nullptr) {
      return 0.0;
    }

    double score = 0.0;
    double size_score = 0.0;

    if (level == 0) {
      size_score = static_cast<double>(tmp_lsm_tree[0]) /
                   static_cast<double>(max_bytes_for_level_base_);
      score = std::max(static_cast<double>(initial_l0_files_n) / 4, size_score);
      // printf(
      //     "  Calculating score for level 0: %lu / %lu = %.4f, max socre: "
      //     "%.4f\n",
      //     tmp_lsm_tree[0], max_bytes_for_level_base_, size_score, score);
    } else if (level == 1) {
      score = static_cast<double>(tmp_lsm_tree[1]) /
              static_cast<double>(max_bytes_for_level_base_);
      // printf("  Calculating score for level 1: %lu / %lu = %.4f\n",
      //        tmp_lsm_tree[1], max_bytes_for_level_base_, score);
    } else {
      uint64_t level_size_limit = GetLevelSizeLimit(level);
      score = static_cast<double>(tmp_lsm_tree[level]) /
              static_cast<double>(level_size_limit);
      // printf("  Calculating score for level %d: %lu / %lu = %.4f\n", level,
      //        tmp_lsm_tree[level], level_size_limit, score);
    }

    return score;
  }
  inline uint64_t GetAllocationScheme() { return allocation_scheme_; }

  uint64_t GetZoneCleaningKickingPoint() {
    if (zc_until_set_) {
      return zc_;
    }
    if (io_zones[0]->max_capacity_ > (1 << 30)) {
      return 30;
    }
    if (io_zones[0]->max_capacity_ == (1 << 30)) {
      return 30;
    }
    return 10;
  }
  //
  explicit ZonedBlockDevice(std::string path, ZbdBackendType backend,
                            std::shared_ptr<Logger> logger,
                            std::shared_ptr<ZenFSMetrics> metrics =
                                std::make_shared<NoZenFSMetrics>());
  virtual ~ZonedBlockDevice();

  IOStatus Open(bool readonly, bool exclusive);

  Zone *GetIOZone(uint64_t offset);

  IOStatus AllocateIOZone(Env::WriteLifeTimeHint file_lifetime, IOType io_type,
                          Zone **out_zone);
  IOStatus AllocateIOZone(std::string fname, bool is_sst, Slice &smallest,
                          Slice &largest, int level,
                          Env::WriteLifeTimeHint file_lifetime, IOType io_type,
                          uint64_t predicted_size, Zone **out_zone,
                          uint64_t min_capacity);
  void SetZoneAllocationFailed() { zone_allocation_state_ = false; }
  bool IsZoneAllocationFailed() { return zone_allocation_state_ == false; }
  IOStatus AllocateMetaZone(Zone **out_meta_zone);

  uint64_t GetFreeSpace();
  uint64_t GetTotalSpace();
  bool PerformZoneCompaction();
  uint64_t GetUsedSpace();
  uint64_t GetReclaimableSpace();
  ////////////////////////////////////
  uint64_t GetFreePercent();
  //////////////////////////////////////
  std::string GetFilename();
  uint32_t GetBlockSize();
  IOStatus RuntimeZoneReset();
  IOStatus RuntimePartialZoneReset(std::vector<bool> &is_reseted);
  IOStatus ResetUnusedIOZones();
  //////////////////////////////////////////////////
  void AddIOBlockedTimeLapse(int s, int e) {
    // std::lock_guard는 범위 기반 잠금을 제공하여, 이 객체가 존재하는 동안
    // 뮤텍스가 잠기고 객체가 소멸될 때 자동으로 잠금이 해제됩니다.
    std::lock_guard<std::mutex> lg_(io_lock_);
    io_block_timelapse_.push_back({gettid(), s, e});
    // zc_io_block_ += (e - s);
    zone_cleaning_io_block_ += (e - s);
  }

  clock_t IOBlockedStartCheckPoint(void) {
    std::lock_guard<std::mutex> lg_(io_lock_);
    clock_t ret = clock();
    io_blocked_thread_n_++;
    return ret;
  }
  // void IOBlockedEndCheckPoint(int start) {
  //   int end = clock();
  //   std::lock_guard<std::mutex> lg_(io_lock_);
  //   io_blocked_thread_n_--;
  //   io_block_timelapse_.push_back({gettid(), start, -1});
  //   if (io_blocked_thread_n_ == 0) {
  //     zc_io_block_ += (end - start);
  //   }
  //   return;
  // }
  void AddZCTimeLapse(int s, int e, long long us, size_t zc_z, size_t copied,
                      bool forced, double zlv) {
    if (forced == true) {
      force_zc_should_triggered_.store(false);
    }
    zc_timelapse_.push_back({zc_z, s, e, us, copied, forced, zlv});
  }
  void AddTimeLapse(int T);
  void AddCumulativeIOBlocking(long ns) {
    cumulative_io_blocking_ += (ns / 1000) / 1000;
  }
  void AddCumulative_1(long ns) { cumulative_1_ += (ns / 1000) / 1000; }
  void AddCumulative_2(long ns) { cumulative_2_ += (ns / 1000) / 1000; }
  void AddCumulative_3(long ns) { cumulative_3_ += (ns / 1000) / 1000; }
  void AddCumulative_4(long ns) { cumulative_4_ += (ns / 1000) / 1000; }
  void AddCumulative_5(long ns) { cumulative_5_ += (ns / 1000) / 1000; }
  void AddCumulative_6(long ns) { cumulative_6_ += (ns / 1000) / 1000; }
  void AddCumulative_7(long ns) { cumulative_7_ += (ns / 1000) / 1000; }

  void AddCalculatelifetimeLapse(long ns) {
    // if (zc_scheme_ == CBZC3) {
    //   // printf("calculate lifetime CBZC3!!\n");
    calculate_lapse += (ns / 1000) / 1000;
    // } else if (zc_scheme_ == CBZC2) {
    //   // std::cout << ns << std::endl;
    //   calculate_lapse += (ns / 1000);
    //   // std::cout << calculate_lapse << std::endl;
    // } else {
    //   calculate_lapse += (ns / 1000) / 1000;
    //   // std::cout << calculate_lapse << std::endl;
    // }
  }

  uint64_t CalculateCapacityRemain() {
    uint64_t ret = 0;
    for (const auto z : io_zones) {
      ret += z->capacity_;
    }
    return ret;
  }
  uint64_t GetEmptyZoneN() {
    uint64_t ret = 0;
    for (auto z : io_zones) {
      if (z->IsEmpty()) {
        ret++;
      }
    }
    return ret;
  }

  uint64_t GetFullZoneN() {
    uint64_t ret = 0;
    for (auto z : io_zones) {
      if (z->IsFull()) {
        ret++;
      }
    }
    return ret;
  }

  bool ShouldZCByEmptyZoneN() {
    if (GetEmptyZoneN() < zc_) {
      return true;
    }
    return false;
  }

  uint64_t CalculateFreePercent(void) {
    // uint64_t device_size = (uint64_t)ZENFS_IO_ZONES * (uint64_t)ZONE_SIZE;
    uint64_t device_size = (uint64_t)io_zones.size() *
                           BYTES_TO_MB(io_zones[0]->max_capacity_);  // MB
    // std::cout << "calculated device_size: " << device_size << std::endl;
    // uint64_t zone_sz = BYTES_TO_MB(zbd_be_->GetZoneSize());  // MB
    // uint64_t device_size = (uint64_t)GetNrZones() * zone_sz;  // MB
    // printf("calcuatefreepercent::io_zones.size() : %ld\n",
    // io_zones.size()); uint64_t device_size = io_zones.size() * zone_sz;  //
    // MB
    // uint64_t device_size = (uint64_t)80 * zone_sz;  // MB
    uint64_t d_free_space = device_size;  // MB
    uint64_t writed = 0;
    // for (const auto z : io_zones) {
    //   // if (z->IsBusy()) {
    //   //   d_free_space -= (uint64_t)ZONE_SIZE;
    //   // } else {
    //   writed += z->wp_ - z->start_;  // BYTE
    //   // }
    // }
    for (const auto z : io_zones) {
      uint64_t tmp = z->wp_ - z->start_;
      if (tmp > z->max_capacity_) {
        tmp = z->max_capacity_;
      }
      writed += tmp;
    }

    // printf("df1 %ld\n", d_free_space);
    // d_free_space -= (writed >> 20);
    d_free_space -= BYTES_TO_MB(writed);
    // printf("df2 %ld\n", d_free_space);
    device_free_space_.store(d_free_space);
    cur_free_percent_ = (d_free_space * 100) / device_size;
    // CalculateResetThreshold();
    // printf("cf %ld\n", cur_free_percent_);
    return cur_free_percent_;
  }

  void LogZoneStats();
  void LogZoneUsage();
  void LogGarbageInfo();

  uint64_t GetZoneSize();
  uint32_t GetNrZones();
  std::vector<Zone *> GetMetaZones() { return meta_zones; }

  void SetFinishTreshold(uint32_t threshold) { finish_threshold_ = threshold; }

  void PutOpenIOZoneToken();
  void PutActiveIOZoneToken();

  void EncodeJson(std::ostream &json_stream);

  void SetZoneDeferredStatus(IOStatus status);

  std::shared_ptr<ZenFSMetrics> GetMetrics() { return metrics_; }

  void GetZoneSnapshot(std::vector<ZoneSnapshot> &snapshot);

  int Read(char *buf, uint64_t offset, int n, bool direct);
  IOStatus InvalidateCache(uint64_t pos, uint64_t size);

  IOStatus ReleaseMigrateZone(Zone *zone);

  IOStatus TakeMigrateZone(Zone **out_zone, Env::WriteLifeTimeHint lifetime,
                           uint32_t min_capacity);
  IOStatus TakeMigrateZone(Slice &smallest, Slice &largest, int level,
                           Zone **out_zone,
                           Env::WriteLifeTimeHint file_lifetime,
                           uint64_t file_size, uint64_t min_capacity,
                           bool *run_gc_worker_, bool is_sst);

  // void AddBytesWritten(uint64_t written) { bytes_written_ += written; };
  // void AddGCBytesWritten(uint64_t written) { gc_bytes_written_ += written;
  // };
  void AddBytesWritten(uint64_t written) { bytes_written_.fetch_add(written); };
  void AddFinishCount(uint64_t count) { finish_count_.fetch_add(count); };
  void AddGCBytesWritten(uint64_t written) {
    gc_bytes_written_.fetch_add(written);
    zc_copied_timelapse_.push_back(written);
  };
  uint64_t GetGCBytesWritten(void) { return gc_bytes_written_.load(); }
  uint64_t GetRC(void) { return reset_count_.load(); }
  uint64_t GetBlocking(void) { return cumulative_io_blocking_; }
  uint64_t GetUserBytesWritten() {
    return bytes_written_.load() - gc_bytes_written_.load();
  };
  uint64_t GetTotalBytesWritten() { return bytes_written_.load(); };
  int GetResetCount() { return reset_count_.load(); }
  uint64_t GetWWP() { return wasted_wp_.load(); }
  // void SetResetScheme(uint32_t r, bool f, uint64_t T) {
  //   std::cout << "zbd_->SetResetScheme: r = " << r << ", f = " << f
  //             << ", T = " << T << std::endl;
  //   reset_scheme_ = r;
  //   reset_at_foreground_ = f;
  //   tuning_point_ = T;
  //   // if(zc!=0){
  //   //   zc_until_set_=true;
  //   //   zc_=zc;
  //   //   until_=until;
  //   // }

  //   // for(uint64_t f=0;f<=100;f++){
  //   //   CalculateResetThreshold(f);
  //   // }
  // }
  void SetResetScheme(uint32_t r, uint32_t partial_reset_scheme, uint64_t T,
                      uint64_t zc, uint64_t until, uint64_t allocation_scheme,
                      uint64_t zc_scheme, double alpha_value,
                      double sigma_value, uint64_t finish_scheme,
                      uint64_t predict_cnt,
                      std::vector<uint64_t> &other_options) {
    reset_scheme_ = r;
    allocation_scheme_ = allocation_scheme;
    zc_scheme_ = zc_scheme;
    alpha_value_ = alpha_value;
    sigma_value_ = sigma_value;
    finish_scheme_ = finish_scheme;
    predict_cnt_ = predict_cnt;
    partial_reset_scheme_ = partial_reset_scheme;
    tuning_point_ = T;
    input_aware_scheme_ = other_options[0];
    cbzc_enabled_ = other_options[1];
    default_extent_size_ = other_options[2];

    std::cout << "zbd_->SetResetScheme: r = " << r << ", T = " << T
              << ", allocation_schme = " << allocation_scheme
              << ", zc_scheme = " << zc_scheme
              << ", finish_scheme = " << finish_scheme
              << ", predict_cnt = " << predict_cnt << std::endl;

    for (auto opt : other_options) {
      printf("other options %lu\n", opt);
    }
    if (zc != 0) {
      zc_until_set_ = true;
      zc_ = zc;
      until_ = until;
    }

    for (uint64_t f = 0; f <= 100; f++) {
      CalculateResetThreshold(f);
    }
    for (uint64_t f = 0; f <= 100; f++) {
      CalculateFinishThreshold(f);
    }
  }
  void SetDBPtr(DB *db_ptr) { db_ptr_ = db_ptr; }

  bool SetSSTFileforZBDNoLock(uint64_t fno, ZoneFile *zoneFile);

  ZoneFile *GetSSTZoneFileInZBDNoLock(uint64_t fno);

  void GiveZenFStoLSMTreeHint(
      std::vector<uint64_t> &compaction_inputs_input_level_fno,
      std::vector<uint64_t> &compaction_inputs_output_level_fno,
      int output_level, bool trivial_move);

  IOStatus RuntimeReset(void);
  double GetMaxInvalidateCompactionScore(std::vector<uint64_t> &file_candidates,
                                         uint64_t *candidate_size, bool stats);
  double GetMaxSameZoneScore(std::vector<uint64_t> &compaction_inputs_fno);
  inline bool RuntimeZoneResetDisabled() {
    return partial_reset_scheme_ == RUNTIME_ZONE_RESET_DISABLED;
  }
  inline bool RuntimeZoneResetOnly() {
    return partial_reset_scheme_ == RUNTIME_ZONE_RESET_ONLY;
  }
  inline bool PartialResetWithZoneReset() {
    return (partial_reset_scheme_ == PARTIAL_RESET_WITH_ZONE_RESET);
  }
  inline bool PartialResetOnly() {
    return partial_reset_scheme_ == PARTIAL_RESET_ONLY;
    //  &&         log2_erase_unit_size_ > 0;
  }
  inline bool PartialResetAtBackground() {
    return partial_reset_scheme_ == PARTIAL_RESET_AT_BACKGROUND;
  }
  inline bool PartialResetAtBackgroundThresholdWithZoneReset() {
    return partial_reset_scheme_ == PARTIAL_RESET_BACKGROUND_T_WITH_ZONE_RESET;
  }
  inline bool ProactiveZoneCleaning() {
    return partial_reset_scheme_ == PROACTIVE_ZONECLEANING;
  }

  uint32_t GetPartialResetScheme() { return partial_reset_scheme_; }

  void PrintZoneToFileStatus(void);

  void SameLevelFileList(int level, std::vector<uint64_t> &fno_list,
                         bool exclude_being_compacted = true);
  void SameLevelFileList(int level, std::vector<uint64_t> &fno_list,
                         std::set<uint64_t> &compacting_files);
          
  void SameLevelFileList(int level, std::unordered_map<uint64_t, uint64_t>& file_map,
                                 bool exclude_being_compacted = true, bool overap = true);  
  void UpperLevelFileList(Slice &smallest, Slice &largest, int level,
                          std::vector<uint64_t> &fno_list);
  bool GetMinMaxKey(uint64_t fno, Slice &smallest, Slice &largest);

  void SetZCRunning(bool v) { zc_running_.store(v); }
  bool GetZCRunning(void) { return zc_running_.load(); }
  void TrivialMoveFiles(int level, std::set<uint64_t> &trivial_set);
  void DownwardAdjacentFileList(Slice &s, Slice &l, int level,
                                std::vector<uint64_t> &fno_list);
  double GetAvgCompressibilityOflevel(int output_level);
  //
 private:
  std::vector<std::pair<uint64_t, uint64_t>> SortedByZoneScore(
      std::vector<uint64_t> &zone_score) {
    std::vector<std::pair<uint64_t, uint64_t>> ret;
    ret.clear();
    for (uint64_t index = 0; index < zone_score.size(); index++) {
      ret.push_back({zone_score[index], index});
    }
    std::sort(ret.rbegin(), ret.rend());
    return ret;
  }

  IOStatus GetZoneDeferredStatus();
  bool GetActiveIOZoneTokenIfAvailable();
  void WaitForOpenIOZoneToken(bool prioritized);
  IOStatus ApplyFinishThreshold();
  // IOStatus FinishCheapestIOZone(bool put_token = true);
  bool FinishProposal(bool put_token);
  bool FinishProposal2(bool put_token);
  bool FinishCheapestIOZone(bool put_token = true);
  IOStatus GetBestOpenZoneMatch(Env::WriteLifeTimeHint file_lifetime,
                                unsigned int *best_diff_out, Zone **zone_out,
                                uint32_t min_capacity = 0);
  IOStatus GetAnyLargestRemainingZone(Zone **zone_out,
                                      uint32_t min_capacity = 0);
  IOStatus AllocateEmptyZone(Zone **zone_out);
  //////////////////////////////////////////////
  IOStatus AllocateAllInvalidZone(Zone **zone_out);
  bool CompactionSimulator(uint64_t predicted_size, int level, Slice &smallest,
                           Slice &largest);
  bool CalculateZoneScore(std::vector<uint64_t> &fno_list,
                          std::vector<uint64_t> &zone_score);
  void AllocateZoneBySortedScore(
      std::vector<std::pair<uint64_t, uint64_t>> &sorted, Zone **allocated_zone,
      uint64_t min_capacity);
  IOStatus AllocateCompactionAwaredZone(Slice &smallest, Slice &largest,
                                        int level,
                                        Env::WriteLifeTimeHint file_lifetime,
                                        uint64_t predicted_size,
                                        Zone **zone_out,
                                        uint64_t min_capacity = 0);
  IOStatus AllocateCompactionAwaredZoneV2(Slice &smallest, Slice &largest,
                                          int level,
                                          Env::WriteLifeTimeHint file_lifetime,
                                          uint64_t predicted_size,
                                          Zone **zone_out,
                                          uint64_t min_capacity = 0);
  IOStatus AllocateMostL0FilesZone(std::vector<uint64_t> &zone_score,
                                   std::vector<uint64_t> &fno_list,
                                   std::vector<bool> &is_input_in_zone,
                                   Zone **zone_out, uint64_t min_capacity);
  void AdjacentFileList(Slice &smallest, Slice &largest, int level,
                        std::vector<uint64_t> &fno_list);

  uint64_t MostLargeUpperAdjacentFile(Slice &s, Slice &l, int level);
  uint64_t MostSmallDownwardAdjacentFile(Slice &s, Slice &l, int level);
  // void SameLevelFileList(int level, std::vector<uint64_t> &fno_list,
  //                        bool exclude_being_compacted = true);
  // int NumLevelFiles(int level);
  IOStatus AllocateSameLevelFilesZone(Slice &smallest, Slice &largest,
                                      const std::vector<uint64_t> &fno_list,
                                      std::vector<bool> &is_input_in_zone,
                                      Zone **zone_out, uint64_t min_capacity);
  IOStatus GetNearestZoneFromZoneFile(ZoneFile *zFile,
                                      std::vector<bool> &is_input_in_zone,
                                      Zone **zone_out, uint64_t min_capacity);
  ////////////////////////////////////////////////////////
  inline uint64_t LazyLog(uint64_t sz, uint64_t fr, uint64_t T);
  inline uint64_t LazyLinear(uint64_t sz, uint64_t fr, uint64_t T);
  inline uint64_t Custom(uint64_t sz, uint64_t fr, uint64_t T);
  inline uint64_t LogLinear(uint64_t sz, uint64_t fr, uint64_t T);
  inline uint64_t LazyExponential(uint64_t sz, uint64_t fr, uint64_t T);
};

}  // namespace ROCKSDB_NAMESPACE

#endif  // !defined(ROCKSDB_LITE) && defined(OS_LINUX)
