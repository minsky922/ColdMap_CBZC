// Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved.
// Copyright (c) 2019-present, Western Digital Corporation
//  This source code is licensed under both the GPLv2 (found in the
//  COPYING file in the root directory) and Apache 2.0 License
//  (found in the LICENSE.Apache file in the root directory).
#if !defined(ROCKSDB_LITE) && !defined(OS_WIN)

#include "zbd_zenfs.h"

#include <assert.h>
#include <errno.h>
#include <fcntl.h>
#include <libzbd/zbd.h>
#include <linux/blkzoned.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ioctl.h>
#include <time.h>
#include <unistd.h>

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <mutex>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

// #include "rocksdb/db/db_impl/db_impl.cc"
#include "rocksdb/env.h"
#include "rocksdb/io_status.h"
#include "snapshot.h"
#include "zbdlib_zenfs.h"
#include "zonefs_zenfs.h"

#define KB (1024)
#define MB (1024 * KB)

/* Number of reserved zones for metadata
 * Two non-offline meta zones are needed to be able
 * to roll the metadata log safely. One extra
 * is allocated to cover for one zone going offline.
 */

// #define ZENFS_SPARE_ZONES (1)

// #define ZENFS_META_ZONES (3)

#define ZENFS_IO_ZONES (100)
// #define ZENFS_IO_ZONES (1843)  // 1.8T

/* Minimum of number of zones that makes sense */
#define ZENFS_MIN_ZONES (32)

namespace ROCKSDB_NAMESPACE {

Zone::Zone(ZonedBlockDevice *zbd, ZonedBlockDeviceBackend *zbd_be,
           std::unique_ptr<ZoneList> &zones, unsigned int idx)
    : zbd_(zbd),        // ZonedBlockDevice 객체의 포인터 초기화
      zbd_be_(zbd_be),  // ZonedBlockDeviceBackend 객체의 포인터 초기화
      busy_(false),     // 초기 상태는 busy 아님
      start_(zbd_be->ZoneStart(zones, idx)),  // 존의 시작 위치 초기화
      max_capacity_(
          zbd_be->ZoneMaxCapacity(zones, idx)),  // 존의 최대 용량 초기화
      wp_(zbd_be->ZoneWp(zones, idx)),
      zidx_(idx) {                // 존의 현재 쓰기 포인터 초기화
  lifetime_ = Env::WLTH_NOT_SET;  // 존의 수명 초기화
  used_capacity_ = 0;             // 사용된 용량 초기화
  capacity_ = 0;                  // 현재 용량 초기화
  zone_sz_ = zbd_be_->GetZoneSize();
  block_sz_ = zbd_be_->GetBlockSize();
  if (zbd_be->ZoneIsWritable(zones, idx))  // 존이 쓰기 가능한 상태인지 확인
    capacity_ =
        max_capacity_ - (wp_ - start_);  // 쓰기 가능한 경우 현재 용량 설정
}

bool Zone::IsUsed() { return (used_capacity_ > 0); }
uint64_t Zone::GetCapacityLeft() { return capacity_; }
bool Zone::IsFull() { return (capacity_ == 0); }
bool Zone::IsEmpty() { return (wp_ == start_); }
bool Zone::IsFinished() { return is_finished_; }
uint64_t Zone::GetZoneNr() { return start_ / zbd_->GetZoneSize(); }  // 존 넘버

void Zone::EncodeJson(std::ostream &json_stream) {
  json_stream << "{";
  json_stream << "\"start\":" << start_ << ",";
  json_stream << "\"capacity\":" << capacity_ << ",";
  json_stream << "\"max_capacity\":" << max_capacity_ << ",";
  json_stream << "\"wp\":" << wp_ << ",";
  json_stream << "\"lifetime\":" << lifetime_ << ",";
  json_stream << "\"used_capacity\":" << used_capacity_;
  json_stream << "}";
}

IOStatus Zone::Reset() {
  bool offline;
  uint64_t max_capacity;

  assert(!IsUsed());
  // assert(IsBusy());
  // is_finished_ = false;
  // struct{
  //   initial_allocated_time_;
  //   reset_time;
  //   std::vector<time> a;

  // }

  if (this_zone_motivation_check_) {
    printf("@@@@@@@@@@@@@@@@@@ reset motivation_lifetime_diffs RESET\n");
    struct timespec motivation_reset_time;
    clock_gettime(CLOCK_MONOTONIC, &motivation_reset_time);
    motivation_zone_lifetime_diff tmp;
    tmp.zone_allocated_time_ = allocated_time_;
    tmp.zone_resetted_time_ = motivation_reset_time;
    tmp.motivation_lifetime_diffs = motivation_lifetime_diffs;
    zbd_->motivation_zone_lifetime_diffs_.push_back(tmp);
    motivation_lifetime_diffs.clear();
    is_allocated_ = false;
  }

  IOStatus ios = zbd_be_->Reset(start_, &offline, &max_capacity);
  if (ios != IOStatus::OK()) return ios;

  if (offline)
    capacity_ = 0;  // 존이 오프라인 상태이면 용량을 0으로 설정
  else
    max_capacity_ = capacity_ = max_capacity;  // 최대 용량과 현재 용량을 새로
                                               // 설정된 최대 용량으로 업데이트

  wp_ = start_;  // 쓰기 포인터(write pointer)를 존의 시작 위치로 재설정
  lifetime_ = Env::WLTH_NOT_SET;  // 존의 수명(lifetime)을 초기화

  //////
  // std::cout << "#####zone Reset" << std::endl;
  ////
  return IOStatus::OK();
}

IOStatus Zone::Finish() {
  // assert(IsBusy());

  IOStatus ios = zbd_be_->Finish(start_);
  if (ios != IOStatus::OK()) return ios;

  capacity_ = 0;
  // is_finished_ = true;
  wp_ = start_ + zbd_->GetZoneSize();
  // zbd_->AddFinishCount(1);
  // finish_count_++;
  //////
  // std::cout << "######zone Finish" << std::endl;
  ////
  return IOStatus::OK();
}

IOStatus Zone::Close() {
  // assert(IsBusy());

  if (!(IsEmpty() || IsFull())) {
    IOStatus ios = zbd_be_->Close(start_);
    if (ios != IOStatus::OK()) return ios;
  }
  //////
  // std::cout << "######Starting close" << std::endl;
  ////

  return IOStatus::OK();
}

IOStatus Zone::Append(char *data, uint32_t size) {
  ZenFSMetricsLatencyGuard guard(zbd_->GetMetrics(), ZENFS_ZONE_WRITE_LATENCY,
                                 Env::Default());
  zbd_->GetMetrics()->ReportThroughput(ZENFS_ZONE_WRITE_THROUGHPUT, size);
  char *ptr = data;
  uint32_t left = size;
  int ret;
  recent_inval_time_ = std::chrono::system_clock::now();
  if (capacity_ < size)
    return IOStatus::NoSpace("Not enough capacity for append");

  assert((size % zbd_->GetBlockSize()) == 0);

  while (left) {
    ret = zbd_be_->Write(ptr, left, wp_);
    if (ret < 0) {
      return IOStatus::IOError(strerror(errno));
    }

    ptr += ret;
    wp_ += ret;
    capacity_ -= ret;
    left -= ret;
    zbd_->AddBytesWritten(ret);
  }

  return IOStatus::OK();
}

inline IOStatus Zone::CheckRelease() {
  if (!Release()) {
    assert(false);
    return IOStatus::Corruption("Failed to unset busy flag of zone " +
                                std::to_string(GetZoneNr()));
  }

  return IOStatus::OK();
}

Zone *ZonedBlockDevice::GetIOZone(uint64_t offset) {
  for (const auto z : io_zones)
    if (z->start_ <= offset && offset < (z->start_ + zbd_be_->GetZoneSize()))
      return z;
  return nullptr;
}

ZonedBlockDevice::ZonedBlockDevice(std::string path, ZbdBackendType backend,
                                   std::shared_ptr<Logger> logger,
                                   std::shared_ptr<ZenFSMetrics> metrics)
    : logger_(logger), metrics_(metrics) {
  if (backend == ZbdBackendType::kBlockDev) {
    zbd_be_ = std::unique_ptr<ZbdlibBackend>(new ZbdlibBackend(path));
    Info(logger_, "New Zoned Block Device: %s", zbd_be_->GetFilename().c_str());
  } else if (backend == ZbdBackendType::kZoneFS) {
    zbd_be_ = std::unique_ptr<ZoneFsBackend>(new ZoneFsBackend(path));
    Info(logger_, "New zonefs backing: %s", zbd_be_->GetFilename().c_str());
  }
  zone_sz_ = zbd_be_->GetZoneSize();
  clock_gettime(CLOCK_MONOTONIC, &motivation_mount_time_);
  sst_file_bitmap_ = new ZoneFile *[1 << 20];
  memset(sst_file_bitmap_, 0, (1 << 20));
  for(int i = 0; i<5000;i++){
    total_deletion_after_copy_seq_distribution_[i]=0;
  }
  for( int i =0;i<10;i++){
    latest_file_operation_sequence_[i]=0;
    CBSC_mispredict_stats_[i]=0;
    CBSC_total_predict_stats_[i]=0;
  }
}

IOStatus ZonedBlockDevice::Open(bool readonly, bool exclusive) {
  std::unique_ptr<ZoneList> zone_rep;
  unsigned int max_nr_active_zones;
  unsigned int max_nr_open_zones;
  Status s;
  uint64_t i = 0;
  uint64_t m = 0;
  // Reserve one zone for metadata and another one for extent migration
  // int reserved_zones = 2;
  int reserved_zones = 1;

  if (!readonly && !exclusive)
    return IOStatus::InvalidArgument("Write opens must be exclusive");

  IOStatus ios = zbd_be_->Open(readonly, exclusive, &max_nr_active_zones,
                               &max_nr_open_zones);
  printf("Debug: max_nr_active_zones=%u, max_nr_open_zones=%u\n",
         max_nr_active_zones, max_nr_open_zones);
  if (ios != IOStatus::OK()) return ios;

  if (zbd_be_->GetNrZones() < ZENFS_MIN_ZONES) {
    return IOStatus::NotSupported("To few zones on zoned backend (" +
                                  std::to_string(ZENFS_MIN_ZONES) +
                                  " required)");
  }

  if (max_nr_active_zones == 0)
    max_nr_active_io_zones_ = zbd_be_->GetNrZones();
  else
    max_nr_active_io_zones_ = max_nr_active_zones - reserved_zones;

  if (max_nr_open_zones == 0)
    max_nr_open_io_zones_ = zbd_be_->GetNrZones();
  else
    max_nr_open_io_zones_ = max_nr_open_zones - reserved_zones;

  max_nr_active_io_zones_ = 12;
  max_nr_open_io_zones_ = 12;
  // open(/home/femu/yc-bccp[/open_zone])
  Info(logger_, "Zone block device nr zones: %u max active: %u max open: %u \n",
       zbd_be_->GetNrZones(), max_nr_active_zones, max_nr_open_zones);
  //
  printf(
      "Zone block device nr zones: %u max active: %u (%d) max open: %u(%d) \n",
      zbd_be_->GetNrZones(), max_nr_active_zones,
      max_nr_active_io_zones_.load(), max_nr_open_zones, max_nr_open_io_zones_);
  //
  zone_rep = zbd_be_->ListZones();
  if (zone_rep == nullptr || zone_rep->ZoneCount() != zbd_be_->GetNrZones()) {
    Error(logger_, "Failed to list zones");
    return IOStatus::IOError("Failed to list zones");
  }

  while (m < ZENFS_META_ZONES && i < zone_rep->ZoneCount()) {
    /* Only use sequential write required zones */
    if (zbd_be_->ZoneIsSwr(zone_rep, i)) {
      if (!zbd_be_->ZoneIsOffline(zone_rep, i)) {
        meta_zones.push_back(new Zone(this, zbd_be_.get(), zone_rep, i));
      }
      m++;
    }
    i++;
  }

  active_io_zones_ = 0;
  open_io_zones_ = 0;
  // uint64_t device_io_capacity = 85899345920;  // 80GB
  // uint64_t device_io_capacity = 10737418240;  // 10GB
  // for (; i < zone_rep->ZoneCount() &&
  //        (io_zones.size() * zbd_be_->GetZoneSize()) < (device_io_capacity);
  //      i++) {
  printf("zone_rep->ZoneCount() %u\n", zone_rep->ZoneCount());
  for (; i < zone_rep->ZoneCount() && i < ZENFS_IO_ZONES + 3; i++) {
    // for (; i < ZENFS_IO_ZONES + 3; i++) {
    /* Only use sequential write required zones */
    /*Sequential Write Required (SWR)*/
    /*Conventional Zone (Non-SWR) :  임의의 위치에 자유롭게
     * 쓰기가 가능*/
    if (zbd_be_->ZoneIsSwr(zone_rep, i)) {
      // printf("ZoneIsSwr : zone_rep: %u, i : %ld\n", zone_rep->ZoneCount(),
      // i);
      if (!zbd_be_->ZoneIsOffline(zone_rep, i)) {
        Zone *newZone = new Zone(this, zbd_be_.get(), zone_rep, i);
        // printf("ZoneIsOffline : zone_rep: %u, i : %ld\n",
        // zone_rep->ZoneCount(),
        //  i);
        if (!newZone->Acquire()) {
          printf("Failed to allocate new Zone at index %ld\n", i);
          assert(false);
          return IOStatus::Corruption("Failed to set busy flag of zone " +
                                      std::to_string(newZone->GetZoneNr()));
        }
        io_zones.push_back(newZone);
        // printf("io zone at %ld\n", i);

        if (zbd_be_->ZoneIsActive(zone_rep, i)) {
          printf("ZoneIsActive i : %ld\n", i);
          printf("active resoruced %lu\n", active_io_zones_.load());
          active_io_zones_++;
          if (zbd_be_->ZoneIsOpen(zone_rep, i)) {
            if (!readonly) {
              newZone->Close();
            }
          }
        }
        IOStatus status = newZone->CheckRelease();
        if (!status.ok()) {
          return status;
        }
      }
    }
  }
  // io_zones[3]->this_zone_motivation_check_=true;

  // uint64_t device_free_space = (ZENFS_IO_ZONES) * (ZONE_SIZE);
  // printf("device free space : %ld\n", device_free_space);
  // device_free_space_.store(device_free_space);
  // for (uint64_t f = 0; f <= 100; f++) {
  //   CalculateResetThreshold(f);
  // }
  // for (uint64_t f = 0; f <= 100; f++) {
  //   CalculateFinishThreshold(f);
  // }
  printf("active io zones :  %lu\n", active_io_zones_.load());
  printf("io_zones.size() : %ld\n", io_zones.size());
  printf("zone sz %lu\n", zone_sz_);
  uint64_t device_free_space = io_zones.size() * zbd_be_->GetZoneSize();
  printf("device free space : %ld\n", BYTES_TO_MB(device_free_space));
  // printf("zone sz %ld\n", zone_sz_);
  // device_free_space_.store(device_free_space);

  start_time_ = time(NULL);

  return IOStatus::OK();
}

uint64_t ZonedBlockDevice::GetFreeSpace() {
  uint64_t free = 0;
  for (const auto z : io_zones) {
    free += z->capacity_;
  }

  // std::cout << "######getFreeSpace" << "\n";

  return free;
}

uint64_t ZonedBlockDevice::GetUsedSpace() {
  uint64_t used = 0;
  for (const auto z : io_zones) {
    used += z->used_capacity_;
  }
  return used;
}

uint64_t ZonedBlockDevice::GetReclaimableSpace() {
  uint64_t reclaimable = 0;
  for (const auto z : io_zones) {
    if (z->IsFull()) reclaimable += (z->max_capacity_ - z->used_capacity_);
  }
  return reclaimable;
}

uint64_t ZonedBlockDevice::GetFreePercent() {
  // uint64_t non_free = zbd_->GetUsedSpace() + zbd_->GetReclaimableSpace();
  uint64_t free = GetFreeSpace();
  return (100 * free) / io_zones.size() * io_zones[0]->max_capacity_;
}

void ZonedBlockDevice::LogZoneStats() {
  uint64_t used_capacity = 0;
  uint64_t reclaimable_capacity = 0;
  uint64_t reclaimables_max_capacity = 0;
  uint64_t active = 0;

  for (const auto z : io_zones) {
    used_capacity += z->used_capacity_;

    if (z->used_capacity_) {
      reclaimable_capacity += z->max_capacity_ - z->used_capacity_;
      reclaimables_max_capacity += z->max_capacity_;
    }

    if (!(z->IsFull() || z->IsEmpty())) active++;
  }

  if (reclaimables_max_capacity == 0) reclaimables_max_capacity = 1;

  Info(logger_,
       "[Zonestats:time(s),used_cap(MB),reclaimable_cap(MB), "
       "avg_reclaimable(%%), active(#), active_zones(#), open_zones(#)] %ld "
       "%lu %lu %lu %lu %ld %ld\n",
       time(NULL) - start_time_, used_capacity / MB, reclaimable_capacity / MB,
       100 * reclaimable_capacity / reclaimables_max_capacity, active,
       active_io_zones_.load(), open_io_zones_.load());
}

void ZonedBlockDevice::LogZoneUsage() {
  int i = 0;
  for (const auto z : io_zones) {
    int64_t used = z->used_capacity_;
    // printf(
    //     "@@@ LogZoneUsage [%d] :  remaining : %lu used : %lu invalid : %lu wp
    //     "
    //     ": %lu\n",
    //     i, z->capacity_, used, z->wp_ - z->start_ - used, z->wp_ -
    //     z->start_);
    if (used > 0) {
      Debug(logger_, "Zone 0x%lX used capacity: %ld bytes (%ld MB)\n",
            z->start_, used, used / MB);
    }
    i++;
  }
}

void ZonedBlockDevice::LogGarbageInfo() {
  int zone_gc_stat[12] = {0};  
  for (auto z : io_zones) {    
    if (!z->Acquire()) {      
      continue;
    }

    if (z->IsEmpty()) {  
      zone_gc_stat[0]++;
      z->Release(); 
      continue;
    }

    double garbage_rate = 0;  
    if (z->IsFull()) {        
      garbage_rate =
          double(z->max_capacity_ - z->used_capacity_) / z->max_capacity_;
    } else {  
      garbage_rate =
          double(z->wp_ - z->start_ - z->used_capacity_) / z->max_capacity_;
    }
    assert(garbage_rate >= 0);                 
    int idx = int((garbage_rate + 0.1) * 10);  
    zone_gc_stat[idx]++;  

    z->Release();  
  }

  std::stringstream ss; 
  ss << "Zone Garbage Stats: [";
  for (int i = 0; i < 12; i++) {  
    ss << zone_gc_stat[i] << " ";
  }
  ss << "]";
  Info(logger_, "%s", ss.str().data());  
}

ZonedBlockDevice::~ZonedBlockDevice() {
  size_t rc = reset_count_.load();
  uint64_t wwp = wasted_wp_.load() / (1 << 20);
  uint64_t zone_sz = BYTES_TO_MB(io_zones[0]->max_capacity_);  // MB
  uint64_t R_wp;
  printf("zone size at ~ %lu\n", zone_sz);
  if (rc == 0) {
    R_wp = 100;
  } else {
    R_wp = (zone_sz * 100 - wwp * 100 / (rc)) / zone_sz;  // MB
  }
  printf("============================================================\n");
  printf("WWP (MB) : %lu, R_wp : %lu\n",
         wasted_wp_.load() / (1 << 20), R_wp);
  printf("NEW WWP(MB) : %lu\n", new_wasted_wp_.load() / (1 << 20));
  if (rc != 0) {
    printf("Runtime zone reset R_wp %lu\n",
           ((zone_sz * 100) - ((wwp * 100) / rc)) / zone_sz);
  }
  printf("ZONE FINISH VALID(MB) %lu\n", finished_valid_data_.load() >> 20);
  printf("ZONE FINISH WWP(MB) : %lu\n", finished_wasted_wp_.load() / (1 << 20));
  // printf("ZC IO Blocking time : %d, Compaction Refused : %lu\n",
  // zc_io_block_,
  //        compaction_blocked_at_amount_.size());
  printf("ZC IO Blocking time : %d\n", zone_cleaning_io_block_);

  printf("============================================================\n");
  uint64_t total_copied = 0;
  size_t rc_zc = 0;
  int io_blocking_sum = 0;
  long long io_blocking_ms_sum = 0;
  double sum_zlv = 0.0;
  // for (size_t i = 0;
  //  i < zc_timelapse_.size() && i < zc_copied_timelapse_.size(); i++) {
  for (size_t i = 0; i < zc_timelapse_.size(); i++) {
    // bool forced = zc_timelapse_[i].forced;
    size_t zc_z = zc_timelapse_[i].zc_z;
    int s = zc_timelapse_[i].s;
    int e = zc_timelapse_[i].e;
    long long us = zc_timelapse_[i].us;
    double zlv = zc_timelapse_[i].zlv;
    io_blocking_sum += e - s + 1;
    io_blocking_ms_sum += us / 1000;
    printf(
        "[%lu] :: %d ~ %d, %llu ms, %ld (MB), Reclaimed Zone : %lu, ZLV: "
        "%f\n",
        i + 1, s, e, us / 1000, (zc_timelapse_[i].copied >> 20), zc_z, zlv);
    // total_copied += zc_copied_timelapse_[i];
    total_copied += zc_timelapse_[i].copied;
    rc_zc += zc_z;
    sum_zlv += zlv;
  }
  printf("Total ZC Copied (MB) :: %lu, Recaimed by ZC :: %lu \n",
         total_copied / (1 << 20), rc_zc);
  printf("Total ZC copied- GC_BYTES_WRITTEN(MB):: %lu \n",
         (gc_bytes_written_.load()) >> 20);
  printf("Reset Count (R+ZC) : %ld+%ld=%ld\n", rc - rc_zc, rc_zc,
         rc);
  printf("Finish Count : %ld\n", finish_count_.load());
  double avg_zlv = 0.0;
  size_t count = zc_timelapse_.size();
  if (count > 0) {
    avg_zlv = sum_zlv / static_cast<double>(count);
  }
  printf("Average ZLV :: %f\n", avg_zlv);
  printf("propagation count : %ld\n",
         propagation_count_.load(std::memory_order_relaxed));

  printf("TOTAL I/O BLOKCING TIME %d\n", io_blocking_sum);
  printf("TOTAL I/O BLOCKING TIME(ms) %llu\n", io_blocking_ms_sum);
  printf("Cumulative I/O Blocking(ms) %lu\n", cumulative_io_blocking_);
  // printf("Cumlative 1(ms) %lu\n", cumulative_1_);
  // printf("Cumlative 2(ms) %lu\n", cumulative_2_);
  // printf("Cumlative 3(ms) %lu\n", cumulative_3_);
  // printf("Cumlative 4(ms) %lu\n", cumulative_4_);
  // printf("Cumlative 5(ms) %lu\n", cumulative_5_);
  // printf("Cumlative 6(ms) %lu\n", cumulative_6_);
  // printf("Cumlative 7(ms) %lu\n", cumulative_7_);

  // if (zbd_->GetZCScheme() == 2) {
  //   printf("Calculate Lifetime Time(us) %lu\n", calculate_lapse);
  // } else {
  printf("Calculate Lifetime Time(ms) %lu\n", calculate_lapse);
  // }
  if (GetUserBytesWritten()) {
    printf("copy/written ratio : %lu/%lu=%lu\n", gc_bytes_written_.load(),
           GetUserBytesWritten(),
           (gc_bytes_written_.load() * 100) / GetUserBytesWritten());
  }

  for (int i = 0; i < 4; i++) {
    uint64_t denom = stats_[i].denominator.load();
    uint64_t numer = stats_[i].numerator.load();
    double ratio = 0.0;
    if (denom != 0) {
      ratio = static_cast<double>(numer) / static_cast<double>(denom);
    }
    printf("level[%d]: 분자(overlapping) %lu / 분모(file size) %lu = %.4f\n", i,
           numer, denom, ratio);
  }

  {
    int i = 4;
    uint64_t denom = stats_[i].denominator.load();
    uint64_t numer = stats_[i].numerator.load();
    double ratio = 0.0;
    if (denom != 0) {
      ratio = static_cast<double>(numer) / static_cast<double>(denom);
    }
    printf("level[%d]: 분자(overlapping) %lu / 분모(file size) %lu = %.4f\n", i,
           numer, denom, ratio);
  }

  if (total_deletion_after_copy_n_.load()) {
    // 누적 시간은 마이크로초, 개수는 n개
    // => 평균 시간도 마이크로초 단위
    auto total_time_us = total_deletion_after_copy_time_.load();
    auto total_count = total_deletion_after_copy_n_.load();
    auto avg_us = total_time_us / total_count;

    printf("avg deletion after time %lu us (total: %lu us, count: %lu)\n",
           avg_us, total_time_us, total_count);
    printf("reset_count_before_full_ %lu reset_size_before_full_ %lu\n",
           reset_count_before_full_.load(), reset_size_before_full_.load());
  }

  for (const auto z : meta_zones) {
    delete z;
  }

  // motivation_mount_time_/

  long mount_mili = (motivation_mount_time_.tv_sec * 1000) +
                    (motivation_mount_time_.tv_nsec / 1000000);
  printf("SEQUENCE DISTRIBUTION\n");
  for(int i = 0; i<SEQ_DIST_MAX;i++){
    printf("%ld ",total_deletion_after_copy_seq_distribution_[i].load());
  }
  printf("FILE SIZE DISTRIBUITION\n");
  for(int i = 0; i<1077;i++){
    printf("%ld ",file_size_distribution_[i].load());
  }
  for (const auto z : io_zones) {
    if (z->this_zone_motivation_check_) {
      // motivation_zone_lifetime_diffs_
      int i = 0;
      printf("z->this_zone_motivation_check_ at index  %lu\n", z->zidx_);
      for (auto res : motivation_zone_lifetime_diffs_) {
        printf("zone alloc/delete time %d\n", i++);
        long start_milliseconds = (res.zone_allocated_time_.tv_sec * 1000) +
                                  (res.zone_allocated_time_.tv_nsec / 1000000);
        long end_milliseconds = (res.zone_resetted_time_.tv_sec * 1000) +
                                (res.zone_resetted_time_.tv_nsec / 1000000);

        printf("%lu\t%lu\n", start_milliseconds - mount_mili,
               end_milliseconds - mount_mili);
        printf("------------------\n");

        // std::vector< std::pair<long,long> > sorted;
        std::vector<start_end_to_be_sorted> sorted;
        sorted.clear();
        for (auto file_time : res.motivation_lifetime_diffs) {
          start_milliseconds = (file_time.file_created_time_.tv_sec * 1000) +
                               (file_time.file_created_time_.tv_nsec / 1000000);
          end_milliseconds = (file_time.file_deleted_time_.tv_sec * 1000) +
                             (file_time.file_deleted_time_.tv_nsec / 1000000);
          bool zcied = file_time.zcied_;
          start_milliseconds -= mount_mili;
          end_milliseconds -= mount_mili;
          sorted.push_back({start_milliseconds, end_milliseconds, zcied});
          // printf("%lu\t%lu\n",start_milliseconds-mount_mili,end_milliseconds-mount_mili);
        }

        std::sort(sorted.begin(), sorted.end(),
                  start_end_to_be_sorted::compare_start_end_to_be_sorted);

        for (auto s : sorted) {
          printf("%lu\t%lu\t%s\n", s.start_milliseconds, s.end_milliseconds,
                 s.zcied ? "zcied" : "");
        }
        printf("------------------\n");
      }
    }
    delete z;
  }
  PrintMisPredictStats();

  std::cout << "[Accumulated Stats] right_vertical_total="
  << predict_right_vertical_.load()
  << ", false_vertical_total=" << predict_false_vertical_.load()
  << ", right_horizontal_total=" << predict_right_horizontal_.load()
  << ", false_horizontal_total=" << predict_false_horizontal_.load()
  << std::endl;

   {
    uint64_t rvert  = predict_right_vertical_.load();
    uint64_t fvert  = predict_false_vertical_.load();
    uint64_t totalv = rvert + fvert;
    if (totalv == 0) {
      std::cout << "[Accumulated] Vertical ratio: N/A" << std::endl;
    } else {
      double ratio_v = double(fvert) / double(totalv);
      std::cout << "[Accumulated] Vertical false/(right+false) = " 
                << ratio_v << " (≈ " << ratio_v * 100.0 << "%)" << std::endl;
    }

    uint64_t rhori  = predict_right_horizontal_.load();
    uint64_t fhori  = predict_false_horizontal_.load();
    uint64_t totalh = rhori + fhori;
    if (totalh == 0) {
      std::cout << "[Accumulated] Horizontal ratio: N/A" << std::endl;
    } else {
      double ratio_h = double(fhori) / double(totalh);
      std::cout << "[Accumulated] Horizontal false/(right+false) = "
                << ratio_h << " (≈ " << ratio_h * 100.0 << "%)" << std::endl;
    }
  }
}

#define LIFETIME_DIFF_NOT_GOOD (100)
#define LIFETIME_DIFF_COULD_BE_WORSE (50)

unsigned int GetLifeTimeDiff(Env::WriteLifeTimeHint zone_lifetime,
                             Env::WriteLifeTimeHint file_lifetime) {
  assert(file_lifetime <= Env::WLTH_EXTREME);

  if ((file_lifetime == Env::WLTH_NOT_SET) ||
      (file_lifetime == Env::WLTH_NONE)) {
    if (file_lifetime == zone_lifetime) {
      return 0;
    } else {
      return LIFETIME_DIFF_NOT_GOOD;
    }
  }

  // if (zone_lifetime > file_lifetime) return zone_lifetime - file_lifetime;
  if (zone_lifetime == file_lifetime) return 0;
  return LIFETIME_DIFF_NOT_GOOD;
  // else return file_lifetime-zone_lifetime;
  // return LIFETIME_DIFF_NOT_GOOD;
  // return zone_lifetime - file_lifetime;
}

IOStatus ZonedBlockDevice::AllocateMetaZone(Zone **out_meta_zone) {
  assert(out_meta_zone);
  *out_meta_zone = nullptr;
  ZenFSMetricsLatencyGuard guard(metrics_, ZENFS_META_ALLOC_LATENCY,
                                 Env::Default());
  metrics_->ReportQPS(ZENFS_META_ALLOC_QPS, 1);

  for (const auto z : meta_zones) {
    /* If the zone is not used, reset and use it */
    if (z->Acquire()) {
      if (!z->IsUsed()) {
        if (!z->IsEmpty() && !z->Reset().ok()) {
          Warn(logger_, "Failed resetting zone!");
          IOStatus status = z->CheckRelease();
          if (!status.ok()) return status;
          continue;
        }
        *out_meta_zone = z;
        return IOStatus::OK();
      }
    }
  }
  assert(true);
  Error(logger_, "Out of metadata zones, we should go to read only now.");
  return IOStatus::NoSpace("Out of metadata zones");
}

void ZonedBlockDevice::GiveZenFStoLSMTreeHint(
    std::vector<uint64_t> &compaction_inputs_input_level_fno,
    std::vector<uint64_t> &compaction_inputs_output_level_fno, int output_level,
    bool trivial_move) {
  ZoneFile *zfile = nullptr;

  if (allocation_scheme_ == LIZA) {
    return;
  }

  LSM_Tree_State cur_state;
  cur_state.set=true;

  for(int l =0 ;l <=cur_max_level_;l++){
    cur_state.lsm_tree_[l] = lsm_tree_[l];
    std::unordered_map<uint64_t, uint64_t> file_map;
    this->SameLevelFileList(l, file_map, false, false);
    cur_state.ssts[l].clear();
    for (const auto& [fno, oscore] : file_map) {
      cur_state.ssts[l].push_back(fno);
      cur_state.oscore_map[l][fno] = oscore;
    }
    // for (const auto& [fno, oscore] : file_map) {
    //   // std::cout << "[GiveZenFStoLSMTreeHint] Level=" << l
    //   //           << ", fno=" << fno
    //   //           << ", o_score=" << oscore << std::endl;
    // }
  }

  int predict_right_vertical= 0;
  int predict_false_vertical=0;
  int predict_right_horizontal= 0;
  int predict_false_horizontal=0;
  if(prev_state_.set==true){
    for (int l = 0; l <= cur_max_level_; l++) {
      if (prev_state_.lsm_tree_[l] == cur_state.lsm_tree_[l]) {
        predict_right_vertical++;
      } else {
        predict_false_vertical++;
      }
    }

    for (int l = 0; l <= cur_max_level_; l++) {    
      std::unordered_map<uint64_t, uint64_t> prev_file_map = prev_state_.oscore_map[l];
      std::unordered_map<uint64_t, uint64_t> cur_file_map = cur_state.oscore_map[l];
    
      for (const auto& [fno, oscore] : cur_file_map) {
        if (prev_file_map.find(fno) == prev_file_map.end()) {
          predict_false_horizontal++; //new
        } else {
          if (prev_file_map[fno] == oscore) {
            predict_right_horizontal++;
          } else {
            predict_false_horizontal++;
          }
        }
      }
    
      // deleted
      for (const auto& [fno, oscore] : prev_file_map) {
        if (cur_file_map.find(fno) == cur_file_map.end()) {
          predict_false_horizontal++;
        }
      }
    }
    
  }

  predict_right_vertical_.fetch_add(predict_right_vertical);
  predict_false_vertical_.fetch_add(predict_false_vertical);
  predict_right_horizontal_.fetch_add(predict_right_horizontal);
  predict_false_horizontal_.fetch_add(predict_false_horizontal);

 //save cur state as prev state
  prev_state_ = cur_state;

  if (trivial_move) {
    if (output_level == 1) {
      // printf("0->1 trivial move??\n");
    }

    for (uint64_t fno : compaction_inputs_input_level_fno) {
      zfile = GetSSTZoneFileInZBDNoLock(fno);

      if (zfile == nullptr) {
        printf("why nullptr? %lu\n", fno);
        continue;
      }
      zfile->level_=output_level;
      uint64_t file_size = zfile->predicted_size_;
      stats_[output_level - 1].denominator.fetch_add(file_size);
      lsm_tree_[output_level - 1].fetch_sub(file_size);
      lsm_tree_[output_level].fetch_add(file_size);
    }
    return;
  }

  //////////////// invaldation compaction


  for (uint64_t fno : compaction_inputs_input_level_fno) {
    zfile = GetSSTZoneFileInZBDNoLock(fno);

    if (zfile == nullptr) {
      printf("why nullptr? %lu\n", fno);
      continue;
    }
    uint64_t file_size = zfile->predicted_size_;
    stats_[output_level - 1].denominator.fetch_add(file_size);  // file_size
    lsm_tree_[output_level - 1].fetch_sub(file_size);
  }
  for (uint64_t fno : compaction_inputs_output_level_fno) {
    zfile = GetSSTZoneFileInZBDNoLock(fno);

    if (zfile == nullptr) {
      printf("why nullptr? %lu\n", fno);
      continue;
    }
    uint64_t file_size = zfile->predicted_size_;
    stats_[output_level - 1].numerator.fetch_add(file_size);  // overlapping
    lsm_tree_[output_level].fetch_sub(file_size);
  }

  std::vector<uint64_t> input_fno = compaction_inputs_input_level_fno;
  input_fno.insert(input_fno.end(), compaction_inputs_output_level_fno.begin(),
                   compaction_inputs_output_level_fno.end());
  if (input_aware_scheme_ == 1) {
    for (auto fno : input_fno) {
      auto zFile = GetSSTZoneFileInZBDNoLock(fno);
      if (!zFile) {
        zFile->selected_as_input_ = true;
      }
    }
  }

  double score = GetMaxSameZoneScore(input_fno);
  uint64_t none;
  double inval_score = GetMaxInvalidateCompactionScore(input_fno, &none, true);
  {
    std::lock_guard<std::mutex> lg(same_zone_score_mutex_);
    same_zone_score_[output_level].push_back(score);
    invalidate_score_[output_level].push_back(inval_score);

    same_zone_score_for_timelapse_[output_level].clear();
    same_zone_score_for_timelapse_[output_level] =
        same_zone_score_[output_level];
  }
}

void ZonedBlockDevice::AddTimeLapse(int T) {
  far_stats_.emplace_back(cur_free_percent_, reset_count_.load(),
                          wasted_wp_.load() / (1 << 20), T, reset_threshold_);
}

inline uint64_t ZonedBlockDevice::LazyLog(uint64_t sz, uint64_t fr,
                                          uint64_t T) {
  T++;
  if (fr >= T) {
    return 0 + (1 << 14);
  }
  return sz * (1 - log(fr + 1) / log(T));
}

inline uint64_t ZonedBlockDevice::LazyLinear(uint64_t sz, uint64_t fr,
                                             uint64_t T) {
  if (fr >= T) {
    return 0 + (1 << 14);
  }
  return sz - (sz * fr) / T;
}
inline uint64_t ZonedBlockDevice::Custom(uint64_t sz, uint64_t fr, uint64_t T) {
  if (fr >= T) {
    return sz - sz * T / 100;
  }
  return sz - (fr * sz) / 100;
}

inline uint64_t ZonedBlockDevice::LogLinear(uint64_t sz, uint64_t fr,
                                            uint64_t T) {
  double ratio;
  if (fr >= T) {
    ratio = (1 - log(fr + 1) / log(101));
    return ratio * sz;
  }
  ratio = (1 - log(T + 1) / log(101)) * 100;
  double slope = (100 - ratio) / T;
  ratio = 100 - slope * fr;

  return ratio * sz / 100;
}
inline uint64_t ZonedBlockDevice::LazyExponential(uint64_t sz, uint64_t fr,
                                                  uint64_t T) {
  if (fr >= T) {
    return 0 + (1 << 14);
  }
  if (fr == 0) {
    return sz;
  }
  double b = pow(100, 1 / (float)T);
  b = pow(b, fr);
  return sz - (b * sz / 100);
}

void ZonedBlockDevice::CalculateFinishThreshold(uint64_t free_percent) {
  uint64_t rt = 0;
  uint64_t max_capacity = io_zones[0]->max_capacity_;
  // uint64_t free_percent = cur_free_percent_;
  switch (finish_scheme_) {
    case FINISH_ENABLE:
      rt = max_capacity;
      break;
    case FINISH_DISABLE:
      rt = 0;
      break;
    case FINISH_PROPOSAL:  // Constant scale
                           // if very high free space ratio, no finish
                           // medium : do finish
                           // if very low free space ratio, no finish
      rt = max_capacity - (max_capacity * free_percent) / 100;
      // if(free_percent>50){
      //   rt=(max_capacity-rt);
      // }
      // rt = LazyLinear(max_capacity, free_percent, tuning_point_);
      break;
    // case kLazy_Log:
    //   rt = LazyLog(max_capacity, free_percent, tuning_point_);
    //   break;
    // case kNoRuntimeLinear:
    // case kLazy_Linear:
    // rt = LazyLinear(max_capacity, free_percent, tuning_point_);
    // break;
    // case kCustom:
    //   rt = Custom(max_capacity, free_percent, tuning_point_);
    //   break;
    // case kLogLinear:
    //   rt = LogLinear(max_capacity, free_percent, tuning_point_);
    //   break;
    // case kLazyExponential:
    //   rt = LazyExponential(max_capacity, free_percent, tuning_point_);
    // break;
    default:
      break;
  }
  // finish_threshold_ = rt;
  printf("%lu : %lu\n", free_percent, rt);
  finish_threshold_arr_[free_percent] = rt;
}

void ZonedBlockDevice::CalculateResetThreshold(uint64_t free_percent) {
  uint64_t rt = 0;
  uint64_t max_capacity = io_zones[0]->max_capacity_;
  // uint64_t free_percent = cur_free_percent_;
  switch (reset_scheme_) {
    case kEager:
      rt = max_capacity;
      break;
    case kLazy:
      rt = 0;
      break;
    case kFAR:  // Constant scale
      rt = max_capacity - (max_capacity * free_percent) / 100;
      break;
    case kLazy_Log:
      rt = LazyLog(max_capacity, free_percent, tuning_point_);
      break;
    case kNoRuntimeLinear:
    case kLazy_Linear:
      rt = LazyLinear(max_capacity, free_percent, tuning_point_);
      break;
    case kCustom:
      rt = Custom(max_capacity, free_percent, tuning_point_);
      break;
    case kLogLinear:
      rt = LogLinear(max_capacity, free_percent, tuning_point_);
      break;
    case kLazyExponential:
      rt = LazyExponential(max_capacity, free_percent, tuning_point_);
      break;
    default:
      break;
  }
  reset_threshold_ = rt;
  reset_threshold_arr_[free_percent] = rt;
}
/* io_zones 벡터의 각 존을 순회하며, 사용되지 않는 IO 존을 재설정합니다.
재설정이 완료되면, 필요에 따라 토큰을 반환합니다.*/
IOStatus ZonedBlockDevice::ResetUnusedIOZones() {
  clock_t reset_latency{0};

  for (size_t i = 0; i < io_zones.size(); i++) {
    const auto z = io_zones[i];
    // if (z->is_finished_) {
    //   printf("resetunused - is finished ? %d\n", z->is_finished_);
    // }

    bool full = z->IsFull();
    if (!full) {
      continue;
    }
    if (z->IsUsed()) {
      continue;
    }
    if (z->IsEmpty()) {
      continue;
    }

    if (z->Acquire()) {
      if (z->IsEmpty()) {
        z->Release();
        continue;
      }
      // if(!ProactiveZoneCleaning() && !full){
      //   z->Release();
      //   continue;
      // }
      if (!full) {
        z->Release();
        continue;
      }

      if (!z->IsUsed()) {
        // bool full=
        if (full) {
          erase_size_zc_.fetch_add(io_zones[i]->max_capacity_);
        } else {
          erase_size_proactive_zc_.fetch_add(io_zones[i]->wp_ -
                                             io_zones[i]->start_);
        }
        // wasted_wp_.fetch_add(io_zones[i]->capacity_);
        clock_t start = clock();
        IOStatus reset_status = z->Reset();
        uint64_t cp = z->GetCapacityLeft();

        wasted_wp_.fetch_add(cp);
        new_wasted_wp_.fetch_add(cp);
        // if (z->is_finished_) {
        //   printf("resetunued - finish++\n");
        //   finished_wasted_wp_.fetch_add(cp);
        //   finish_count_.fetch_add(1);
        //   z->finish_count_++;
        // }
        clock_t end = clock();
        reset_latency += (end - start);
        runtime_reset_reset_latency_.fetch_add(reset_latency);
        if (!reset_status.ok()) {
          z->Release();
          return reset_status;
        }
        reset_count_.fetch_add(1);
        // z->is_finished_ = false;
        if (!full) {
          PutActiveIOZoneToken();
          // PutOpenIOZoneToken();
        }
      }

      z->Release();
    }
  }
  return IOStatus::OK();  // 모든 존에 대한 작업이 성공적으로 완료되면
                          // IOStatus::OK()를 반환
}

IOStatus ZonedBlockDevice::RuntimeZoneReset() {
  size_t total_invalid = 0;

  // uint64_t zeu_size = 1 << 30;
  uint64_t zeu_size = io_zones[0]->max_capacity_;
  (void)(zeu_size);
  IOStatus reset_status = IOStatus::OK();
  for (size_t i = 0; i < io_zones.size(); i++) {
    const auto z = io_zones[i];
    // if (z->is_finished_) {
    //   printf("resetunused - is finished ? %d\n", z->is_finished_);
    // }
    // if (is_reseted[i]) {
    //   continue;
    // }
    if (z->IsEmpty()) {
      continue;
    }
    if (z->IsUsed()) {
      continue;
    }
    if (z->Acquire()) {
      if (z->IsEmpty()) {
        z->Release();
        continue;
      }
      if (z->IsUsed()) {
        z->Release();
        continue;
      }

      bool full = z->IsFull();
      uint64_t cp = z->GetCapacityLeft();

      total_invalid = z->wp_ - z->start_ < z->max_capacity_
                          ? (z->wp_ - z->start_)
                          : z->max_capacity_;

      if ((z->max_capacity_ - total_invalid) >
          reset_threshold_arr_[cur_free_percent_]) {
        goto no_reset;
      }
      erase_size_.fetch_add(total_invalid);
      if (total_invalid % zeu_size) {
        new_wasted_wp_.fetch_add(zeu_size - (total_invalid % zeu_size));
      }
      // if (z->is_finished_) {
      //   if (total_invalid % zeu_size) {
      //     printf("runtime - finish++\n");
      //     finished_wasted_wp_.fetch_add(zeu_size - (total_invalid %
      //     zeu_size)); finish_count_.fetch_add(1); z->finish_count_++;
      //   }
      // }

      reset_status = z->Reset();

      wasted_wp_.fetch_add(cp);
      if (cp) {
        reset_count_before_full_.fetch_add(1);
        reset_size_before_full_.fetch_add(cp);
      }
      if (!reset_status.ok()) return reset_status;
      // is_reseted[i] = true;
      // z->is_fin ished_ = false;
      reset_count_.fetch_add(1);
      z->reset_count_++;

    no_reset:
      if (!full && z->IsEmpty()) {  // not full -> empty
        PutActiveIOZoneToken();
      }
      z->Release();
    }
  }

  return IOStatus::OK();
}

void ZonedBlockDevice::WaitForOpenIOZoneToken(bool prioritized) {
  long allocator_open_limit;

  /* Avoid non-priortized allocators from starving prioritized ones */
  if (prioritized) {
    allocator_open_limit = max_nr_open_io_zones_;
  } else {
    allocator_open_limit = max_nr_open_io_zones_ - 1;
  }

  /* Wait for an open IO Zone token - after this function returns
   * the caller is allowed to write to a closed zone. The callee
   * is responsible for calling a PutOpenIOZoneToken to return the resource
   */
  std::unique_lock<std::mutex> lk(zone_resources_mtx_);
  zone_resources_.wait(lk, [this, allocator_open_limit] {
    if (open_io_zones_.load() < allocator_open_limit) {
      open_io_zones_++;
      return true;
    } else {
      return false;
    }
  });
}

IOStatus ZonedBlockDevice::RuntimeReset(void) {
  IOStatus s = IOStatus::OK();
  if (RuntimeZoneResetDisabled()) {
    printf("RuntimeZoneResetDisabled!!\n");
    return s;
  }
  if (ProactiveZoneCleaning()) {
    return s;
  }
  // std::vector<bool> is_reseted;
  // auto start_chrono = std::chrono::high_resolution_clock::now();
  // is_reseted.assign(io_zones.size(), false);
  switch (GetPartialResetScheme()) {
    case RUNTIME_ZONE_RESET_ONLY:
      // printf("RUNTIME_ZONE_RESET_ONLY!!\n");
      s = RuntimeZoneReset();

      break;
    default:
      break;
  }
  return s;
}

bool ZonedBlockDevice::GetActiveIOZoneTokenIfAvailable() {
  /* Grap an active IO Zone token if available - after this function returns
   * the caller is allowed to write to a closed zone. The callee
   * is responsible for calling a PutActiveIOZoneToken to return the resource
   */
  std::unique_lock<std::mutex> lk(zone_resources_mtx_);
  if (active_io_zones_.load() < max_nr_active_io_zones_) {
    active_io_zones_++;
    return true;
  }
  return false;
}

void ZonedBlockDevice::PutOpenIOZoneToken() {
  {
    std::unique_lock<std::mutex> lk(zone_resources_mtx_);
    open_io_zones_--;
  }
  zone_resources_.notify_one();
}

void ZonedBlockDevice::PutActiveIOZoneToken() {
  {
    std::unique_lock<std::mutex> lk(zone_resources_mtx_);
    active_io_zones_--;
  }
  zone_resources_.notify_one();
}

IOStatus ZonedBlockDevice::ApplyFinishThreshold() {
  IOStatus s;

  if (finish_threshold_ == 0) return IOStatus::OK();

  for (const auto z : io_zones) {
    if (z->Acquire()) {
      bool within_finish_threshold =
          z->capacity_ < (z->max_capacity_ * finish_threshold_ / 100);
      if (!(z->IsEmpty() || z->IsFull()) && within_finish_threshold) {
        /* If there is less than finish_threshold_% remaining capacity in a
         * non-open-zone, finish the zone */
        s = z->Finish();
        if (!s.ok()) {
          z->Release();
          Debug(logger_, "Failed finishing zone");
          return s;
        }
        s = z->CheckRelease();
        if (!s.ok()) return s;
        PutActiveIOZoneToken();
      } else {
        s = z->CheckRelease();
        if (!s.ok()) return s;
      }
    }
  }

  return IOStatus::OK();
}

// get mininal valid data with upper threshold
bool ZonedBlockDevice::FinishProposal(bool put_token) {
  IOStatus s;
  Zone *finish_victim = nullptr;
  uint64_t finish_threshold_now = finish_threshold_arr_[cur_free_percent_];
  // uint64_t finish_score;

  for (const auto z : io_zones) {
    if (z->Acquire()) {
      if (z->IsEmpty() || z->IsFull()) {
        s = z->CheckRelease();
        // if (!s.ok()) return s;
        if (!s.ok()) return false;
        continue;
      }
      if ((z->wp_ - z->start_) < finish_threshold_now) {
        z->Release();
        continue;
      }

      if (finish_victim == nullptr) {
        finish_victim = z;
        continue;
      }
      if (finish_victim->used_capacity_ > z->used_capacity_) {
        s = finish_victim->CheckRelease();
        // if (!s.ok()) return s;
        if (!s.ok()) return false;
        finish_victim = z;
      } else {
        s = z->CheckRelease();
        // if (!s.ok()) return s;
        if (!s.ok()) return false;
      }
    }
  }

  // If all non-busy zones are empty or full, we should return success.
  if (finish_victim == nullptr) {
    Info(logger_, "All non-busy zones are empty or full, skip.");
    // return IOStatus::OK();
    return false;
  }
  uint64_t valid_data = finish_victim->used_capacity_;
  uint64_t cp = finish_victim->capacity_;

  s = finish_victim->Finish();
  // IOStatus release_status =
  finish_victim->CheckRelease();

  // if (s.ok()) {
  //   PutActiveIOZoneToken();
  // }
  if (put_token) {
    PutActiveIOZoneToken();
  }

  finished_valid_data_.fetch_add(valid_data);
  finished_wasted_wp_.fetch_add(cp);
  finish_count_.fetch_add(1);

  return true;
}

bool ZonedBlockDevice::FinishProposal2(bool put_token) {
  IOStatus s;
  Zone *finish_victim = nullptr;

  for (const auto z : io_zones) {
    if (z->Acquire()) {
      if (z->IsEmpty() || z->IsFull()) {
        s = z->CheckRelease();
        // if (!s.ok()) return s;
        if (!s.ok()) return false;
        continue;
      }
      if (finish_victim == nullptr) {
        finish_victim = z;
        continue;
      }
      if (finish_victim->capacity_ > z->capacity_) {
        s = finish_victim->CheckRelease();
        // if (!s.ok()) return s;
        if (!s.ok()) return false;
        finish_victim = z;
      } else {
        s = z->CheckRelease();
        // if (!s.ok()) return s;
        if (!s.ok()) return false;
      }
    }
  }

  // If all non-busy zones are empty or full, we should return success.
  if (finish_victim == nullptr) {
    Info(logger_, "All non-busy zones are empty or full, skip.");
    // return IOStatus::OK();
    return false;
  }
  uint64_t valid_data = finish_victim->used_capacity_;
  uint64_t cp = finish_victim->capacity_;
  uint64_t written = finish_victim->wp_ - finish_victim->start_;
  // printf("written %lu finish_threshold_arr_[cur_free_percent_] %lu free
  // %lu\n",
  //   written>>20,finish_threshold_arr_[cur_free_percent_]>>20,cur_free_percent_);

  if (written < finish_threshold_arr_[cur_free_percent_]) {
    // printf("NO FINISH\n");
    finish_victim->CheckRelease();
    return false;
  }
  // printf("1 finish_victim->capacity_: %lu\n",
  //        finish_victim->capacity_ / (1 << 20));
  s = finish_victim->Finish();
  IOStatus release_status = finish_victim->CheckRelease();

  // if (s.ok()) {
  //   PutActiveIOZoneToken();
  // }
  if (put_token) {
    PutActiveIOZoneToken();
  }

  // if (!release_status.ok()) {
  //   return release_status;
  // }
  // uint64_t cp = finish_victim->capacity_;

  // printf("2 finish_victim->capacity_: %lu\n", cp / (1 << 20));
  // printf("After finish_victim->capacity_: %lu\n",
  //  finish_victim->capacity_ / (1 << 20));
  // finish_victim->is_finished_ = true;
  // printf("FINISH OK, %lu \n",cp>>20);
  finished_valid_data_.fetch_add(valid_data);
  finished_wasted_wp_.fetch_add(cp);
  finish_count_.fetch_add(1);
  // printf("Zone Finish!!! \n");
  // printf(
  //     "Finish complete: Zone start: 0x%lx, capacity left: %lu,
  //     open_io_zones_: "
  //     "%ld\n",
  //     finish_victim->start_, cp, open_io_zones_.load());

  return true;
}

bool ZonedBlockDevice::FinishCheapestIOZone(bool put_token) {
  IOStatus s;
  Zone *finish_victim = nullptr;

  for (const auto z : io_zones) {
    if (z->Acquire()) {
      if (z->IsEmpty() || z->IsFull()) {
        s = z->CheckRelease();
        // if (!s.ok()) return s;
        if (!s.ok()) return false;
        continue;
      }
      if (finish_victim == nullptr) {
        finish_victim = z;
        continue;
      }
      if (finish_victim->capacity_ > z->capacity_) {
        s = finish_victim->CheckRelease();
        // if (!s.ok()) return s;
        if (!s.ok()) return false;
        finish_victim = z;
      } else {
        s = z->CheckRelease();
        // if (!s.ok()) return s;
        if (!s.ok()) return false;
      }
    }
  }

  // If all non-busy zones are empty or full, we should return success.
  if (finish_victim == nullptr) {
    Info(logger_, "All non-busy zones are empty or full, skip.");
    // return IOStatus::OK();
    return false;
  }
  uint64_t valid_data = finish_victim->used_capacity_;
  uint64_t cp = finish_victim->GetCapacityLeft();
  // printf("1 finish_victim->capacity_: %lu\n",
  //        finish_victim->capacity_ / (1 << 20));
  s = finish_victim->Finish();
  IOStatus release_status = finish_victim->CheckRelease();

  // if (s.ok()) {
  //   PutActiveIOZoneToken();
  // }
  if (put_token) {
    PutActiveIOZoneToken();
  }

  // if (!release_status.ok()) {
  //   return release_status;
  // }
  // uint64_t cp = finish_victim->capacity_;

  // printf("2 finish_victim->capacity_: %lu\n", cp / (1 << 20));
  // printf("After finish_victim->capacity_: %lu\n",
  //  finish_victim->capacity_ / (1 << 20));
  // finish_victim->is_finished_ = true;
  finished_valid_data_.fetch_add(valid_data);
  finished_wasted_wp_.fetch_add(cp);
  finish_count_.fetch_add(1);

  // printf("Zone Finish!!! \n");
  // printf(
  //     "Finish complete: Zone start: 0x%lx, capacity left: %lu,
  //     open_io_zones_: "
  //     "%ld\n",
  //     finish_victim->start_, cp, open_io_zones_.load());

  return true;
}

IOStatus ZonedBlockDevice::GetAnyLargestRemainingZone(Zone **zone_out,
                                                      uint32_t min_capacity) {
  IOStatus s = IOStatus::OK();
  Zone *allocated_zone = nullptr;

  for (const auto z : io_zones) {
    if (!z->Acquire()) {
      continue;
    }
    if (z->IsEmpty()) {
      z->Release();
      continue;
    }
    if (z->capacity_ > min_capacity) {
      if (allocated_zone) {
        s = allocated_zone->CheckRelease();
        if (!s.ok()) {
          return s;
        }
      }
      allocated_zone = z;
      min_capacity = z->capacity_;
      continue;
    }

    s = z->CheckRelease();
    if (!s.ok()) {
      return s;
    }
  }

  *zone_out = allocated_zone;
  return s;
}

IOStatus ZonedBlockDevice::GetBestOpenZoneMatch(
    Env::WriteLifeTimeHint file_lifetime, unsigned int *best_diff_out,
    Zone **zone_out, uint32_t min_capacity) {
  unsigned int best_diff = LIFETIME_DIFF_NOT_GOOD;
  Zone *allocated_zone = nullptr;
  IOStatus s;

  for (const auto z : io_zones) {
    if (z->Acquire()) {
      if ((z->used_capacity_ > 0) && !z->IsFull() &&
          z->capacity_ >= min_capacity) {
        unsigned int diff = GetLifeTimeDiff(z->lifetime_, file_lifetime);

        if (diff == 0) {
          allocated_zone = z;
          best_diff = diff;

          break;
        }
        // else{
        z->Release();
        continue;
        // }

        // if (diff < best_diff) {
        //   if (allocated_zone != nullptr) {
        //     s = allocated_zone->CheckRelease();
        //     if (!s.ok()) {
        //       IOStatus s_ = z->CheckRelease();
        //       if (!s_.ok()) return s_;
        //       return s;
        //     }
        //   }
        //   allocated_zone = z;
        //   best_diff = diff;
        // } else {
        //   s = z->CheckRelease();
        //   if (!s.ok()) return s;
        // }
      }
      // else {

      s = z->CheckRelease();
      if (!s.ok()) return s;
      // }
    }
  }

  *best_diff_out = best_diff;
  *zone_out = allocated_zone;

  return IOStatus::OK();
}

// IOStatus ZonedBlockDevice::GetBestOpenZoneMatch(
//     Env::WriteLifeTimeHint file_lifetime, unsigned int *best_diff_out,
//     Zone **zone_out, uint32_t min_capacity) {
//   unsigned int best_diff = LIFETIME_DIFF_NOT_GOOD;
//   Zone *allocated_zone = nullptr;
//   IOStatus s;

//   for (const auto z : io_zones) {
//     if (z->Acquire()) {
//       if ((z->used_capacity_ > 0) && !z->IsFull() &&
//           z->capacity_ >= min_capacity) {
//         unsigned int diff = GetLifeTimeDiff(z->lifetime_, file_lifetime);
//         if (diff < best_diff) {
//           if (allocated_zone != nullptr) {
//             s = allocated_zone->CheckRelease();
//             if (!s.ok()) {
//               IOStatus s_ = z->CheckRelease();
//               if (!s_.ok()) return s_;
//               return s;
//             }
//           }
//           allocated_zone = z;
//           best_diff = diff;
//         } else {
//           s = z->CheckRelease();
//           if (!s.ok()) return s;
//         }
//       } else {
//         s = z->CheckRelease();
//         if (!s.ok()) return s;
//       }
//     }
//   }

//   *best_diff_out = best_diff;
//   *zone_out = allocated_zone;

//   return IOStatus::OK();
// }

IOStatus ZonedBlockDevice::AllocateEmptyZone(Zone **zone_out) {
  IOStatus s;
  Zone *allocated_zone = nullptr;
  for (const auto z : io_zones) {
    if (z->Acquire()) {
      if (z->IsEmpty()) {
        allocated_zone = z;
        break;
      } else {
        s = z->CheckRelease();
        if (!s.ok()) return s;
      }
    }
  }
  *zone_out = allocated_zone;
  return IOStatus::OK();
}

IOStatus ZonedBlockDevice::InvalidateCache(uint64_t pos, uint64_t size) {
  int ret = zbd_be_->InvalidateCache(pos, size);

  if (ret) {
    return IOStatus::IOError("Failed to invalidate cache");
  }
  return IOStatus::OK();
}

int ZonedBlockDevice::Read(char *buf, uint64_t offset, int n, bool direct) {
  int ret = 0;
  int left = n;
  int r = -1;

  while (left) {
    r = zbd_be_->Read(buf, left, offset, direct);
    if (r <= 0) {
      if (r == -1 && errno == EINTR) {
        continue;
      }
      break;
    }
    ret += r;
    buf += r;
    left -= r;
    offset += r;
  }

  if (r < 0) return r;
  return ret;
}

IOStatus ZonedBlockDevice::ReleaseMigrateZone(Zone *zone) {
  IOStatus s = IOStatus::OK();
  {
    std::unique_lock<std::mutex> lock(migrate_zone_mtx_);
    migrating_ = false;
    if (zone != nullptr) {
      bool full = zone->IsFull();
      if (zone->this_zone_motivation_check_) {
        struct timespec timespec;
        clock_gettime(CLOCK_MONOTONIC, &timespec);
        if (zone->is_allocated_ == false) {
          zone->is_allocated_ = true;
          zone->allocated_time_ = timespec;
        }
      }
      s = zone->CheckRelease();
      // PutOpenIOZoneToken();
      if (full) {
        PutActiveIOZoneToken();
      }
      Info(logger_, "ReleaseMigrateZone: %lu", zone->start_);
    }
  }
  migrate_resource_.notify_one();
  return s;
}

IOStatus ZonedBlockDevice::TakeMigrateZone(Zone **out_zone,
                                           Env::WriteLifeTimeHint file_lifetime,
                                           uint32_t min_capacity) {
  std::unique_lock<std::mutex> lock(migrate_zone_mtx_);
  migrate_resource_.wait(lock, [this] { return !migrating_; });
  IOStatus s;
  migrating_ = true;

  unsigned int best_diff = LIFETIME_DIFF_NOT_GOOD;

  while (CalculateCapacityRemain() > min_capacity) {
    // for (int i = 0; i < 2; i++) {
    s = GetBestOpenZoneMatch(file_lifetime, &best_diff, out_zone, min_capacity);

    if (s.ok() && (*out_zone) != nullptr) {
      Info(logger_, "TakeMigrateZone: %lu", (*out_zone)->start_);
      break;
    } else {
      s = GetAnyLargestRemainingZone(out_zone, min_capacity);
    }

    if (s.ok() && (*out_zone) != nullptr) {
      Info(logger_, "TakeMigrateZone: %lu", (*out_zone)->start_);
      break;
    }
    if (GetActiveIOZoneTokenIfAvailable()) {
      s = AllocateEmptyZone(out_zone);
      if (s.ok() && (*out_zone) != nullptr) {
        Info(logger_, "TakeMigrateZone: %lu", (*out_zone)->start_);
        break;
      }
      if ((*out_zone) == nullptr) {
        PutActiveIOZoneToken();
      }
    }

    s = GetAnyLargestRemainingZone(out_zone, min_capacity);
    if (s.ok() && (*out_zone) != nullptr) {
      Info(logger_, "TakeMigrateZone: %lu", (*out_zone)->start_);

      break;
    }

    s = ResetUnusedIOZones();
    if (!s.ok()) {
      return s;
    }
  }

  if (s.ok() && (*out_zone) != nullptr) {
    Info(logger_, "TakeMigrateZone: %lu", (*out_zone)->start_);
  } else {
    migrating_ = false;
  }

  return s;
}

IOStatus ZonedBlockDevice::TakeMigrateZone(Slice &smallest, Slice &largest,
                                           int level, Zone **out_zone,
                                           Env::WriteLifeTimeHint file_lifetime,
                                           uint64_t file_size,
                                           uint64_t min_capacity,
                                           bool *run_gc_worker_, bool is_sst) {
  std::unique_lock<std::mutex> lock(migrate_zone_mtx_);
  migrate_resource_.wait(lock, [this] { return !migrating_; });
  IOStatus s;
  migrating_ = true;

  unsigned int best_diff = LIFETIME_DIFF_NOT_GOOD;

  (void)(run_gc_worker_);

  if (db_ptr_ == nullptr) {
    return s;
  }

  while (CalculateCapacityRemain() > min_capacity) {
    if (is_sst) {
      if (allocation_scheme_ == CAZA) {
        // printf("AllocateCompactionAwaredZone\n");
        s = AllocateCompactionAwaredZone(smallest, largest, level,
                                         file_lifetime, file_size, out_zone,
                                         min_capacity);
      } else if (allocation_scheme_ == CAZA_ADV) {
        // printf("AllocateCompactionAwaredZoneV2\n");
        s = AllocateCompactionAwaredZoneV2(smallest, largest, level,
                                           file_lifetime, file_size, out_zone,
                                           min_capacity);
        // printf("AllocateCompactionAwaredZoneV2 - finished!!\n");
        if (s.ok() && (*out_zone) != nullptr) {
          // printf("AllocateCompactionAwaredZoneV2 - successed!!\n");
          break;
        }
        s = AllocateCompactionAwaredZone(smallest, largest, level,
                                         file_lifetime, file_size, out_zone,
                                         min_capacity);
        // if (GetActiveIOZoneTokenIfAvailable()) {
        //   s = AllocateEmptyZone(out_zone);
        //   if (s.ok() && (*out_zone) != nullptr) {
        //     printf("CAZA2-AllocateEmptyZone\n");
        //     Info(logger_, "TakeMigrateZone: %lu", (*out_zone)->start_);
        //     (*out_zone)->lifetime_ = file_lifetime;
        //     break;
        //   } else {
        //     PutActiveIOZoneToken();
        //     printf("CAZA2-PutActiveIOZoneToken\n");
        //   }
        // } else {
        //   AllocateAllInvalidZone(out_zone);
        //   printf("CAZA2-AllocateAllInvalidZone\n");
        // }

      } else {
        // printf("I am LIZA!\n");
      }

      if (s.ok() && (*out_zone) != nullptr) {
        // printf(
        //     "AllocateCompactionAwaredZone - successed!! : min_capacity :
        //     %lu\n", min_capacity);
        break;
      }
    }

    s = GetBestOpenZoneMatch(file_lifetime, &best_diff, out_zone, min_capacity);

    if (s.ok() && (*out_zone) != nullptr) {
      Info(logger_, "TakeMigrateZone: %lu", (*out_zone)->start_);
      // printf("GetBest : min_capacity : %lu\n", min_capacity);
      break;
    }

    if (finish_scheme_ != FINISH_ENABLE) {
      AllocateAllInvalidZone(out_zone);
      if (*out_zone) {
        break;
      }
      GetAnyLargestRemainingZone(out_zone, min_capacity);
      if (*out_zone) {
        break;
      }

      // if(GetFullZoneN()>io_zones.size()-zc_-4){
      if (GetActiveIOZoneTokenIfAvailable()) {
        AllocateEmptyZone(out_zone);  // 빈 영역 할당
        if (*out_zone != nullptr) {
          (*out_zone)->lifetime_ = file_lifetime;
          break;
        }
        PutActiveIOZoneToken();
      }
      // }

    } else if (finish_scheme_ == FINISH_ENABLE) {
      while (true) {
        if (GetActiveIOZoneTokenIfAvailable()) {
          break;
        }
        if (FinishCheapestIOZone(false)) {
          break;
        }
      }
      s = AllocateEmptyZone(out_zone);  // 빈 영역 할당
      if (*out_zone != nullptr) {
        // PutActiveIOZoneToken();
        (*out_zone)->lifetime_ = file_lifetime;
        break;
      }
      PutActiveIOZoneToken();

      AllocateAllInvalidZone(out_zone);
      if (*out_zone) {
        break;
      }
      GetAnyLargestRemainingZone(out_zone, min_capacity);
      if (*out_zone) {
        break;
      }
    }

    // std::cout << "finish_scheme_: " << finish_scheme_ << std::endl;
    // if (!finish_scheme_) {
    //   if (!GetActiveIOZoneTokenIfAvailable()) {
    //     printf("Takemigrate - finish!!\n");
    //     s = FinishCheapestIOZone(false);
    //     // s = FinishCheapestIOZone(true);
    //     if (!s.ok()) {
    //       PutOpenIOZoneToken();
    //     }
    //     s = AllocateEmptyZone(out_zone);
    //     if (s.ok() && (*out_zone) != nullptr) {
    //       Info(logger_, "TakeMigrateZone: %lu", (*out_zone)->start_);
    //       (*out_zone)->lifetime_ = file_lifetime;
    //       break;
    //     } else {
    //       printf("Before PutActiveIOZoneToken: open_io_zones_=%ld\n",
    //              open_io_zones_.load());
    //       PutActiveIOZoneToken();
    //       printf("After PutActiveIOZoneToken: open_io_zones_=%ld\n",
    //              open_io_zones_.load());
    //     }
    //   }
    //   s = AllocateEmptyZone(out_zone);

    // } else {
    //   s = AllocateAllInvalidZone(out_zone);
    //   if (s.ok() && (*out_zone) != nullptr) {
    //     printf("takemigrate - allocateallinvalidzone!!\n");
    //     break;
    //   }
    // }
    // if (s.ok() && (*out_zone) != nullptr) {
    //   Info(logger_, "TakeMigrateZone: %lu", (*out_zone)->start_);
    //   // printf("Empty: min_capacity : %lu\n", min_capacity);
    //   break;
    // }

    // s = GetAnyLargestRemainingZone(out_zone, min_capacity);
    // if (s.ok() && (*out_zone) != nullptr) {
    //   Info(logger_, "TakeMigrateZone: %lu", (*out_zone)->start_);
    //   // printf("2Getany: min_capacity : %lu\n", min_capacity);
    //   break;
    // }
    // reset:
    s = ResetUnusedIOZones();
    if (!s.ok()) {
      return s;
    }
  }

  if (s.ok() && (*out_zone) != nullptr) {
    // printf("LAST: min_capacity : %lu\n", min_capacity);
    Info(logger_, "TakeMigrateZone: %lu", (*out_zone)->start_);
  } else {
    migrating_ = false;
  }

  return s;
}

// IOStatus ZonedBlockDevice::TakeMigrateZone(Slice &smallest, Slice &largest,
//                                            int level, Zone **out_zone,
//                                            Env::WriteLifeTimeHint
//                                            file_lifetime, uint64_t file_size,
//                                            uint64_t min_capacity,
//                                            bool *run_gc_worker_, bool is_sst)
//                                            {
//   std::unique_lock<std::mutex> lock(migrate_zone_mtx_);
//   migrate_resource_.wait(lock, [this] { return !migrating_; });
//   // printf("#############TakeMigateZone!!!\n");

//   IOStatus s;
//   migrating_ = true;

//   int blocking_time = 0;
//   unsigned int best_diff = LIFETIME_DIFF_NOT_GOOD;

//   (void)(smallest);
//   (void)(largest);
//   (void)(level);
//   std::vector<uint64_t> none;

//   WaitForOpenIOZoneToken(false);

//   while (CalculateCapacityRemain() > min_capacity) {
//     if ((*run_gc_worker_) == false) {
//       migrating_ = false;
//       printf("migrating_ false\n");
//       return IOStatus::OK();
//     }
//     if (is_sst) {
//       if (allocation_scheme_ == CAZA) {
//         AllocateCompactionAwaredZone(smallest, largest, level, file_lifetime,
//                                      std::vector<uint64_t>(0), file_size,
//                                      out_zone, min_capacity);
//       } else if (allocation_scheme_ == CAZA_ADV) {
//         printf("AllocateCompactionAwaredZoneV2\n");
//         AllocateCompactionAwaredZoneV2(smallest, largest, level,
//         file_lifetime,
//                                        std::vector<uint64_t>(0), file_size,
//                                        out_zone, min_capacity);
//         if (s.ok() && (*out_zone) != nullptr) {
//           break;
//         }

//         if (GetActiveIOZoneTokenIfAvailable()) {
//           s = AllocateEmptyZone(out_zone);
//           if (s.ok() && (*out_zone) != nullptr) {
//             Info(logger_, "TakeMigrateZone: %lu", (*out_zone)->start_);
//             (*out_zone)->lifetime_ = file_lifetime;
//             break;
//           } else {
//             PutActiveIOZoneToken();
//           }
//         } else {
//           AllocateAllInvalidZone(out_zone);
//         }

//       } else {
//         // printf("I am LIZA!\n");
//       }

//       if (s.ok() && (*out_zone) != nullptr) {
//         break;
//       }
//     }

//     s = GetBestOpenZoneMatch(file_lifetime, &best_diff, out_zone,
//     min_capacity);

//     if (s.ok() && (*out_zone) != nullptr) {
//       Info(logger_, "TakeMigrateZone: %lu", (*out_zone)->start_);
//       break;
//     }

//     s = GetAnyLargestRemainingZone(out_zone, min_capacity);
//     if (s.ok() && (*out_zone) != nullptr) {
//       Info(logger_, "TakeMigrateZone: %lu", (*out_zone)->start_);
//       break;
//     }

//     if (GetActiveIOZoneTokenIfAvailable()) {
//       s = AllocateEmptyZone(out_zone);
//       if (s.ok() && (*out_zone) != nullptr) {
//         Info(logger_, "TakeMigrateZone: %lu", (*out_zone)->start_);
//         (*out_zone)->lifetime_ = file_lifetime;
//         break;
//       } else {
//         PutActiveIOZoneToken();
//       }
//     }

//     s = ResetUnusedIOZones();

//     blocking_time++;

//     if (!s.ok()) {
//       return s;
//     }
//   }

//   if (s.ok() && (*out_zone) != nullptr) {
//     Info(logger_, "TakeMigrateZone: %lu", (*out_zone)->start_);
//   } else {
//     migrating_ = false;
//   }
//   // (*out_zone)->state_ = Zone::State::OPEN;
//   return s;
// }

IOStatus ZonedBlockDevice::AllocateIOZone(Env::WriteLifeTimeHint file_lifetime,
                                          IOType io_type, Zone **out_zone) {
  // IOStatus s;
  // s = ResetUnusedIOZones();
  Zone *allocated_zone = nullptr;
  unsigned int best_diff = LIFETIME_DIFF_NOT_GOOD;  // 수명차이
  int new_zone = 0;                                 // 새로운 영역인지 여부
  IOStatus s;

  // std::cout << "@@@ zbd::AllocateIOZone - life_time: " << file_lifetime
  //           << "// out_zone: " << out_zone << "\n";

  // I/O 유형에 따라 적절한 추적 태그 설정
  auto tag = ZENFS_WAL_IO_ALLOC_LATENCY;
  if (io_type != IOType::kWAL) {
    // L0 플러시는 중간 수명을 갖음
    if (file_lifetime == Env::WLTH_MEDIUM) {
      tag = ZENFS_L0_IO_ALLOC_LATENCY;
    } else {
      tag = ZENFS_NON_WAL_IO_ALLOC_LATENCY;
    }
  }
  // Latency Guard 생성
  ZenFSMetricsLatencyGuard guard(metrics_, tag, Env::Default());
  metrics_->ReportQPS(ZENFS_IO_ALLOC_QPS, 1);  // QPS(초당 작업 수) 리포트

  // Deferred I/O 에러가 있는지 확인
  s = GetZoneDeferredStatus();
  if (!s.ok()) {
    return s;
  }
  // WAL이 아닌 경우 Finish Threshold 적용
  // if (io_type != IOType::kWAL) {
  //   s = ApplyFinishThreshold();
  //   if (!s.ok()) {
  //     return s;
  //   }
  // }

  WaitForOpenIOZoneToken(io_type == IOType::kWAL);  // I/O 토큰 대기

  // allocatecompactionawaredzone //

  /* Try to fill an already open zone(with the best life time diff) */
  // 열린 영역에서 수명과 가장 잘 맞는 영역을 찾음
  s = GetBestOpenZoneMatch(file_lifetime, &best_diff, &allocated_zone);
  if (!s.ok()) {
    PutOpenIOZoneToken();  // 오류 발생 시 토큰 반환
    return s;
  }

  // Holding allocated_zone if != nullptr
  // 최적 영역이 없으면 새로 할당
  if (best_diff >= LIFETIME_DIFF_COULD_BE_WORSE) {
    bool got_token = GetActiveIOZoneTokenIfAvailable();  // 새로운 영역 열기

    /* If we did not get a token, try to use the best match, even if the life
     * time diff not good but a better choice than to finish an existing zone
     * and open a new one
     */
    // 이미 할당된 영역이 없으면 새 영역을 할당
    if (allocated_zone != nullptr) {
      if (!got_token && best_diff == LIFETIME_DIFF_COULD_BE_WORSE) {
        Debug(logger_,
              "Allocator: avoided a finish by relaxing lifetime diff "
              "requirement\n");
      } else {
        s = allocated_zone->CheckRelease();
        if (!s.ok()) {
          PutOpenIOZoneToken();
          if (got_token) PutActiveIOZoneToken();
          return s;
        }
        allocated_zone = nullptr;
      }
    }

    /* If we haven't found an open zone to fill, open a new zone */
    if (allocated_zone == nullptr) {
      /* We have to make sure we can open an empty zone */
      while (!got_token && !GetActiveIOZoneTokenIfAvailable()) {
        // s = FinishCheapestIOZone();
        if (!s.ok()) {
          PutOpenIOZoneToken();
          return s;
        }
      }

      s = AllocateEmptyZone(&allocated_zone);  // 빈 영역 할당
      //
      if (s.ok() && allocated_zone == nullptr) {
        s = GetAnyLargestRemainingZone(&allocated_zone);
      }
      //
      if (!s.ok()) {
        PutActiveIOZoneToken();
        PutOpenIOZoneToken();
        return s;
      }

      if (allocated_zone != nullptr) {
        // assert(allocated_zone->IsBusy());
        allocated_zone->lifetime_ = file_lifetime;  // 새 영역에 수명 설정
        new_zone = true;
      } else {
        PutActiveIOZoneToken();
      }
    }
  }
  // 할당된 영역 정보 출력
  if (allocated_zone) {
    // assert(allocated_zone->IsBusy());
    Debug(logger_,
          "Allocating zone(new=%d) start: 0x%lx wp: 0x%lx lt: %d file lt: %d\n",
          new_zone, allocated_zone->start_, allocated_zone->wp_,
          allocated_zone->lifetime_, file_lifetime);
  } else {
    PutOpenIOZoneToken();  // 할당 실패 시 토큰 반환
  }

  if (io_type != IOType::kWAL) {
    LogZoneStats();  // WAL이 아닐 경우 영역 통계 로깅
  }

  *out_zone = allocated_zone;

  metrics_->ReportGeneral(ZENFS_OPEN_ZONES_COUNT, open_io_zones_);
  metrics_->ReportGeneral(ZENFS_ACTIVE_ZONES_COUNT, active_io_zones_);

  return IOStatus::OK();
}

IOStatus ZonedBlockDevice::AllocateIOZone(
    std::string fname, bool is_sst, Slice &smallest, Slice &largest, int level,
    Env::WriteLifeTimeHint file_lifetime, IOType io_type,
    uint64_t predicted_size, Zone **out_zone, uint64_t min_capacity) {
  Zone *allocated_zone = nullptr;
  unsigned int best_diff = LIFETIME_DIFF_NOT_GOOD;  // 수명차이
  // int new_zone = 0;  // 새로운 영역인지 여부
  IOStatus s;

  // (void)(input_fno);
  (void)(fname);
  // printf("predicted_size : %lu\n",predicted_size);
  // I/O 유형에 따라 적절한 추적 태그 설정
  auto tag = ZENFS_WAL_IO_ALLOC_LATENCY;
  if (io_type != IOType::kWAL) {
    // L0 플러시는 중간 수명을 갖음
    if (file_lifetime == Env::WLTH_MEDIUM) {
      tag = ZENFS_L0_IO_ALLOC_LATENCY;
    } else {
      tag = ZENFS_NON_WAL_IO_ALLOC_LATENCY;
    }
  }
  // Latency Guard 생성
  ZenFSMetricsLatencyGuard guard(metrics_, tag, Env::Default());
  metrics_->ReportQPS(ZENFS_IO_ALLOC_QPS, 1);  // QPS(초당 작업 수) 리포트

  // Deferred I/O 에러가 있는지 확인
  s = GetZoneDeferredStatus();
  if (!s.ok()) {
    return s;
  }
  // WAL이 아닌 경우 Finish Threshold 적용
  // if (io_type != IOType::kWAL) {
  //   s = ApplyFinishThreshold();
  //   if (!s.ok()) {
  //     return s;
  //   }
  // }

  WaitForOpenIOZoneToken(io_type == IOType::kWAL);  // I/O 토큰 대기

  //////////////////// allocatecompactionawaredzone ///////////////////
  if (is_sst && level >= 0 && allocation_scheme_ == CAZA) {
    s = AllocateCompactionAwaredZone(smallest, largest, level, file_lifetime,
                                     predicted_size, &allocated_zone,
                                     min_capacity);
    // printf("allocateiozone-AllocateCompactionAwaredZone\n");

    if (!s.ok()) {
      PutOpenIOZoneToken();
      return s;
    }
  } else if (is_sst && level >= 0 && allocation_scheme_ == CAZA_ADV) {
    s = AllocateCompactionAwaredZoneV2(smallest, largest, level, file_lifetime,
                                       predicted_size, &allocated_zone,
                                       min_capacity);
    // // printf("allocateiozone-AllocateCompactionAwaredZoneV2\n");
    if (!s.ok()) {
      // printf("allocateiozone-putopen!!\n");
      PutOpenIOZoneToken();
      return s;
    }
    if (allocated_zone == nullptr) {
      s = AllocateCompactionAwaredZone(smallest, largest, level, file_lifetime,
                                       predicted_size, &allocated_zone,
                                       min_capacity);
    }
    // if (allocated_zone != nullptr) {
    //   // printf("allocateiozone-go to end\n");
    //   goto end;
    // }

    // if (!GetActiveIOZoneTokenIfAvailable()) {
    //   // printf("allocateiozone-finishchepest!!\n");
    //   // FinishCheapestIOZone(false);
    //   FinishCheapestIOZone(false);
    // }
    // s = AllocateEmptyZone(&allocated_zone);
    // // printf("allocateiozone-allocate emptyzone!!\n");
    // if (allocated_zone == nullptr) {
    //   // printf("allocateiozone-allocate emptyzone - putactive!!\n");
    //   PutActiveIOZoneToken();
    // }
  }

  if (allocated_zone != nullptr) {
    goto end;
  }
  //////////////////////////////////////////////////////////////////////

  /* Try to fill an already open zone(with the best life time diff) */
  // 열린 영역에서 수명과 가장 잘 맞는 영역을 찾음
  s = GetBestOpenZoneMatch(file_lifetime, &best_diff, &allocated_zone);
  if (!s.ok()) {
    PutOpenIOZoneToken();  // 오류 발생 시 토큰 반환
    return s;
  }
  // defauLt

  if (allocated_zone != nullptr) {
    if (best_diff >= LIFETIME_DIFF_COULD_BE_WORSE) {
      allocated_zone->CheckRelease();
      allocated_zone = nullptr;
    }
  }

  if (allocated_zone == nullptr) {
    if (finish_scheme_ == FINISH_DISABLE) {
      // printf("FINISH_DISABLE\n");
      if (GetActiveIOZoneTokenIfAvailable()) {
        AllocateEmptyZone(&allocated_zone);  // 빈 영역 할당
        if (allocated_zone != nullptr) {
          allocated_zone->lifetime_ = file_lifetime;
          goto end;
        }
        PutActiveIOZoneToken();
      }

      AllocateAllInvalidZone(&allocated_zone);
      if (allocated_zone) {
        goto end;
      }
      GetAnyLargestRemainingZone(&allocated_zone);
      if (allocated_zone) {
        goto end;
      }
      PutOpenIOZoneToken();

      return IOStatus::OK();

    } else if (finish_scheme_ == FINISH_ENABLE) {
      // printf("FINISH_ENABLE\n");
      while (true) {
        if (GetActiveIOZoneTokenIfAvailable()) {
          break;
        }
        if (FinishCheapestIOZone(false)) {
          break;
        }
      }
      // 무조건 active 있음
      s = AllocateEmptyZone(&allocated_zone);  // 빈 영역 할당
      if (allocated_zone != nullptr) {
        // PutActiveIOZoneToken();
        allocated_zone->lifetime_ = file_lifetime;
        goto end;
      }
      PutActiveIOZoneToken();

      AllocateAllInvalidZone(&allocated_zone);
      if (allocated_zone) {
        goto end;
      }
      GetAnyLargestRemainingZone(&allocated_zone);
      if (allocated_zone) {
        goto end;
      }
      PutOpenIOZoneToken();
      return IOStatus::OK();
    } else {  // finish_scheme_ == FINISH_PROPOSAL 2 or 3
              // printf("FINISH_PROPOSAL\n");
      AllocateEmptyZone(&allocated_zone);
      if (allocated_zone != nullptr) {
        if (GetActiveIOZoneTokenIfAvailable()) {
          goto end;
        }
        // if(level>1){
        // if(PredictCompactionScore(level)>1.0){ // is cold level
        // if(finish_scheme_==FINISH_PROPOSAL){
        if (FinishProposal(false)) {
          goto end;
        }
        // }else if(finish_scheme_==FINISH_PROPOSAL2){

        // }
        // if(FinishProposal2(false)){
        //     goto end;
        // }
        // }
        // }
        allocated_zone->Release();
        allocated_zone = nullptr;
      }
      // {
      AllocateAllInvalidZone(&allocated_zone);
      if (allocated_zone) {
        goto end;
      }
      GetAnyLargestRemainingZone(&allocated_zone);
      if (allocated_zone) {
        goto end;
      }
      PutOpenIOZoneToken();
      return IOStatus::OK();
      // }
    }
  }

  /////////////////////////////////////
  // // if (best_diff >= LIFETIME_DIFF_COULD_BE_WORSE) {
  //   bool got_token = GetActiveIOZoneTokenIfAvailable();  // 새로운 영역 열기

  //   /* If we did not get a token, try to use the best match, even if the life
  //    * time diff not good but a better choice than to finish an existing zone
  //    * and open a new one
  //    */
  //   // 이미 할당된 영역이 없으면 새 영역을 할당
  //   if (allocated_zone != nullptr) {
  //     if (!got_token && best_diff == LIFETIME_DIFF_COULD_BE_WORSE) {
  //       Debug(logger_,
  //             "Allocator: avoided a finish by relaxing lifetime diff "
  //             "requirement\n");
  //     } else {
  //       // s = allocated_zone->CheckRelease();
  //       // if (!s.ok()) {
  //       //   PutOpenIOZoneToken();
  //       //   if (got_token) PutActiveIOZoneToken();
  //       //   return s;
  //       // }
  //       allocated_zone = nullptr;
  //     }
  //   }
  //   if (allocated_zone == nullptr) {
  //     /* We have to make sure we can open an empty zone */
  //     // while (!got_token && !GetActiveIOZoneTokenIfAvailable()) {
  //     //   printf("allocateiozone - finish!!\n");
  //     //   s = FinishCheapestIOZone(false);
  //     //   if (!s.ok()) {
  //     //     PutOpenIOZoneToken();
  //     //     return s;
  //     //   }
  //     // }
  //     if (!finish_scheme_) {
  //       if (!GetActiveIOZoneTokenIfAvailable()) {
  //         s = FinishCheapestIOZone(false);
  //         // s = FinishCheapestIOZone(true);
  //         if (!s.ok()) {
  //           PutOpenIOZoneToken();
  //           return s;
  //         }
  //       }

  //       s = AllocateEmptyZone(&allocated_zone);  // 빈 영역 할당
  //       if ()
  //         if allocate
  //           zone fail, put active token
  //     } else {
  //       s = AllocateAllInvalidZone(&allocated_zone);
  //       if (s.ok() && allocated_zone != nullptr) {
  //         printf("allocateiozone - allocateallinvalidzone!!\n");
  //       }
  //     }
  //     //
  //     if (s.ok() && allocated_zone == nullptr) {
  //       s = GetAnyLargestRemainingZone(&allocated_zone);
  //     }

  //     if (allocated_zone != nullptr) {
  //       // assert(allocated_zone->IsBusy());
  //       // allocated_zone->lifetime_ = file_lifetime;  // 새 영역에 수명 설정
  //       // new_zone = true;
  //     } else {
  //       PutActiveIOZoneToken();
  //     }
  //   }
  // }

  // // if (allocated_zone == nullptr && GetActiveIOZoneTokenIfAvailable()) {
  // //   s = AllocateEmptyZone(&allocated_zone);
  // //   if (allocated_zone != nullptr && s.ok()) {
  // //     // assert(allocated_zone->IsBusy());
  // //     allocated_zone->lifetime_ = file_lifetime;
  // //     new_zone = true;
  // //   } else {
  // //     PutActiveIOZoneToken();
  // //   }
  // // }
  // // if (s.ok() && allocated_zone == nullptr) {
  // //   s = GetAnyLargestRemainingZone(&allocated_zone, min_capacity);
  // //   if (allocated_zone) {
  // //     allocated_zone->lifetime_ = file_lifetime;
  // //   }
  // // }

  // if (allocated_zone) {
  //   Debug(logger_,
  //         "Allocating zone(new=%d) start: 0x%lx wp: 0x%lx lt: %d file lt:
  //         %d\n", new_zone, allocated_zone->start_, allocated_zone->wp_,
  //         allocated_zone->lifetime_, file_lifetime);
  // } else {
  //   PutOpenIOZoneToken();
  // }

  // if (io_type != IOType::kWAL) {
  //   LogZoneStats();
  // }
  /////////////////////////////////////////////////////

end:
  // printf("allocateiozone - end!!\n");
  *out_zone = allocated_zone;

  // printf("ZENFS_ACTIVE_ZONES_COUNT: %ld\n", active_io_zones_.load());
  // printf("ZENFS_OPEN_ZONES_COUNT: %ld\n", open_io_zones_.load());

  metrics_->ReportGeneral(ZENFS_OPEN_ZONES_COUNT, open_io_zones_);
  metrics_->ReportGeneral(ZENFS_ACTIVE_ZONES_COUNT, active_io_zones_);

  return IOStatus::OK();
}

std::string ZonedBlockDevice::GetFilename() { return zbd_be_->GetFilename(); }

uint32_t ZonedBlockDevice::GetBlockSize() { return zbd_be_->GetBlockSize(); }

uint64_t ZonedBlockDevice::GetZoneSize() { return zbd_be_->GetZoneSize(); }

uint32_t ZonedBlockDevice::GetNrZones() { return zbd_be_->GetNrZones(); }

void ZonedBlockDevice::EncodeJsonZone(std::ostream &json_stream,
                                      const std::vector<Zone *> zones) {
  bool first_element = true;
  json_stream << "[";
  for (Zone *zone : zones) {
    if (first_element) {
      first_element = false;
    } else {
      json_stream << ",";
    }
    zone->EncodeJson(json_stream);
  }

  json_stream << "]";
}

void ZonedBlockDevice::EncodeJson(std::ostream &json_stream) {
  json_stream << "{";
  json_stream << "\"meta\":";
  EncodeJsonZone(json_stream, meta_zones);
  json_stream << ",\"io\":";
  EncodeJsonZone(json_stream, io_zones);
  json_stream << "}";
}

IOStatus ZonedBlockDevice::GetZoneDeferredStatus() {
  std::lock_guard<std::mutex> lock(zone_deferred_status_mutex_);
  return zone_deferred_status_;
}

void ZonedBlockDevice::SetZoneDeferredStatus(IOStatus status) {
  std::lock_guard<std::mutex> lk(zone_deferred_status_mutex_);
  if (!zone_deferred_status_.ok()) {
    zone_deferred_status_ = status;
  }
}

void ZonedBlockDevice::GetZoneSnapshot(std::vector<ZoneSnapshot> &snapshot) {
  for (auto *zone : io_zones) {
    snapshot.emplace_back(*zone);
  }
}

bool ZonedBlockDevice::SetSSTFileforZBDNoLock(uint64_t fno,
                                              ZoneFile *zoneFile) {
  auto sst = sst_file_bitmap_[fno];
  if (sst != nullptr) {  // already set
    return false;
  }
  // printf("zonedBlockDevice::SetSSTFileforZBDNoLock->fno: %lu\n", fno);
  zoneFile->fno_ = fno;
  sst_file_bitmap_[fno] = zoneFile;
  return true;
}

ZoneFile *ZonedBlockDevice::GetSSTZoneFileInZBDNoLock(uint64_t fno) {
  auto ret = sst_file_bitmap_[fno];
  // printf("zonedBlockDevice::GetSSTFileforZBDNoLock->fno: %lu\n", fno);
  if (ret == nullptr) {
    return nullptr;
  }
  if (ret != nullptr && ret->IsDeleted()) {
    return nullptr;
  }
  return ret;
}

bool ZonedBlockDevice::GetMinMaxKey(uint64_t fno, Slice &smallest,
                                    Slice &largest) {
  ZoneFile *zone_file = GetSSTZoneFileInZBDNoLock(fno);
  if (zone_file == nullptr) {
    return false;
  }
  smallest = zone_file->smallest_;
  largest = zone_file->largest_;
  return true;
}

void ZonedBlockDevice::SameLevelFileList(int level,
                                         std::vector<uint64_t> &fno_list,
                                         std::set<uint64_t> &compacting_files) {
  assert(db_ptr_ != nullptr);
  fno_list.clear();
  if (db_ptr_ == nullptr) {
    // printf("ZonedBlockDevice::SameLevelFileList!!!!\n");
    return;
  }
  // printf("zbd::samelevelfilelist->level: %d", level);
  db_ptr_->SameLevelFileList(level, fno_list, compacting_files);
}

void ZonedBlockDevice::SameLevelFileList(int level,
                                         std::vector<uint64_t> &fno_list,
                                         bool exclude_being_compacted) {
  assert(db_ptr_ != nullptr);
  fno_list.clear();
  // printf("level %d",level);
  db_ptr_->SameLevelFileList(level, fno_list, exclude_being_compacted);
}

void ZonedBlockDevice::SameLevelFileList(int level,
  std::unordered_map<uint64_t, uint64_t> &file_map,
  bool exclude_being_compacted, bool overap) {
  assert(db_ptr_ != nullptr);
  
  file_map.clear();
  // printf("level %d",level);
  db_ptr_->SameLevelFileList(level, file_map, exclude_being_compacted, overap);
}

void ZonedBlockDevice::UpperLevelFileList(Slice &smallest, Slice &largest,
                                          int level,
                                          std::vector<uint64_t> &fno_list) {
  assert(db_ptr_ != nullptr);
  fno_list.clear();
  if (db_ptr_ == nullptr) {
    return;
  }

  db_ptr_->UpperLevelFileList(smallest, largest, level, fno_list);
}

IOStatus ZonedBlockDevice::AllocateCompactionAwaredZoneV2(
    Slice &smallest, Slice &largest, int level,
    Env::WriteLifeTimeHint file_lifetime, uint64_t predicted_size,
    Zone **zone_out, uint64_t min_capacity) {
  if (allocation_scheme_ == LIZA) {
    return IOStatus::OK();
  }

  (void)(file_lifetime);
  IOStatus s;
  Zone *allocated_zone = nullptr;
  std::vector<uint64_t> fno_list;
  uint64_t max_invalid_data = 0;
  std::vector<uint64_t> zone_score(io_zones.size(), 0);
  std::vector<std::pair<uint64_t, uint64_t>> sorted;
  std::vector<bool> is_input_in_zone(io_zones.size(), false);
  // (void)(input_fno);
  (void)(predicted_size);
  (void)(max_invalid_data);

  if (db_ptr_ == nullptr) {
    // printf("AllocateCompactionAwaredZoneV2 - db_ptr is nullptr!!");
    return IOStatus::OK();
  }

  if (level == 0) {
    goto l0;
  }

  if (IS_BIG_SSTABLE(predicted_size)) {
    double upper_level_score = PredictCompactionScore(level - 1);
    double this_level_score = PredictCompactionScore(level);
    double l0_score = PredictCompactionScore(0);

    if (level == 1) {
      // append to downward
      // printf("CAZA - Level 1!!\n");
      fno_list.clear();
      zone_score.clear();
      zone_score.assign(io_zones.size(), 0);
      DownwardAdjacentFileList(smallest, largest, level, fno_list);
      if (CalculateZoneScore(fno_list, zone_score)) {
        // printf("CAZA - Level 1 - calculated!!\n");
        sorted = SortedByZoneScore(zone_score);
        AllocateZoneBySortedScore(sorted, &allocated_zone, min_capacity);
      }
    } else if (level == 2 && (l0_score > this_level_score) &&
               (l0_score > upper_level_score)) {
      // printf("CAZA - Level 2!!\n");
      fno_list.clear();
      zone_score.clear();
      zone_score.assign(io_zones.size(), 0);
      DownwardAdjacentFileList(smallest, largest, level, fno_list);

      if (CalculateZoneScore(fno_list, zone_score)) {
        sorted = SortedByZoneScore(zone_score);
        AllocateZoneBySortedScore(sorted, &allocated_zone, min_capacity);
      }
    }
    // if upper level occur first == false
    else if (this_level_score > upper_level_score) {
      // printf("this_level_score > upper_level_score!!\n");
      fno_list.clear();
      zone_score.clear();
      zone_score.assign(io_zones.size(), 0);
      DownwardAdjacentFileList(smallest, largest, level, fno_list);

      if (CalculateZoneScore(fno_list, zone_score)) {
        sorted = SortedByZoneScore(zone_score);
        AllocateZoneBySortedScore(sorted, &allocated_zone, min_capacity);
      }

    } else {  // upper_level_score>this_level_score
      // printf("upper_level_score>this_level_score!!\n");
      uint64_t upper_level_sst_fno =
          MostLargeUpperAdjacentFile(smallest, largest, level);

      ZoneFile *zfile = GetSSTZoneFileInZBDNoLock(upper_level_sst_fno);

      if (zfile && (IS_BIG_SSTABLE(zfile->GetFileSize()) || level == 1)) {
        // append to upper,this zfile
        GetNearestZoneFromZoneFile(zfile, is_input_in_zone, &allocated_zone,
                                   min_capacity);
      } else {
        // append to downward
        fno_list.clear();
        zone_score.clear();
        zone_score.assign(io_zones.size(), 0);
        DownwardAdjacentFileList(smallest, largest, level, fno_list);
        if (CalculateZoneScore(fno_list, zone_score)) {
          sorted = SortedByZoneScore(zone_score);
          AllocateZoneBySortedScore(sorted, &allocated_zone, min_capacity);
        }
      }
    }

  } else {
    uint64_t upper_level_sst_fno =
        MostLargeUpperAdjacentFile(smallest, largest, level);
    // printf("MostLargeUpperAdjacentFile!!\n");
    ZoneFile *zfile = GetSSTZoneFileInZBDNoLock(upper_level_sst_fno);
    if (level == 1) {
      fno_list.clear();
      // zone_score.assign(0,zone_score.size());
      zone_score.clear();
      zone_score.assign(io_zones.size(), 0);
      SameLevelFileList(0, fno_list);
      SameLevelFileList(1, fno_list);
      s = AllocateMostL0FilesZone(zone_score, fno_list, is_input_in_zone,
                                  &allocated_zone, min_capacity);
    } else if (zfile && IS_BIG_SSTABLE(zfile->predicted_size_)) {
      // append to upper zfile
      GetNearestZoneFromZoneFile(zfile, is_input_in_zone, &allocated_zone,
                                 min_capacity);
    } else {
      // append to same levels
      fno_list.clear();
      zone_score.clear();
      zone_score.assign(io_zones.size(), 0);
      SameLevelFileList(level, fno_list);
      s = AllocateSameLevelFilesZone(smallest, largest, fno_list,
                                     is_input_in_zone, &allocated_zone,
                                     min_capacity);
    }
  }
  if (allocated_zone != nullptr) {
    *zone_out = allocated_zone;
    return IOStatus::OK();
  }

/////////////////////////////
l0:
  // if level 0, most level 0 zone
  if (level == 0 || level == 100) {
    // printf("CAZA - Level 0!!\n");
    fno_list.clear();
    zone_score.clear();
    zone_score.assign(io_zones.size(), 0);
    SameLevelFileList(0, fno_list);
    SameLevelFileList(1, fno_list);
    s = AllocateMostL0FilesZone(zone_score, fno_list, is_input_in_zone,
                                &allocated_zone, min_capacity);
  }

  if (!s.ok()) {
    // printf("AllocateMostL0FilesZone - !s.ok()\n");
    return s;
  }
  if (allocated_zone != nullptr) {
    // printf("CAZA 3\n");
    *zone_out = allocated_zone;
    return s;
  }
  return s;
}

IOStatus ZonedBlockDevice::AllocateCompactionAwaredZone(
    Slice &smallest, Slice &largest, int level,
    Env::WriteLifeTimeHint file_lifetime, uint64_t predicted_size,
    Zone **zone_out, uint64_t min_capacity) {
  ////////CAZA/////////
  if (allocation_scheme_ == LIZA) {
    return IOStatus::OK();
  }
  if (db_ptr_ == nullptr) {
    // printf("AllocateCompactionAwaredZone - db_ptr is nullptr!!");
    return IOStatus::OK();
  }

  // printf("AllocateCompactionAwaredZone!!\n");
  (void)(file_lifetime);
  IOStatus s;
  uint64_t cur_score;
  uint64_t cur_invalid_data;
  bool no_near_level_files = true;
  Zone *allocated_zone = nullptr;
  Zone *target_zone;
  std::vector<uint64_t> fno_list;
  uint64_t max_score = 0;
  uint64_t max_invalid_data = 0;
  std::vector<bool> is_input_in_zone(io_zones.size(), false);
  // (void)(input_fno);
  (void)(predicted_size);
  (void)(cur_invalid_data);
  (void)(max_invalid_data);

  // zone valid overlapping capacity
  // 1. find UPPER/LOWER OVERLAPP RANGE zone
  std::vector<uint64_t> zone_score(io_zones.size(), 0);
  std::vector<std::pair<uint64_t, uint64_t>> sorted;
  if (level == 0) {
    goto l0;
  }
  {
    fno_list.clear();
    zone_score.clear();
    zone_score.assign(io_zones.size(), 0);
    AdjacentFileList(smallest, largest, level, fno_list);

    for (auto fno : fno_list) {
      ZoneFile *zFile = GetSSTZoneFileInZBDNoLock(fno);
      if (zFile == nullptr) {
        continue;
      }
      // if (zFile->selected_as_input_) {
      //   continue;
      // }
      auto extents = zFile->GetExtents();
      for (auto extent : extents) {
        if (!extent->zone_->IsFull()) {
          // zone_->index do not skip meta,spare zone
          // 존 스코어는 각 존이 얼마나 많은 데이터를 보유하고 있는지에 대한
          // 지표로, 이를 통해 특정 존이 얼마나 가득 찼는지 알 수 있다.
          zone_score[extent->zone_->zidx_ - ZENFS_META_ZONES -
                     ZENFS_SPARE_ZONES] += extent->length_;
          no_near_level_files = false;  // 상위 또는 하위 레벨에 겹치는 Key
                                        // 범위를 가진 파일이 존재
        }
      }
    }
  }

  sorted = SortedByZoneScore(zone_score);

  if (!no_near_level_files) {
    for (auto zidx : sorted) {
      cur_score = zidx.first;
      target_zone = io_zones[zidx.second];

      if (cur_score == 0 || target_zone->IsFull()) {
        continue;
      }

      if (cur_score < max_score) {
        continue;
      }

      if (!target_zone->Acquire()) {
        continue;
      }
      if (target_zone->IsEmpty()) {
        target_zone->Release();
        continue;
      }

      if (target_zone->capacity_ <= min_capacity) {
        target_zone->Release();
        continue;
      }
      allocated_zone = target_zone;
      break;
    }
  }

  if (allocated_zone != nullptr) {
    // printf("CAZA l1 <= \n");
    *zone_out = allocated_zone;
    return IOStatus::OK();
  }

l0:
  if (level == 0 || level == 1 || level == 100) {
    fno_list.clear();
    // zone_score.assign(0, zone_score.size());
    zone_score.clear();
    zone_score.assign(io_zones.size(), 0);
    SameLevelFileList(0, fno_list);
    SameLevelFileList(1, fno_list);
    s = AllocateMostL0FilesZone(zone_score, fno_list, is_input_in_zone,
                                &allocated_zone, min_capacity);
    if (allocated_zone != nullptr) {
      // printf("CAZA 2.1\n");
    }
  } else {  // if other level, same level but near key-sstfile zone
    fno_list.clear();
    // zone_score.assign(0,zone_score.size());
    zone_score.clear();
    zone_score.assign(io_zones.size() - ZENFS_META_ZONES - ZENFS_SPARE_ZONES,
                      0);
    SameLevelFileList(level, fno_list);
    s = AllocateSameLevelFilesZone(smallest, largest, fno_list,
                                   is_input_in_zone, &allocated_zone,
                                   min_capacity);
    if (allocated_zone != nullptr) {
      // printf("CAZA 2.2\n");
    }
  }

  if (!s.ok()) {
    return s;
  }
  if (allocated_zone != nullptr) {
    // printf("CAZA 3\n");
    *zone_out = allocated_zone;
    return s;
  }

  if (allocated_zone != nullptr) {
    *zone_out = allocated_zone;
    allocated_zone->lifetime_ = file_lifetime;
  }

  return s;
}

IOStatus ZonedBlockDevice::AllocateMostL0FilesZone(
    std::vector<uint64_t> &zone_score, std::vector<uint64_t> &fno_list,
    std::vector<bool> &is_input_in_zone, Zone **zone_out,
    uint64_t min_capacity) {
  Zone *allocated_zone = nullptr;
  Zone *target_zone = nullptr;
  IOStatus s;
  uint64_t max_score = 0;
  uint64_t cur_score;
  bool no_same_level_files = true;
  (void)(is_input_in_zone);

  // printf("AllocateMostL0FilesZone start!!\n");

  {
    // for (const auto &fno : fno_list) {
    //   std::cout << "fno: " << fno << std::endl;
    // }

    // std::lock_guard<std::mutex> lg(sst_file_map_lock_);
    for (auto fno : fno_list) {
      ZoneFile *zFile = GetSSTZoneFileInZBDNoLock(fno);
      if (zFile == nullptr) {
        continue;
      }
      // if (zFile->selected_as_input_) {
      //   continue;
      // }
      auto extents = zFile->GetExtents();
      for (auto e : extents) {
        if (!e->zone_->IsFull()) {
          // std::cout << "zidx : "
          //           << e->zone_->zidx_ - ZENFS_META_ZONES - ZENFS_SPARE_ZONES
          //           << std::endl;
          zone_score[e->zone_->zidx_ - ZENFS_META_ZONES - ZENFS_SPARE_ZONES] +=
              e->length_;
          no_same_level_files = false;
        }
      }
    }
  }
  if (no_same_level_files) {
    return IOStatus::OK();
  }

  //////////////////////////////
  auto sorted = SortedByZoneScore(zone_score);
  // printf("AllocateMostL0FilesZone sorted!!\n");

  for (auto zidx : sorted) {
    cur_score = zidx.first;
    target_zone = io_zones[zidx.second];

    if (cur_score == 0) {
      continue;
    }
    if (cur_score < max_score) {
      continue;
    }
    if (!target_zone->Acquire()) {
      continue;
    }
    if (target_zone->capacity_ <= min_capacity || target_zone->IsFull() ||
        target_zone->IsEmpty()) {
      target_zone->Release();
      continue;
    }
    allocated_zone = target_zone;
    break;
  }
  // printf("AllocateMostL0FilesZone\n");
  *zone_out = allocated_zone;
  return IOStatus::OK();
}

IOStatus ZonedBlockDevice::AllocateSameLevelFilesZone(
    Slice &smallest, Slice &largest, const std::vector<uint64_t> &fno_list,
    std::vector<bool> &is_input_in_zone, Zone **zone_out,
    uint64_t min_capacity) {
  Zone *allocated_zone = nullptr;
  IOStatus s;
  const Comparator *icmp = db_ptr_->GetDefaultICMP();
  ZoneFile *zFile;
  size_t idx;
  size_t l_idx;
  size_t r_idx;
  size_t fno_list_sz = fno_list.size();
  if (fno_list.empty()) {
    return IOStatus::OK();
  }

  {
    // std::lock_guard<std::mutex> lg(sst_file_map_lock_);
    if (fno_list_sz == 1) {
      zFile = GetSSTZoneFileInZBDNoLock(fno_list[0]);
      if (zFile != nullptr) {
        if (!zFile->selected_as_input_) {
          s = GetNearestZoneFromZoneFile(zFile, is_input_in_zone,
                                         &allocated_zone, min_capacity);
          if (!s.ok()) {
            return s;
          }
          *zone_out = allocated_zone;
          return s;
        }
      }
    }
    // fno_list is increasing order : db/version_set.h line 580
    for (idx = 0; idx < fno_list_sz; idx++) {
      zFile = GetSSTZoneFileInZBDNoLock(fno_list[idx]);
      if (zFile == nullptr) {
        continue;
      }
      // if (zFile->selected_as_input_) {
      //   continue;
      // }
      int res = icmp->Compare(largest, zFile->smallest_);
      if (res <= 0) {
        res = icmp->Compare(smallest, zFile->largest_);
        assert(res <= 0);
        break;
      }
    }

    l_idx = idx - 1;
    r_idx = idx;

    // it is most smallest key file
    if (idx == 0) {
      for (auto it = fno_list.begin(); it != fno_list.end(); it++) {
        zFile = GetSSTZoneFileInZBDNoLock(*it);
        if (zFile == nullptr) {
          continue;
        }
        // if (zFile->selected_as_input_) {
        //   continue;
        // }
        s = GetNearestZoneFromZoneFile(zFile, is_input_in_zone, &allocated_zone,
                                       min_capacity);
        if (!s.ok()) {
          return s;
        }
        if (allocated_zone != nullptr) {
          break;
        }
      }
    }
    // it is most largest key file
    else if (idx == fno_list_sz) {
      for (auto it = fno_list.rbegin(); it != fno_list.rend(); it++) {
        zFile = GetSSTZoneFileInZBDNoLock(*it);
        if (zFile == nullptr) {
          continue;
        }
        // if (zFile->selected_as_input_) {
        //   continue;
        // }
        s = GetNearestZoneFromZoneFile(zFile, is_input_in_zone, &allocated_zone,
                                       min_capacity);
        if (!s.ok()) {
          return s;
        }
        if (allocated_zone != nullptr) {
          break;
        }
      }
    }
    // it is middle key file
    else {
      for (bool flip = true; ((l_idx < r_idx) || (r_idx < fno_list_sz));
           flip = !flip) {
        if (flip) {
          if (l_idx < r_idx) {
            zFile = GetSSTZoneFileInZBDNoLock(fno_list[l_idx]);
            if (zFile == nullptr) {
              l_idx--;
              continue;
            }
            // if (zFile->selected_as_input_) {
            //   l_idx--;
            //   continue;
            // }
            s = GetNearestZoneFromZoneFile(zFile, is_input_in_zone,
                                           &allocated_zone, min_capacity);
            if (!s.ok()) {
              return s;
            }
            if (allocated_zone != nullptr) {
              break;
            }
            l_idx--;
          }
        } else {
          if (r_idx < fno_list_sz) {
            zFile = GetSSTZoneFileInZBDNoLock(fno_list[r_idx]);
            if (zFile == nullptr) {
              r_idx++;
              continue;
            }
            // if (zFile->selected_as_input_) {
            //   r_idx--;
            //   continue;
            // }
            s = GetNearestZoneFromZoneFile(zFile, is_input_in_zone,
                                           &allocated_zone, min_capacity);
            if (!s.ok()) {
              return s;
            }
            if (allocated_zone != nullptr) {
              break;
            }
            r_idx++;
          }
        }
      }
    }
  }

  *zone_out = allocated_zone;
  return s;
}

double ZonedBlockDevice::GetAvgCompressibilityOflevel(int output_level) {
  if (db_ptr_ == nullptr) {
    return 0;
  }
  return db_ptr_->GetAvgCompressibilityOflevel(output_level);
}

void ZonedBlockDevice::AdjacentFileList(Slice &smallest, Slice &largest,
                                        int level,
                                        std::vector<uint64_t> &fno_list) {
  assert(db_ptr_ != nullptr);
  fno_list.clear();
  if (db_ptr_ == nullptr) {
    return;
  }
  db_ptr_->AdjacentFileList(smallest, largest, level, fno_list);
}

uint64_t ZonedBlockDevice::MostSmallDownwardAdjacentFile(Slice &s, Slice &l,
                                                         int level) {
  if (db_ptr_ == nullptr) {
    return 0;
  }
  return db_ptr_->MostSmallDownwardAdjacentFile(s, l, level);
}
uint64_t ZonedBlockDevice::MostLargeUpperAdjacentFile(Slice &s, Slice &l,
                                                      int level) {
  if (db_ptr_ == nullptr) {
    return 0;
  }
  return db_ptr_->MostLargeUpperAdjacentFile(s, l, level);
}
void ZonedBlockDevice::DownwardAdjacentFileList(
    Slice &s, Slice &l, int level, std::vector<uint64_t> &fno_list) {
  if (db_ptr_ == nullptr) {
    return;
  }
  db_ptr_->DownwardAdjacentFileList(s, l, level, fno_list);
}
void ZonedBlockDevice::TrivialMoveFiles(int level,
                                        std::set<uint64_t> &trivial_set) {
  if (db_ptr_ == nullptr) {
    return;
  }
  db_ptr_->TrivialMoveFiles(level, trivial_set);
}

// return most large one
// uint64_t ZonedBlockDevice::MostLargeUpperAdjacentFile(Slice& smallest ,Slice&
// largest, int level){
//   assert(db_ptr_!=nullptr);

//   return db_ptr_->MostLargeUpperAdjacentFile(smallest,largest,level);
// }

IOStatus ZonedBlockDevice::GetNearestZoneFromZoneFile(
    ZoneFile *zFile, std::vector<bool> &is_input_in_zone, Zone **zone_out,
    uint64_t min_capacity) {
  // IOStatus s;

  // Zone* allocated_zone=nullptr;
  std::vector<std::pair<uint64_t, uint64_t>> zone_score(io_zones.size(),
                                                        {0, 0});
  auto extents = zFile->GetExtents();
  (void)(is_input_in_zone);

  for (auto e : extents) {
    uint64_t zidx = e->zone_->zidx_ - ZENFS_META_ZONES - ZENFS_SPARE_ZONES;
    zone_score[zidx].second = zidx;
    zone_score[zidx].first += e->length_;
  }

  std::sort(zone_score.rbegin(), zone_score.rend());

  for (auto zscore : zone_score) {
    uint64_t score = zscore.first;
    uint64_t zidx = zscore.second;
    // printf("zscore : %lu zidx %lu\n",score>>20,zidx);

    if (score == 0) {
      break;
    }
    Zone *z = io_zones[zidx];

    // if(z->IsEmpty()){
    //   continue;
    // }
    if (!z->Acquire()) {
      continue;
    }
    if (z->capacity_ <= min_capacity || z->IsFull() || z->IsEmpty()) {
      z->Release();
      continue;
    }
    // printf("return %lu\n",zidx);
    *zone_out = io_zones[zidx];
    return IOStatus::OK();
  }

  return IOStatus::OK();
}

IOStatus ZonedBlockDevice::AllocateAllInvalidZone(Zone **zone_out) {
  IOStatus s;
  Zone *allocated_zone = nullptr;

  for (const auto z : io_zones) {
    if (!z->Acquire()) {
      continue;
    }
    if (z->IsEmpty()) {
      z->Release();
      continue;
    }
    if (z->IsUsed()) {
      z->Release();
      continue;
    }

    allocated_zone = z;
    break;
  }

  *zone_out = allocated_zone;
  return IOStatus::OK();
}

bool ZonedBlockDevice::CalculateZoneScore(std::vector<uint64_t> &fno_list,
                                          std::vector<uint64_t> &zone_score) {
  bool there_is_near_level_files = false;

  // for (const auto &fno : fno_list) {
  //   std::cout << "calculatezs - fno: " << fno << std::endl;
  // }

  for (auto fno : fno_list) {
    ZoneFile *zFile = GetSSTZoneFileInZBDNoLock(fno);
    if (zFile == nullptr) {
      continue;
    }
    // if (zFile->selected_as_input_) {
    //   continue;
    // }
    auto extents = zFile->GetExtents();
    for (auto extent : extents) {
      if (!extent->zone_->IsFull()) {
        // zone_->index do not skip meta,spare zone
        zone_score[extent->zone_->zidx_ - ZENFS_META_ZONES -
                   ZENFS_SPARE_ZONES] += extent->length_;
        there_is_near_level_files = true;
      }
    }
  }

  return there_is_near_level_files;
}

double ZonedBlockDevice::GetMaxSameZoneScore(
    std::vector<uint64_t> &compaction_inputs_fno) {
  uint64_t total_size = 0, initial_total_size;
  uint64_t sst_in_zone_square = 0;
  double score = 0.0;
  std::vector<uint64_t> sst_in_zone(io_zones.size(), 0);
  for (uint64_t fno : compaction_inputs_fno) {
    ZoneFile *zFile = GetSSTZoneFileInZBDNoLock(fno);
    auto extents = zFile->GetExtents();
    for (ZoneExtent *extent : extents) {
      uint64_t zidx =
          extent->zone_->zidx_ - ZENFS_META_ZONES - ZENFS_SPARE_ZONES;
      sst_in_zone[zidx] += extent->length_;
      total_size += extent->length_;
    }
  }
  initial_total_size = total_size;
  /*
  2048
  1024 -> 0.5
  512
  512

  (512*512 + 512*512) / (1024*1024) * 0.5


  1024 -> 0.5
  1024 -> 0.5

  1024 -> 0.5
  768
  256
  (768^2 + 256 ^2) ÷ (1024^2) × 0.5= 0.3125
  */
  for (size_t i = 0; i < io_zones.size(); i++) {
    if (sst_in_zone[i] > io_zones[i]->max_capacity_ - (1 << 24)) {
      score += ((double)sst_in_zone[i] / (double)initial_total_size);
      total_size -= sst_in_zone[i];
      continue;
    }

    sst_in_zone_square += (sst_in_zone[i] * sst_in_zone[i]);
  }

  if (total_size > 0 && initial_total_size > 0) {
    score +=
        (double(sst_in_zone_square / total_size) / (double)initial_total_size);
  }
  // printf("\t\t\tscore : %lf     \n",score);
  return score;
}

double ZonedBlockDevice::GetMaxInvalidateCompactionScore(
    std::vector<uint64_t> &file_candidates, uint64_t *candidate_size,
    bool stats) {
  // printf("io zones n %lu\n", io_zones.size());
  std::vector<bool> is_sst_in_zone(io_zones.size(), false);
  std::vector<uint64_t> sst_in_zone(io_zones.size(), 0);
  std::vector<uint64_t> zone_score(io_zones.size(), 0);
  // printf("io zones n %lu 22222222\n", io_zones.size());
  uint64_t total_size = 0;
  uint64_t zidx;
  uint64_t zone_size = io_zones[0]->max_capacity_;
  uint64_t zone_score_sum = 0;
  uint64_t sst_in_zone_n = 0;
  uint64_t zone_score_max = 0;

  for (uint64_t fno : file_candidates) {
    ZoneFile *zFile = GetSSTZoneFileInZBDNoLock(fno);
    if (zFile == nullptr) {
      continue;
    }
    auto extents = zFile->GetExtents();
    for (ZoneExtent *extent : extents) {
      zidx = extent->zone_->zidx_ - ZENFS_META_ZONES - ZENFS_SPARE_ZONES;
      is_sst_in_zone[zidx] = true;
      sst_in_zone[zidx] += extent->length_;
    }
  }

  for (size_t i = 0; i < io_zones.size(); i++) {
    if (is_sst_in_zone[i] == false) {
      continue;
    }
    total_size += sst_in_zone[i];
    uint64_t relative_wp = io_zones[i]->wp_ - io_zones[i]->start_;
    relative_wp = relative_wp > io_zones[i]->max_capacity_
                      ? io_zones[i]->max_capacity_
                      : relative_wp;
    uint64_t after_valid_capacity =
        io_zones[i]->used_capacity_ - sst_in_zone[i];
    uint64_t after_invalid_capacity = relative_wp - after_valid_capacity;

    if (after_invalid_capacity > zone_size) {  // error point
      printf("after_invalid_capacity %lu > zone_size %lu??\n",
             after_invalid_capacity, zone_size);
      zone_score[i] = 100;
    } else {
      if (relative_wp != 0) {
        zone_score[i] = (after_invalid_capacity) * 100 / (relative_wp);
      } else {
        zone_score[i] = 0;
      }
    }

    if (zone_score[i] > zone_score_max) {
      zone_score_max = zone_score[i];
    }

    if (stats) {
      zone_score_sum += zone_score[i];
    } else {
      zone_score_sum += zone_score[i] * zone_score[i];
    }
    sst_in_zone_n++;
  }

  *candidate_size = total_size;

  if (sst_in_zone_n == 0) {
    return 0;
  }

  // if free space high, capacity first
  // if free space low, invalid score first

  (void)(zone_score_sum);
  (void)(zone_score_max);
  return (double)zone_score_sum / (double)sst_in_zone_n;
}

void ZonedBlockDevice::AllocateZoneBySortedScore(
    std::vector<std::pair<uint64_t, uint64_t>> &sorted, Zone **allocated_zone,
    uint64_t min_capacity) {
  uint64_t cur_score;
  Zone *target_zone;
  for (auto zidx : sorted) {
    cur_score = zidx.first;
    target_zone = io_zones[zidx.second];

    if (cur_score == 0) {
      break;
    }
    if (target_zone->IsFull()) {
      continue;
    }

    if (!target_zone->Acquire()) {
      continue;
    }

    if (target_zone->capacity_ <= min_capacity || target_zone->IsEmpty()) {
      target_zone->Release();
      continue;
    }
    (*allocated_zone) = target_zone;
    break;
  }
}

std::array<uint64_t, 10> ZonedBlockDevice::GetCurrentLSMTree() {
  std::array<uint64_t, 10> lsm_tree_snapshot;

  for (int i = 0; i < 10; ++i) {
    lsm_tree_snapshot[i] = lsm_tree_[i].load();
  }

  return lsm_tree_snapshot;
}

}  // namespace ROCKSDB_NAMESPACE

#endif  // !defined(ROCKSDB_LITE) && !defined(OS_WIN)
