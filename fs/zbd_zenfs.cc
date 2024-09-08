// Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved.
// Copyright (c) 2019-present, Western Digital Corporation
//  This source code is licensed under both the GPLv2 (found in the
//  COPYING file in the root directory) and Apache 2.0 License
//  (found in the LICENSE.Apache file in the root directory).
/* Zoned Block Device(ZBD)와 관련된 기능을 구현하고 정의하는 파일입니다. ZBD의
 * 초기화, 읽기, 쓰기, 리셋 등의 작업을 포함하며, ZBD와의 인터페이스를
 * 제공합니다.*/
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
#define ZENFS_META_ZONES (3)

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
      wp_(zbd_be->ZoneWp(zones, idx)) {  // 존의 현재 쓰기 포인터 초기화
  lifetime_ = Env::WLTH_NOT_SET;         // 존의 수명 초기화
  used_capacity_ = 0;                    // 사용된 용량 초기화
  capacity_ = 0;                         // 현재 용량 초기화
  if (zbd_be->ZoneIsWritable(zones, idx))  // 존이 쓰기 가능한 상태인지 확인
    capacity_ =
        max_capacity_ - (wp_ - start_);  // 쓰기 가능한 경우 현재 용량 설정
}

bool Zone::IsUsed() { return (used_capacity_ > 0); }
uint64_t Zone::GetCapacityLeft() { return capacity_; }
bool Zone::IsFull() { return (capacity_ == 0); }
bool Zone::IsEmpty() { return (wp_ == start_); }
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

/* 존이 재설정되면 용량, 쓰기 포인터, 수명 등이 초기화 */
IOStatus Zone::Reset() {
  bool offline;
  uint64_t max_capacity;

  assert(!IsUsed());
  assert(IsBusy());

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
/* 존을 마무리하고 용량을 0으로 설정하며, 쓰기 포인터를 존의 끝으로
이동시킵니다. 존이 더 이상 쓰기 작업을 받지 않도록 설정합니다. 일반적으로 데이터
쓰기가 완료된 후 존을 마무리할 때 사용됩니다.*/
IOStatus Zone::Finish() {
  assert(IsBusy());

  IOStatus ios = zbd_be_->Finish(start_);
  if (ios != IOStatus::OK()) return ios;

  capacity_ = 0;
  wp_ = start_ + zbd_->GetZoneSize();
  //////
  // std::cout << "######zone Finish" << std::endl;
  ////
  return IOStatus::OK();
}

/* 존을 닫습니다. 존이 비어 있지 않거나 가득 차 있지 않을 때만 닫기 작업을
수행합니다. 일반적으로 존을 더 이상 사용하지 않을 때 호출됩니다. 이는 존이 가득
찼거나 비어 있지 않은 경우에만 수행됩니다.*/
IOStatus Zone::Close() {
  assert(IsBusy());

  if (!(IsEmpty() || IsFull())) {
    IOStatus ios = zbd_be_->Close(start_);
    if (ios != IOStatus::OK()) return ios;
  }
  //////
  // std::cout << "######Starting close" << std::endl;
  ////

  return IOStatus::OK();
}

/* 주어진 데이터를 존에 추가하는 작업을 수행합니다. 이 함수는 데이터를 쓰기 전에
 * 존의 용량이 충분한지 확인하고, 데이터 크기가 블록 크기의 배수인지 검증합니다.
 * 쓰기 작업이 완료되면 쓰기 포인터와 용량을 업데이트하고, 쓰여진 바이트 수를
 * 기록합니다. 또한, 메트릭 시스템에 쓰기 지연 시간과 처리량을 보고합니다. 모든
 * 데이터가 성공적으로 쓰여지면 IOStatus::OK()를 반환합니다*/
IOStatus Zone::Append(char *data, uint32_t size) {
  ZenFSMetricsLatencyGuard guard(zbd_->GetMetrics(), ZENFS_ZONE_WRITE_LATENCY,
                                 Env::Default());
  zbd_->GetMetrics()->ReportThroughput(ZENFS_ZONE_WRITE_THROUGHPUT, size);
  char *ptr = data;
  uint32_t left = size;
  int ret;

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

/* 존의 사용 플래그(busy flag)를 해제하려고 시도하며, 해제가 실패하면 오류
 * 상태를 반환합니다. 성공적으로 해제되면 IOStatus::OK()를 반환*/
inline IOStatus Zone::CheckRelease() {
  if (!Release()) {
    assert(false);
    return IOStatus::Corruption("Failed to unset busy flag of zone " +
                                std::to_string(GetZoneNr()));
  }

  return IOStatus::OK();
}

/* 주어진 오프셋이 포함된 존을 io_zones 리스트에서 찾아 반환합니다*/
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
}

/* ZonedBlockDevice 객체를 초기화하고 열기 위한 작업을 수행합니다. 함수는 다음
작업을 수행합니다:

쓰기 열기가 독점적인지 확인합니다.
디바이스를 엽니다.
디바이스의 존 수가 최소 요구 사항을 충족하는지 확인합니다.
메타데이터 존과 I/O 존을 설정합니다.
시작 시간을 설정합니다.

이 작업을 통해 디바이스가 준비되고, 사용할 수 있는 존이 올바르게 설정됩니다.*/
IOStatus ZonedBlockDevice::Open(bool readonly, bool exclusive) {
  std::unique_ptr<ZoneList> zone_rep;
  unsigned int max_nr_active_zones;
  unsigned int max_nr_open_zones;
  Status s;
  uint64_t i = 0;
  uint64_t m = 0;
  // Reserve one zone for metadata and another one for extent migration
  int reserved_zones = 2;

  if (!readonly && !exclusive)
    return IOStatus::InvalidArgument("Write opens must be exclusive");

  IOStatus ios = zbd_be_->Open(readonly, exclusive, &max_nr_active_zones,
                               &max_nr_open_zones);
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

  Info(logger_, "Zone block device nr zones: %u max active: %u max open: %u \n",
       zbd_be_->GetNrZones(), max_nr_active_zones, max_nr_open_zones);
  //
  printf("Zone block device nr zones: %u max active: %u max open: %u \n",
         zbd_be_->GetNrZones(), max_nr_active_zones, max_nr_open_zones);
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
  uint64_t device_io_capacity = 85899345920;  // 80GB
  // uint64_t device_io_capacity = 10737418240;  // 10GB
  for (; i < zone_rep->ZoneCount() &&
         (io_zones.size() * zbd_be_->GetZoneSize()) < (device_io_capacity);
       i++) {
    // for (; i < zone_rep->ZoneCount(); i++) {
    /* Only use sequential write required zones */
    if (zbd_be_->ZoneIsSwr(zone_rep, i)) {
      if (!zbd_be_->ZoneIsOffline(zone_rep, i)) {
        Zone *newZone = new Zone(this, zbd_be_.get(), zone_rep, i);
        if (!newZone->Acquire()) {
          assert(false);
          return IOStatus::Corruption("Failed to set busy flag of zone " +
                                      std::to_string(newZone->GetZoneNr()));
        }
        io_zones.push_back(newZone);
        printf("io zone at %ld\n", i);
        if (zbd_be_->ZoneIsActive(zone_rep, i)) {
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
  // uint64_t device_free_space = (ZENFS_IO_ZONES) * (ZONE_SIZE);
  // printf("device free space : %ld\n", device_free_space);
  // device_free_space_.store(device_free_space);

  printf("io_zones.size() : %ld\n", io_zones.size());
  printf("zbd_be_->GetZoneSize(): %ld\n", zbd_be_->GetZoneSize());
  uint64_t device_free_space = io_zones.size() * zbd_be_->GetZoneSize();
  printf("device free space : %ld\n", BYTES_TO_MB(device_free_space));
  // printf("zone sz %ld\n", zone_sz_);
  device_free_space_.store(device_free_space);

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
  // 존 쓰레기 통계 벡터를 로그로 기록합니다.
  //
  // 벡터의 각 값은 특정 쓰레기 비율을 가진 존의 수를 나타냅니다.
  // 각 인덱스의 쓰레기 비율: [0%, <10%, <20%, ... <100%, 100%]
  // 예를 들어, `[100, 1, 2, 3....]`는 100개의 존이 비어 있고,
  // 1개의 존은 쓰레기 비율이 10% 미만, 2개의 존은 쓰레기 비율이 10% ~ 20%
  // 사이임을 의미합니다.
  //
  // 데이터를 읽기만 하므로 io_zones를 잠글 필요가 없으며, 결과의 정확성을 높일
  // 필요도 없습니다.
  int zone_gc_stat[12] = {0};  // 각 쓰레기 비율 구간별 존의 수를 저장하는 배열
  for (auto z : io_zones) {  // 모든 I/O 존에 대해 반복
    if (!z->Acquire()) {     // 존이 사용 가능하지 않으면 건너뜀
      continue;
    }

    if (z->IsEmpty()) {  // 존이 비어 있는 경우
      zone_gc_stat[0]++;
      z->Release();  // 존을 해제
      continue;
    }

    double garbage_rate = 0;  // 쓰레기 비율 초기화
    if (z->IsFull()) {        // 존이 가득 찬 경우
      garbage_rate =
          double(z->max_capacity_ - z->used_capacity_) / z->max_capacity_;
    } else {  // 존이 가득 차지 않은 경우
      garbage_rate =
          double(z->wp_ - z->start_ - z->used_capacity_) / z->max_capacity_;
    }
    assert(garbage_rate >= 0);  // 쓰레기 비율이 0 이상인지 확인
    int idx = int((garbage_rate + 0.1) * 10);  // 쓰레기 비율을 인덱스로 변환
    zone_gc_stat[idx]++;  // 해당 쓰레기 비율 구간의 존 수 증가

    z->Release();  // 존을 해제
  }

  std::stringstream ss;  // 로그 메시지를 저장할 문자열 스트림
  ss << "Zone Garbage Stats: [";
  for (int i = 0; i < 12; i++) {  // 각 쓰레기 비율 구간의 존 수를 스트림에 추가
    ss << zone_gc_stat[i] << " ";
  }
  ss << "]";
  Info(logger_, "%s", ss.str().data());  // 로그 메시지를 기록
}

ZonedBlockDevice::~ZonedBlockDevice() {
  size_t rc = reset_count_.load();
  uint64_t wwp = wasted_wp_.load() / (1 << 20);
  uint64_t R_wp;
  if (rc == 0) {
    R_wp = 100;
  } else {
    R_wp = (ZONE_SIZE * 100 - wwp * 100 / (rc)) / ZONE_SIZE;
  }
  printf("============================================================\n");
  printf("FAR STAT 1 :: WWP (MB) : %lu, R_wp : %lu\n",
         wasted_wp_.load() / (1 << 20), R_wp);
  // printf("ZC IO Blocking time : %d, Compaction Refused : %lu\n",
  // zc_io_block_,
  //        compaction_blocked_at_amount_.size());
  printf("ZC IO Blocking time : %d\n", zc_io_block_);

  printf("============================================================\n");
  uint64_t total_copied = 0;
  size_t rc_zc = 0;
  int io_blocking_sum = 0;
  // long long io_blocking_ms_sum = 0;
  for (size_t i = 0;
       i < zc_timelapse_.size() && i < zc_copied_timelapse_.size(); i++) {
    bool forced = zc_timelapse_[i].forced;
    size_t zc_z = zc_timelapse_[i].zc_z;
    int s = zc_timelapse_[i].s;
    int e = zc_timelapse_[i].e;
    // long long us = zc_timelapse_[i].us;
    io_blocking_sum += e - s + 1;
    // io_blocking_ms_sum += us;
    printf("[%lu] :: %d ~ %d, %ld (MB), Reclaimed Zone : %lu [%s]\n", i + 1, s,
           e, (zc_copied_timelapse_[i]) / (1 << 20), zc_z,
           forced ? "FORCED" : " ");
    total_copied += zc_copied_timelapse_[i];
    rc_zc += zc_z;
  }
  printf("Total ZC Copied (MB) :: %lu, Recaimed by ZC :: %lu \n",
         total_copied / (1 << 20), rc_zc);
  printf("FAR STAT  :: Reset Count (R+ZC) : %ld+%ld=%ld\n", rc - rc_zc, rc_zc,
         rc);
  // for (size_t i = 0; i < io_block_timelapse_.size(); i++) {
  //   int s = io_block_timelapse_[i].s;
  //   int e = io_block_timelapse_[i].e;
  //   pid_t tid = io_block_timelapse_[i].tid;
  //   printf("[%lu] :: (%d) %d ~ %d\n", i + 1, tid % 100, s, e);
  // }

  printf("TOTAL I/O BLOKCING TIME %d\n", io_blocking_sum);
  // printf("TOTAL I/O BLOCKING TIME(ms) %llu\n", io_blocking_ms_sum);

  for (const auto z : meta_zones) {
    delete z;
  }

  for (const auto z : io_zones) {
    delete z;
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

  if (zone_lifetime > file_lifetime) return zone_lifetime - file_lifetime;
  if (zone_lifetime == file_lifetime) return LIFETIME_DIFF_COULD_BE_WORSE;

  return LIFETIME_DIFF_NOT_GOOD;
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

void ZonedBlockDevice::CalculateResetThreshold() {
  uint64_t rt = 0;
  uint64_t max_capacity = io_zones[0]->max_capacity_;
  uint64_t free_percent = cur_free_percent_;
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
}
/* io_zones 벡터의 각 존을 순회하며, 사용되지 않는 IO 존을 재설정합니다.
재설정이 완료되면, 필요에 따라 토큰을 반환합니다.*/
IOStatus ZonedBlockDevice::ResetUnusedIOZones() {
  clock_t reset_latency{0};
  for (const auto z : io_zones) {
    // Acquire 함수는 존을 사용할 수 있는지 확인하고 가능하면 락을 거는 역할
    if (z->Acquire()) {
      // 현재 존이 비어있지 않고 사용중이지 않은 경우를 확인
      if (!z->IsEmpty() && !z->IsUsed()) {
        bool full = z->IsFull();  // 현재 존이 가득 찬 상태인지 확인하여 full
                                  // 변수에 저장
        uint64_t cp = z->GetCapacityLeft();
        //(cp > reset_threshold_) {
        //   IOStatus release_status = z->CheckRelease();
        //   if(!release_status.ok()){
        //     return release_status;
        //   }
        //   continue;
        // }
        clock_t start = clock();
        IOStatus reset_status = z->Reset();  // z->Reset(): 현재 존을 재설정
        reset_count_.fetch_add(1);
        wasted_wp_.fetch_add(cp);

        clock_t end = clock();
        reset_latency += (end - start);
        runtime_reset_reset_latency_.fetch_add(end - start);
        IOStatus release_status =
            z->CheckRelease();  // z->CheckRelease(): 재설정 후 존을 해제
        if (!reset_status.ok())
          return reset_status;  // 두 함수의 상태를 확인하여 오류가 발생하면
                                // 해당 상태를 반환
        if (!release_status.ok()) return release_status;  //
        if (!full)
          PutActiveIOZoneToken();  // 현재 존이 가득 차지 않았다면,
                                   // PutActiveIOZoneToken 함수를 호출하여
                                   // 토큰을 반환
      } else {  // 현재 존이 비어있거나 사용 중인 경우, CheckRelease 함수로 해제
                // 상태를 확인하고, 오류가 발생하면 해당 상태를 반환
        IOStatus release_status = z->CheckRelease();
        if (!release_status.ok()) return release_status;
      }
    }
  }
  return IOStatus::OK();  // 모든 존에 대한 작업이 성공적으로 완료되면
                          // IOStatus::OK()를 반환
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

IOStatus ZonedBlockDevice::FinishCheapestIOZone() {
  IOStatus s;
  Zone *finish_victim = nullptr;

  for (const auto z : io_zones) {
    if (z->Acquire()) {
      if (z->IsEmpty() || z->IsFull()) {
        s = z->CheckRelease();
        if (!s.ok()) return s;
        continue;
      }
      if (finish_victim == nullptr) {
        finish_victim = z;
        continue;
      }
      if (finish_victim->capacity_ > z->capacity_) {
        s = finish_victim->CheckRelease();
        if (!s.ok()) return s;
        finish_victim = z;
      } else {
        s = z->CheckRelease();
        if (!s.ok()) return s;
      }
    }
  }

  // If all non-busy zones are empty or full, we should return success.
  if (finish_victim == nullptr) {
    Info(logger_, "All non-busy zones are empty or full, skip.");
    return IOStatus::OK();
  }

  s = finish_victim->Finish();
  IOStatus release_status = finish_victim->CheckRelease();

  if (s.ok()) {
    PutActiveIOZoneToken();
  }

  if (!release_status.ok()) {
    return release_status;
  }

  return s;
}

IOStatus ZonedBlockDevice::GetAnyLargestRemainingZone(Zone **zone_out,
                                                      uint32_t min_capacity) {
  IOStatus s = IOStatus::OK();
  Zone *allocated_zone = nullptr;

  for (const auto z : io_zones) {
    if (!z->Acquire()) {
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
        if (diff <= best_diff) {
          if (allocated_zone != nullptr) {
            s = allocated_zone->CheckRelease();
            if (!s.ok()) {
              IOStatus s_ = z->CheckRelease();
              if (!s_.ok()) return s_;
              return s;
            }
          }
          allocated_zone = z;
          best_diff = diff;
        } else {
          s = z->CheckRelease();
          if (!s.ok()) return s;
        }
      } else {
        s = z->CheckRelease();
        if (!s.ok()) return s;
      }
    }
  }

  *best_diff_out = best_diff;
  *zone_out = allocated_zone;

  return IOStatus::OK();
}

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
      s = zone->CheckRelease();
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
  // auto s =
  //     GetBestOpenZoneMatch(file_lifetime, &best_diff, out_zone,
  //     min_capacity);
  // if (s.ok() && (*out_zone) != nullptr) {
  //   Info(logger_, "TakeMigrateZone: %lu", (*out_zone)->start_);
  // } else {
  //   migrating_ = false;
  // }

  // return s;
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
      printf("TakeMigrateZone: %lu", (*out_zone)->start_);
      break;
    }
    s = AllocateEmptyZone(out_zone);
    if (s.ok() && (*out_zone) != nullptr) {
      Info(logger_, "TakeMigrateZone: %lu", (*out_zone)->start_);
      printf("TakeMigrateZone - emptyzone : %lu", (*out_zone)->start_);
      break;
    }

    // If AllocateEmptyZone fails, try to allocate any zone that meets the
    // min_capacity requirement
    s = GetAnyLargestRemainingZone(out_zone, min_capacity);
    if (s.ok() && (*out_zone) != nullptr) {
      Info(logger_, "TakeMigrateZone: %lu", (*out_zone)->start_);
      printf("TakeMigrateZone - emptyzone failed : %lu", (*out_zone)->start_);
      break;
    }

    s = ResetUnusedIOZones();
    sleep(1);
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

IOStatus ZonedBlockDevice::AllocateIOZone(Env::WriteLifeTimeHint file_lifetime,
                                          IOType io_type, Zone **out_zone) {
  // IOStatus s;
  // s = ResetUnusedIOZones();
  Zone *allocated_zone = nullptr;
  unsigned int best_diff = LIFETIME_DIFF_NOT_GOOD;
  int new_zone = 0;
  IOStatus s;

  // std::cout << "@@@ zbd::AllocateIOZone - life_time: " << file_lifetime
  //           << "// out_zone: " << out_zone << "\n";

  auto tag = ZENFS_WAL_IO_ALLOC_LATENCY;
  if (io_type != IOType::kWAL) {
    // L0 flushes have lifetime MEDIUM
    if (file_lifetime == Env::WLTH_MEDIUM) {
      tag = ZENFS_L0_IO_ALLOC_LATENCY;
    } else {
      tag = ZENFS_NON_WAL_IO_ALLOC_LATENCY;
    }
  }

  ZenFSMetricsLatencyGuard guard(metrics_, tag, Env::Default());
  metrics_->ReportQPS(ZENFS_IO_ALLOC_QPS, 1);

  // Check if a deferred IO error was set
  s = GetZoneDeferredStatus();
  if (!s.ok()) {
    return s;
  }

  if (io_type != IOType::kWAL) {
    s = ApplyFinishThreshold();
    if (!s.ok()) {
      return s;
    }
  }

  WaitForOpenIOZoneToken(io_type == IOType::kWAL);

  /* Try to fill an already open zone(with the best life time diff) */
  s = GetBestOpenZoneMatch(file_lifetime, &best_diff, &allocated_zone);
  if (!s.ok()) {
    PutOpenIOZoneToken();
    return s;
  }
  if (allocated_zone == nullptr) {
    // printf("GetBestOpenZone return nullptr\n");
  }

  // Holding allocated_zone if != nullptr

  if (best_diff >= LIFETIME_DIFF_COULD_BE_WORSE) {
    bool got_token = GetActiveIOZoneTokenIfAvailable();

    /* If we did not get a token, try to use the best match, even if the life
     * time diff not good but a better choice than to finish an existing zone
     * and open a new one
     */
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
        s = FinishCheapestIOZone();
        if (!s.ok()) {
          PutOpenIOZoneToken();
          return s;
        }
      }

      s = AllocateEmptyZone(&allocated_zone);
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
        assert(allocated_zone->IsBusy());
        allocated_zone->lifetime_ = file_lifetime;
        new_zone = true;
      } else {
        PutActiveIOZoneToken();
      }
    }
  }

  if (allocated_zone) {
    assert(allocated_zone->IsBusy());
    Debug(logger_,
          "Allocating zone(new=%d) start: 0x%lx wp: 0x%lx lt: %d file lt: %d\n",
          new_zone, allocated_zone->start_, allocated_zone->wp_,
          allocated_zone->lifetime_, file_lifetime);
  } else {
    PutOpenIOZoneToken();
  }

  if (io_type != IOType::kWAL) {
    LogZoneStats();
  }

  *out_zone = allocated_zone;

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

}  // namespace ROCKSDB_NAMESPACE

#endif  // !defined(ROCKSDB_LITE) && !defined(OS_WIN)
