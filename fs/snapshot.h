// Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved.
// Copyright (c) 2019-present, Western Digital Corporation
//  This source code is licensed under both the GPLv2 (found in the
//  COPYING file in the root directory) and Apache 2.0 License
//  (found in the LICENSE.Apache file in the root directory).

/* ZenFS에서의 스냅샷은 Zoned Block Device(ZBD)와 관련된 다양한 통계 정보와
 * 상태를 캡처하여 저장 */
/* 스냅샷 기능을 구현하는 파일로, 파일 시스템의 상태를 특정 시점에서 저장하고
 * 복구하는 기능을 담당*/
#pragma once

#include <string>
#include <vector>

#include "io_zenfs.h"
#include "zbd_zenfs.h"

namespace ROCKSDB_NAMESPACE {

// Indicate what stats info we want.
struct ZenFSSnapshotOptions {
  // Global zoned device stats info
  bool zbd_ = 0;  // 전역 Zoned Block Device 통계 정보
  // Per zone stats info
  bool zone_ = 0;  // 개별 존 통계 정보
  // Get all file->extents & extent->file mappings
  bool zone_file_ = 0;  // 파일->익스텐트 및 익스텐트->파일 매핑 정보
  bool trigger_report_ = 0;            // 보고서 트리거
  bool log_garbage_ = 0;               // 로그 가비지 정보
  bool as_lock_free_as_possible_ = 1;  // 가능한 락프리로 처리
};

/* ZBD의 전역적인 공간 정보(사용된 공간, 남은 공간, 회수 가능한 공간 등)를 캡처
 */
class ZBDSnapshot {
 public:
  uint64_t free_space;
  uint64_t used_space;
  uint64_t reclaimable_space;

 public:
  ZBDSnapshot() =
      default;  // 아무런 초기화도 수행하지 않고, 멤버 변수를 기본 값으로 설정
  ZBDSnapshot(const ZBDSnapshot&) =
      default;  // ZBDSnapshot 객체의 값을 그대로 복사
  ZBDSnapshot(ZonedBlockDevice& zbd)  // ZonedBlockDevice 객체를 매개변수로
                                      // 받아와 ZBDSnapshot 객체를 초기화
      : free_space(
            zbd.GetFreeSpace()),  // free_space 멤버 변수를 zbd.GetFreeSpace()의
                                  // 반환 값으로 초기화
        used_space(zbd.GetUsedSpace()),
        reclaimable_space(zbd.GetReclaimableSpace()) {}
};

/* 각 존의 시작 위치, 용량, 사용된 용량 등을 포함한 정보를 캡처 */
class ZoneSnapshot {
 public:
  uint64_t start;  // 존의 시작 위치
  uint64_t wp;     // 존의 쓰기 포인터

  uint64_t capacity;       // 존의 총 용량
  uint64_t used_capacity;  // 존에서 사용된 용량
  uint64_t max_capacity;   // 존의 최대 용량

 public:
  ZoneSnapshot(const Zone& zone)  // ZoneSnapshot의 각 멤버 변수는 Zone 객체의
                                  // 해당 멤버 변수 값으로 설정
      : start(zone.start_),
        wp(zone.wp_),
        capacity(zone.capacity_),
        used_capacity(zone.used_capacity_),
        max_capacity(zone.max_capacity_) {}
};

/* 각 존의 익스텐트(연속된 데이터 블록)의 시작 위치, 길이, 존의 시작 위치 등을
 * 캡처("파일"의 특정 부분(익스텐트)을 기록) */
class ZoneExtentSnapshot {
 public:
  uint64_t start;
  uint64_t length;
  uint64_t zone_start;
  std::string filename;

 public:
  ZoneExtentSnapshot(const ZoneExtent& extent, const std::string fname)
      : start(extent.start_),
        length(extent.length_),
        zone_start(extent.zone_->start_),
        filename(fname) {}
};

/* 각 파일의 ID, 이름, 익스텐트 목록을 캡처 */
class ZoneFileSnapshot {
 public:
  uint64_t file_id;
  std::string filename;
  std::vector<ZoneExtentSnapshot> extents;

 public:
  ZoneFileSnapshot(ZoneFile& file)
      : file_id(file.GetID()),
        filename(file.GetFilename()) {  // 파일id 초기화, 파일이름 초기화
    for (const auto* extent : file.GetExtents()) {  // 익스텐트 목록 초기화
      extents.emplace_back(*extent, filename);
    }
  }
};

/* ZenFS 파일 시스템의 스냅샷 */
class ZenFSSnapshot {
 public:
  ZenFSSnapshot() {}
  /* 이동 할당 연산자 - 다른 ZenFSSnapshot 객체로부터 자원을 이동하여 현재
   * 객체를 초기화 */
  ZenFSSnapshot& operator=(ZenFSSnapshot&& snapshot) {
    zbd_ = snapshot.zbd_;
    zones_ = std::move(snapshot.zones_);
    zone_files_ = std::move(snapshot.zone_files_);
    extents_ = std::move(snapshot.extents_);
    return *this;
  }

 public:
  ZBDSnapshot zbd_;  // 사용된 공간, 남은 공간, 회수 가능한 공간
  std::vector<ZoneSnapshot>
      zones_;  // 존의 시작 위치, 쓰기 포인터, 총 용량, 사용된 용량
  std::vector<ZoneFileSnapshot> zone_files_;  // 파일의 ID, 이름, 익스텐트 목록
  std::vector<ZoneExtentSnapshot>
      extents_;  // 익스텐트의 시작 위치, 길이, 존의 시작 위치, 파일 이름
};

}  // namespace ROCKSDB_NAMESPACE
