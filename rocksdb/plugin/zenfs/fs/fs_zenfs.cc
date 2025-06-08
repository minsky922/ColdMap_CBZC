// Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved.
// Copyright (c) 2019-present, Western Digital Corporation
//  This source code is licensed under both the GPLv2 (found in the
//  COPYING file in the root directory) and Apache 2.0 License
//  (found in the LICENSE.Apache file in the root directory).
#if !defined(ROCKSDB_LITE) && defined(OS_LINUX)
//
#include "fs_zenfs.h"
///
#include <chrono>
#include <cmath>
#include <ctime>
#include <iostream>
#include <numeric>

/// cur_ops_
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

IOStatus ZenMetaLog::Read(Slice* slice) {
  char* data = (char*)slice->data(); 
  size_t read = 0;                   
  size_t to_read = slice->size();     
  int ret;                            

  if (read_pos_ >= zone_->wp_) {
    slice->clear();
    return IOStatus::OK();
  }

  if ((read_pos_ + to_read) > (zone_->start_ + zone_->max_capacity_)) {
    return IOStatus::IOError("Read across zone");
  }

  while (read < to_read) {
    ret = zbd_->Read(data + read, read_pos_, to_read - read, false);

    if (ret == -1 && errno == EINTR) continue;  
    if (ret < 0)
      return IOStatus::IOError("Read failed"); 

    read += ret;       
    read_pos_ += ret; 
  }

  return IOStatus::OK();
}

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

ZenFS::ZenFS(ZonedBlockDevice* zbd, std::shared_ptr<FileSystem> aux_fs,
             std::shared_ptr<Logger> logger)
    : FileSystemWrapper(aux_fs), zbd_(zbd), logger_(logger) {
  Info(logger_, "ZenFS initializing");

  Info(logger_, "ZenFS parameters: block device: %s, aux filesystem: %s",
       zbd_->GetFilename().c_str(), target()->Name());

  next_file_id_ = 1;

  metadata_writer_.zenFS = this;

  memset(file_size_dist, 0, sizeof(file_size_dist));
}

ZenFS::~ZenFS() {
  Status s;                             
  Info(logger_, "ZenFS shutting down");  
  zbd_->LogZoneUsage();                  
  LogFiles();                           

  run_gc_worker_ = false;  
  run_bg_reset_worker_ = false;

  if (gc_worker_) {      
    gc_worker_->join();  
  }

  if (bg_reset_worker_) {
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

  meta_log_.reset(nullptr);  
  uint64_t non_free = zbd_->GetUsedSpace() + zbd_->GetReclaimableSpace();
  uint64_t free = zbd_->GetFreeSpace();
  uint64_t free_percent = (100 * free) / (free + non_free);

  printf("@@~Zenfs Last Free percent freepercent %ld \n", free_percent);
  ClearFiles();  
  delete zbd_;   
}

void ZenFS::BackgroundStatTimeLapse() {
  while (run_bg_reset_worker_) {
    free_percent_ = zbd_->CalculateFreePercent();
    zbd_->AddTimeLapse(mount_time_);
    /*
    stressed test
    1. cur ops
    2. 10 second throughput
    3. COPY
    4. BLOCKING
    5. RESET COUNT

    */
    sleep(1);
    int cur_time = mount_time_.fetch_add(1);
    // uint64_t total = 0;
    if (cur_time % 50 == 0|| run_bg_reset_worker_==false) {
      struct timespec timespec;
      clock_gettime(CLOCK_MONOTONIC, &timespec);
      uint64_t gc_bytes_written = zbd_->GetGCBytesWritten();
      uint64_t rc = zbd_->GetRC();
      uint64_t cur_ops = cur_ops_.load();
      uint64_t cur_get_ops = cur_get_ops_.load();
      uint64_t cur_scan_ops = cur_scan_ops_.load();
      uint64_t cur_io_blocking = zbd_->GetBlocking();

      // uint64_t total_time_ms =
      //     zbd_->total_deletion_after_copy_time_.load() / 1000;
      // uint64_t total_count = zbd_->total_deletion_after_copy_n_.load();
      // uint64_t total_deletion_after_copy_size =zbd_->total_deletion_after_copy_size_;
      // uint64_t avg_ms = total_count == 0 ? 0 : (total_time_ms / total_count);
      // uint64_t actual_cost_benefit_score =
      //     zbd_->actual_cost_benefit_score_.load();
      // uint64_t average_actual_cost_benefit_score;

      

      ////
 
      // uint64_t cost_benefit_score_sum_sequence_mb = zbd_->cost_benefit_score_sum_sequence_mb_.load();
      // uint64_t total_deletion_after_copy_seq = zbd_->total_deletion_after_copy_seq_.load()/100;
      // uint64_t copy_but_not_deleted_size = 0;
      // if(cur_time % 100 == 0|| run_bg_reset_worker_==false){
      //   //  get copied but not deleted
      //   int cur_fops_sequence = zbd_->file_operation_sequence_.load();

      //   std::lock_guard<std::mutex> lg_(files_mtx_);
      //   for (auto f : files_) {
      //     int tmp_n=0;
      //     int tmp_seq = 0;
      //     size_t tmp_size=0;
      //     long diff_ns_sum = 0;
      //     for (auto ext : f.second->GetExtents()) {
      //       if (ext->is_zc_copied_ && ext->ZC_COPIED_STATE==ZC_COPIED) {
      //         diff_ns_sum +=
      //             (timespec.tv_sec - ext->zc_copied_ts_.tv_sec) * 1000000000L +
      //             (timespec.tv_nsec - ext->zc_copied_ts_.tv_nsec);
              
      //         copy_but_not_deleted_size += ext->length_;
      //         tmp_size+=ext->length_;
      //         // actual_cost_benefit_score +=
      //         //     ((ext->length_ >> 20) * (diff_ns / 1000 / 1000));
      //             // is_copied=true;
      //             tmp_n++;
      //             tmp_seq+=(cur_fops_sequence-ext->zc_copied_sequence_);
      //         }

      //     }
      //     if(tmp_n){
      //       actual_cost_benefit_score +=
      //       ((tmp_size >> 20) * (diff_ns_sum / 1000 / 1000));
      //       // if(run_bg_reset_worker_==false){
      //       //   printf("SEQ NOT DELETED UNTIL END %d\n",(tmp_seq/tmp_n));
      //       // }
      //       total_deletion_after_copy_seq+=(tmp_seq/tmp_n);
      //       total_count++;
      //       cost_benefit_score_sum_sequence_mb += ((tmp_size>>20)*tmp_seq)/tmp_n;
      //     }
      //   }
      //   total_deletion_after_copy_size+=copy_but_not_deleted_size;
      // }

      
      
      // uint64_t A_CD_sequence =  total_count == 0 ? 0 : 
      //              ( total_deletion_after_copy_seq) /  total_count;
      // uint64_t WA_CD_sequence = (total_deletion_after_copy_size>>20) == 0 ? 0 :
      // cost_benefit_score_sum_sequence_mb / (total_deletion_after_copy_size >> 20);

      // average_actual_cost_benefit_score = zbd_->GetGCBytesWritten() == 0 ? 0 :
      //     actual_cost_benefit_score / (zbd_->GetGCBytesWritten() >> 20);
      // average_actual_cost_benefit_score = (total_deletion_after_copy_size>>20 == 0) ? 0 :
      // actual_cost_benefit_score / (total_deletion_after_copy_size >> 20);
      // if (cur_get_ops || cur_scan_ops) {
      //   printf("%d\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t",
      //          cur_time, zbd_->CalculateFreePercent(), cur_ops,cur_get_ops,cur_scan_ops,
      //          gc_bytes_written, cur_io_blocking, rc,
      //          zbd_->open_io_zones_.load(), zbd_->active_io_zones_.load(),
      //          total_time_ms, total_count, avg_ms, actual_cost_benefit_score,average_actual_cost_benefit_score,total_deletion_after_copy_size);

      // } else {
      //   printf("%d\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t",
      //          cur_time, zbd_->CalculateFreePercent(), cur_ops,
      //          gc_bytes_written, cur_io_blocking, rc,
      //          zbd_->open_io_zones_.load(), zbd_->active_io_zones_.load(),
      //          total_time_ms, total_count, avg_ms, actual_cost_benefit_score,average_actual_cost_benefit_score,total_deletion_after_copy_size);
      // }
      // printf("%d\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t",
      //          cur_time, zbd_->CalculateFreePercent(), cur_ops,cur_get_ops,cur_scan_ops,
      //          (gc_bytes_written>>20), cur_io_blocking, rc);


      printf("%d\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t",
               cur_time, zbd_->CalculateFreePercent(), cur_ops,cur_get_ops,cur_scan_ops,
               (gc_bytes_written>>20), cur_io_blocking, rc);
      printf("\n");
      // printf("A_CD_sequence = %lu / %lu = %lu\n", ( zbd_->total_deletion_after_copy_seq_.load()/100),
      // (total_count),A_CD_sequence);
      // printf("WA_CD_sequence = %lu / %lu = %lu\n", ( cost_benefit_score_sum_sequence_mb),
      // (total_deletion_after_copy_size >> 20),WA_CD_sequence);
      zbd_->PrintMisPredictStats();
      // printf("\n");
      fflush(stdout);
      fflush(stderr);

      // for (int i = 0; i < 5; i++) {
      //   total += file_size_dist[i];
      // }

      // if (total == 0) {
      //   memset(file_size_dist, 0, sizeof(file_size_dist));
      //   printf("[FileSizeDist] No files counted yet.\n");
      // }

      // double perc0 = (file_size_dist[0] * 100.0) / total;
      // double perc1 = (file_size_dist[1] * 100.0) / total;
      // double perc2 = (file_size_dist[2] * 100.0) / total;
      // double perc3 = (file_size_dist[3] * 100.0) / total;
      // double perc4 = (file_size_dist[4] * 100.0) / total;

      // printf("[FileSizeDist] Total files: %llu\n", (unsigned long
      // long)total); printf("  0 ~ 32MB :   %llu (%.2f%%)\n",
      //        (unsigned long long)file_size_dist[0], perc0);
      // printf("  33 ~ 63MB:   %llu (%.2f%%)\n",
      //        (unsigned long long)file_size_dist[1], perc1);
      // printf("  64 ~ 128MB:  %llu (%.2f%%)\n",
      //        (unsigned long long)file_size_dist[2], perc2);
      // printf(" 129 ~ 256MB:  %llu (%.2f%%)\n",
      //        (unsigned long long)file_size_dist[3], perc3);
      // printf(" Over 256MB:   %llu (%.2f%%)\n",
      //        (unsigned long long)file_size_dist[4], perc4);
      // // puts("");
      // fflush(stdout);
      // fflush(stderr);
    }
  }


//    { 
//     int cur_fops_sequence = zbd_->file_operation_sequence_.load();

//     std::lock_guard<std::mutex> lg_(files_mtx_);
//     for (auto f : files_) {
//       int tmp_n=0;
//       int tmp_seq = 0;

//       for (auto ext : f.second->GetExtents()) {
//         if (ext->is_zc_copied_ && ext->ZC_COPIED_STATE==ZC_COPIED) {

//               tmp_n++;
//               tmp_seq+=(cur_fops_sequence-ext->zc_copied_sequence_);
//           }

//       }
//       if(tmp_n){

//         printf("SEQ NOT DELETED UNTIL END %d\n",(tmp_seq/tmp_n));

//       }
//     }

// }
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
      return 1;
  }
}

//===================================================================//
// void ZenFS::NormalizeZoneLifetimes() {
//   double min_lifetime = std::numeric_limits<double>::max();
//   double max_lifetime = 0.0;

//   // 최소값과 최대값 계산
//   for (const auto& [zone_start, entry] : zone_lifetime_map_) {
//     double average_lifetime = entry.first / entry.second;
//     min_lifetime = std::min(min_lifetime, average_lifetime);
//     max_lifetime = std::max(max_lifetime, average_lifetime);
//   }

//   // 0~100 정규화
//   double range = max_lifetime - min_lifetime;
//   for (auto& [zone_start, entry] : zone_lifetime_map_) {
//     double average_lifetime = entry.first / entry.second;

//     uint64_t normalized_lifetime;
//     if (range > 0) {
//       normalized_lifetime = static_cast<uint64_t>(
//           ((average_lifetime - min_lifetime) * 100) / range);
//     } else {
//       normalized_lifetime = 50;
//     }

//     entry.first = normalized_lifetime;  // 정규화된 Lifetime 저장
//   }
// }

// uint64_t ZenFS::CalculateZoneLifetimeVariance() {
//   // 모든 존의 Lifetime 값을 모아 분산 계산
//   std::vector<uint64_t> lifetimes;
//   for (const auto& [zone_start, entry] : zone_lifetime_map_) {
//     lifetimes.push_back(entry.first);  // 정규화된 Lifetime 값
//   }

//   // 평균 계산
//   uint64_t sum = std::accumulate(lifetimes.begin(), lifetimes.end(), 0ULL);
//   uint64_t mean = sum / lifetimes.size();

//   // 분산 계산
//   uint64_t variance = 0;
//   for (const auto& lifetime : lifetimes) {
//     variance += (lifetime - mean) * (lifetime - mean);
//   }
//   return variance / lifetimes.size();
// }

// double ZenFS::CalculateZoneLifetimeVariance() {
//   size_t n = zone_lifetime_map_.size();
//   if (n == 0) {
//     return 0;
//   }

//   double sum = 0;
//   double sum_of_squares = 0;

//   for (const auto& [zone_start, entry] : zone_lifetime_map_) {
//     // double mean_lifetime = static_cast<double>(
//     //     static_cast<double>(entry.first) / entry.second * 100);
//     double total_lifetime = std::get<0>(entry);  // 존의 lifetime 총합
//     int file_count = std::get<1>(entry);         // 존에 포함된 파일 수
//     double mean_lifetime = total_lifetime / file_count;
//     sum += mean_lifetime;
//     sum_of_squares += mean_lifetime * mean_lifetime;
//   }

//   double mean = sum / n;

//   double variance = (sum_of_squares / n) - (mean * mean);

//   return variance * 1000;
// }

// void ZenFS::CalculateHorizontalLifetimes(
//     std::map<int, std::vector<std::pair<uint64_t, double>>>& level_file_map)
//     {
//   for (int level = 0; level < 6; level++) {
//     std::vector<uint64_t> fno_list;
//     std::set<uint64_t> compacting_files;

//     if (db_ptr_ != nullptr) {
//       zbd_->SameLevelFileList(level, fno_list, compacting_files);
//     }

//     std::vector<std::pair<uint64_t, double>> file_with_normalized_index;

//     // 컴팩션 중이 아닌 파일 수 계산
//     size_t num_non_compacting_files = fno_list.size() -
//     compacting_files.size();

//     for (size_t i = 0, non_compacting_index = 0; i < fno_list.size(); ++i) {
//       uint64_t fno = fno_list[i];
//       double normalized_index;
//       // fno = trial_move -> 0.0 (cold)

//       // 컴팩션 중인 파일은 인덱스를 1.0으로 설정
//       if (compacting_files.find(fno) != compacting_files.end()) {
//         normalized_index = 1.0;
//       } else {
//         // 정규화된 인덱스 계산
//         normalized_index =
//             level == 0
//                 ? 1.0
//                 : 1.0 - static_cast<double>(non_compacting_index) /
//                             static_cast<double>(num_non_compacting_files -
//                             1);
//         non_compacting_index++;
//       }

// 상위 레벨과 겹치는 파일은 최대 Lifetime으로 계산 if (level > 0) {
//   Slice smallest, largest;
//   if (zbd_->GetMinMaxKey(fno, smallest, largest)) {
//     std::vector<uint64_t> upper_fno_list;
//     zbd_->UpperLevelFileList(smallest, largest, level, upper_fno_list);

//     for (uint64_t upper_fno : upper_fno_list) {
//       auto it = std::find_if(level_file_map[level - 1].begin(),
//                              level_file_map[level - 1].end(),
//                              [&](const std::pair<uint64_t, double>& pair) {
//                                return pair.first == upper_fno;
//                              });
//       if (it != level_file_map[level - 1].end()) {
//         // std::cout << "Level : " << level
//         //           << "Comparing Lifetime: Current Max: "
//         //           << max_upper_lifetime << ", Upper File: " <<
//         //           upper_fno
//         //           << ",fno : " << fno << ", Lifetime: " <<
//         // it->second
//         //           << std::endl;
//         normalized_index = std::max(normalized_index, it->second);
//       }
//     }
//   }
// }

//       // 결과 저장
//       file_with_normalized_index.emplace_back(fno, normalized_index);
//     }
//     // 최종적으로 map에 저장
//     level_file_map[level] = std::move(file_with_normalized_index);
//   }

//   // for (const auto& level : level_file_map) {
//   //   std::cout << "Level " << level.first << ": [";
//   //   for (size_t i = 0; i < level.second.size(); ++i) {
//   //     std::cout << "(" << level.second[i].first << ", "
//   //               << level.second[i].second << ")";
//   //     if (i < level.second.size() - 1) {
//   //       std::cout << ", ";
//   //     }
//   //   }
//   //   std::cout << "]" << std::endl;
//   // }
// }

// void ZenFS::ReCalculateLifetimes() {
//   std::map<int, std::vector<std::pair<uint64_t, double>>> level_file_map;

//   CalculateHorizontalLifetimes(level_file_map);

//   zone_lifetime_map_.clear();

//   // 모든 level의 vertical_lifetime을 구하여 최소값과 최대값을 찾음
//   std::vector<double> vertical_lifetimes(6);
//   double max_vertical_lifetime = 0.0;
//   double min_vertical_lifetime = std::numeric_limits<double>::max();

//   for (int level = 0; level < 6; level++) {
//     double vertical_lifetime = zbd_->PredictCompactionScore(level);
//     vertical_lifetimes[level] = vertical_lifetime;
//     max_vertical_lifetime = std::max(max_vertical_lifetime,
//     vertical_lifetime); min_vertical_lifetime =
//     std::min(min_vertical_lifetime, vertical_lifetime);
//   }

//   std::vector<double> normalized_vertical_lifetimes(6);
//   double range = max_vertical_lifetime - min_vertical_lifetime;
//   // std::cout << "==================================================="
//   // << std::endl;
//   for (int level = 0; level < 6; level++) {
//     normalized_vertical_lifetimes[level] =
//         (range > 0)
//             ? (vertical_lifetimes[level] - min_vertical_lifetime) / range
//             : 0;  // 모든 값이 같으면 0으로 설정

//     // std::cout << "Level: " << level
//     //           << ", Original: " << vertical_lifetimes[level]
//     //           << ", Normalized: " << normalized_vertical_lifetimes[level]
//     //           << std::endl;
//   }
//   // std::cout << "==================================================="
//   // << std::endl;

//   double alpha_value = zbd_->GetAlphaValue();
//   double alpha_ = alpha_value;
//   double beta_ = 1 - alpha_;

//   // double alpha_ = alpha_value * 2;
//   // double beta_ = 1;

//   // 해당 레벨의 파일들에 대해 수평 및 수직 lifetime 계산
//   // for (const auto& file_pair : level_file_map[level]) {
//   //   uint64_t fno = file_pair.first;
//   //   double horizontal_lifetime = file_pair.second;
//   for (const auto& [level, file_lifetimes] : level_file_map) {
//     double vertical_lifetime_ = normalized_vertical_lifetimes[level];
//     for (const auto& [fno, horizontal_lifetime] : file_lifetimes) {
//       double sst_lifetime_value =
//           alpha_ * (1 - horizontal_lifetime) + beta_ * (1 -
//           vertical_lifetime_);

//       // std::cout << "Level: " << level
//       //           << ", vertical Lifetime: " << vertical_lifetime_
//       //           << ", horizontal_lifetime: " << horizontal_lifetime
//       //           << ", sst_lifetime: " << sst_lifetime_value << std::endl;

//       ZoneFile* zone_file = zbd_->GetSSTZoneFileInZBDNoLock(fno);
//       if (zone_file == nullptr) {
//         continue;
//       }
//       if (zone_file->IsDeleted()) {
//         continue;
//       }

//       // 각 ZoneFile의 extents를 순회하여 파일이 속한 존을 찾음
//       for (const auto* extent : zone_file->GetExtents()) {
//         uint64_t zone_start = extent->zone_->start_;  // 존의 시작 위치 사용

//         // 해당 존의 lifetime 합산 및 파일 수를 업데이트
//         if (zone_lifetime_map_.find(zone_start) == zone_lifetime_map_.end())
//         {
//           // zone_lifetime_map_[zone_start] = {sst_lifetime_value, 1};

//           zone_lifetime_map_[zone_start] = {
//               sst_lifetime_value, 1, {sst_lifetime_value}};

//         } else {
//           // zone_lifetime_map_[zone_start].first +=
//           //     sst_lifetime_value;                      // 수명 추가
//           // zone_lifetime_map_[zone_start].second += 1;  // 파일 수 추가

//           auto& entry = zone_lifetime_map_[zone_start];
//           std::get<0>(entry) += sst_lifetime_value;  // 총합에 추가
//           std::get<1>(entry) += 1;                   // 파일 수 증가
//           std::get<2>(entry).push_back(
//               sst_lifetime_value);  // 각 파일의 lifetime 저장
//         }
//       }
//     }
//   }
//   // }

//   // 각 존의 평균 lifetime 계산 및 출력
//   // for (const auto& zone_entry : zone_lifetime_map_) {
//   //   uint64_t zone_start = zone_entry.first;
//   //   double total_lifetime = zone_entry.second.first;
//   //   int file_count = zone_entry.second.second;

//   //   double average_lifetime =
//   //       total_lifetime / file_count;  // 존의 평균 lifetime 계산

//   //   std::cout << "Zone starting at " << zone_start
//   //             << " has average lifetime: " << average_lifetime <<
//   //             std::endl;
//   // }

//   /* wal -> hottest */
//   for (auto& kv : files_) {
//     const std::string& fname = kv.first;
//     std::shared_ptr<ZoneFile> zoneFile = kv.second;
//     (void)fname;

//     if (!zoneFile) {
//       continue;
//     }
//     if (zoneFile->IsDeleted()) {
//       continue;
//     }

//     if (zoneFile->is_wal_) {
//       // printf("RecalculateLifetimes - WAL\n");
//       double wal_lifetime_value = 0.0;

//       for (const auto* extent : zoneFile->GetExtents()) {
//         uint64_t zone_start = extent->zone_->start_;

//         // 처음 본 zone이면 초기화
//         if (zone_lifetime_map_.find(zone_start) == zone_lifetime_map_.end())
//         {
//           // printf("ReCalculateLifetimes-WAL : not found\n");
//           zone_lifetime_map_[zone_start] = {
//               0.0,  // 총합
//               1,    // 파일 수
//               {}    // 파일별 lifetime 리스트
//           };
//         }

//         auto& entry = zone_lifetime_map_[zone_start];
//         std::get<0>(entry) += wal_lifetime_value;
//         std::get<1>(entry) += 1;
//         std::get<2>(entry).push_back(wal_lifetime_value);
//       }
//     }
//   }
// }

//==========================================================================//

void ZenFS::CalculateHorizontalLifetimes(
    std::map<int, std::vector<FileInfo_>>& level_file_map) {
  for (int level = 0; level < 6; level++) {
    std::vector<uint64_t> fno_list;
    std::set<uint64_t> compacting_files;
    if (db_ptr_ != nullptr) {
      zbd_->SameLevelFileList(level, fno_list, compacting_files);
    }

    std::set<uint64_t> trivial_set;
    zbd_->TrivialMoveFiles(level, trivial_set);

    // for (auto fno : trivial_set) {
    //   // std::cout << "Trivial move candidate => level " << level
    //   //           << ", fno = " << fno << std::endl;
    // }
    const size_t total_files = fno_list.size();
    const size_t compacting_count = compacting_files.size();
    const size_t trivial_count = trivial_set.size();

    const size_t num_normal_files =
        total_files - compacting_count - trivial_count;

    std::vector<FileInfo_> file_info_vec;
    file_info_vec.reserve(total_files);

    size_t normal_index = 0;

    for (uint64_t fno : fno_list) {
      bool is_compacting =
          (compacting_files.find(fno) != compacting_files.end());
      bool is_trivial = (trivial_set.find(fno) != trivial_set.end());

      double normalized_index = 1.0;

      if (!is_compacting && !is_trivial) {
        if (level == 0) {
          // Level 0이면 무조건 1.0
          normalized_index = 1.0;
        } else {
          // Level > 0
          if (num_normal_files > 1) {
            normalized_index =
                1.0 -
                static_cast<double>(normal_index) /
                    static_cast<double>(num_normal_files - 1);  //(kminoverlap)
            // static_cast<double>(normal_index) /
            // static_cast<double>(num_normal_files - 1); //(OAZA)
          } else {
            normalized_index = 1.0;
          }
        }
        normal_index++;
      }

      FileInfo_ sst;
      sst.fno = fno;
      sst.horizontal_lifetime = normalized_index;
      sst.is_compacting = is_compacting;
      sst.is_trivial = is_trivial;
      sst.sst_lifetime_value_ = 0;

      file_info_vec.push_back(std::move(sst));
    }

    level_file_map[level] = std::move(file_info_vec);
  }
}

void ZenFS::ReCalculateLifetimes() {
  // std::map<int, std::vector<FileInfo_>> level_file_map;
  level_file_map_.clear();
  CalculateHorizontalLifetimes(level_file_map_);

  zone_lifetime_map_.clear();

  std::vector<double> vertical_lifetimes(6);
  double max_vertical_lifetime = 0.0;
  double min_vertical_lifetime = std::numeric_limits<double>::max();

  for (int level = 0; level < 6; level++) {
    double vertical_lifetime = zbd_->PredictCompactionScore(level);
    vertical_lifetimes[level] = vertical_lifetime;
    max_vertical_lifetime = std::max(max_vertical_lifetime, vertical_lifetime);
    min_vertical_lifetime = std::min(min_vertical_lifetime, vertical_lifetime);
  }

  std::vector<double> normalized_vertical_lifetimes(6);
  double range = max_vertical_lifetime - min_vertical_lifetime;

  for (int level = 0; level < 6; level++) {
    normalized_vertical_lifetimes[level] =
        (range > 0)
            ? (vertical_lifetimes[level] - min_vertical_lifetime) / range
            : 0;  // 모든 값이 같으면 0으로 설정

    // std::cout << "Level: " << level
    //           << ", Original: " << vertical_lifetimes[level]
    //           << ", Normalized: " << normalized_vertical_lifetimes[level]
    //           << std::endl;
  }

  double alpha_value = zbd_->GetAlphaValue();
  double alpha_ = alpha_value;
  double beta_ = 1.0 - alpha_;

  for (auto& [level, file_infos] : level_file_map_) {
    double vertical_lifetime = normalized_vertical_lifetimes[level];

    for (auto& sst : file_infos) {
      double sst_lifetime_value = 0.0;
      if (sst.is_compacting) {
        sst_lifetime_value = 0.0;
      } else if (sst.is_trivial) {
        sst_lifetime_value = 1.0;
      } else {
        sst_lifetime_value = alpha_ * (1 - sst.horizontal_lifetime) +
                             beta_ * (1 - vertical_lifetime);
      }
      sst.sst_lifetime_value_ = sst_lifetime_value;
    }
  }

  uint64_t predict_cnt = zbd_->GetPredictCnt();

  // std::cout << "predict cnt : " << predict_cnt << std::endl;

  PredictCompaction(predict_cnt);

  // printf("prdeictcompaction completed\n");

  for (auto& [level, file_infos] : level_file_map_) {
    for (auto& sst : file_infos) {
      ZoneFile* zone_file = zbd_->GetSSTZoneFileInZBDNoLock(sst.fno);
      if (!zone_file || zone_file->IsDeleted()) {
        continue;
      }

      double sst_lifetime_value = sst.sst_lifetime_value_;
      for (const auto* extent : zone_file->GetExtents()) {
        uint64_t zone_start = extent->zone_->start_;  

        auto& entry = zone_lifetime_map_[zone_start];

        entry.total_lifetime += sst_lifetime_value;
        entry.file_count += 1;
        entry.file_lifetimes[sst.fno] = sst_lifetime_value;
      }
    }
  }

  // ====================================================================
  // for (auto& [level, file_infos] : level_file_map_) {
  //   double vertical_lifetime = normalized_vertical_lifetimes[level];
  //   for (auto& sst : file_infos) {
  //     // 컴팩션 중은 0.0, trivial move는 1.0, 그 외는 기존 계산
  //     double sst_lifetime_value = 0.0;
  //     if (sst.is_compacting) {
  //       sst_lifetime_value = 0.0;
  //     } else if (sst.is_trivial) {
  //       sst_lifetime_value = 1.0;
  //     } else {
  //       sst_lifetime_value = alpha_ * (1 - sst.horizontal_lifetime) +
  //                            beta_ * (1 - vertical_lifetime);
  //     }

  //     ZoneFile* zone_file = zbd_->GetSSTZoneFileInZBDNoLock(sst.fno);
  //     if (zone_file == nullptr) {
  //       continue;
  //     }
  //     if (zone_file->IsDeleted()) {
  //       continue;
  //     }

  //     sst.sst_lifetime_value_ = sst_lifetime_value;

  //     PredictCompaction(4);

  //     // 각 ZoneFile의 extents를 순회하여 파일이 속한 존을 찾음
  //     for (const auto* extent : zone_file->GetExtents()) {
  //       uint64_t zone_start = extent->zone_->start_;  // 존의 시작 위치 사용

  //       // 해당 존의 lifetime 합산 및 파일 수를 업데이트
  //       if (zone_lifetime_map_.find(zone_start) == zone_lifetime_map_.end())
  //       {
  //         zone_lifetime_map_[zone_start] = {
  //             sst_lifetime_value, 1, {sst_lifetime_value}};

  //       } else {
  //         auto& entry = zone_lifetime_map_[zone_start];
  //         std::get<0>(entry) += sst_lifetime_value;  // 총합에 추가
  //         std::get<1>(entry) += 1;                   // 파일 수 증가
  //         std::get<2>(entry).push_back(
  //             sst_lifetime_value);  // 각 파일의 lifetime 저장
  //       }
  //     }
  //   }
  // }
  // }===============================================================

  // 각 존의 평균 lifetime 계산 및 출력
  // for (const auto& zone_entry : zone_lifetime_map_) {
  //   uint64_t zone_start = zone_entry.first;
  //   double total_lifetime = zone_entry.second.first;
  //   int file_count = zone_entry.second.second;

  //   double average_lifetime =
  //       total_lifetime / file_count;  // 존의 평균 lifetime 계산

  //   std::cout << "Zone starting at " << zone_start
  //             << " has average lifetime: " << average_lifetime <<
  //             std::endl;
  // }

  /* wal -> hottest */
  for (auto& kv : files_) {
    const std::string& fname = kv.first;
    std::shared_ptr<ZoneFile> zoneFile = kv.second;
    (void)fname;

    if (!zoneFile) {
      continue;
    }
    if (zoneFile->IsDeleted()) {
      continue;
    }

    if (zoneFile->is_wal_) {
      // printf("RecalculateLifetimes - WAL\n");
      double wal_lifetime_value = 0.0;

      for (const auto* extent : zoneFile->GetExtents()) {
        uint64_t zone_start = extent->zone_->start_;

        auto& entry = zone_lifetime_map_[zone_start];
        entry.total_lifetime += wal_lifetime_value;
        entry.file_count++;
      }
    }
  }
}

// trivial move -> iteration 포함 안함

// void ZenFS::PredictCompaction(int step) {
//   std::array<uint64_t, 10> tmp_lsm_tree = zbd_->GetCurrentLSMTree();
//   // std::set<uint64_t> fno_already_propagated;
//   std::unordered_set<int> excluded_levels;
//   excluded_levels.clear();
//   int initial_l0_files_n = 0;
//   {
//     auto it_l0 = level_file_map_.find(0);
//     if (it_l0 != level_file_map_.end()) {
//       for (const auto& file : it_l0->second) {
//         if (!file.is_compacting) {
//           initial_l0_files_n++;
//         }
//       }
//     }
//   }
//   // printf("Initial L0 non-compacting files: %d\n", initial_l0_files_n);

//   fno_already_propagated.clear();
//   fno_not_should_selected_as_pivot_again.clear();
//   while (step > 0) {
//     std::vector<uint64_t> unpivot_fno_list;
//     unpivot_fno_list.clear();

//     // PredictCompactionImpl(pivot_level, tmp_lsm_tree, pivot_fno,
//     //                       unpivot_fno_list, initial_l0_files_n);
//     uint64_t pivot_level = GetMaxLevelScoreLevel(
//         tmp_lsm_tree, initial_l0_files_n, excluded_levels);
//     if (pivot_level == (uint64_t)-1) {
//       // printf("[PredictCompaction] No available level!\n");
//       break;
//     }

//     uint64_t pivot_fno = GetMaxHorizontalFno(pivot_level);
//     // printf("[Impl] pivot_level=%lu, pivot_fno=%lu\n", pivot_level,
//     // pivot_fno);

//     if (pivot_fno == 0) {
//       excluded_levels.insert(pivot_level);

//       continue;

//       // pivot_level = GetMaxLevelScoreLevel(tmp_lsm_tree, initial_l0_files_n,
//       //                                     excluded_levels);

//       // pivot_fno = GetMaxHorizontalFno(pivot_level);
//       // // printf("[Fix] pivot_fno was 0, so reselect => level=%lu, fno=%lu\n",
//       // //        pivot_level, pivot_fno);

//       // if (pivot_fno != 0) {
//       //   GetOverlappingFno(pivot_fno, pivot_level, unpivot_fno_list);
//       // } else {
//       //   // printf("[Fix] Still no pivot fno => reduce step, go next loop\n");
//       //   step--;
//       //   continue;
//       // }

//     }else{
//       GetOverlappingFno(pivot_fno, pivot_level, unpivot_fno_list);

//     }
    


//     // if (unpivot_fno_list.empty()) {
//     //   // if trivial move, return
//     //   // printf(
//     //   //     "[PredictCompaction] unpivot_fno_list is empty. pivot_fno=%lu, "
//     //   //     "level=%lu\n",
//     //   // //     pivot_fno, pivot_level);

//     //   // fno_already_propagated.insert(pivot_fno);

//     //   // continue;
//     //   // return;
//     // }

//     if (fno_already_propagated.find(pivot_fno) !=
//         fno_already_propagated.end()) {
//       // fno_already_propagated.insert(pivot_fno);
//       continue;
//     }



//     // for (auto it = unpivot_fno_list.begin(); it != unpivot_fno_list.end();) {
//     //   if (fno_already_propagated.find(*it) != fno_already_propagated.end()) {
//     //     it = unpivot_fno_list.erase(it);
//     //   } else {
//     //     ++it;
//     //   }
//     // }
//     // }  
//     bool should_not_selected_again = false;

//     for (auto it = unpivot_fno_list.begin(); it != unpivot_fno_list.end();) {
//       if (fno_already_propagated.find(*it) != fno_already_propagated.end()) {
//         // it = unpivot_fno_list.erase(it);
//         // printf("Removed an already propagated fno from unpivot_fno_list.\n");
//         // fno_not_should_selected_as_pivot_again.insert(pivot_fno);

//         should_not_selected_again = true;

//         // continue;
//         break;
//       } else {
//         ++it;
//       }
//     }

//     if (should_not_selected_again == true) {
//       fno_not_should_selected_as_pivot_again.insert(pivot_fno);
//       continue;
//     }

//     ZoneFile* pivot_file = zbd_->GetSSTZoneFileInZBDNoLock(pivot_fno);
//     if (pivot_file == nullptr) {
//       // printf("[PredictCompaction] no pivot file (fno=%lu)\n", pivot_fno);
//       continue;
//     }
//     if (pivot_file->IsDeleted()) {
//       // printf("[PredictCompaction] pivot_file is deleted (fno=%lu)\n",
//       //  pivot_fno);
//       continue;
//     }

//     if (pivot_level == 0) {
//       uint64_t total_l0_size = 0;
//       std::vector<uint64_t> l0_files;

//       auto it = level_file_map_.find(0);
//       if (it != level_file_map_.end()) {
//         auto& file_infos = it->second;
//         for (auto& sst : file_infos) {
//           uint64_t sst_size = 0;
//           ZoneFile* zf = zbd_->GetSSTZoneFileInZBDNoLock(sst.fno);
//           if (zf == nullptr) {
//             continue;
//           }
//           if (zf->IsDeleted()) {
//             continue;
//           }
//           sst_size = zf->GetFileSize();

//           total_l0_size += sst_size;
//           l0_files.push_back(sst.fno);
//         }
//       }
//       // std::cout << "total_l0_size: " << total_l0_size << std::endl;
//       // std::cout << "tmplsmtree[0]" << tmp_lsm_tree[0] << std::endl;

//       uint64_t unpivot_total = 0;
//       for (auto fno : unpivot_fno_list) {
//         ZoneFile* zf = zbd_->GetSSTZoneFileInZBDNoLock(fno);
//         if (!zf || zf->IsDeleted()) {
//           continue;
//         }
//         unpivot_total += zf->GetFileSize();
//       }

//       uint64_t total_input_size = total_l0_size + unpivot_total;

//       double compressibility = 0.0;

//       compressibility = zbd_->GetAvgCompressibilityOflevel(pivot_level + 1);

//       if (compressibility <= 0.0) {
//         compressibility = 1.0;
//       }

//       uint64_t compressed_size = static_cast<uint64_t>(
//           static_cast<double>(total_input_size) * compressibility);

//       initial_l0_files_n = 0;

//       if (tmp_lsm_tree[0] < total_l0_size) {
//         tmp_lsm_tree[0] = 0;
//       } else {
//         tmp_lsm_tree[0] -= total_l0_size;
//       }
//       // std::cout << "After tmplsmtree[0]" << tmp_lsm_tree[0] << std::endl;
//       // std::cout << "tmplsmtree[1]" << tmp_lsm_tree[1] << std::endl;

//       if (tmp_lsm_tree[1] < unpivot_total) {
//         tmp_lsm_tree[1] = 0;
//       } else {
//         tmp_lsm_tree[1] -= unpivot_total;
//       }

//       // tmp_lsm_tree[1] += total_l0_size;
//       tmp_lsm_tree[1] += compressed_size;

//       // std::cout << "After tmplsmtree[1]" << tmp_lsm_tree[1] << std::endl;

//       for (auto fno : l0_files) {
//         fno_already_propagated.insert(fno);
//       }

//       for (auto& f : unpivot_fno_list) {
//         fno_already_propagated.insert(f);
//       }

//       Propagation(pivot_fno, unpivot_fno_list);

//       step--;
//       // printf("l0 step!\n");
//       continue;  // L0는 전부 끝났으니 다음 루프
//     }
    
//     // !L0
//     uint64_t file_size = pivot_file->GetFileSize();
//     uint64_t unpivot_total = 0;
//     for (auto fno : unpivot_fno_list) {
//       ZoneFile* zf = zbd_->GetSSTZoneFileInZBDNoLock(fno);
//       if (!zf || zf->IsDeleted()) {
//         continue;
//       }
//       unpivot_total += zf->GetFileSize();
//     }
//     uint64_t total_input_size = file_size + unpivot_total;

//     // std::cout << "pivot_file size: " << file_size << std::endl;
//     // tmp_lsm_tree[pivot_level] -= pivot_file->GetFileSize();

//     double compressibility = 0.0;

//     compressibility = zbd_->GetAvgCompressibilityOflevel(pivot_level + 1);
//     if (unpivot_total == 0) {
//       compressibility = 1.0;
//     }
//     // if (compressibility == 1.0) {
//     //   compressibility = 1.0;
//     // }

//     uint64_t compressed_size = static_cast<uint64_t>(
//         static_cast<double>(total_input_size) * compressibility);

//     if (tmp_lsm_tree[pivot_level] < file_size) {
//       tmp_lsm_tree[pivot_level] = 0;
//     } else {
//       tmp_lsm_tree[pivot_level] -= file_size;
//     }
//     if (tmp_lsm_tree[pivot_level + 1] < unpivot_total) {
//       tmp_lsm_tree[pivot_level + 1] = 0;
//     } else {
//       tmp_lsm_tree[pivot_level + 1] -= unpivot_total;
//     }
//     // tmp_lsm_tree[pivot_level + 1] += file_size;
//     tmp_lsm_tree[pivot_level + 1] += compressed_size;

//     fno_already_propagated.insert(pivot_fno);
//     for (auto& f : unpivot_fno_list) {
//       fno_already_propagated.insert(f);
//     }
//     Propagation(pivot_fno, unpivot_fno_list);
//     step--;
//     // printf("step!\n");
//   }
// }


void ZenFS::PredictCompaction(int step) {
  std::array<uint64_t, 10> tmp_lsm_tree = zbd_->GetCurrentLSMTree();
  // std::set<uint64_t> fno_already_propagated;
  std::unordered_set<int> excluded_levels;
  excluded_levels.clear();
  int initial_l0_files_n = 0;
  {
    auto it_l0 = level_file_map_.find(0);
    if (it_l0 != level_file_map_.end()) {
      for (const auto& file : it_l0->second) {
        if (!file.is_compacting) {
          initial_l0_files_n++;
        }
      }
    }
  }
  // printf("Initial L0 non-compacting files: %d\n", initial_l0_files_n);

  fno_already_propagated.clear();
  while (step > 0) {
    // uint64_t pivot_level = 0;
    // uint64_t pivot_fno = 0;
    std::vector<uint64_t> unpivot_fno_list;
    unpivot_fno_list.clear();

    // PredictCompactionImpl(pivot_level, tmp_lsm_tree, pivot_fno,
    //                       unpivot_fno_list, initial_l0_files_n);
    uint64_t pivot_level = GetMaxLevelScoreLevel(
        tmp_lsm_tree, initial_l0_files_n, excluded_levels);
    if (pivot_level == (uint64_t)-1) {
      printf("[PredictCompaction] No available level!\n");
      break;
    }

    uint64_t pivot_fno = GetMaxHorizontalFno(pivot_level);
    // printf("[Impl] pivot_level=%lu, pivot_fno=%lu\n", pivot_level,
    // pivot_fno);

    if (pivot_fno == 0) {
      excluded_levels.insert(pivot_level);
      pivot_level = GetMaxLevelScoreLevel(tmp_lsm_tree, initial_l0_files_n,
                                          excluded_levels);

      pivot_fno = GetMaxHorizontalFno(pivot_level);
      // printf("[Fix] pivot_fno was 0, so reselect => level=%lu, fno=%lu\n",
      //        pivot_level, pivot_fno);

      // if (pivot_fno != 0) {
      //   GetOverlappingFno(pivot_fno, pivot_level, unpivot_fno_list);
      // } else {
      //   printf("[Fix] Still no pivot fno => reduce step, go next loop\n");
      //   step--;
      //   continue;
      // }
      if(pivot_fno==0){
        step--;
        continue;
      }
    }

    GetOverlappingFno(pivot_fno, pivot_level, unpivot_fno_list);

    if (fno_already_propagated.find(pivot_fno) !=
        fno_already_propagated.end()) {
      // printf("fno_already_propagated\n");
      // continue;
      return;
    }

    for (auto it = unpivot_fno_list.begin(); it != unpivot_fno_list.end();) {
      if (fno_already_propagated.find(*it) != fno_already_propagated.end()) {
        // 제거
        it = unpivot_fno_list.erase(it);
        // printf("Removed an already propagated fno from unpivot_fno_list.\n");
      } else {
        ++it;
      }
      
      // return;
      
    }

    ZoneFile* pivot_file = zbd_->GetSSTZoneFileInZBDNoLock(pivot_fno);
    if (pivot_file == nullptr) {
      printf("[PredictCompaction] no pivot file (fno=%lu)\n", pivot_fno);
      continue;
    }
    if (pivot_file->IsDeleted()) {
      printf("[PredictCompaction] pivot_file is deleted (fno=%lu)\n",
             pivot_fno);
      continue;
    }

    if (pivot_level == 0) {
      uint64_t total_l0_size = 0;
      std::vector<uint64_t> l0_files;

      auto it = level_file_map_.find(0);
      if (it != level_file_map_.end()) {
        auto& file_infos = it->second;
        for (auto& sst : file_infos) {
          uint64_t sst_size = 0;
          ZoneFile* zf = zbd_->GetSSTZoneFileInZBDNoLock(sst.fno);
          if (zf == nullptr) {
            continue;
          }
          if (zf->IsDeleted()) {
            continue;
          }
          sst_size = zf->GetFileSize();

          total_l0_size += sst_size;
          l0_files.push_back(sst.fno);
        }
      }
      // std::cout << "total_l0_size: " << total_l0_size << std::endl;
      // std::cout << "tmplsmtree[0]" << tmp_lsm_tree[0] << std::endl;
      initial_l0_files_n = 0;
      // tmp_lsm_tree[0] -= total_l0_size;
      if (tmp_lsm_tree[0] < total_l0_size) {
        tmp_lsm_tree[0] = 0;
      } else {
        tmp_lsm_tree[0] -= total_l0_size;
      }
      // std::cout << "After tmplsmtree[0]" << tmp_lsm_tree[0] << std::endl;
      // std::cout << "tmplsmtree[1]" << tmp_lsm_tree[1] << std::endl;
      tmp_lsm_tree[1] += total_l0_size;
      // std::cout << "After tmplsmtree[1]" << tmp_lsm_tree[1] << std::endl;
      for (auto fno : l0_files) {
        fno_already_propagated.insert(fno);
      }

      Propagation(pivot_fno, unpivot_fno_list);

      step--;
      // printf("l0 step!\n");
      continue;  // L0는 전부 끝났으니 다음 루프
    }
    // !L0
    uint64_t file_size = pivot_file->GetFileSize();
    // std::cout << "pivot_file size: " << file_size << std::endl;
    // tmp_lsm_tree[pivot_level] -= pivot_file->GetFileSize();
    if (tmp_lsm_tree[pivot_level] < file_size) {
      tmp_lsm_tree[pivot_level] = 0;
    } else {
      tmp_lsm_tree[pivot_level] -= file_size;
    }
    tmp_lsm_tree[pivot_level + 1] += file_size;

    fno_already_propagated.insert(pivot_fno);
    for (auto& f : unpivot_fno_list) {
      fno_already_propagated.insert(f);
    }
    Propagation(pivot_fno, unpivot_fno_list);
    step--;
    // printf("step!\n");
  }
}

void ZenFS::GetOverlappingFno(uint64_t pivot_fno, uint64_t pivot_level,
                              std::vector<uint64_t>& unpivot_fno_list) {
  Slice smallest, largest;
  if (zbd_->GetMinMaxKey(pivot_fno, smallest, largest)) {
    zbd_->DownwardAdjacentFileList(smallest, largest, pivot_level,
                                   unpivot_fno_list);
  }
}

void ZenFS::Propagation(uint64_t pivot_fno,
                        std::vector<uint64_t> unpivot_fno_list) {
  double pivot_lifetime = 0.0;
  bool found_pivot = false;

  for (auto& [level, file_infos] : level_file_map_) {
    for (auto& sst : file_infos) {
      if (sst.fno == pivot_fno) {
        pivot_lifetime = sst.sst_lifetime_value_;
        found_pivot = true;
        break;
      }
    }
    if (found_pivot) break;
  }

  if (!found_pivot) {
    // std::cout << "[Propagation] pivot_fno(" << pivot_fno << ") not found.\n";
    return;
  }

  // unpivot_fno_list에 있는 모든 파일의 sst_lifetime_value_를
  // pivot_lifetime으로 설정
  for (uint64_t unpivot_fno : unpivot_fno_list) {
    bool found_unpivot = false;
    for (auto& [level, file_infos] : level_file_map_) {
      for (auto& sst : file_infos) {
        if (sst.fno == unpivot_fno) {
          sst.sst_lifetime_value_ = pivot_lifetime;
          found_unpivot = true;
          zbd_->AddPropagationCount(1);
          break;
        }
      }
      if (found_unpivot) break;
    }
  }
}

uint64_t ZenFS::GetMaxLevelScoreLevel(
    std::array<uint64_t, 10>& tmp_lsm_tree, int initial_l0_files_n,
    std::unordered_set<int>& excluded_levels) {
  int max_level = -1;
  double max_score = -1;

  for (int level = 0; level < 6; ++level) {
    if (excluded_levels.find(level) != excluded_levels.end()) {
      // printf("Skipping level %d\n", level);
      continue;
    }
    double score = zbd_->PredictCompactionScoreTmp(level, tmp_lsm_tree,
                                                   initial_l0_files_n);

    if (score > max_score) {
      max_score = score;
      max_level = level;
    }
  }

  if (max_score <= 0.0) {
    return (uint64_t)-1;
  }

  return max_level;
}

uint64_t ZenFS::GetMaxHorizontalFno(int pivot_level) {
  auto it = level_file_map_.find(pivot_level);
  if (it == level_file_map_.end() || it->second.empty()) {
    // printf("no level || no files");
  }

  const auto& files = it->second;

  uint64_t max_fno = 0;
  double max_horizontal_lifetime = -1.0;

  // printf("Files in level %d: %zu\n", pivot_level, files.size());

  for (const auto& file : files) {
    if (fno_already_propagated.find(file.fno) != fno_already_propagated.end()) {
      // printf("getmaxhorizontalfno skip\n");
      continue;
    }

    if (fno_not_should_selected_as_pivot_again.find(file.fno) !=
        fno_not_should_selected_as_pivot_again.end()) {
      // printf("getmaxhorizontalfno skip\n");
      continue;
    }
    if (file.horizontal_lifetime > max_horizontal_lifetime) {
      max_horizontal_lifetime = file.horizontal_lifetime;
      max_fno = file.fno;
    }
  }

  return max_fno;
}

//====================================================================================//

// void ZenFS::Adv_ReCalculateLifetimes() {
//   std::map<int, std::vector<std::pair<uint64_t, double>>> level_file_map;

//   CalculateHorizontalLifetimes(level_file_map);

//   zone_lifetime_map_.clear();

//   std::vector<double> vertical_lifetimes(6);
//   double max_vertical_lifetime = 0.0;
//   double min_vertical_lifetime = std::numeric_limits<double>::max();

//   for (int level = 0; level < 6; ++level) {
//     double vertical_lifetime = zbd_->PredictCompactionScore(level);
//     vertical_lifetimes[level] = vertical_lifetime;
//     max_vertical_lifetime = std::max(max_vertical_lifetime,
//     vertical_lifetime); min_vertical_lifetime =
//     std::min(min_vertical_lifetime, vertical_lifetime);
//   }

//   std::vector<double> normalized_vertical_lifetimes(6);
//   double range = max_vertical_lifetime - min_vertical_lifetime;
//   for (int level = 0; level < 6; ++level) {
//     normalized_vertical_lifetimes[level] =
//         (range > 0)
//             ? (vertical_lifetimes[level] - min_vertical_lifetime) / range
//             : 0;
//   }

//   double alpha_value = zbd_->GetAlphaValue();
//   double alpha_ = alpha_value;
//   double beta_ = 1 - alpha_;

//   // previous level에 original lifetime 저장
//   std::map<uint64_t, double> previous_level_map;

//   for (const auto& [level, file_lifetimes] : level_file_map) {
//     double vertical_lifetime_ = normalized_vertical_lifetimes[level];

//     for (const auto& [fno, horizontal_lifetime] : file_lifetimes) {
//       /* sst lifetime value는 높을수록 cold, 낮을수록 hot */
//       double sst_lifetime_value =
//           alpha_ * (1 - horizontal_lifetime) + beta_ * (1 -
//           vertical_lifetime_);

//       // std::cout << "Initial Lifetime (Level " << level << ", SSTable "
//       << fno
//           //           << "): " << sst_lifetime_value << std::endl;

//           /*첫순회는 level0이라 upper 로직 패스 후 zone_lifetime_map에
//           Initial
//            lifetime저장하고 previous_level_map에 Initial lifetime들 저장*/

//           /*이후 level i부터 upper level과 겹치는 파일들을
//           previous_level_map에
//            저장된 level i-1의 Initial lifetime과 level i의 changable
//            lifetime을 비교해 Hottest value로 치환, 전파된 값은
//            zone_lifetime_map에 저장.*/

//           /* Initial lifetime 저장구조를 둠으로써 독립적인 전파를 시행*/

//           if (level > 0) {
//         Slice smallest, largest;
//         if (zbd_->GetMinMaxKey(fno, smallest, largest)) {
//           std::vector<uint64_t> upper_fno_list;
//           zbd_->UpperLevelFileList(smallest, largest, level, upper_fno_list);

//           for (uint64_t upper_fno : upper_fno_list) {
//             auto it = previous_level_map.find(upper_fno);
//             if (it != previous_level_map.end()) {
//               /* hottest value 선택 */
//               sst_lifetime_value = std::min(sst_lifetime_value, it->second);
//             }
//           }
//         }
//       }

//       // std::cout << "Final Lifetime (Level " << level << ", SSTable " <<
//       fno
//           //           << "): " << sst_lifetime_value << std::endl;

//           ZoneFile* zone_file = zbd_->GetSSTZoneFileInZBDNoLock(fno);
//       if (zone_file != nullptr) {
//         for (const auto* extent : zone_file->GetExtents()) {
//           uint64_t zone_start = extent->zone_->start_;

//           if (zone_lifetime_map_.find(zone_start) ==
//           zone_lifetime_map_.end()) {
//             // zone_lifetime_map_[zone_start] = {sst_lifetime_value, 1};
//             zone_lifetime_map_[zone_start] = {
//                 sst_lifetime_value, 1, {sst_lifetime_value}};
//           } else {
//             // zone_lifetime_map_[zone_start].first += sst_lifetime_value;
//             // zone_lifetime_map_[zone_start].second += 1;
//             auto& entry = zone_lifetime_map_[zone_start];
//             std::get<0>(entry) += sst_lifetime_value;  // 총합에 추가
//             std::get<1>(entry) += 1;                   // 파일 수 증가
//             std::get<2>(entry).push_back(
//                 sst_lifetime_value);  // 각 파일의 lifetime 저장
//           }
//         }
//       } else {
//         std::cerr << "Error: Cannot find ZoneFile for fno: " << fno
//                   << std::endl;
//       }
//     }

//     previous_level_map.clear();
//     for (const auto& [fno, horizontal_lifetime] : file_lifetimes) {
//       previous_level_map[fno] =
//           alpha_ * (1 - horizontal_lifetime) +
//           beta_ * (1 - normalized_vertical_lifetimes[level]);
//     }
//   }
// }

void ZenFS::ZoneCleaning(bool forced) {
  uint64_t zc_scheme = zbd_->GetZCScheme();
  if (db_ptr_ == nullptr) {
    // printf("ZenFS::ZoneCleaning - db_ptr is nullptr!!");
    return;
  }
  if(zbd_->coldest_type_set_==false){
    for(int i =0 ;i<10;i++){
      zbd_->check_coldest_[i]=false;
    }
    // coldest_type_ = zbd_->GetCBSCColdestType();
    zbd_->SetCBSCColdestType();
    zbd_->coldest_type_set_= true;
  }

  struct timespec start_ts, end_ts;
  clock_gettime(CLOCK_MONOTONIC, &start_ts);

  ReCalculateLifetimes();

  // NormalizeZoneLifetimes();
  // double cur_variance = CalculateZoneLifetimeVariance();
  // std::cout << "Zone Lifetime Variance: " << cur_variance << std::endl;
  clock_gettime(CLOCK_MONOTONIC, &end_ts);
  long elapsed_ns_ts = (end_ts.tv_sec - start_ts.tv_sec) * 1000000000 +
                       (end_ts.tv_nsec - start_ts.tv_nsec);
  zbd_->AddCalculatelifetimeLapse(elapsed_ns_ts);
  // }
  // printf("zonecleaning->zc_scheme : %lu\n", zc_scheme);
  // uint64_t zone_size = zbd_->GetZoneSize();
  size_t should_be_copied = 0;
  double ZLV = 0;
  int start = GetMountTime();
  struct timespec start_timespec, end_timespec;

  ZenFSSnapshot snapshot;
  ZenFSSnapshotOptions options;
  options.zone_ = 1;
  options.zone_file_ = 1;
  options.log_garbage_ = 1;

  GetZenFSSnapshot(snapshot, options);
  size_t all_inval_zone_n = 0;
  // std::vector<std::pair<uint64_t, uint64_t>> victim_candidate;
  struct ZoneInfo {
    double score;
    uint64_t zone_start;
    uint64_t garbage_percent_approx;
    double ZoneLifetimeValue;
  };
  std::vector<ZoneInfo> victim_candidate;
  std::set<uint64_t> migrate_zones_start;

  for (const auto& zone : snapshot.zones_) {
    // zone.capacity == 0 -> full-zone
    if (zone.capacity != 0) {  // skip not-full zone
      continue;
    }
    if (zone.used_capacity == zone.max_capacity) {  // skip all valid zone
      continue;
    }
    uint64_t garbage_percent_approx =
        100 - 100 * zone.used_capacity / zone.max_capacity;  
    if (zone.used_capacity > 0) {  // 유효 데이터(valid data)가 있는 경우

      if (zc_scheme == CBZC1 || zc_scheme == CBZC2) {
        struct timespec start_age_ts, end_age_ts;
        clock_gettime(CLOCK_MONOTONIC, &start_age_ts);
        auto current_time = std::chrono::system_clock::now();
        uint64_t total_age = 0;
        rocksdb::IOOptions io_options;
        auto lifetime_hints = GetWriteLifeTimeHints();
        for (const auto& zone_file : snapshot.zone_files_) {
          for (const auto& extent : zone_file.extents) {
            if (extent.zone_start == zone.start) {
              if (zc_scheme == CBZC1) {
                // printf("CBZC1!!!\n");
                // CBZC1 - file creation time based
                uint64_t file_mod_time = 0;
                IOStatus s = GetFileModificationTime(
                    extent.filename, io_options, &file_mod_time, nullptr);

                if (s.ok()) {
                  uint64_t file_age =
                      std::chrono::duration_cast<std::chrono::seconds>(
                          current_time.time_since_epoch())
                          .count() -
                      file_mod_time;
                  total_age += file_age;
                }
                // std::cout << "Total_age: " << total_age << std::endl;
              } else {
                // CBZC2 - LIZC
                // printf("CBZC2!!!\n");
                auto it = lifetime_hints.find(extent.filename);
                if (it != lifetime_hints.end()) {
                  Env::WriteLifeTimeHint hint = it->second;
                  // uint64_t file_age = EstimateFileAge(hint);
                  uint64_t file_age = hint;
                  total_age += file_age;
                }
              }
            }
          }
        }
        clock_gettime(CLOCK_MONOTONIC, &end_age_ts);
        long elapsed_ns_age =
            (end_age_ts.tv_sec - start_age_ts.tv_sec) * 1000000000 +
            (end_age_ts.tv_nsec - start_age_ts.tv_nsec);
        zbd_->AddCalculatelifetimeLapse(elapsed_ns_age);

        /* cost-benefit */
        /* benefit/cost = (free space generated * age of data) / cost
                        = (1-u) * age / (1+u) */
        /* cost = u(유효데이터(zone.used_capacity)(비율)를 읽는 비용)
                + u(valid data copy) = 2u */
        /* age = 세그먼트 내 가장 최근에 수정된 블록의 시간을 공간이 여유
          상태를 유지할 시간의 추정치로 사용(즉, 가장 최근에 수정된 블록의
          나이)
         */
        /* u = valid * garbage_percent_approx(free space+invalid)
             = 1 - valid = 1-u  ex) 80 %
         cost = 2u = (100 - gpa) * 2  */

        uint64_t cost = (100 - garbage_percent_approx) * 2;
        uint64_t benefit = garbage_percent_approx * total_age;
        if (cost != 0) {
          double cost_benefit_score = benefit / cost;
          victim_candidate.push_back(
              {cost_benefit_score, zone.start, garbage_percent_approx, 0.0});
        }
      } else if (zc_scheme == CBZC5) {
        auto now = std::chrono::system_clock::now();
        auto age = std::chrono::duration_cast<std::chrono::milliseconds>(
                       now - zone.recent_inval_time)
                       .count();

        // std::cout << "Zone age (ms): " << age << std::endl;

        double cost = 2 * (static_cast<double>(zone.used_capacity) /
                           static_cast<double>(zone.max_capacity));

        double freeSpace =
            static_cast<double>(zone.max_capacity - zone.used_capacity) /
            static_cast<double>(zone.max_capacity);

        double benefit = freeSpace * age;

        if (cost != 0) {
          double cost_benefit_score =
              static_cast<double>(benefit) / static_cast<double>(cost);
          victim_candidate.push_back(
              {cost_benefit_score, zone.start, garbage_percent_approx, 0.0});
          // std::cout << "cost_benefit_score : " << cost_benefit_score
          //           << std::endl;
        }

        // std::cout << "cost : " << cost << std::endl;
        // std::cout << "freeSpace : " << freeSpace << std::endl;
        // std::cout << "benefit : " << benefit << std::endl;
        // std::cout << "garbage : " << garbage_percent_approx << std::endl;
        // printf("============================================\n");
      } else if(zc_scheme==CBZC6){


        auto now = std::chrono::system_clock::now();
        uint64_t timestamp_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                                    now.time_since_epoch()).count();
        uint64_t total_age=0;
        

        // for(uint64_t i = 0; i<(zone.max_capacity>>12);i++){
        //   if(zone.i_bitmap[i]==0){
        //     total_age  += timestamp_ms-zone.v_bitmap[i];
        //     (void)(timestamp_ms);
        //     continue;
        //   }
        //   total_age  += (zone.i_bitmap[i]-zone.v_bitmap[i]);
        // }
        // total_age>>=30;
      uint64_t extent_n=0;

      // for (auto& ext : snapshot.extents_) {
      //   // if (migrate_zones_start.find(ext.zone_start) != migrate_zones_start.end()) {
      //   if(zone.start==ext.zone_start){
      //     // total_age+=timestamp_ms-ext.create_time;
      //     extent_n++;
      //   }
      // }
      // if(extent_n==0){
      //   continue;
      // }
      // extent_n=0;

      for(auto& file : snapshot.zone_files_){
        if(file.extents.size()){
          // if (migrate_zones_start.find(file.extents[0].start) != migrate_zones_start.end()) {
          if(zone.start==file.extents[0].start){
            total_age+=timestamp_ms-file.extents[0].create_time;
            extent_n++;
          }
        }
      }


      if(extent_n){
        total_age/=extent_n;
      }else{
        total_age=0;
        // continue;
      }


        // uint64_t cost = (100 - ((garbage_percent_approx*garbage_percent_approx)/100) ) * 2;
        // uint64_t benefit =  ((garbage_percent_approx*garbage_percent_approx)/100)  * total_age;

        uint64_t cost = (100 - garbage_percent_approx) * 2;
        uint64_t benefit = (garbage_percent_approx * total_age);


        if (cost != 0) {
          double cost_benefit_score = benefit / cost;
          victim_candidate.push_back(
              {cost_benefit_score, zone.start, garbage_percent_approx, 0.0});
          //     printf("garbage_percent_approx %lu total_age %lu cost %lu benefit %lu cost_benefit_score %f\n",
          // garbage_percent_approx,total_age,cost,benefit,cost_benefit_score);
        }
        // else{
        //   victim_candidate.push_back(
        //       {DBL_MAX, zone.start, garbage_percent_approx, 0.0});
        // }

      }else {
        // printf("CBZC3!!");
        // std::cout << "zc_scheme: " << zc_scheme << std::endl;
        uint64_t zone_start = zone.start;

        double average_lifetime = 0;

        if (zone_lifetime_map_.find(zone_start) != zone_lifetime_map_.end()) {
          auto& entry = zone_lifetime_map_[zone_start];

          double total_lifetime = entry.total_lifetime;  // 존의 lifetime 총합
          int file_count = entry.file_count;             // 존에 포함된 파일 수

          average_lifetime = total_lifetime / file_count;

          // std::map<double, int> value_counts;
          // for (const auto& value : lifetime_values) {
          //   value_counts[value]++;
          // }

          /* benefit = sigma^ZLV * (1-simga)^(free space)
            0<ZLV<1, 0<free space<1
            0<sigma<1 */
          // uint64_t ZoneLifetimeValue = 0;

          // uint64_t min_variance = 20;
          // uint64_t max_variance = 700;
          // double sigma_min = 0.5;
          // double sigma_max = 2.0;

          // double variance_weight =
          //     (static_cast<double>(cur_variance) - min_variance) /
          //     (max_variance - min_variance);
          // double sigma = sigma_min + (sigma_max - sigma_min) *
          // variance_weight; double weighted_age = pow(ZoneLifetimeValue,
          // sigma);
          // uint64_t u = 100 * zone.used_capacity / zone.max_capacity;
          // uint64_t freeSpace = 100 * (zone.max_capacity -
          // zone.used_capacity)
          // /
          //                      zone.max_capacity;
          // uint64_t cost = 100 + u;
          // uint64_t benefit = freeSpace * ZoneLifetimeValue;
          // uint64_t benefit = freeSpace * weighted_age;

          // double sigma = cur_variance;

          double u = 2 * (static_cast<double>(zone.used_capacity) /
                          static_cast<double>(zone.max_capacity));

          double freeSpace =
              static_cast<double>(zone.max_capacity - zone.used_capacity) /
              static_cast<double>(zone.max_capacity);

          // double weighted_age = pow(sigma, average_lifetime);
          // double weighted_age = pow(average_lifetime * 100, sigma);
          double weighted_age = average_lifetime * 100;
          // double weighted_freeSpace = pow(1-sigma, freeSpace);
          double weighted_freeSpace = freeSpace * 100;

          double cost = 2 * u * 100;
          if (weighted_age == 0) {
            weighted_age = 0.1;
          }
          double benefit = weighted_freeSpace * weighted_age;

          // double cost = 2 * (static_cast<double>(zone.used_capacity) /
          //                    static_cast<double>(zone.max_capacity));
          // double benefit = (static_cast<double>(zone.max_capacity) -
          //                   static_cast<double>(zone.used_capacity)) *
          //                  average_lifetime;  // free space * lifetime

          // std::cout << "cost : " << cost << std::endl;
          // std::cout << "freespace : " << freeSpace << std::endl;
          // std::cout << "weighted freespace : " << weighted_freeSpace <<
          // std::endl; std::cout << "age : " << average_lifetime <<
          // std::endl; std::cout << "weighted age : " << weighted_age <<
          // std::endl; std::cout << "sigma : " << sigma << std::endl;
          // std::cout << "benefit : " << benefit << std::endl;

          if (cost != 0) {
            double cost_benefit_score =
                static_cast<double>(benefit) / static_cast<double>(cost);
            victim_candidate.push_back({cost_benefit_score, zone.start,
                                        garbage_percent_approx,
                                        average_lifetime});

            // std::cout << "cost-benefit score: " << cost_benefit_score
            //           << ", zone start: " << zone_start
            //           << ", Garbage Percentage: " << garbage_percent_approx
            //           << "%" << std::endl;

            // std::cout << "  Lifetime Values(" << file_count << "): [";
            // bool first = true;
            // for (const auto& [value, count] : value_counts) {
            //   if (!first) {
            //     std::cout << ", ";
            //   }
            //   first = false;
            //   std::cout << value * 100 << " : " << count;
            // }
            // std::cout << "]" << ",freespace: " << weighted_freeSpace
            //           << ",age: " << weighted_age << std::endl;
          }
        }
      }
    } else {  // 유효 데이터가 없는 경우
      all_inval_zone_n++;
      // std::cout << "all_inal_zone..." << std::endl;
    }
  }

  if (zc_scheme == GREEDY) {
    // GREEDY에서는 garbage_percent_approx 를 기준으로 내림차순 정렬
    std::sort(victim_candidate.begin(), victim_candidate.end(),
              [](const ZoneInfo& a, const ZoneInfo& b) {
                return a.garbage_percent_approx > b.garbage_percent_approx;
              });
  } else {
    // CBZC에서는 cost_benefit_score (score)를 기준으로 내림차순 정렬
    std::sort(
        victim_candidate.begin(), victim_candidate.end(),
        [](const ZoneInfo& a, const ZoneInfo& b) { return a.score > b.score; });
  }

  // std::cout
  //     <<
  //     "-------------------------------------------------------------------"
  //     << std::endl;
  // for (const auto& candidate : victim_candidate) {
  //   std::cout << "cost-benefit score: " << candidate.score
  //             << ", zone start: " << candidate.zone_start
  //             << ", Garbage Percentage: " <<
  //             candidate.garbage_percent_approx
  //             << "%" << std::endl;
  // }

  if (!victim_candidate.empty()) {
    for (const auto& candidate : victim_candidate) {
      migrate_zones_start.emplace(candidate.zone_start);
      ZLV = candidate.ZoneLifetimeValue;
      // std::cout << "ZLV: " << ZLV << std::endl;

      // std::cout << "[Picked] cost-benefit score: " << candidate.score
      //           << ", zone start: " << candidate.zone_start
      //           << ", Garbage Percentage: " <<
      //           candidate.garbage_percent_approx
      //           << "%" << std::endl;

      break;
    }
  }

  // uint64_t threshold = 0;
  // uint64_t reclaimed_zone_n = 1;

  // // 청소할 존 수 계산
  // reclaimed_zone_n = reclaimed_zone_n > victim_candidate.size()
  //                        ? victim_candidate.size()
  //                        : reclaimed_zone_n;

  // // 청소 대상 존 선택
  // for (size_t i = 0;
  //      (i < reclaimed_zone_n && migrate_zones_start.size() <
  //      reclaimed_zone_n); i++) {
  //   migrate_zones_start.emplace(victim_candidate[i].zone_start);
  // }

  std::vector<ZoneExtentSnapshot*> migrate_exts;
  for (auto& ext : snapshot.extents_) {
    if (migrate_zones_start.find(ext.zone_start) != migrate_zones_start.end()) {
      should_be_copied += ext.length;
      migrate_exts.push_back(&ext);
    }
  }

  if (migrate_zones_start.size() > 0) {
    IOStatus s;
    Info(logger_, "Garbage collecting %d extents \n", (int)migrate_exts.size());
    // printf("Garbage collecting %d extents \n", (int)migrate_exts.size());
    clock_gettime(CLOCK_MONOTONIC, &start_timespec);
    s = MigrateExtents(migrate_exts);
    clock_gettime(CLOCK_MONOTONIC, &end_timespec);
    // printf("GetGCBytesWritten %lu\n",zbd_->GetGCBytesWritten());
    if (!s.ok()) {
      printf("Garbage collection failed\n");
      Error(logger_, "Garbage collection failed");
    }
    // 종료 시간 기록
    int end = GetMountTime();
    if (should_be_copied > 0) {
      long elapsed_ns_timespec =
          (end_timespec.tv_sec - start_timespec.tv_sec) * 1000000000 +
          (end_timespec.tv_nsec - start_timespec.tv_nsec);
      // long long microseconds =
      //     std::chrono::duration_cast<std::chrono::microseconds>(elapsed)
      //         .count();
      zbd_->AddCumulativeIOBlocking(elapsed_ns_timespec);
      zbd_->AddZCTimeLapse(start, end, (elapsed_ns_timespec / 1000),
                           migrate_zones_start.size(), should_be_copied, forced,
                           ZLV);
    }
    zc_triggerd_count_.fetch_add(1);
  }
}

void ZenFS::GCWorker() {
  while (run_gc_worker_) {
    // free_percent_ = zbd_->CalculateFreePercent();
    // if (free_percent_ > 20) {
    //   usleep(100 * 1000);
    // }
    // zbd_->SetZCRunning(false);
    // // std::cout << "GCWorker : free_percent_ : " << free_percent_ << "\n";
    // while (zbd_->CalculateFreePercent() < zbd_->zc_) {
    //   zbd_->SetZCRunning(true);
    //   ZoneCleaning(true);
    // }

    // zbd_->SetZCRunning(false);

    // free_percent_ = zbd_->CalculateFreePercent();
    // if (!zbd_->ShouldZCByEmptyZoneN()) {
    //   usleep(100 * 1000);
    //   continue;
    // }
    bool shoudl_zc = zbd_->zc_ > zbd_->CalculateFreePercent() ||
                     zbd_->ShouldZCByEmptyZoneN();
    if (!shoudl_zc) {
      usleep(100 * 1000);
      continue;
    }

    zbd_->SetZCRunning(false);
    // std::cout << "GCWorker : free_percent_ : " << free_percent_ << "\n";
    int try_n = 0;
    while (zbd_->ShouldZCByEmptyZoneN()) {
      zbd_->SetZCRunning(true);
      ZoneCleaning(true);
      try_n++;
      if (try_n > 8) {
        break;
      }
    }
    try_n = 0;
    while (zbd_->CalculateFreePercent() < zbd_->zc_) {
      zbd_->SetZCRunning(true);
      ZoneCleaning(true);
      try_n++;
      if (try_n > 8) {
        break;
      }
    }

    zbd_->SetZCRunning(false);

    //   usleep(1000 * 1000 * 10);
    //   uint64_t non_free = zbd_->GetUsedSpace() +
    //   zbd_->GetReclaimableSpace();
    //   uint64_t free = zbd_->GetFreeSpace();
    //   uint64_t free_percent = (100 * free) / (free + non_free);
    //   ZenFSSnapshot snapshot;
    //   ZenFSSnapshotOptions options;
    //   if (free_percent > GC_START_LEVEL) continue;
    //   options.zone_ = 1; options.zone_file_ = 1; options.log_garbage_ =
    //   1;
    //   GetZenFSSnapshot(snapshot, options);
    //   uint64_t threshold = (100 - GC_SLOPE * (GC_START_LEVEL -
    //   free_percent));
    //   std::set<uint64_t> migrate_zones_start; 
    //   for (const auto& zone : snapshot.zones_) {  
    //     if (zone.capacity == 0) {  
    //       uint64_t garbage_percent_approx =  
    //           100 - 100 * zone.used_capacity / zone.max_capacity;
    //       if (garbage_percent_approx >
    //               threshold &&  
    //           garbage_percent_approx < 100) {  
    //         migrate_zones_start.emplace(
    //             zone.start); 
    //       }
    //     }
    //   }
    //   std::vector<ZoneExtentSnapshot*> migrate_exts;
    //   for (auto& ext : snapshot.extents_) {
    //     if (migrate_zones_start.find(ext.zone_start) !=
    //         migrate_zones_start.end()) {
    //       migrate_exts.push_back(&ext);
    //     }
    //   }
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

  // auto cur_zc_time = std::chrono::system_clock::now();
  struct timespec cur_zc_ts;
  int cur_fops_sequence = zbd_->file_operation_sequence_.load();
  clock_gettime(CLOCK_MONOTONIC, &cur_zc_ts);

  if (!s.ok()) {
    return s;
  }
  // uint64_t sum_diff_ns = 0;
  // size_t deleted_after_copy_extents_n = 0;
  // size_t deleted_after_copy_size = 0;
  for (size_t i = 0; i < new_extents.size(); ++i) {
    ZoneExtent* old_ext = old_extents[i];
    if (old_ext->start_ != new_extents[i]->start_) {

      if (old_ext->is_zc_copied_ == false) {
        new_extents[i]->zc_copied_sequence_= cur_fops_sequence;
        new_extents[i]->zc_copied_ts_ = cur_zc_ts;
      }else{
        new_extents[i]->zc_copied_sequence_ = old_ext->zc_copied_sequence_;
        new_extents[i]->zc_copied_ts_ = old_ext->zc_copied_ts_;
      }
      new_extents[i]->is_zc_copied_ = true;
      // if (old_ext->is_zc_copied_ == true) {
      //   long diff_ns =
      //       (cur_zc_ts.tv_sec - old_ext->zc_copied_ts_.tv_sec) * 1000000000L +
      //       (cur_zc_ts.tv_nsec - old_ext->zc_copied_ts_.tv_nsec);
      //   sum_diff_ns += diff_ns;
      //   deleted_after_copy_extents_n++;
      //   deleted_after_copy_size += old_ext->length_;
      // }else{
      //   // new_extents[i]->zc_copied_sequence_= cur_fops_sequence;

      // }
      
      // new_extents[i]->zc_copied_sequence_= cur_fops_sequence;

      if(old_ext->ZC_COPIED_STATE==NO_COPIED){
        new_extents[i]->ZC_COPIED_STATE=ZC_COPIED;
      }
      // else if(old_ext->ZC_COPIED_STATE == ZC_COPIED){
      //   new_extents[i]->ZC_COPIED_STATE=COPIED_COPIED;
      // }
      // else {
      //   // new_extents[i]->is_zc_copied_ = true;
      //   // new_extents[i]->zc_copied_ts_ = cur_zc_ts;
      // }

      if (old_ext->zone_->this_zone_motivation_check_) {
        std::lock_guard<std::mutex> lg(
            old_ext->zone_->motivation_lifetime_diffs_lock_);

        old_ext->zone_->motivation_lifetime_diffs.push_back(
            {zoneFile->created_time_, cur_zc_ts, true});
      }
      if (new_extents[i]->zone_->this_zone_motivation_check_) {
        zoneFile->created_time_ = cur_zc_ts;
      }

      old_ext->zone_->used_capacity_ -= old_ext->length_;
    }
    delete old_ext;
  }

  // uint64_t sum_diff_us = sum_diff_ns / 1000;
  // if (deleted_after_copy_extents_n != 0) {
  //   zbd_->total_deletion_after_copy_time_.fetch_add(
  //       (sum_diff_us / deleted_after_copy_extents_n));
  //   zbd_->total_deletion_after_copy_n_.fetch_add(1);
  //       zbd_->total_deletion_after_copy_size_.fetch_add(deleted_after_copy_size);
  //   uint64_t actual_cost_benefit_score =
  //       (deleted_after_copy_size >> 20) *
  //       ((sum_diff_us / 1000) / deleted_after_copy_extents_n);
  //   zbd_->actual_cost_benefit_score_.fetch_add(actual_cost_benefit_score);
  // }

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
inline bool ends_with(std::string const& value, std::string const& ending) {
  if (ending.size() > value.size()) return false;
  return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

/* Must hold files_mtx_ */
IOStatus ZenFS::DeleteFileNoLock(std::string fname, const IOOptions& options,
                                 IODebugContext* dbg) {
  std::shared_ptr<ZoneFile> zoneFile(nullptr);
  IOStatus s;

  fname = FormatPathLexically(fname);
  zoneFile = GetFileNoLock(fname);
  if (zoneFile != nullptr) {
    // uint64_t file_size_bytes = zoneFile->GetFileSize();
    // uint64_t file_size_mb = file_size_bytes >> 20;

    // if (file_size_mb <= 32)
    //   file_size_dist[0]++;
    // else if (file_size_mb <= 63)
    //   file_size_dist[1]++;
    // else if (file_size_mb <= 128)
    //   file_size_dist[2]++;
    // else if (file_size_mb <= 256)
    //   file_size_dist[3]++;
    // else
    //   file_size_dist[4]++;
    if (ends_with(fname, ".log")) {
      if(zoneFile->GetFileSize()>>20 < 1024){
        zbd_->file_size_distribution_[zoneFile->GetFileSize()>>20].fetch_add(1);
      }else{
        printf("DeleteFileNoLock fsize %lu\n",zoneFile->GetFileSize()>>20);
      }

    }
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
  // if (s.ok()) s = zbd_->ResetUnusedIOZones();
  if (s.ok()) {
    s = zbd_->RuntimeReset();
  }
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
    zoneFile->AddLinkName(fname);  // 파일 링크 추가

    /* RocksDB does not set the right io type(!)*/
    zoneFile->is_sst_ = ends_with(fname, ".sst");
    // wal
    if (ends_with(fname, ".log")) {
      zoneFile->SetIOType(IOType::kWAL);
      zoneFile->is_wal_ = true;
      zoneFile->SetSparse(!file_opts.use_direct_writes);
      int seq = zbd_->file_operation_sequence_.fetch_add(1);
      zbd_->latest_file_operation_sequence_[SeqWAL] = seq;
      {
        
        // if(coldest_type_set_== true){
        //   // todo
        //   if(coldest_type_!=SeqWAL){
        //     zbd_->CBSC_mispredict_stats_[coldest_type_].fetch_add(1);
        //   }
        //   zbd_->CBSC_total_predict_stats_[coldest_type_].fetch_add(1);
        //   coldest_type_set_=false;
        // }
        
      }
      if(zbd_->coldest_type_set_== true){
        // todo
        std::lock_guard<std::mutex> lg(zbd_->coldest_type_lock_);
        bool ok = true;
        zbd_->check_coldest_[SeqWAL]=true;
        
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
          if(zbd_->coldest_type_!=SeqWAL){
            zbd_->CBSC_mispredict_stats_[zbd_->coldest_type_].fetch_add(1);
          }
          zbd_->CBSC_total_predict_stats_[zbd_->coldest_type_].fetch_add(1);
          zbd_->coldest_type_set_=false;
        }


        else{
          if(zbd_->coldest_type_==SeqWAL){
            zbd_->CBSC_mispredict_stats_[zbd_->coldest_type_].fetch_add(1);
            zbd_->CBSC_total_predict_stats_[zbd_->coldest_type_].fetch_add(1);
            zbd_->coldest_type_set_=false;
          }
        }
      }
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

  // if (resetIOZones) s = zbd_->ResetUnusedIOZones();
  if (resetIOZones) {
    s = zbd_->RuntimeReset();
  }

  return s;
}

IOStatus ZenFS::DeleteFile(const std::string& fname, const IOOptions& options,
                           IODebugContext* dbg) {
  IOStatus s;

  Debug(logger_, "DeleteFile: %s \n", fname.c_str());

  files_mtx_.lock();
  if(ends_with(fname, ".log")){
    int seq= zbd_->file_operation_sequence_.fetch_add(1);
    zbd_->latest_file_operation_sequence_[SeqWAL] = seq;
    {
      // std::lock_guard<std::mutex> lg(coldest_type_lock_);
      if(zbd_->coldest_type_set_== true){
        // todo
        std::lock_guard<std::mutex> lg(zbd_->coldest_type_lock_);
        bool ok = true;
        zbd_->check_coldest_[SeqWAL]=true;
        
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
          if(zbd_->coldest_type_!=SeqWAL){
            zbd_->CBSC_mispredict_stats_[zbd_->coldest_type_].fetch_add(1);
          }
          zbd_->CBSC_total_predict_stats_[zbd_->coldest_type_].fetch_add(1);
          zbd_->coldest_type_set_=false;
        }       
        else{
          if(zbd_->coldest_type_==SeqWAL){
            zbd_->CBSC_mispredict_stats_[zbd_->coldest_type_].fetch_add(1);
            zbd_->CBSC_total_predict_stats_[zbd_->coldest_type_].fetch_add(1);
            zbd_->coldest_type_set_=false;
          }
        }
      }

      
    }
  }
  s = DeleteFileNoLock(fname, options, dbg);

  files_mtx_.unlock();
  // if (s.ok()) s = zbd_->ResetUnusedIOZones();
  if (s.ok()) {
    s = zbd_->RuntimeReset();
  }
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

// void ZenFS::SetResetScheme(uint32_t r, bool f, uint64_t T) {
//   std::cout << "ZenFS::SetResetScheme: r = " << r << ", f = " << f
//             << ", T = " << T << std::endl;
//   zbd_->SetResetScheme(r, f, T);
//   run_bg_reset_worker_ = true;
//   if (gc_worker_ != nullptr) {
//     if (bg_reset_worker_ == nullptr) {
//       bg_reset_worker_.reset(
//           new std::thread(&ZenFS::BackgroundStatTimeLapse, this));
//     }
//   }
// }

void ZenFS::SetResetScheme(uint32_t r, uint32_t partial_reset_scheme,
                           uint64_t T, uint64_t zc, uint64_t until,
                           uint64_t allocation_scheme, uint64_t zc_scheme,
                           double alpha_value, double sigma_value,
                           uint64_t finish_scheme, uint64_t predict_cnt,
                           std::vector<uint64_t>& other_options) {
  std::cout << "ZenFS::SetResetScheme: r = " << r << ", T = " << T
            << ", allocation_schme = " << allocation_scheme
            << ", zc_scheme = " << zc_scheme
            << ", finish_scheme = " << finish_scheme
            << ", predict_cnt = " << predict_cnt << std::endl;
  zbd_->SetResetScheme(r, partial_reset_scheme, T, zc, until, allocation_scheme,
                       zc_scheme, alpha_value, sigma_value, finish_scheme,
                       predict_cnt, other_options);
  run_bg_reset_worker_ = true;
  if (gc_worker_ != nullptr) {
    if (bg_reset_worker_ == nullptr) {
      bg_reset_worker_.reset(
          new std::thread(&ZenFS::BackgroundStatTimeLapse, this));
    }
  }

  // if (partial_reset_scheme == PARTIAL_RESET_AT_BACKGROUND ||
  //     partial_reset_scheme == PARTIAL_RESET_BACKGROUND_T_WITH_ZONE_RESET) {
  //   printf("PARTIAL RESET AT BG\n");
  //   if (bg_partial_reset_worker_ == nullptr) {
  //     run_bg_partial_reset_worker_ = true;
  //     bg_partial_reset_worker_.reset(
  //         new std::thread(&ZenFS::PartialResetWorker, this, T));
  //   }
  //   return;
  // }
  // switch (partial_reset_scheme) {
  //   case RUNTIME_ZONE_RESET_DISABLED:
  //     printf("RUNTIME_ZONE_RESET_DISABLED\n");
  //     return;
  //   case RUNTIME_ZONE_RESET_ONLY:
  //     printf("RUNTIME_ZONE_RESET_ONLY\n");
  //     return;
  //   case PARTIAL_RESET_WITH_ZONE_RESET:
  //     printf("PARTIAL_RESET_WITH_ZONE_RESET\n");
  //     return;
  //   case PARTIAL_RESET_ONLY:
  //     printf("PARTIAL_RESET_ONLY\n");
  //     return;
  //   case PROACTIVE_ZONECLEANING:
  //     printf("PROACTIVE_ZONECLEANING\n");
  //     return;
  //   default:
  //     printf("UNKNOWN SCHEME!\n");
  //     // exit(-1);
  //     break;
  // }
}

double ZenFS::GetMaxInvalidateCompactionScore(
    std::vector<uint64_t>& file_candidates, uint64_t* candidate_size) {
  return zbd_->GetMaxInvalidateCompactionScore(file_candidates, candidate_size,
                                               false);
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
  // if (s.ok()) s = zbd_->ResetUnusedIOZones();
  if (s.ok()) {
    s = zbd_->RuntimeReset();
  }
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

  if (!readonly) {
    Info(logger_, "Resetting unused IO Zones..");
    IOStatus status = zbd_->ResetUnusedIOZones();
    if (!status.ok()) return status;
    Info(logger_, "  Done");

    if (superblock_->IsGCEnabled()) {
      Info(logger_, "Starting garbage collection worker");
      run_gc_worker_ = true;
      gc_worker_.reset(new std::thread(&ZenFS::GCWorker, this));
      run_bg_reset_worker_ = true;
      if (bg_reset_worker_ == nullptr) {
        printf("starting bg_reset_worker");
        bg_reset_worker_.reset(
            new std::thread(&ZenFS::BackgroundStatTimeLapse, this));
      }
    }
  }

  LogFiles();

  return Status::OK();
}

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
    // std::unordered_map<uint64_t, std::vector<ZoneFileSnapshot*>>
    // zone_to_files;

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

IOStatus ZenFS::MigrateExtents(
    const std::vector<ZoneExtentSnapshot*>& extents) {
  // struct timespec start1, end1;
  // struct timespec start2, end2;
  // struct timespec start3, end3;
  // clock_gettime(CLOCK_MONOTONIC, &start1);
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
  // clock_gettime(CLOCK_MONOTONIC, &end1);
  // long elapsed1 = (end1.tv_sec - start1.tv_sec) * 1000000000 +
  //                 (end1.tv_nsec - start1.tv_nsec);
  // zbd_->AddCumulative_1(elapsed1);

  // clock_gettime(CLOCK_MONOTONIC, &start2);
  // 파일 익스텐트 마이그레이션 및 존 재설정
  for (const auto& it : file_extents) {
    s = MigrateFileExtents(it.first, it.second);
    if (!s.ok()) break;  // 마이 그레이션 실패하면 break
  }
  // clock_gettime(CLOCK_MONOTONIC, &end2);
  // long elapsed2 = (end2.tv_sec - start2.tv_sec) * 1000000000 +
  //                 (end2.tv_nsec - start2.tv_nsec);
  // zbd_->AddCumulative_2(elapsed2);

  // clock_gettime(CLOCK_MONOTONIC, &start3);
  s = zbd_->ResetUnusedIOZones();
  // clock_gettime(CLOCK_MONOTONIC, &end3);
  // long elapsed3 = (end3.tv_sec - start3.tv_sec) * 1000000000 +
  //                 (end3.tv_nsec - start3.tv_nsec);
  // zbd_->AddCumulative_3(elapsed3);
  return s;
}

IOStatus ZenFS::MigrateFileExtents(
    const std::string& fname,
    const std::vector<ZoneExtentSnapshot*>& migrate_exts) {
  IOStatus s = IOStatus::OK();
  uint64_t copied = 0;
 
  Info(logger_, "MigrateFileExtents, fname: %s, extent count: %lu",
       fname.data(), migrate_exts.size());

  // The file may be deleted by other threads, better double check.

  // struct timespec start1, end1;
  // struct timespec start2, end2;
  // struct timespec start3, end3;
  // struct timespec start4, end4;
  // struct timespec start5, end5;
  // struct timespec start6, end6;
  // struct timespec start7, end7;
  // clock_gettime(CLOCK_MONOTONIC, &start4);
  auto zfile = GetFile(fname);
  if (zfile == nullptr) {
    return IOStatus::OK();
  }

  // Don't migrate open for write files and prevent write reopens while we
  // migrate
  if (!zfile->TryAcquireWRLock()) {
    return IOStatus::OK();
  }

  std::vector<ZoneExtent*> new_extent_list;
  std::vector<ZoneExtent*> extents = zfile->GetExtents();
  for (const auto* ext : extents) {
    new_extent_list.push_back(
        new ZoneExtent(ext->start_, ext->length_, ext->zone_));
  }
  // clock_gettime(CLOCK_MONOTONIC, &end4);
  // long elapsed4 = (end4.tv_sec - start4.tv_sec) * 1000000000 +
  //                 (end4.tv_nsec - start4.tv_nsec);
  // zbd_->AddCumulative_4(elapsed4);

  // Modify the new extent list
  // clock_gettime(CLOCK_MONOTONIC, &start5);
  for (ZoneExtent* ext : new_extent_list) {
    // Check if current extent need to be migrated
    // 유효 데이터(즉, 마이그레이션할 필요가 있는 익스텐트)를 확인
    // clock_gettime(CLOCK_MONOTONIC, &start1);
    auto it = std::find_if(migrate_exts.begin(), migrate_exts.end(),
                           [&](const ZoneExtentSnapshot* ext_snapshot) {
                             return ext_snapshot->start == ext->start_ &&
                                    ext_snapshot->length == ext->length_;
                           });
    // clock_gettime(CLOCK_MONOTONIC, &end1);
    // long elapsed1 = (end1.tv_sec - start1.tv_sec) * 1000000000 +
    // (end1.tv_nsec - start1.tv_nsec);
    // zbd_->AddCumulative_1(elapsed1);

    if (it == migrate_exts.end()) {
      Info(logger_, "Migrate extent not found, ext_start: %lu", ext->start_);
      continue;
    }

    Zone* target_zone = nullptr;

    // Allocate a new migration zone.
    // s = zbd_->TakeMigrateZone(&target_zone, zfile->GetWriteLifeTimeHint(),
    //                           ext->length_);
    // clock_gettime(CLOCK_MONOTONIC, &start7);
    s = zbd_->TakeMigrateZone(zfile->smallest_, zfile->largest_, zfile->level_,
                              &target_zone, zfile->GetWriteLifeTimeHint(),
                              zfile->predicted_size_, ext->length_,
                              &run_gc_worker_, zfile->IsSST());
    // clock_gettime(CLOCK_MONOTONIC, &end7);
    // long elapsed7 = (end7.tv_sec - start7.tv_sec) * 1000000000 +
    // (end7.tv_nsec - start7.tv_nsec);
    // zbd_->AddCumulative_7(elapsed7);
    if (!run_gc_worker_) {
      printf("MigrateFileExtents - !run_gc_worker\n");
      return IOStatus::OK();
    }
    if (!s.ok()) {
      printf("MigrateFileExtents - !s.ok\n");
      continue;
    }
    // clock_gettime(CLOCK_MONOTONIC, &start2);
    if (target_zone == nullptr) {
      zbd_->ReleaseMigrateZone(target_zone);
      printf("MigrateFileExtents - Migrate Zone Acquire Failed, Ignore Task\n");
      Info(logger_, "Migrate Zone Acquire Failed, Ignore Task.");
      continue;
    }
    // clock_gettime(CLOCK_MONOTONIC, &end2);
    // long elapsed2 = (end2.tv_sec - start2.tv_sec) * 1000000000 +
    // (end2.tv_nsec - start2.tv_nsec);
    // zbd_->AddCumulative_2(elapsed2);

    // clock_gettime(CLOCK_MONOTONIC, &start3);
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
    // clock_gettime(CLOCK_MONOTONIC, &end3);
    // long elapsed3 = (end3.tv_sec - start3.tv_sec) * 1000000000 +
    // (end3.tv_nsec - start3.tv_nsec);
    // zbd_->AddCumulative_3(elapsed3);

    // clock_gettime(CLOCK_MONOTONIC, &start4);
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
    // clock_gettime(CLOCK_MONOTONIC, &end4);
    // long elapsed4 = (end4.tv_sec - start4.tv_sec) * 1000000000 +
    // (end4.tv_nsec - start4.tv_nsec);
    // zbd_->AddCumulative_4(elapsed4);
  }
  // clock_gettime(CLOCK_MONOTONIC, &end5);
  // long elapsed5 = (end5.tv_sec - start5.tv_sec) * 1000000000 +
  // (end5.tv_nsec - start5.tv_nsec);
  // zbd_->AddCumulative_5(elapsed5);

  // clock_gettime(CLOCK_MONOTONIC, &start6);
  zbd_->AddGCBytesWritten(copied);
  SyncFileExtents(zfile.get(), new_extent_list);
  zfile->ReleaseWRLock();
  // clock_gettime(CLOCK_MONOTONIC, &end6);
  // long elapsed6 = (end6.tv_sec - start6.tv_sec) * 1000000000 +
  // (end6.tv_nsec - start6.tv_nsec);
  // zbd_->AddCumulative_6(elapsed6);
  // printf(
  //     "MigrateFileExtents Finished, fname: %s, extent count: %lu, copied :
  //     "
  //     "%lu\n",
  //     fname.data(), migrate_exts.size(), (copied) >> 20);
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
