// Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved.
// Copyright (c) 2019-present, Western Digital Corporation
//  This source code is licensed under both the GPLv2 (found in the
//  COPYING file in the root directory) and Apache 2.0 License
//  (found in the LICENSE.Apache file in the root directory).
/*  I/O 관련 기능을 구현하고 정의하는 파일입니다. 데이터 입출력 처리, 버퍼 관리,
 * I/O 작업의 스케줄링 등 ZenFS의 핵심적인 I/O 작업을 처리합니다*/
#if !defined(ROCKSDB_LITE) && !defined(OS_WIN)

#include "io_zenfs.h"

#include <assert.h>
#include <errno.h>
#include <fcntl.h>
#include <libzbd/zbd.h>
#include <linux/blkzoned.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ioctl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <chrono>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "rocksdb/env.h"
#include "rocksdb/file_system.h"
#include "rocksdb/status.h"
#include "util/coding.h"

#include "zbd_zenfs.h"  


namespace ROCKSDB_NAMESPACE {

ZoneExtent::ZoneExtent(uint64_t start, uint64_t length, Zone* zone)
    : start_(start), length_(length), zone_(zone), is_zc_copied_(false),ZC_COPIED_STATE(NO_COPIED) {
    
        auto now = std::chrono::system_clock::now();
        create_time_ = std::chrono::duration_cast<std::chrono::milliseconds>(
                                    now.time_since_epoch()).count();
    }

ZoneExtent::ZoneExtent(uint64_t start, uint64_t length, Zone* zone,
                       std::string fname, ZoneFile* zfile)
    : start_(start),
      length_(length),
      zone_(zone),
      is_invalid_(false),
      fname_(fname),
      header_size_(0),
      zfile_(zfile),
      is_zc_copied_(false),
      ZC_COPIED_STATE(NO_COPIED) {
  zc_copied_ts_.tv_sec = 0;
  zc_copied_ts_.tv_nsec = 0;

        auto now = std::chrono::system_clock::now();
        create_time_ = std::chrono::duration_cast<std::chrono::milliseconds>(
                                    now.time_since_epoch()).count();
  if (zone == nullptr) {
    return;
  }

  // zone->PushExtent(this);

  // uint64_t align = (length_ + header_size_) % block_sz;
  // if (align) {
  //   pad_size_ = block_sz - align;
  // }
}

Status ZoneExtent::DecodeFrom(Slice* input) {
  if (input->size() != (sizeof(start_) + sizeof(length_)))
    return Status::Corruption("ZoneExtent", "Error: length missmatch");

  GetFixed64(input, &start_);
  GetFixed64(input, &length_);
  return Status::OK();
}

void ZoneExtent::EncodeTo(std::string* output) {
  PutFixed64(output, start_);
  PutFixed64(output, length_);
}

void ZoneExtent::EncodeJson(std::ostream& json_stream) {
  json_stream << "{";
  json_stream << "\"start\":" << start_ << ",";
  json_stream << "\"length\":" << length_;
  json_stream << "}";
}

enum ZoneFileTag : uint32_t {
  kFileID = 1,
  kFileNameDeprecated = 2,
  kFileSize = 3,
  kWriteLifeTimeHint = 4,
  kExtent = 5,
  kModificationTime = 6,
  kActiveExtentStart = 7,
  kIsSparse = 8,
  kLinkedFilename = 9,
};

void ZoneFile::EncodeTo(std::string* output, uint32_t extent_start) {
  PutFixed32(output, kFileID);
  PutFixed64(output, file_id_);

  PutFixed32(output, kFileSize);
  PutFixed64(output, file_size_);

  PutFixed32(output, kWriteLifeTimeHint);
  PutFixed32(output, (uint32_t)lifetime_);

  for (uint32_t i = extent_start; i < extents_.size(); i++) {
    std::string extent_str;

    PutFixed32(output, kExtent);
    extents_[i]->EncodeTo(&extent_str);
    PutLengthPrefixedSlice(output, Slice(extent_str));
  }

  PutFixed32(output, kModificationTime);
  PutFixed64(output, (uint64_t)m_time_);

  /* We store the current extent start - if there is a crash
   * we know that this file wrote the data starting from
   * active extent start up to the zone write pointer.
   * We don't need to store the active zone as we can look it up
   * from extent_start_ */
  PutFixed32(output, kActiveExtentStart);
  PutFixed64(output, extent_start_);

  if (is_sparse_) {
    PutFixed32(output, kIsSparse);
  }

  for (uint32_t i = 0; i < linkfiles_.size(); i++) {
    PutFixed32(output, kLinkedFilename);
    PutLengthPrefixedSlice(output, Slice(linkfiles_[i]));
  }
}

void ZoneFile::EncodeJson(std::ostream& json_stream) {
  json_stream << "{";
  json_stream << "\"id\":" << file_id_ << ",";
  json_stream << "\"size\":" << file_size_ << ",";
  json_stream << "\"hint\":" << lifetime_ << ",";
  json_stream << "\"filename\":[";

  bool first_element = true;
  for (const auto& name : GetLinkFiles()) {
    if (first_element) {
      first_element = false;
    } else {
      json_stream << ",";
    }
    json_stream << "\"" << name << "\"";
  }
  json_stream << "],";

  json_stream << "\"extents\":[";

  first_element = true;
  for (ZoneExtent* extent : extents_) {
    if (first_element) {
      first_element = false;
    } else {
      json_stream << ",";
    }
    extent->EncodeJson(json_stream);
  }
  json_stream << "]}";
}

Status ZoneFile::DecodeFrom(Slice* input) {
  uint32_t tag = 0;

  GetFixed32(input, &tag);
  if (tag != kFileID || !GetFixed64(input, &file_id_))
    return Status::Corruption("ZoneFile", "File ID missing");

  while (true) {
    Slice slice;
    ZoneExtent* extent;
    Status s;

    if (!GetFixed32(input, &tag)) break;

    switch (tag) {
      case kFileSize:
        if (!GetFixed64(input, &file_size_))
          return Status::Corruption("ZoneFile", "Missing file size");
        break;
      case kWriteLifeTimeHint:
        uint32_t lt;
        if (!GetFixed32(input, &lt))
          return Status::Corruption("ZoneFile", "Missing life time hint");
        lifetime_ = (Env::WriteLifeTimeHint)lt;
        break;
      case kExtent:
        extent = new ZoneExtent(0, 0, nullptr);
        GetLengthPrefixedSlice(input, &slice);
        s = extent->DecodeFrom(&slice);
        if (!s.ok()) {
          delete extent;
          return s;
        }
        extent->zone_ = zbd_->GetIOZone(extent->start_);
        if (!extent->zone_)
          return Status::Corruption("ZoneFile", "Invalid zone extent");
        extent->zone_->used_capacity_ += extent->length_;
        extents_.push_back(extent);
        break;
      case kModificationTime:
        uint64_t ct;
        if (!GetFixed64(input, &ct))
          return Status::Corruption("ZoneFile", "Missing creation time");
        m_time_ = (time_t)ct;
        break;
      case kActiveExtentStart:
        uint64_t es;
        if (!GetFixed64(input, &es))
          return Status::Corruption("ZoneFile", "Active extent start");
        extent_start_ = es;
        break;
      case kIsSparse:
        is_sparse_ = true;
        break;
      case kLinkedFilename:
        if (!GetLengthPrefixedSlice(input, &slice))
          return Status::Corruption("ZoneFile", "LinkFilename missing");

        if (slice.ToString().length() == 0)
          return Status::Corruption("ZoneFile", "Zero length Linkfilename");

        linkfiles_.push_back(slice.ToString());
        break;
      default:
        return Status::Corruption("ZoneFile", "Unexpected tag");
    }
  }

  MetadataSynced();
  return Status::OK();
}

Status ZoneFile::MergeUpdate(std::shared_ptr<ZoneFile> update, bool replace) {
  if (file_id_ != update->GetID())
    return Status::Corruption("ZoneFile update", "ID missmatch");

  SetFileSize(update->GetFileSize());
  SetWriteLifeTimeHint(update->GetWriteLifeTimeHint());
  SetFileModificationTime(update->GetFileModificationTime());

  if (replace) {
    ClearExtents();
  }

  std::vector<ZoneExtent*> update_extents = update->GetExtents();
  for (long unsigned int i = 0; i < update_extents.size(); i++) {
    ZoneExtent* extent = update_extents[i];
    Zone* zone = extent->zone_;
    zone->used_capacity_ += extent->length_;
    extents_.push_back(new ZoneExtent(extent->start_, extent->length_, zone));
  }
  extent_start_ = update->GetExtentStart();
  is_sparse_ = update->IsSparse();
  MetadataSynced();

  linkfiles_.clear();
  for (const auto& name : update->GetLinkFiles()) linkfiles_.push_back(name);

  return Status::OK();
}

ZoneFile::ZoneFile(ZonedBlockDevice* zbd, uint64_t file_id,
                   MetadataWriter* metadata_writer, FileSystemWrapper* zenfs)
    : zbd_(zbd),
      active_zone_(NULL),
      extent_start_(NO_EXTENT),
      extent_filepos_(0),
      lifetime_(Env::WLTH_NOT_SET),
      io_type_(IOType::kUnknown),
      file_size_(0),
      file_id_(file_id),
      nr_synced_extents_(0),
      m_time_(0),
      metadata_writer_(metadata_writer),
      zenfs_(zenfs) {}  ///

std::string ZoneFile::GetFilename() { return linkfiles_[0]; }
time_t ZoneFile::GetFileModificationTime() { return m_time_; }

uint64_t ZoneFile::GetFileSize() { return file_size_; }
void ZoneFile::SetFileSize(uint64_t sz) { file_size_ = sz; }
void ZoneFile::SetFileModificationTime(time_t mt) { m_time_ = mt; }
void ZoneFile::SetIOType(IOType io_type) { io_type_ = io_type; }

ZoneFile::~ZoneFile() { ClearExtents(); }

void ZoneFile::ClearExtents() {
  uint64_t zc_scheme = zbd_->GetZCScheme();
  // auto cur_deletion_time = std::chrono::system_clock::now();
  // std::chrono::microseconds sum_diff_time(0);
  // uint64_t sum_diff_time_uint64t = 0;
  // size_t deleted_after_copy_extents_n = 0;

  struct timespec cur_deletion_ts;
  clock_gettime(CLOCK_MONOTONIC, &cur_deletion_ts);
  int cur_fops_sequence = zbd_->file_operation_sequence_.load();
  uint64_t sum_diff_ns = 0;
  size_t deleted_after_copy_extents_n = 0;
  size_t deleted_after_copy_size = 0;
  // bool motivation_check_=false;
  size_t deleted_after_copy_sequence_sum=0;

  Zone* motivation_check_zone = nullptr;
  for (auto e = std::begin(extents_); e != std::end(extents_); ++e) {
    Zone* zone = (*e)->zone_;

    assert(zone && zone->used_capacity_ >= (*e)->length_);
    zone->used_capacity_ -= (*e)->length_;
    if (zone->this_zone_motivation_check_) {
      motivation_check_zone = zone;
    }
    if (zc_scheme == CBZC5) {
      zone->recent_inval_time_ = std::chrono::system_clock::now();
    }

    if ((*e)->is_zc_copied_ && ((*e)->ZC_COPIED_STATE==ZC_COPIED)) {
      // (현재삭제시각 - 복사시각) -> 나노초로 계산
      long diff_ns =
          (cur_deletion_ts.tv_sec - (*e)->zc_copied_ts_.tv_sec) * 1000000000L +
          (cur_deletion_ts.tv_nsec - (*e)->zc_copied_ts_.tv_nsec);

      if (diff_ns < 0) {
        diff_ns = 0;
      }
      sum_diff_ns += diff_ns;
      deleted_after_copy_extents_n++;
      deleted_after_copy_size += (*e)->length_;
      deleted_after_copy_sequence_sum+=(cur_fops_sequence-(*e)->zc_copied_sequence_);
    }
    delete *e;
  }
  extents_.clear();

  uint64_t sum_diff_us = sum_diff_ns / 1000;

  if (motivation_check_zone) {
    printf("@@@@@@@@@@@@@@@@@@ reset motivation_lifetime_diffs CLEAREXTENTS\n");
    std::lock_guard<std::mutex> lg(
        motivation_check_zone->motivation_lifetime_diffs_lock_);

    motivation_check_zone->motivation_lifetime_diffs.push_back(
        {created_time_, cur_deletion_ts, false});
  }

  // zbd_->total_deletion_after_copy_time_.fetch_add(sum_diff_us);
  // zbd_->total_deletion_after_copy_n_.fetch_add(deleted_after_copy_extents_n);
  if (deleted_after_copy_extents_n) {
    zbd_->total_deletion_after_copy_time_.fetch_add(
        sum_diff_us / deleted_after_copy_extents_n);
    zbd_->total_deletion_after_copy_n_.fetch_add(1);
    zbd_->total_deletion_after_copy_size_.fetch_add(deleted_after_copy_size);
    uint64_t actual_cost_benefit_score =
        (deleted_after_copy_size >> 20) *
        ((sum_diff_us / 1000) / deleted_after_copy_extents_n);
    zbd_->actual_cost_benefit_score_.fetch_add(actual_cost_benefit_score);
    
    zbd_->total_deletion_after_copy_seq_.fetch_add((deleted_after_copy_sequence_sum*100)/deleted_after_copy_extents_n);
    if(deleted_after_copy_sequence_sum/deleted_after_copy_extents_n>SEQ_DIST_MAX){
      printf("SEQUENCE OVER SEQ_DIST_MAX : %lu",(deleted_after_copy_sequence_sum/deleted_after_copy_extents_n));
    }else{
      zbd_->total_deletion_after_copy_seq_distribution_[(deleted_after_copy_sequence_sum/deleted_after_copy_extents_n)]++;
    }
   
    uint64_t tmp = ((deleted_after_copy_size >> 20) * deleted_after_copy_sequence_sum)/(deleted_after_copy_extents_n);
    // tmp = tmp * 
    zbd_->cost_benefit_score_sum_sequence_mb_.fetch_add(tmp);

  }
}

IOStatus ZoneFile::CloseActiveZone() {
  struct timespec timespec;
  clock_gettime(CLOCK_MONOTONIC, &timespec);
  IOStatus s = IOStatus::OK();
  if (active_zone_) {
    bool full = active_zone_->IsFull();
    s = active_zone_->Close();

    if (active_zone_->this_zone_motivation_check_) {
      printf("@@@@@@@@@@@@@@@@@@ reset motivation_lifetime_diffs ALLOC %s\n",
             linkfiles_[0].c_str());
      if (active_zone_->is_allocated_ == false) {
        active_zone_->is_allocated_ = true;
        active_zone_->allocated_time_ = timespec;
      }
      created_time_ = timespec;
    }

    ReleaseActiveZone();
    if (!s.ok()) {
      return s;
    }
    zbd_->PutOpenIOZoneToken();
    if (full) {
      zbd_->PutActiveIOZoneToken();
    }
  }
  return s;
}

void ZoneFile::AcquireWRLock() {
  open_for_wr_mtx_.lock();
  open_for_wr_ = true;
}

bool ZoneFile::TryAcquireWRLock() {
  if (!open_for_wr_mtx_.try_lock()) return false;
  open_for_wr_ = true;
  return true;
}

void ZoneFile::ReleaseWRLock() {
  assert(open_for_wr_);
  open_for_wr_ = false;
  open_for_wr_mtx_.unlock();
}

bool ZoneFile::IsOpenForWR() { return open_for_wr_; }

IOStatus ZoneFile::CloseWR() {
  IOStatus s;
  /* Mark up the file as being closed */
  extent_start_ = NO_EXTENT;
  s = PersistMetadata();
  if (!s.ok()) return s;
  ReleaseWRLock();
  return CloseActiveZone();
}

IOStatus ZoneFile::PersistMetadata() {
  assert(metadata_writer_ != NULL);
  return metadata_writer_->Persist(this);
}

ZoneExtent* ZoneFile::GetExtent(uint64_t file_offset, uint64_t* dev_offset) {
  for (unsigned int i = 0; i < extents_.size(); i++) {
    if (file_offset < extents_[i]->length_) {
      *dev_offset = extents_[i]->start_ + file_offset;
      return extents_[i];
    } else {
      file_offset -= extents_[i]->length_;
    }
  }
  return NULL;
}

IOStatus ZoneFile::InvalidateCache(uint64_t pos, uint64_t size) {
  ReadLock lck(this);
  uint64_t offset = pos;
  uint64_t left = size;
  IOStatus s = IOStatus::OK();

  if (left == 0) {
    left = GetFileSize();
  }

  while (left) {
    uint64_t dev_offset;
    ZoneExtent* extent = GetExtent(offset, &dev_offset);

    if (!extent) {
      s = IOStatus::IOError("Extent not found while invalidating cache");
      break;
    }

    uint64_t extent_end = extent->start_ + extent->length_;
    uint64_t invalidate_size = std::min(left, extent_end - dev_offset);

    s = zbd_->InvalidateCache(dev_offset, invalidate_size);
    if (!s.ok()) break;

    left -= invalidate_size;
    offset += invalidate_size;
  }

  return s;
}

IOStatus ZoneFile::PositionedRead(uint64_t offset, size_t n, Slice* result,
                                  char* scratch, bool direct) {
  ZenFSMetricsLatencyGuard guard(zbd_->GetMetrics(), ZENFS_READ_LATENCY,
                                 Env::Default());
  zbd_->GetMetrics()->ReportQPS(ZENFS_READ_QPS, 1);

  ReadLock lck(this);

  char* ptr;
  uint64_t r_off;
  size_t r_sz;
  ssize_t r = 0;
  size_t read = 0;
  ZoneExtent* extent;
  uint64_t extent_end;
  IOStatus s;

  if (offset >= file_size_) {
    *result = Slice(scratch, 0);
    return IOStatus::OK();
  }

  r_off = 0;
  extent = GetExtent(offset, &r_off);
  if (!extent) {
    /* read start beyond end of (synced) file data*/
    *result = Slice(scratch, 0);
    return s;
  }
  extent_end = extent->start_ + extent->length_;

  /* Limit read size to end of file */
  if ((offset + n) > file_size_)
    r_sz = file_size_ - offset;
  else
    r_sz = n;

  ptr = scratch;

  while (read != r_sz) {
    size_t pread_sz = r_sz - read;

    if ((pread_sz + r_off) > extent_end) pread_sz = extent_end - r_off;

    /* We may get some unaligned direct reads due to non-aligned extent lengths,
     * so increase read request size to be aligned to next blocksize boundary.
     */
    bool aligned = (pread_sz % zbd_->GetBlockSize() == 0);

    size_t bytes_to_align = 0;
    if (direct && !aligned) {
      bytes_to_align = zbd_->GetBlockSize() - (pread_sz % zbd_->GetBlockSize());
      pread_sz += bytes_to_align;
      aligned = true;
    }

    r = zbd_->Read(ptr, r_off, pread_sz, direct && aligned);
    if (r <= 0) break;

    /* Verify and update the the bytes read count (if read size was incremented,
     * for alignment purposes).
     */
    if ((size_t)r <= pread_sz - bytes_to_align)
      pread_sz = (size_t)r;
    else
      pread_sz -= bytes_to_align;

    ptr += pread_sz;
    read += pread_sz;
    r_off += pread_sz;

    if (read != r_sz && r_off == extent_end) {
      extent = GetExtent(offset + read, &r_off);
      if (!extent) {
        /* read beyond end of (synced) file data */
        break;
      }
      r_off = extent->start_;
      extent_end = extent->start_ + extent->length_;
    }
  }

  if (r < 0) {
    s = IOStatus::IOError("pread error\n");
    read = 0;
  }

  *result = Slice((char*)scratch, read);
  return s;
}

void ZoneFile::PushExtent() {
  uint64_t length;

  assert(file_size_ >= extent_filepos_);

  if (!active_zone_) return;

  length = file_size_ - extent_filepos_;
  if (length == 0) return;

  assert(length <= (active_zone_->wp_ - extent_start_));
  extents_.push_back(new ZoneExtent(extent_start_, length, active_zone_));
  // extents_.push_back(
  //     new ZoneExtent(extent_start_, length, active_zone_, filename, this));
  // if(active_zone_->wp_==)
  active_zone_->used_capacity_ += length;
  extent_start_ = active_zone_->wp_;
  extent_filepos_ = file_size_;
  // if (zbd_->GetZCScheme() == CBZC5) {
  active_zone_->recent_inval_time_ = std::chrono::system_clock::now();
  // }
}

IOStatus ZoneFile::AllocateNewZone(uint64_t min_capacity) {
  Zone* zone = nullptr;
  // IOStatus s = zbd_->AllocateIOZone(lifetime_, io_type_, &zone);
  int try_n = 0;
  // (void)(min_capacity);  // 활성화할땐 지우기
  IOStatus s = zbd_->AllocateIOZone(linkfiles_[0], IsSST(), smallest_, largest_,
                                    level_, lifetime_, io_type_,
                                    predicted_size_, &zone, min_capacity);
  // input_fno_.clear();
  if (zone == nullptr) {
    int start = zenfs_->GetMountTime();
    while (zbd_->CalculateCapacityRemain() > (1 << 20) * 128) {
      // s = zbd_->AllocateIOZone(lifetime_, io_type_, &zone);
      s = zbd_->AllocateIOZone(linkfiles_[0], IsSST(), smallest_, largest_,
                               level_, lifetime_, io_type_, predicted_size_,
                               &zone, min_capacity);
      try_n++;
      if (zone != nullptr) {
        break;
      }
      // zbd_->ResetUnusedIOZones();
      zbd_->RuntimeReset();
    }
    int end = zenfs_->GetMountTime();

    zbd_->AddIOBlockedTimeLapse(start, end);
  }
  if (!s.ok()) {
    // std::cout << "zf::allocatenewzone=>!s.ok(): " << s;
    return s;
  }

  if (!zone) {
    printf("Error :: zone allocation fail\n");
    zbd_->SetZoneAllocationFailed();
    return IOStatus::NoSpace("Zone allocation failure :: AllocateNewZone\n");
  }

  SetActiveZone(zone);
  extent_start_ = active_zone_->wp_;
  extent_filepos_ = file_size_;

  /* Persist metadata so we can recover the active extent using
     the zone write pointer in case there is a crash before syncing */
  return PersistMetadata();
}

/* Byte-aligned writes without a sparse header */
IOStatus ZoneFile::BufferedAppend(char* buffer, uint32_t data_size) {
  uint32_t left = data_size;
  uint32_t wr_size;
  uint32_t block_sz = GetBlockSize();
  IOStatus s;

  if (active_zone_ == NULL) {
    s = AllocateNewZone();
    if (!s.ok()) return s;
  }

  std::string filename;
  if (linkfiles_.size()) {
    filename = linkfiles_[0];
  } else {
    filename = "NONE(BufferedAppend)";
  }

  while (left) {
    wr_size = left;
    if (wr_size > active_zone_->capacity_) wr_size = active_zone_->capacity_;

    /* Pad to the next block boundary if needed */
    uint32_t align = wr_size % block_sz;
    uint32_t pad_sz = 0;

    if (align) pad_sz = block_sz - align;

    /* the buffer size s aligned on block size, so this is ok*/
    if (pad_sz) memset(buffer + wr_size, 0x0, pad_sz);

    uint64_t extent_length = wr_size;

    s = active_zone_->Append(buffer, wr_size + pad_sz);
    if (!s.ok()) return s;

    // ZoneExtent* new_ext = new ZoneExtent(extent_start_, extent_length,
    //                                      active_zone_, filename, this);

    // extents_.push_back(new_ext);

    extents_.push_back(
        new ZoneExtent(extent_start_, extent_length, active_zone_));

    extent_start_ = active_zone_->wp_;
    active_zone_->used_capacity_ += extent_length;
    file_size_ += extent_length;
    left -= extent_length;

    if (active_zone_->capacity_ == 0) {
      s = CloseActiveZone();
      if (!s.ok()) {
        return s;
      }
      if (left) {
        memmove((void*)(buffer), (void*)(buffer + wr_size), left);
      }
      s = AllocateNewZone();
      if (!s.ok()) return s;
    }
  }

  return IOStatus::OK();
}

/* Byte-aligned, sparse writes with inline metadata
   the caller reserves 8 bytes of data for a size header */
IOStatus ZoneFile::SparseAppend(char* sparse_buffer, uint32_t data_size) {
  uint32_t left = data_size;
  uint32_t wr_size;
  uint32_t block_sz = GetBlockSize();
  IOStatus s;

  if (active_zone_ == NULL) {
    s = AllocateNewZone();
    if (!s.ok()) return s;
  }

  while (left) {
    wr_size = left + ZoneFile::SPARSE_HEADER_SIZE;
    if (wr_size > active_zone_->capacity_) wr_size = active_zone_->capacity_;

    /* Pad to the next block boundary if needed */
    uint32_t align = wr_size % block_sz;
    uint32_t pad_sz = 0;

    if (align) pad_sz = block_sz - align;

    /* the sparse buffer has block_sz extra bytes tail allocated for padding, so
     * this is safe */
    if (pad_sz) memset(sparse_buffer + wr_size, 0x0, pad_sz);

    uint64_t extent_length = wr_size - ZoneFile::SPARSE_HEADER_SIZE;
    EncodeFixed64(sparse_buffer, extent_length);

    s = active_zone_->Append(sparse_buffer, wr_size + pad_sz);
    if (!s.ok()) return s;

    extents_.push_back(
        new ZoneExtent(extent_start_ + ZoneFile::SPARSE_HEADER_SIZE,
                       extent_length, active_zone_));

    extent_start_ = active_zone_->wp_;
    active_zone_->used_capacity_ += extent_length;
    file_size_ += extent_length;
    left -= extent_length;

    if (active_zone_->capacity_ == 0) {
      s = CloseActiveZone();
      if (!s.ok()) {
        return s;
      }
      if (left) {
        memmove((void*)(sparse_buffer + ZoneFile::SPARSE_HEADER_SIZE),
                (void*)(sparse_buffer + wr_size), left);
      }
      s = AllocateNewZone();
      if (!s.ok()) return s;
    }
  }

  return IOStatus::OK();
}

IOStatus ZoneFile::CAZAAppend(const char* data, uint32_t size, bool positioned,
                              uint64_t offset) {
  IOStatus s = IOStatus::OK();
  // printf("@@@ CAZASstBufferedAppend called : %ld %u\n",fno_,size);
  char* buf;
  if (!IsSST()) {
    return IOStatus::IOError("CAZASstBufferedAppend only apply to sst file");
  }
  if (is_sparse_) {
    return IOStatus::IOError("sst file should not be sparse");
  }

  // if(is_sst_ended_){
  //   return IOStatus::IOError("SST file could not Append after file ended");
  // }
  // buf=new char[size];
  int r = posix_memalign((void**)&buf, sysconf(_SC_PAGE_SIZE), size);
  if (r < 0) {
    printf("posix memalign error here@@@@@@@@@@@2\n");
    return IOStatus::IOError("CAZASstBufferedAppend fail to allocate memory");
  }

  memcpy(buf, data, size);
  SSTBuffer* sst_buffer = new SSTBuffer(buf, size, positioned, offset);

  if (sst_buffer == nullptr) {
    return IOStatus::IOError("CAZASstBufferedAppend fail to allocate memory");
  }
  sst_buffers_.push_back(sst_buffer);

  return s;
}

IOStatus ZonedWritableFile::CAZAFlushSST() {
  IOStatus s;

  if (!zoneFile_->IsSST()) {
    return IOStatus::OK();
  }

  zoneFile_->fno_ = fno_;
  zoneFile_->GetZbd()->SetSSTFileforZBDNoLock(fno_, zoneFile_.get());

  if (zoneFile_->GetZbd()->GetAllocationScheme() == LIZA) {
    return IOStatus::OK();
  }
  // printf("CAZAFlushSST - CAZA!!\n");
  std::vector<SSTBuffer*>* sst_buffers = zoneFile_->GetSSTBuffers();
  zoneFile_->predicted_size_ = 0;
  for (auto it : *sst_buffers) {
    zoneFile_->predicted_size_ += it->size_;
  }

  // printf("predcited size :\t%lu\n", (zoneFile_->predicted_size_ >> 20));

  for (auto it : *sst_buffers) {
    if (buffered) {
      buffer_mtx_.lock();
      s = BufferedWrite(it->content_, it->size_);
      buffer_mtx_.unlock();
      if (!s.ok()) {
        return s;
      }
    } else {
      s = zoneFile_->Append(it->content_, it->size_);
      if (!s.ok()) {
        return s;
      }
      wp += it->size_;
    }
    delete it;
  }
  // printf("%lu
  // %d\n",zoneFile_->predicted_size_,IS_BIG_SSTABLE(zoneFile_->predicted_size_));
  if(zoneFile_->level_ ==0){
    int seq = zoneFile_->GetZbd()->file_operation_sequence_.fetch_add(1);
    zoneFile_->GetZbd()->latest_file_operation_sequence_[SeqL0L1andFlush] = seq;

    if(zoneFile_->GetZbd()->coldest_type_set_== true){
      // todo
      std::lock_guard<std::mutex> lg(zoneFile_->GetZbd()->coldest_type_lock_);
      bool ok = true;
      zoneFile_->GetZbd()->check_coldest_[SeqL0L1andFlush]=true;
      
      for(int i =0;i<10;i++){
        if(zoneFile_->GetZbd()->latest_file_operation_sequence_[i]==0){
          continue;
        }
        if(zoneFile_->GetZbd()->check_coldest_[i]==true){
          continue;
        }
        ok=false;
        break;
      }

      if(ok==true){
        if(zoneFile_->GetZbd()->coldest_type_!=SeqL0L1andFlush){
          zoneFile_->GetZbd()->CBSC_mispredict_stats_[zoneFile_->GetZbd()->coldest_type_].fetch_add(1);
        }
        zoneFile_->GetZbd()->CBSC_total_predict_stats_[zoneFile_->GetZbd()->coldest_type_].fetch_add(1);
        zoneFile_->GetZbd()->coldest_type_set_=false;
      }else{
        // if(zoneFile_->GetZbd()->coldest_type_==)
      }
    }


  }
  //   zoneFile_->GetZbd()->lsm_tree_[0].fetch_add(1);
  // }else{
  zoneFile_->GetZbd()->lsm_tree_[zoneFile_->level_].fetch_add(
      zoneFile_->predicted_size_);
  if(zoneFile_->predicted_size_>>20 < 1024){
    zoneFile_->GetZbd()->file_size_distribution_[zoneFile_->predicted_size_>>20].fetch_add(1);
  
  }else{
    printf("CAZAFlushSST fsize %lu\n",zoneFile_->predicted_size_>>20);
  }

  // }
  // input_fno_.clear();
  return s;
}


/* Assumes that data and size are block aligned */
IOStatus ZoneFile::Append(void* data, int data_size) {
  uint32_t left = data_size;  
  uint32_t wr_size,
      offset = 0;  
  IOStatus s = IOStatus::OK(); 


  if (!active_zone_) {
    s = AllocateNewZone();
    if (!s.ok()) return s;
  }

  while (left) {
    if (active_zone_->capacity_ == 0) {
      PushExtent();

      s = CloseActiveZone();
      if (!s.ok()) {
        return s;
      }

      s = AllocateNewZone();
      if (!s.ok()) return s;
    }

    wr_size = left;
    if (wr_size > active_zone_->capacity_) wr_size = active_zone_->capacity_;

    s = active_zone_->Append((char*)data + offset, wr_size);
    if (!s.ok()) return s;

    file_size_ += wr_size;
    left -= wr_size;
    offset += wr_size;
  }

  return IOStatus::OK(); 
}

IOStatus ZoneFile::RecoverSparseExtents(uint64_t start, uint64_t end,
                                        Zone* zone) {
  /* Sparse writes, we need to recover each individual segment */
  IOStatus s;
  uint32_t block_sz = GetBlockSize();
  uint64_t next_extent_start = start;
  char* buffer;
  int recovered_segments = 0;
  int ret;

  ret = posix_memalign((void**)&buffer, sysconf(_SC_PAGESIZE), block_sz);
  if (ret) {
    return IOStatus::IOError("Out of memory while recovering");
  }

  while (next_extent_start < end) {
    uint64_t extent_length;

    ret = zbd_->Read(buffer, next_extent_start, block_sz, false);
    if (ret != (int)block_sz) {
      s = IOStatus::IOError("Unexpected read error while recovering");
      break;
    }

    extent_length = DecodeFixed64(buffer);
    if (extent_length == 0) {
      s = IOStatus::IOError("Unexpected extent length while recovering");
      break;
    }
    recovered_segments++;

    zone->used_capacity_ += extent_length;
    extents_.push_back(new ZoneExtent(next_extent_start + SPARSE_HEADER_SIZE,
                                      extent_length, zone));

    uint64_t extent_blocks = (extent_length + SPARSE_HEADER_SIZE) / block_sz;
    if ((extent_length + SPARSE_HEADER_SIZE) % block_sz) {
      extent_blocks++;
    }
    next_extent_start += extent_blocks * block_sz;
  }

  free(buffer);
  return s;
}

IOStatus ZoneFile::Recover() {
  /* If there is no active extent, the file was either closed gracefully
     or there were no writes prior to a crash. All good.*/
  if (!HasActiveExtent()) return IOStatus::OK();

  /* Figure out which zone we were writing to */
  Zone* zone = zbd_->GetIOZone(extent_start_);

  if (zone == nullptr) {
    return IOStatus::IOError(
        "Could not find zone for extent start while recovering");
  }

  if (zone->wp_ < extent_start_) {
    return IOStatus::IOError("Zone wp is smaller than active extent start");
  }

  /* How much data do we need to recover? */
  uint64_t to_recover = zone->wp_ - extent_start_;

  /* Do we actually have any data to recover? */
  if (to_recover == 0) {
    /* Mark up the file as having no missing extents */
    extent_start_ = NO_EXTENT;
    return IOStatus::OK();
  }

  /* Is the data sparse or was it writted direct? */
  if (is_sparse_) {
    IOStatus s = RecoverSparseExtents(extent_start_, zone->wp_, zone);
    if (!s.ok()) return s;
  } else {
    /* For non-sparse files, the data is contigous and we can recover directly
       any missing data using the WP */
    zone->used_capacity_ += to_recover;
    extents_.push_back(new ZoneExtent(extent_start_, to_recover, zone));
  }

  /* Mark up the file as having no missing extents */
  extent_start_ = NO_EXTENT;

  /* Recalculate file size */
  file_size_ = 0;
  for (uint32_t i = 0; i < extents_.size(); i++) {
    file_size_ += extents_[i]->length_;
  }

  return IOStatus::OK();
}

void ZoneFile::ReplaceExtentList(std::vector<ZoneExtent*> new_list) {
  assert(IsOpenForWR() && new_list.size() > 0);
  assert(new_list.size() == extents_.size());

  WriteLock lck(this);
  extents_ = new_list;
}

void ZoneFile::AddLinkName(const std::string& linkf) {
  linkfiles_.push_back(linkf);
}

IOStatus ZoneFile::RenameLink(const std::string& src, const std::string& dest) {
  auto itr = std::find(linkfiles_.begin(), linkfiles_.end(), src);
  if (itr != linkfiles_.end()) {
    linkfiles_.erase(itr);
    linkfiles_.push_back(dest);
  } else {
    return IOStatus::IOError("RenameLink: Failed to find the linked file");
  }
  return IOStatus::OK();
}

IOStatus ZoneFile::RemoveLinkName(const std::string& linkf) {
  assert(GetNrLinks());
  auto itr = std::find(linkfiles_.begin(), linkfiles_.end(), linkf);
  if (itr != linkfiles_.end()) {
    linkfiles_.erase(itr);
  } else {
    return IOStatus::IOError("RemoveLinkInfo: Failed to find the link file");
  }
  return IOStatus::OK();
}

IOStatus ZoneFile::SetWriteLifeTimeHint(Env::WriteLifeTimeHint lifetime,
                                        int level) {
  // lifetime_ = lifetime;
  // printf("SetWriteLifeTimeHint : %s %d %d\n", linkfiles_[0].c_str(),
  // lifetime,
  //        level);
  (void)(lifetime);
  // wal-log 파일은 수명이 SHORT(2)
  // sst 파일은 3부터시작
  // L0,L1- SST MEDIUM(3), L2-LONG(4), L3-EXTREME(5)
  // SHORT는 SSTable이 아닌,
  // Write-Ahead Log (WAL)과 LSM-TREE의 요약 정보를 포함하는 Manifest와 같은
  // 데이터에 할당된다.
  // lifetime에 .log는 SHORT라 정수변환수 2로 찍히는데
  // level로 수명을 정수로 변환수에 -1씩해서 SHORT(1), MEDIUM(2), LONG(3),
  // EXTREME(4)으로 변환
  if (is_wal_) {
    lifetime_ = Env::WLTH_SHORT;
    return IOStatus::OK();
  }
  switch (level) {
    case 0:
      /* fall through */
    case 1:
      lifetime_ = Env::WLTH_SHORT;
      break;
    case 2:
      lifetime_ = Env::WLTH_MEDIUM;
      break;
    case 3:
      lifetime_ = Env::WLTH_LONG;
      break;
    default:
      lifetime_ = Env::WLTH_EXTREME;
      break;
  }
  // printf("%d -> %d\n", level, lifetime_);

  return IOStatus::OK();
}
void ZonedWritableFile::SetMinMaxKeyAndLevel(const Slice& s, const Slice& l,
                                             const int output_level) {
  if (output_level < 0) {
    printf(
        "@@@ ZonedWritableFile::SetMinMaxAndLEvel :: failed , level should be "
        "> 0\n");
    return;
  }
  // printf("set min max : fno :%ld at %d\n", zoneFile_->fno_, output_level);

  zoneFile_->smallest_ = s;
  zoneFile_->largest_ = l;
  zoneFile_->level_ = output_level;

  zoneFile_->GetZbd()->cur_max_level_ = std::max(zoneFile_->GetZbd()->cur_max_level_,output_level);

  return;
}

void ZoneFile::ReleaseActiveZone() {
  assert(active_zone_ != nullptr);
  bool ok = active_zone_->Release();
  assert(ok);
  (void)ok;
  active_zone_ = nullptr;
}

void ZoneFile::SetActiveZone(Zone* zone) {
  assert(active_zone_ == nullptr);
  assert(zone->IsBusy());
  active_zone_ = zone;
}

// **sparse_buffer**는 데이터가 연속되지 않고 중간에 비어 있는(데이터가 없는)
// 부분을 포함한 파일, 또한 파일의 데이터를 메모리에서 처리할 때, 중간에 비어
// 있는 공간을 효율적으로 관리하고, 필요한 경우 패딩(block padding)을 추가하여
// 데이터를 저장하거나 처리하는 역할
ZonedWritableFile::ZonedWritableFile(ZonedBlockDevice* zbd, bool _buffered,
                                     std::shared_ptr<ZoneFile> zoneFile) {
  assert(zoneFile->IsOpenForWR());
  wp = zoneFile->GetFileSize();

  buffered = _buffered;
  block_sz = zbd->GetBlockSize();
  zoneFile_ = zoneFile;
  buffer_pos = 0;
  sparse_buffer = nullptr;
  buffer = nullptr;

  if (buffered) {
    if (zoneFile->IsSparse()) {
      size_t sparse_buffer_sz;

      sparse_buffer_sz =
          1024 * 1024 + block_sz; /* one extra block size for padding */
      // 메모리를 페이지 크기에 맞춰 정렬하여 할당
      int ret = posix_memalign((void**)&sparse_buffer, sysconf(_SC_PAGESIZE),
                               sparse_buffer_sz);

      if (ret) sparse_buffer = nullptr;

      assert(sparse_buffer != nullptr);

      buffer_sz = sparse_buffer_sz - ZoneFile::SPARSE_HEADER_SIZE - block_sz;
      buffer = sparse_buffer + ZoneFile::SPARSE_HEADER_SIZE;
    } else {
      buffer_sz = 1024 * 1024;
      int ret =
          posix_memalign((void**)&buffer, sysconf(_SC_PAGESIZE), buffer_sz);

      if (ret) buffer = nullptr;
      assert(buffer != nullptr);
    }
  }

  open = true;
}

ZonedWritableFile::~ZonedWritableFile() {
  IOStatus s = CloseInternal();
  if (buffered) {
    if (sparse_buffer != nullptr) {
      free(sparse_buffer);
    } else {
      free(buffer);
    }
  }

  if (!s.ok()) {
    zoneFile_->GetZbd()->SetZoneDeferredStatus(s);
  }
}

MetadataWriter::~MetadataWriter() {}

IOStatus ZonedWritableFile::Truncate(uint64_t size,
                                     const IOOptions& /*options*/,
                                     IODebugContext* /*dbg*/) {
  zoneFile_->SetFileSize(size);
  return IOStatus::OK();
}
// 데이터를 디스크에 동기화
IOStatus ZonedWritableFile::DataSync() {
  if (buffered) {
    IOStatus s;
    buffer_mtx_.lock();
    /* Flushing the buffer will result in a new extent added to the list*/
    s = FlushBuffer();
    buffer_mtx_.unlock();
    if (!s.ok()) {
      return s;
    }

    /* We need to persist the new extent, if the file is not sparse,
     * as we can't use the active zone WP, which is block-aligned, to recover
     * the file size */
    if (!zoneFile_->IsSparse()) return zoneFile_->PersistMetadata();
  } else {
    /* For direct writes, there is no buffer to flush, we just need to push
       an extent for the latest written data */
    zoneFile_->PushExtent();
  }

  return IOStatus::OK();
}

IOStatus ZonedWritableFile::Fsync(const IOOptions& /*options*/,
                                  IODebugContext* /*dbg*/) {
  IOStatus s;
  // 동기화 시간 측정
  ZenFSMetricsLatencyGuard guard(zoneFile_->GetZBDMetrics(),
                                 zoneFile_->GetIOType() == IOType::kWAL
                                     ? ZENFS_WAL_SYNC_LATENCY
                                     : ZENFS_NON_WAL_SYNC_LATENCY,
                                 Env::Default());
  // QPS(Queries Per Second) 성능 지표를 기록하는 함수 호출
  zoneFile_->GetZBDMetrics()->ReportQPS(ZENFS_SYNC_QPS, 1);
  // 데이터를 먼저 동기화
  s = DataSync();
  if (!s.ok()) return s;

  /* As we've already synced the metadata in DataSync, no need to do it again */

  if (buffered && !zoneFile_->IsSparse()) return IOStatus::OK();

  return zoneFile_->PersistMetadata();
}

// Fsync와 달리 파일 메타데이터까지 추가로 동기화하지는 않는다

IOStatus ZonedWritableFile::Sync(const IOOptions& /*options*/,
                                 IODebugContext* /*dbg*/) {
  return DataSync();
}

IOStatus ZonedWritableFile::Flush(const IOOptions& /*options*/,
                                  IODebugContext* /*dbg*/) {
  return IOStatus::OK();
}

IOStatus ZonedWritableFile::RangeSync(uint64_t offset, uint64_t nbytes,
                                      const IOOptions& /*options*/,
                                      IODebugContext* /*dbg*/) {
  if (wp < offset + nbytes) return DataSync();

  return IOStatus::OK();
}

IOStatus ZonedWritableFile::Close(const IOOptions& /*options*/,
                                  IODebugContext* /*dbg*/) {
  return CloseInternal();
}

IOStatus ZonedWritableFile::CloseInternal() {
  if (!open) {
    return IOStatus::OK();
  }

  IOStatus s = DataSync();
  if (!s.ok()) return s;

  s = zoneFile_->CloseWR();
  if (!s.ok()) return s;

  open = false;
  return s;
}

IOStatus ZonedWritableFile::FlushBuffer() {
  IOStatus s;

  if (buffer_pos == 0) return IOStatus::OK();

  if (zoneFile_->IsSparse()) {
    s = zoneFile_->SparseAppend(sparse_buffer, buffer_pos);
  } else {
    s = zoneFile_->BufferedAppend(buffer, buffer_pos);
  }

  if (!s.ok()) {
    return s;
  }

  wp += buffer_pos;
  buffer_pos = 0;

  return IOStatus::OK();
}

IOStatus ZonedWritableFile::BufferedWrite(const Slice& slice) {
  uint32_t data_left = slice.size();
  char* data = (char*)slice.data();
  IOStatus s;

  while (data_left) {
    uint32_t buffer_left = buffer_sz - buffer_pos;
    uint32_t to_buffer;

    if (!buffer_left) {
      s = FlushBuffer();
      if (!s.ok()) return s;
      buffer_left = buffer_sz;
    }

    to_buffer = data_left;
    if (to_buffer > buffer_left) {
      to_buffer = buffer_left;
    }

    memcpy(buffer + buffer_pos, data, to_buffer);
    buffer_pos += to_buffer;
    data_left -= to_buffer;
    data += to_buffer;
  }

  return IOStatus::OK();
}

IOStatus ZonedWritableFile::BufferedWrite(char* data, uint32_t data_left) {
  // uint32_t data_left = slice.size();
  IOStatus s;

  while (data_left) {
    uint32_t buffer_left = buffer_sz - buffer_pos;
    uint32_t to_buffer;

    if (!buffer_left) {
      s = FlushBuffer();
      if (!s.ok()) return s;
      buffer_left = buffer_sz;
    }

    to_buffer = data_left;
    if (to_buffer > buffer_left) {
      to_buffer = buffer_left;
    }

    memcpy(buffer + buffer_pos, data, to_buffer);
    buffer_pos += to_buffer;
    data_left -= to_buffer;
    data += to_buffer;
  }

  return IOStatus::OK();
}

IOStatus ZonedWritableFile::Append(const Slice& data,
                                   const IOOptions& /*options*/,
                                   IODebugContext* /*dbg*/) {
  IOStatus s;
  ZenFSMetricsLatencyGuard guard(zoneFile_->GetZBDMetrics(),
                                 zoneFile_->GetIOType() == IOType::kWAL
                                     ? ZENFS_WAL_WRITE_LATENCY
                                     : ZENFS_NON_WAL_WRITE_LATENCY,
                                 Env::Default());
  zoneFile_->GetZBDMetrics()->ReportQPS(ZENFS_WRITE_QPS, 1);
  zoneFile_->GetZBDMetrics()->ReportThroughput(ZENFS_WRITE_THROUGHPUT,
                                               data.size());
  //
  if (zoneFile_->IsSST() && zoneFile_->GetAllocationScheme() != LIZA) {
    // printf("append->CAZAAppend!!\n");
    return zoneFile_->CAZAAppend(data.data(), data.size(), true, 0);
  }
  //
  if (buffered) {
    buffer_mtx_.lock();
    s = BufferedWrite(data);
    buffer_mtx_.unlock();
  } else {
    s = zoneFile_->Append((void*)data.data(), data.size());
    if (s.ok()) wp += data.size();
  }

  return s;
}

IOStatus ZonedWritableFile::PositionedAppend(const Slice& data, uint64_t offset,
                                             const IOOptions& /*options*/,
                                             IODebugContext* /*dbg*/) {
  IOStatus s;
  ZenFSMetricsLatencyGuard guard(zoneFile_->GetZBDMetrics(),
                                 zoneFile_->GetIOType() == IOType::kWAL
                                     ? ZENFS_WAL_WRITE_LATENCY
                                     : ZENFS_NON_WAL_WRITE_LATENCY,
                                 Env::Default());
  zoneFile_->GetZBDMetrics()->ReportQPS(ZENFS_WRITE_QPS, 1);
  zoneFile_->GetZBDMetrics()->ReportThroughput(ZENFS_WRITE_THROUGHPUT,
                                               data.size());
  //
  // printf("positionedAppend!!\n");
  // if (zoneFile_->is_wal_) {
  //   uint64_t lifetime = zoneFile_->GetWriteLifeTimeHint();
  //   // printf("WAL : %llu\n", (unsigned long long)lifetime);
  // }
  // if (zoneFile_->is_sst_) {
  //   uint64_t lifetime = zoneFile_->GetWriteLifeTimeHint();
  //   if (lifetime == 2) {
  //     // printf("SST : %llu\n", (unsigned long long)lifetime);
  //   }
  // }

  if ((zoneFile_->is_wal_ && zoneFile_->GetZCRunning_()) ||
      (zoneFile_->is_sst_ && zoneFile_->GetWriteLifeTimeHint() == 2 &&
       zoneFile_->GetZCRunning_())) {
    // printf("WAL || FLUSH\n");
    while (zoneFile_->GetZCRunning_()) {
      //
    }
  }

  if (zoneFile_->IsSST() && zoneFile_->GetAllocationScheme() != LIZA) {
    // printf("append->CAZAAppend!!\n");
    s = zoneFile_->CAZAAppend(data.data(), data.size(), true, offset);
    return s;
  }

  //
  if (offset != wp) {
    assert(false);
    return IOStatus::IOError("positioned append not at write pointer");
  }

  if (buffered) {
    buffer_mtx_.lock();
    s = BufferedWrite(data);
    buffer_mtx_.unlock();
  } else {
    s = zoneFile_->Append((void*)data.data(), data.size());
    if (s.ok()) wp += data.size();
  }

  return s;
}

void ZonedWritableFile::SetWriteLifeTimeHint(Env::WriteLifeTimeHint hint) {
  // printf("ZonedWritableFile::SetWriteLifeTimeHint: level_ = %d\n", level_);
  if (zoneFile_->is_sst_) {
    zoneFile_->fno_ = fno_;
    // zoneFile_->input_fno_ = input_fno_;
    zoneFile_->GetZbd()->SetSSTFileforZBDNoLock(fno_, zoneFile_.get());
    zoneFile_->level_ = level_;
  }
  zoneFile_->SetWriteLifeTimeHint(hint, level_);
  // auto lifetime_ = zoneFile_->GetWriteLifeTimeHint();
  // printf("level_: %d, hint: %d -> lifetime_ : %d\n", level_, hint,
  // lifetime_);
  // std::vector<uint64_t> fno_list;
  // zoneFile_->GetZbd()->SameLevelFileList(level_, fno_list);
  // for (auto fno : fno_list) {
  //   ZoneFile* zFile = zoneFile_->GetZbd()->GetSSTZoneFileInZBDNoLock(fno);
  //   if (zFile != nullptr) {
  //   }
  // }
}

IOStatus ZonedSequentialFile::Read(size_t n, const IOOptions& /*options*/,
                                   Slice* result, char* scratch,
                                   IODebugContext* /*dbg*/) {
  IOStatus s;

  s = zoneFile_->PositionedRead(rp, n, result, scratch, direct_);
  if (s.ok()) rp += result->size();

  return s;
}

IOStatus ZonedSequentialFile::Skip(uint64_t n) {
  if (rp + n >= zoneFile_->GetFileSize())
    return IOStatus::InvalidArgument("Skip beyond end of file");
  rp += n;
  return IOStatus::OK();
}

IOStatus ZonedSequentialFile::PositionedRead(uint64_t offset, size_t n,
                                             const IOOptions& /*options*/,
                                             Slice* result, char* scratch,
                                             IODebugContext* /*dbg*/) {
  return zoneFile_->PositionedRead(offset, n, result, scratch, direct_);
}

IOStatus ZonedRandomAccessFile::Read(uint64_t offset, size_t n,
                                     const IOOptions& /*options*/,
                                     Slice* result, char* scratch,
                                     IODebugContext* /*dbg*/) const {
  return zoneFile_->PositionedRead(offset, n, result, scratch, direct_);
}

IOStatus ZoneFile::MigrateData(uint64_t offset, uint32_t length,
                               Zone* target_zone) {
  uint32_t step = 128 << 10;
  // uint32_t step = length;
  uint32_t read_sz = step;
  int block_sz = zbd_->GetBlockSize();
  // uint64_t align = step % block_sz;

  // if (align) {
  //   step += (block_sz - align);
  // }

  assert(offset % block_sz == 0);
  if (offset % block_sz != 0) {
    return IOStatus::IOError("MigrateData offset is not aligned!\n");
  }

  char* buf;
  int ret = posix_memalign((void**)&buf, block_sz, step);
  if (ret) {
    return IOStatus::IOError("failed allocating alignment write buffer\n");
  }

  int pad_sz = 0;
  while (length > 0) {
    read_sz = length > read_sz ? read_sz : length;
    pad_sz = read_sz % block_sz == 0 ? 0 : (block_sz - (read_sz % block_sz));

    int r = zbd_->Read(buf, offset, read_sz + pad_sz, true);
    if (r < 0) {
      free(buf);
      return IOStatus::IOError(strerror(errno));
    }
    target_zone->Append(buf, r);
    length -= read_sz;
    offset += r;
  }

  free(buf);

  return IOStatus::OK();
}

}  // namespace ROCKSDB_NAMESPACE

#endif  // !defined(ROCKSDB_LITE) && !defined(OS_WIN)
