/*
 *  Copyright (c) 2012 Daisuke Okanohara
 *
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 *
 *   1. Redistributions of source code must retain the above Copyright
 *      notice, this list of conditions and the following disclaimer.
 *
 *   2. Redistributions in binary form must reproduce the above Copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *
 *   3. Neither the name of the authors nor the names of its contributors
 *      may be used to endorse or promote products derived from this
 *      software without specific prior written permission.
 */

#include <cassert>
#include "Const.hpp"
#include "Util.hpp"
#include "EnumCoder.hpp"
#include "RSDic.hpp"

using namespace std;

namespace rsdic{

  RSDic::RSDic() : bits_(),
		   pointer_blocks_(),
		   rank_blocks_(),
		   select_one_inds_(),
		   select_zero_inds_(),
		   rank_small_blocks_(),
		   num_(0), 
		   one_num_(0),
		   zero_num_(0),
		   prev_one_num_(0),
		   offset_(0),
		   buf_(0),
		   last_encoded_pos_(0){}

  RSDic::~RSDic(){
  }

  void RSDic::Clear() {
    bits_.clear();
    pointer_blocks_.clear();
    rank_blocks_.clear();
    select_one_inds_.clear();
    select_zero_inds_.clear();
    rank_small_blocks_.clear();
    num_ = 0;
    one_num_ = 0;
    zero_num_ = 0;
    prev_one_num_ = 0;
    offset_ = 0;
    buf_ = 0;
    last_encoded_pos_ = 0;
  }

  void RSDic::Delete()
  {
    bits_.clear();
    bits_.shrink_to_fit();
    pointer_blocks_.clear();
    pointer_blocks_.shrink_to_fit();
    rank_blocks_.clear();
    rank_blocks_.shrink_to_fit();
    select_one_inds_.clear();
    select_one_inds_.shrink_to_fit();
    rank_small_blocks_.clear();
    rank_small_blocks_.shrink_to_fit();
    num_ = 0;
    one_num_ = 0;
    zero_num_ = 0;
    prev_one_num_ = 0;
    offset_ = 0;
    buf_ = 0;
    last_encoded_pos_ = 0;
  }

  void RSDic::PushBack(const bool bit){
    if (num_ % kSmallBlockSize == 0){
      WriteBlock();
    }
    if (bit){
      buf_ |= (1LLU << (num_ % kSmallBlockSize));
      if ((one_num_ % kSelectBlockSize) == 0){
	select_one_inds_.push_back(num_ / kLargeBlockSize);
      }
      ++one_num_;
    } else {
      if ((zero_num_ % kSelectBlockSize) == 0){
	select_zero_inds_.push_back(num_ / kLargeBlockSize);
      }
      ++zero_num_;
    }
    ++num_;  
  }

  void RSDic::LastProcess(){
    WriteBlock();
    bits_.resize(bits_.size());
    vector<uint64_t> (bits_).swap(bits_);
    pointer_blocks_.resize(pointer_blocks_.size());
    vector<rsdic_uint> (pointer_blocks_).swap(pointer_blocks_);
    rank_blocks_.resize(rank_blocks_.size());
    vector<rsdic_uint> (rank_blocks_).swap(rank_blocks_);
    select_one_inds_.resize(select_one_inds_.size());
    vector<rsdic_uint> (select_one_inds_).swap(select_one_inds_);
    select_zero_inds_.resize(select_zero_inds_.size());
    vector<rsdic_uint> (select_zero_inds_).swap(select_zero_inds_);
    rank_small_blocks_.resize(rank_small_blocks_.size());
    vector<uint8_t> (rank_small_blocks_).swap(rank_small_blocks_);
  }
  

  void RSDic::WriteBlock(){
    if (num_ > 0) {
      if(num_ != rank_small_blocks_.size()*64){
	uint64_t rank_sb = one_num_ - prev_one_num_;
	rank_small_blocks_.push_back(rank_sb);
	prev_one_num_ = one_num_;
	uint64_t len = EnumCoder::Len(rank_sb);
	uint64_t code = 0;
	if (len == kSmallBlockSize){
	  code = buf_; // use raw 
	} else {
	  code = EnumCoder::Encode(buf_, rank_sb);
	}
	uint64_t new_size =  Util::Floor(offset_ + len, kSmallBlockSize);
	if (new_size > bits_.size()) {
	  bits_.push_back(0);
	}
	Util::SetSlice(bits_, offset_, len, code);
	buf_ = 0;
	offset_ += len;
	last_encoded_pos_ = num_;
      }
    }
    if ((num_ % kLargeBlockSize) == 0){
      rank_blocks_.push_back(one_num_);
      pointer_blocks_.push_back(offset_);
    }
  }


  bool RSDic::GetBit(uint64_t pos) const{
    if(pos < last_encoded_pos_){
      uint64_t lblock = pos / kLargeBlockSize;
      uint64_t pointer = pointer_blocks_[lblock];
      uint64_t sblock = pos / kSmallBlockSize;
      for (uint64_t i = lblock * kSmallBlockPerLargeBlock; i < sblock; ++i){
	pointer += EnumCoder::Len(rank_small_blocks_[i]);
      }
      uint64_t rank_sb = rank_small_blocks_[sblock];
      uint64_t code = Util::GetSlice(bits_, pointer, EnumCoder::Len(rank_sb));
      return EnumCoder::GetBit(code, rank_sb, pos % kSmallBlockSize);
    }
    else{
      return (buf_ >> (pos - last_encoded_pos_)) & 1LLU;
    }
  }

  uint64_t RSDic::Rank(uint64_t pos, bool bit) const{
    uint64_t lblock = pos / kLargeBlockSize;
    uint64_t pointer = pointer_blocks_[lblock];
    uint64_t sblock = pos / kSmallBlockSize;
    uint64_t rank = rank_blocks_[lblock];
    for (uint64_t i = lblock * kSmallBlockPerLargeBlock; i < sblock; ++i){
      uint64_t rank_sb = rank_small_blocks_[i];
      rank += rank_sb;
      pointer += EnumCoder::Len(rank_sb);
    }
    if (pos % kSmallBlockSize == 0){
      return Util::GetNum(bit, rank, pos);
    }
    if(pos > last_encoded_pos_){
      rank += (EnumCoder::PopCount(buf_ & ((1LLU << (pos - last_encoded_pos_)) - 1)));
      return Util::GetNum(bit, rank, pos); 
    }
    uint64_t rank_sb = rank_small_blocks_[sblock];
    uint64_t code = Util::GetSlice(bits_, pointer, EnumCoder::Len(rank_sb));
    rank += EnumCoder::Rank(code, rank_sb, pos % kSmallBlockSize);
    return Util::GetNum(bit, rank, pos);
  }

  pair<uint64_t, uint64_t> RSDic::GetBitAndRank(uint64_t pos) const{
    uint64_t lblock = pos / kLargeBlockSize;
    uint64_t pointer = pointer_blocks_[lblock];
    uint64_t sblock = pos / kSmallBlockSize;
    uint64_t rank = rank_blocks_[lblock];
    for (uint64_t i = lblock * kSmallBlockPerLargeBlock; i < sblock; ++i){
      uint64_t rank_sb = rank_small_blocks_[i];
      rank += rank_sb;
      pointer += EnumCoder::Len(rank_sb);
    }
    uint64_t rank_sb = rank_small_blocks_[sblock];
    uint64_t code = Util::GetSlice(bits_, pointer, EnumCoder::Len(rank_sb));
    rank += EnumCoder::Rank(code, rank_sb, pos % kSmallBlockSize);
    uint64_t ret_bit = EnumCoder::GetBit(code, rank_sb, pos % kSmallBlockSize);
    return make_pair(ret_bit, Util::GetNum(ret_bit, rank, pos));
  }


  uint64_t RSDic::Select1(uint64_t ind) const{
    uint64_t select_ind = ind / kSelectBlockSize;
    uint64_t lblock = select_one_inds_[select_ind];
    for (; lblock < rank_blocks_.size(); ++lblock){
      if (ind < rank_blocks_[lblock]) break;
    }
    --lblock;
    uint64_t sblock = lblock * kSmallBlockPerLargeBlock;
    uint64_t pointer = pointer_blocks_[lblock];
    uint64_t remain = ind - rank_blocks_[lblock] + 1;

    for (; sblock < rank_small_blocks_.size(); ++sblock){
      const uint64_t rank_sb = rank_small_blocks_[sblock];
      if (remain <= rank_sb) break;
      remain -= rank_sb;
      pointer += EnumCoder::Len(rank_sb);
    }
    if(sblock == rank_small_blocks_.size()){
      return sblock * kSmallBlockSize + EnumCoder::SelectRaw(buf_, remain);
    }
    uint64_t rank_sb = rank_small_blocks_[sblock];
    uint64_t code = Util::GetSlice(bits_, pointer, EnumCoder::Len(rank_sb));
    return sblock * kSmallBlockSize + EnumCoder::Select1(code, rank_sb, remain);
  }

  uint64_t RSDic::Select0(uint64_t ind) const{
    uint64_t select_ind = ind / kSelectBlockSize;
    uint64_t lblock = select_zero_inds_[select_ind];
    for (; lblock < rank_blocks_.size(); ++lblock){
      if (lblock * kLargeBlockSize - rank_blocks_[lblock] > ind) break;
    }
    --lblock;

    uint64_t sblock = lblock * kSmallBlockPerLargeBlock;
    uint64_t pointer = pointer_blocks_[lblock];
    uint64_t remain = ind - lblock * kLargeBlockSize + rank_blocks_[lblock] + 1;

    for (; sblock < rank_small_blocks_.size(); ++sblock){
      const uint64_t rank_sb = kSmallBlockSize - rank_small_blocks_[sblock];
      if (remain <= rank_sb) break;
      remain -= rank_sb;
      pointer += EnumCoder::Len(rank_sb);
    }
  
    if(sblock == rank_small_blocks_.size()){
      return sblock * kSmallBlockSize + EnumCoder::SelectRaw(~buf_, remain);
    }
    uint64_t rank_sb = rank_small_blocks_[sblock];
    uint64_t code = Util::GetSlice(bits_, pointer, EnumCoder::Len(rank_sb));
    return sblock * kSmallBlockSize + EnumCoder::Select0(code, rank_sb, remain);
  }

  uint64_t RSDic::Select(uint64_t ind, bool bit) const{
    if (bit) return Select1(ind);
    else return Select0(ind);
  }

  void RSDic::Save(ostream& os){
    LastProcess();
    Save(os, bits_);
    Save(os, pointer_blocks_);
    Save(os, rank_blocks_);
    Save(os, select_one_inds_);
    Save(os, select_zero_inds_);
    Save(os, rank_small_blocks_);
    os.write((const char*)&num_, sizeof(num_));
    os.write((const char*)&one_num_, sizeof(one_num_));
    os.write((const char*)&zero_num_, sizeof(zero_num_));
    os.write((const char*)&prev_one_num_, sizeof(prev_one_num_));
    os.write((const char*)&offset_, sizeof(offset_));
    os.write((const char*)&buf_, sizeof(buf_));
    os.write((const char*)&last_encoded_pos_, sizeof(last_encoded_pos_));
  }

  void RSDic::Load(istream& is){
    Load(is, bits_);
    Load(is, pointer_blocks_);
    Load(is, rank_blocks_);
    Load(is, select_one_inds_);
    Load(is, select_zero_inds_);
    Load(is, rank_small_blocks_);
    is.read((char*)&num_, sizeof(num_));
    is.read((char*)&one_num_, sizeof(one_num_));
    is.read((char*)&zero_num_, sizeof(zero_num_));
    is.read((char*)&prev_one_num_, sizeof(prev_one_num_));
    is.read((char*)&offset_, sizeof(offset_));
    is.read((char*)&buf_, sizeof(buf_));
    is.read((char*)&last_encoded_pos_, sizeof(last_encoded_pos_));
  }

  uint64_t RSDic::GetUsageBytes() const{
    /*
      cout << endl
      << "bits:" << bits_.size() * sizeof(bits_[0]) << endl
      << " ptb:" << pointer_blocks_.size() * sizeof(pointer_blocks_[0]) << endl
      << "  rb:" << rank_blocks_.size() * sizeof(rank_blocks_[0]) << endl
      << " soi:" << select_one_inds_.size() * sizeof(select_one_inds_[0]) << endl
      << " soz:" <<     select_zero_inds_.size() * sizeof(select_zero_inds_[0]) << endl
      << " rsb:" <<     rank_small_blocks_.size() * sizeof(rank_small_blocks_[0]) << endl;
    */
    return
      bits_.size() * sizeof(bits_[0]) +
      pointer_blocks_.size() * sizeof(pointer_blocks_[0]) +
      rank_blocks_.size() * sizeof(rank_blocks_[0]) +
      select_one_inds_.size() * sizeof(select_one_inds_[0]) +
      select_zero_inds_.size() * sizeof(select_zero_inds_[0]) +
      rank_small_blocks_.size() * sizeof(rank_small_blocks_[0]) +
      sizeof(num_) +
      sizeof(one_num_) +
      sizeof(zero_num_) +
      sizeof(prev_one_num_)+
      sizeof(offset_) +
      sizeof(buf_);
  }

  bool RSDic::operator == (const RSDic& bv) const{
    return
      bits_ == bv.bits_ &&
      pointer_blocks_ == bv.pointer_blocks_ &&
      rank_blocks_ == bv.rank_blocks_ &&
      select_one_inds_ == bv.select_one_inds_ &&
      select_zero_inds_ == bv.select_zero_inds_ &&
      rank_small_blocks_ == bv.rank_small_blocks_ &&
      num_ == bv.num_ &&
      one_num_ == bv.one_num_;
  }


}
