/* hash_table_rp.cpp
MIT License

Copyright (c) 2018 Yoshimasa Takabatake

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "hash_table_rp.hpp"

using namespace std;

namespace solca_comp{
void HashTableRP::Init(){

  hash_index_ = 7;
  hash_size_ = kPrimes[hash_index_];
  num_frules_ = 0;
  insert_ht_pos_ = 0;
  num_not_used_ = hash_size_;
  max_pair_ = 0;
  max_occ_ = 0;
  is_first_ = true;
  CFunc::ResizeVec(datas_,
                   hash_size_);
  CFunc::ResizeVec(hash_table_,
                   hash_size_);
  CFunc::ResizeVec(next_,
                   hash_size_);
  CFunc::ResizeVec(not_used_,
                   hash_size_);

  for(size_t i = 0; i < hash_size_; i++){
    datas_[i] = make_tuple(hash_size_,
                           hash_size_,
                           0);
    hash_table_[i] = hash_size_;
    next_[i] = hash_size_;
    not_used_[i] = hash_size_ - i - 1;
  }
}

void HashTableRP::Clear(){
  if(!is_first_ && num_frules_ < hash_size_/2 && hash_index_){
    hash_size_ = kPrimes[--hash_index_];
    CFunc::ResizeVec(datas_,
                     hash_size_);
    CFunc::ResizeVec(hash_table_,
                     hash_size_);
    CFunc::ResizeVec(next_,
                     hash_size_);
    CFunc::ResizeVec(not_used_,
                     hash_size_);
  }
  is_first_ = false;
  num_frules_ = 0;
  num_not_used_ = hash_size_;
  max_pair_ = 0;
  max_occ_ = 0;
  for(size_t i = 0;
      i < hash_size_;
      i++){
    get<0>(datas_[i]) = hash_size_;
    get<1>(datas_[i]) = hash_size_;
    get<2>(datas_[i]) = 0;
    hash_table_[i] = hash_size_;
    next_[i] = hash_size_;
    not_used_[i] = hash_size_ - i - 1;
  }
}

void HashTableRP::Add(const uint64_t kLeft,
                      const uint64_t kRight,
                      const uint64_t kOcc){
  // std::cout << __func__ << ": hash_index_ = " << hash_index_ << ", hash_size_ = " << hash_size_ << ", num_frules_ = " << num_frules_ << ", kLeft = " << kLeft << ", kRight = " << kRight << ", kOcc = " << kOcc << std::endl;
  insert_ht_pos_ = (((uint64_t)(kLeft % hash_size_) * kPrimes[27]) + kRight) % hash_size_;
  uint32_t pos = hash_table_[insert_ht_pos_];
  while(pos != hash_size_){
    if(get<0>(datas_[pos]) == kLeft &&
       get<1>(datas_[pos]) == kRight){
      get<2>(datas_[pos]) += kOcc;
      if(get<2>(datas_[pos]) > max_occ_){
        max_occ_ = get<2>(datas_[pos]);
        max_pair_ = pos;
      }
      return ;
    }
    pos = next_[pos];
  }
  uint32_t insert_pos = not_used_[--num_not_used_];
  InsertToHash(insert_pos);

  get<0>(datas_[insert_pos]) = kLeft;
  get<1>(datas_[insert_pos]) = kRight;
  get<2>(datas_[insert_pos]) = kOcc;

  if(kOcc > max_occ_){
    max_occ_ = kOcc;
    max_pair_ = insert_pos;
  }
  num_frules_++;
  if(num_frules_ == hash_size_){
    ReHash();
  }
}

void HashTableRP::Subtract(const uint64_t kLeft,
                           const uint64_t kRight,
                           const uint64_t kOcc){
  uint64_t hash_pos = (((kLeft % hash_size_) * kPrimes[27]) + kRight) % hash_size_;
  uint64_t pos = hash_table_[hash_pos];
  if(kLeft == get<0>(datas_[pos]) &&
     kRight == get<1>(datas_[pos])){
    get<2>(datas_[pos]) -= kOcc;
    if(get<2>(datas_[pos]) == 0){
      hash_table_[hash_pos] = next_[pos];
      next_[pos] = hash_size_;
      get<0>(datas_[pos]) = hash_size_;
      get<1>(datas_[pos]) = hash_size_;
      get<2>(datas_[pos]) = 0;
      not_used_[num_not_used_++] = pos;
      num_frules_--;
    }
  }
  else{
    uint64_t prev_pos = pos;
    pos = next_[pos];
    while(pos != hash_size_){
      if(kLeft == get<0>(datas_[pos]) &&
         kRight == get<1>(datas_[pos])){
        get<2>(datas_[pos]) -= kOcc;
        if(get<2>(datas_[pos]) == 0){
          hash_table_[hash_pos] = next_[pos];
          next_[pos] = hash_size_;
          get<0>(datas_[pos]) = hash_size_;
          get<1>(datas_[pos]) = hash_size_;
          get<2>(datas_[pos]) = 0;
          not_used_[num_not_used_++] = pos;
          num_frules_--;
        }
        break;
      }
      prev_pos = pos;
      pos = next_[pos];
    }
  }
}

uint32_t HashTableRP::GetMaxOcc(){
  return max_occ_;
}

pair<uint32_t, uint32_t> HashTableRP::GetMaxPair(){
  return make_pair(get<0>(datas_[max_pair_]),
                   get<1>(datas_[max_pair_]));
}

void HashTableRP::ReHash(){
  // std::cout << __func__ << ": hash_index_ = " << hash_index_ << ", hash_size_ = " << hash_size_ << ", num_frules_ = " << num_frules_ << std::endl;
  hash_size_ = kPrimes[++hash_index_];
  num_not_used_ = hash_size_;
  CFunc::ResizeVec(datas_,
                   hash_size_);
  CFunc::ResizeVec(hash_table_,
                   hash_size_);
  CFunc::ResizeVec(next_,
                   hash_size_);
  CFunc::ResizeVec(not_used_,
                   hash_size_);
  for(size_t i = 0;
      i < hash_size_;
      i++){
    hash_table_[i] = hash_size_;
    next_[i] = hash_size_;
    not_used_[i] = hash_size_ - i - 1;
  }
  uint64_t tmp_num_frules = num_frules_;
  num_frules_ = 0;
  for(size_t i = 0;
      i < tmp_num_frules;
      i++){
    Add(get<0>(datas_[i]),
        get<1>(datas_[i]),
        get<2>(datas_[i]));
  }
}
  /*
void HashTableRP::PrintCRDict(){
  for(size_t i = 0; i < hash_size_; i++){
    std::cout << i << " " << var_[i] << " " << left_[i] << " " << right_[i] << " " << next_[i] << endl;
  }

  cout << "hash" << endl;
  for(size_t i = 0; i < hash_size_; i++){
    cout << i << " " << hash_table_[i] << endl;
  }
}
  */
void HashTableRP::InsertToHash(const uint64_t kInsertPos){
  uint64_t tmp = hash_table_[insert_ht_pos_];
  hash_table_[insert_ht_pos_] = kInsertPos;
  next_[kInsertPos] = tmp;
}

uint64_t HashTableRP::ByteSize() const {
  return sizeof(HashTableRP) +
    6 * hash_size_ * sizeof(uint32_t);
}

}
