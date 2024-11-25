/* hash_table_rp.hpp
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

#ifndef HASH_TABLE_RP_HPP_
#define HASH_TABLE_RP_HPP_

#include <iostream>
#include <cstdint>
#include <vector>
#include <tuple>
#include "constant_numbers.hpp"
#include "common_functions.hpp"


namespace solca_comp{
using Data = std::tuple<uint32_t, //left
                        uint32_t, //right
                        uint32_t>; // occ

class HashTableRP{//hash_table for Repair
private:
  std::vector<Data> datas_;
  std::vector<uint32_t> hash_table_;
  std::vector<uint32_t> next_;
  std::vector<uint32_t> not_used_;
  uint32_t hash_size_;
  uint32_t hash_index_;
  uint32_t num_frules_;
  uint64_t insert_ht_pos_;
  uint32_t num_not_used_;
  uint64_t max_pair_;
  uint64_t max_occ_;
  bool is_first_;
public:
  HashTableRP():datas_(),
                hash_table_(0),
                next_(0),
                not_used_(),
                hash_size_(0),
                hash_index_(0),
                num_frules_(0),
                insert_ht_pos_(0),
                num_not_used_(0),
                max_pair_(),
                max_occ_(),
                is_first_(){};
  ~HashTableRP(){};

  void Init();
  void Clear();
  void Add(const uint64_t kLeft,
           const uint64_t kRight,
           const uint64_t kOcc);
  void Subtract(const uint64_t kLeft,
                const uint64_t kRight,
                const uint64_t kOcc);
  uint32_t GetMaxOcc();
  std::pair<uint32_t, uint32_t> GetMaxPair();
  uint64_t ByteSize() const;
private:
  void InsertToHash(const uint64_t kInsertPos);
  // void PrintCRDict();
  void ReHash();

};
}
#endif // HASH_TABLE_RP_HPP_
