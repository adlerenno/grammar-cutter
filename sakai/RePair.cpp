#include "RePair.hpp"

using namespace std;

namespace slp_repair
{
  void NaiveReplace
  (
   std::vector<uint32_t> & v,
   const std::pair<uint32_t, uint32_t> & bigram,
   const uint32_t newId
   ) {
    uint32_t n = 0;
    for (uint32_t i = 0; i < v.size(); ++i) {
      if (i + 1 >= v.size()) {
        v[n++] = v[i];
      } else if (v[i] == bigram.first && v[i+1] == bigram.second) {
        v[n++] = newId;
        ++i;
      } else {
        v[n++] = v[i];
      }
    }
    v.resize(n);
  }


  void RePair::CalcVocc() {
    // std::cout << __func__ << std::endl;
    uint32_t bitIdx = skelton_.size() - 1;
    uint32_t leafIdx = leafVar_.size();
    uint32_t varIdx = GetNumNonnull() - 1; // idx of variable [0..numNonnull)
    vocc_[varIdx] = vocc_[varIdx+1] = vocc_[varIdx+2] = 1;

    stack_.push(varIdx);
    while (varIdx > 0) {
      auto parent = stack_.top();
      stack_.pop();
      if (!(parent & (UINT32_C(1) << 31))) {
        stack_.push(parent | (UINT32_C(1) << 31)); // mark that next time left child is processed
      } else {
        parent &= itmmti::ctcbits::UINTW_MAX(31); // remove the MSB bit that is used for flag
      }

      if (IsOpenBit(--bitIdx)) { // leaf node
        const uint32_t child = leafVar_[--leafIdx];
        if (child >= 256) {
          vocc_[child - 256] += vocc_[parent];
        }
      } else {
        const auto child = --varIdx;
        stack_.push(child);
        vocc_[child] += vocc_[parent];
      }
    }
    ClearStack();
  } 


  void RePair::KickOffRecomp()
  {
    // std::cout << __func__ << std::endl;
    CalcVocc();

    auto & vbeg = vstrBeg_[0];
    auto & vstr = vstr_[0];
    const uint32_t length = skelton_.size();
    uint32_t leaf_ri = 0; // read idx
    uint32_t leaf_wi = 0; // write idx
    uint32_t leafBit_wi = 0; // write idx
    PushVal(vbeg, 0);
    for (uint32_t i = 0; i < length; ++i) {
      if (IsOpenBit(i)) { // leaf
        const uint32_t var = leafVar_[leaf_ri++];
        if (var >= 256) {
          leafBit_.writeBit(true, leafBit_wi++);
          leafVar_[leaf_wi++] = var - 256;
        } else {
          leafBit_.writeBit(false, leafBit_wi++);
        }
        stack_.push(var);
      } else { // internal node
        const auto right = stack_.top();
        stack_.pop();
        const auto left = stack_.top();
        stack_.pop();
        stack_.push(256);
        if (left < 256) {
          if (left == right) {
            PushBlock(vstr, std::make_pair(left, 2));
          } else {
            PushVal(vstr, GetCharVal(left));
            if (right < 256) {
              PushVal(vstr, GetCharVal(right));
            }
          }
        } else if (right < 256) {
          PushVal(vstr, GetCharVal(right));
        }
        PushVal(vbeg, vstr.size());
      }
    }
    ClearStack();
    PushVal(vbeg, vstr.size()); // for string located at left of starting variable
    PushVal(vbeg, vstr.size()); // for string located at right of starting variable
    leafBit_.resize(leafBit_wi);
    leafVar_.resize(leaf_wi);

    BuildLRTBL(vbeg, vstr);
    ComputeFreq(vbeg, vstr);
  }


  void RePair::BuildLRTBL
  (
   const BlockVecT & vbeg,
   const BlockVecT & vstr
   ) {
    // std::cout << __func__  << std::endl;

    uint32_t inner_i = 0;
    uint32_t leafBit_i = 0;
    uint32_t leaf_i = 0;
    uint32_t length = skelton_.size();
    lrTableBeg_.resize(0);
    lrTable_.resize(0);

    PushVal(lrTableBeg_, 0);
    for (uint32_t i = 0; i < length; ++i) {
      if (IsOpenBit(i)) { // leaf
        if (leafBit_.readBit(leafBit_i++)) {
          stack_.push(leafVar_[leaf_i++]);
        } else {
          stack_.push(UINT32_MAX);
        }
      } else {
        auto right = stack_.top();
        stack_.pop();
        auto left = stack_.top();
        stack_.pop();
        stack_.push(inner_i);
        BuildLRTBL(vbeg, vstr, inner_i, left, right);
        PushVal(lrTableBeg_, lrTable_.size());
        // {//debug
        //   // lrTableBeg_.printStatistics(std::cout, true);
        //   // lrTable_.printStatistics(std::cout, true);
        //   const auto lch = LeftMC(inner_i);
        //   const auto rch = RightMC(inner_i);
        //   const auto lmblock = LeftMBlock(inner_i);
        //   const auto rmblock = RightMBlock(inner_i);
        //   std::cout << "inner_i = " << inner_i << ", lch = " << lch << ", rch = " << rch
        //             << ", beg = " << lrTableBeg_[inner_i] << ", end = " << lrTableBeg_[inner_i + 1] << std::endl;
        //   std::cout << "inner_i = " << inner_i << ", lmblock = (" << lmblock.first << ", " << lmblock.second
        //             << "), rmblock = (" << rmblock.first << ", " << rmblock.second << ")" << std::endl;
        // }
        ++inner_i;
      }
    }
    ClearStack();

    // free unused blocks
    ShrinkBV(lrTable_);
    ShrinkBV(lrTableBeg_);
  }


  void RePair::BuildLRTBL
  (
   const BlockVecT & vbeg,
   const BlockVecT & vstr,
   const uint32_t idx,
   const uint32_t left,
   const uint32_t right
   ) {
    // {// debug
    //   std::cout << __func__  << ": idx = " << idx << ", left = " << left << ", right = " << right << std::endl;
    //   if (left != UINT32_MAX) {
    //     auto lblock = LeftMBlock(left);
    //     std::cout << __func__ << ": lblock = (" << lblock.first << ", " << lblock.second << ")" << std::endl;
    //   }
    //   if (right != UINT32_MAX) {
    //     auto rblock = RightMBlock(right);
    //     std::cout << __func__ << ": rblock = (" << rblock.first << ", " << rblock.second << ")" << std::endl;
    //   }
    // }

    { // computing left most block
      std::pair<uint32_t, uint32_t> lblock = {UINT32_MAX, 0};
      if (left != UINT32_MAX) {
        lblock = LeftMBlock(left);
      }
      if (!(lblock.second & 1)) {
        const auto temp = LeftMBlockInStr(vbeg, vstr, idx);
        if (lblock.first == UINT32_MAX) {
          lblock.first = temp.first;
        }
        if (lblock.first == temp.first) {
          lblock.second += temp.second;
        } else if (temp.first != UINT32_MAX) {
          lblock.second |= 1;
        }
      }
      if (!(lblock.second & 1) && right != UINT32_MAX) {
        const auto temp = LeftMBlock(right);
        if (lblock.first == UINT32_MAX) {
          lblock.first = temp.first;
        }
        if (lblock.first == temp.first) {
          lblock.second += temp.second;
        } else {
          lblock.second |= 1;
        }
      }

      if (!(lblock.second & 1)) { // variable is repeating
        lblock.second >>= 1;
        lrTableBeg_[idx] = SetIsSB(lrTableBeg_[idx]);
        PushBlock(lrTable_, lblock);
        return;
      }
      lblock.second >>= 1;
      PushBlock(lrTable_, lblock);
    }
    
    { // computing right most block
      std::pair<uint32_t, uint32_t> rblock = {UINT32_MAX, 0};
      if (right != UINT32_MAX) {
        rblock = RightMBlock(right);
      }
      if (!(rblock.second & 1)) { // right variable is null or it is a single block
        const auto temp = RightMBlockInStr(vbeg, vstr, idx);
        if (rblock.first == UINT32_MAX) {
          rblock.first = temp.first;
        }
        if (rblock.first == temp.first) {
          rblock.second += temp.second;
        } else if (temp.first != UINT32_MAX) {
          rblock.second |= 1;
        }
      }
      if (!(rblock.second & 1) && left != UINT32_MAX) {
        const auto temp = RightMBlock(left);
        if (rblock.first == UINT32_MAX) {
          rblock.first = temp.first;
        }
        if (rblock.first == temp.first) {
          rblock.second += temp.second;
        }
      }

      rblock.second >>= 1;
      PushBlock(lrTable_, rblock);
    }
  }


  void RePair::ComputeFreq
  (
   const BlockVecT & vbeg,
   const BlockVecT & vstr
   ) {
    // std::cout << __func__ << std::endl;
    uint32_t inner_i = 0;
    uint32_t leafBit_i = 0;
    uint32_t leaf_i = 0;
    uint32_t length = skelton_.size();
    const uint32_t numNonnull = GetNumNonnull();

    if (numNonnull) {
      for (uint32_t i = 0; i < length; ++i) {
        if (IsOpenBit(i)) { // leaf
          if (leafBit_.readBit(leafBit_i++)) {
            stack_.push(leafVar_[leaf_i++]);
          } else {
            stack_.push(UINT32_MAX);
          }
        } else {
          auto right = stack_.top();
          stack_.pop();
          auto left = stack_.top();
          stack_.pop();
          stack_.push(inner_i);
          ComputeFreq(vbeg, vstr, inner_i, left, right, false, false);
          ++inner_i;
        }
      }
      vocc_[numNonnull] = vocc_[numNonnull+1] = 1;
      // top string
      ComputeFreq(vbeg, vstr, inner_i, UINT32_MAX, numNonnull - 1, true, false);
      ComputeFreq(vbeg, vstr, ++inner_i, numNonnull - 1, UINT32_MAX, false, true);
      const uint32_t val = lrTableBeg_[numNonnull - 1];
      if (IsSB(val)) { // special case: a repeating block spanning the starting variable
        auto block = ReadBlock(lrTable_, UnsetIsSB(val));
        auto temp = RightMBlockInStr(vbeg, vstr, numNonnull);
        if (block.first == temp.first) {
          block.second += (temp.second >> 1);
        }
        temp = LeftMBlockInStr(vbeg, vstr, numNonnull + 1);
        if (block.first == temp.first) {
          block.second += (temp.second >> 1);
        }
        if (block.second >= 2) {
          freqTable_.Add(block.first, block.first, block.second / 2);
        }
      }
    } else {
      ComputeFreq(vbeg, vstr, 0, UINT32_MAX, UINT32_MAX, true, true);
    }
    ClearStack();
  }


  void RePair::ComputeFreq
  (
   const BlockVecT & vbeg,
   const BlockVecT & vstr,
   const uint32_t idx,
   const uint32_t left,
   const uint32_t right,
   bool lCountFlag, // when it is 0, skip counting repeating bigram
   bool rCountFlag // when it is 0, skip counting repeating bigram
   ) {
    // std::cout << __func__ << ": idx = " << idx << ", left = " << left << ", right = " << right << ", lCountFlag = " << lCountFlag << ", rCountFlag = " << rCountFlag << std::endl;

    std::pair<uint32_t, uint32_t> block;
    const uint32_t vocc = vocc_[idx];
    uint32_t i = vbeg[idx];
    const uint32_t end = vbeg[idx+1];
    bool isBlockRead = false;
    if (left != UINT32_MAX) {
      isBlockRead = true;
      block = RightMBlock(left);
      if (block.second & 1) {
        block.second >>= 1;
        lCountFlag = true;
      }
      if (i < end) {
        const auto temp = ReadBlock(vstr, i);
        i += 1 + (temp.second > 1);
        if (block.first == temp.first) {
          block.second += temp.second;
        } else {
          if (lCountFlag && block.second >= 2) {
            freqTable_.Add(block.first, block.first, vocc * (block.second / 2));
          }
          freqTable_.Add(block.first, temp.first, vocc);
          block = temp;
          lCountFlag = true;
        }
      }
    } else if (i < end) {
      isBlockRead = true;
      block = ReadBlock(vstr, i);
      i += 1 + (block.second > 1);
    }

    while (i < end) {
      const auto temp = ReadBlock(vstr, i);
      i += 1 + (temp.second > 1);
      if (lCountFlag && block.second >= 2) {
        freqTable_.Add(block.first, block.first, vocc * (block.second / 2));
      }
      freqTable_.Add(block.first, temp.first, vocc);
      // std::cout << "block = (" << block.first << ", " << block.second << ")" << std::endl;
      // std::cout << "temp = (" << temp.first << ", " << temp.second << ")" << std::endl;
      // std::cout << "vocc = " << vocc << std::endl;
      block = temp;
      lCountFlag = true;
    }

    if (right != UINT32_MAX) {
      auto temp = LeftMBlock(right);
      if (temp.second & 1) {
        temp.second >>= 1;
        rCountFlag = true;
      }
      if (isBlockRead) {
        if (block.first == temp.first) {
          block.second += temp.second;
        } else {
          if (lCountFlag && block.second >= 2) {
            freqTable_.Add(block.first, block.first, vocc * (block.second / 2));
          }
          freqTable_.Add(block.first, temp.first, vocc);
          block = temp;
          lCountFlag = true;
        }
      } else {
        block = temp;
      }
    }
    if (lCountFlag && rCountFlag && block.second >= 2) {
      freqTable_.Add(block.first, block.first, vocc * (block.second / 2));
    }
  }


  void RePair::ReplaceNonRepeating
  (
   const uint8_t turn,
   const std::pair<uint32_t, uint32_t> & bigram
   ) {
    const uint32_t numNonnull = GetNumNonnull();
    auto & vbeg2 = vstrBeg_[1 - turn];
    auto & vstr2 = vstr_[1 - turn];
    PushVal(vbeg2, 0);

    if (numNonnull) {
      nonnullRank_.clear();
      uint32_t inner_i = 0;
      uint32_t leafBit_i = 0;
      uint32_t leaf_i = 0;
      uint32_t leafBit_wi = 0;
      uint32_t leaf_wi = 0;
      uint32_t skelton_wi = 0;
      uint32_t numNonnull = 0;
      const uint32_t length = skelton_.size();
      for (uint32_t i = 0; i < length; ++i) {
        skelton_.writeBit(skelton_.readBit(i), skelton_wi++);
        if (IsOpenBit(i)) { // leaf
          const bool lb = leafBit_.readBit(leafBit_i++);
          if (lb) { // The leaf is an internal node stored in leafVar_.
            const uint32_t var = leafVar_[leaf_i++];
            // leafVar_.printStatistics(true);
            // std::cout << __func__ << ": var = " << var << ", leaf_i = " << leaf_i << std::endl;
            if (nonnullRank_.readBit(var)) {
              leafVar_[leaf_wi++] = nonnullRank_.rank_1(var) - 1;
              leafBit_.writeBit(true, leafBit_wi++);
            } else {
              leafBit_.writeBit(false, leafBit_wi++);
            }
            stack_.push(var);
          } else {
            leafBit_.writeBit(false, leafBit_wi++);
            stack_.push(UINT32_MAX);
          }
        } else {
          const auto right = stack_.top();
          stack_.pop();
          const auto left = stack_.top();
          stack_.pop();
          stack_.push(inner_i);
          const auto beg = vstr2.size();
          ReplaceNonRepeating(turn, bigram, inner_i, left, right);
          if ((left != UINT32_MAX && nonnullRank_.readBit(left)) ||
              (right != UINT32_MAX && nonnullRank_.readBit(right)) ||
              beg < vstr2.size()) { // nonnull
            PushVal(vbeg2, vstr2.size());
            vocc_[numNonnull++] = vocc_[inner_i];
            nonnullRank_.appendBit(true);
          } else { // variable becomes null
            nonnullRank_.appendBit(false);
            skelton_wi -= 2;
            leafBit_wi -= 2;
            leafBit_.writeBit(false, leafBit_wi++);
            // std::cout << __func__ << ": skelton_wi = " << skelton_wi << ", leafBit_wi = " << leafBit_wi << std::endl;
          }
          ++inner_i;
        }
      }
      skelton_.resize(skelton_wi);
      leafVar_.resize(leaf_wi);
      leafBit_.resize(leafBit_wi);
      // top string
      const auto start = inner_i - 1;
      if (numNonnull) {
        ReplaceNonRepeating(turn, bigram, inner_i, UINT32_MAX - 1, start);
        PushVal(vbeg2, vstr2.size());
        ReplaceNonRepeating(turn, bigram, inner_i+1, start, UINT32_MAX - 1);
        PushVal(vbeg2, vstr2.size());
      } else {
        auto & vbeg1 = vstrBeg_[turn];
        auto & vstr1 = vstr_[turn];
        uint32_t llen = 0;
        uint32_t nlen = 0;
        for (uint32_t j = 0; j < 2; ++j) {
          for (uint32_t i = vbeg1[inner_i + j]; i < vbeg1[inner_i + j + 1]; ++i) {
            auto temp = ReadBlock(vstr1, i);
            i += 1 + (temp.second > 1);
            if (temp.first == bigram.first) {
              llen += temp.second;
            } else if (llen == 1) {
              if (temp.first == bigram.second) {
                ++nlen;
                if (temp.second > 1) {
                  PushBlock(vstr2, std::make_pair(repairVarId_, nlen));
                  PushBlock(vstr2, std::make_pair(bigram.second, temp.second - 1));
                  nlen = 0;
                }
              } else {
                if (nlen) {
                  PushBlock(vstr2, std::make_pair(repairVarId_, nlen));
                  nlen = 0;
                }
                PushVal(vstr2, GetCharVal(bigram.first));
                PushBlock(vstr2, temp);
              }
              llen = 0;
            } else {
              if (nlen) {
                PushBlock(vstr2, std::make_pair(repairVarId_, nlen));
                nlen = 0;
              }
              if (llen) { // llen > 1 in this case
                PushBlock(vstr2, std::make_pair(bigram.first, llen - (temp.first == bigram.second)));
                if (temp.first == bigram.second) {
                  if (temp.second > 1) {
                    PushBlock(vstr2, std::make_pair(repairVarId_, 1));
                    --temp.second;
                    nlen = 0;
                  } else {
                    nlen = 1;
                  }
                }
                llen = 0;
              }
              PushBlock(vstr2, temp);
            }
          }

          if (j == 0) {
            const uint32_t lc = LeftMC(start);
            const uint32_t rc = RightMC(start);
            if (lc == bigram.first) {
              ++llen;
            } else { // lc must be bigram.second
              if (llen == 1) {
                ++nlen;
              } else {
                if (nlen) {
                  PushBlock(vstr2, std::make_pair(repairVarId_, nlen));
                  nlen = 0;
                }
                if (llen) {
                  nlen = 1;
                  PushBlock(vstr2, std::make_pair(bigram.first, llen - 1));
                }
              }
              llen = (rc == bigram.first);
            }
          }
        }
        if (nlen) {
          PushBlock(vstr2, std::make_pair(repairVarId_, nlen));
        }
        if (llen) {
          PushBlock(vstr2, std::make_pair(bigram.first, llen));
        }
        PushVal(vbeg2, vstr2.size());
        PoolBlocks(vstr1, vbeg1[inner_i], vbeg1[inner_i + 2]);
        PoolBlocks(vbeg1, inner_i, inner_i+2);
      }
    } else {
      ReplaceNonRepeating(turn, bigram, 0, UINT32_MAX - 1, UINT32_MAX - 1);
    }
    PoolBlocks(vstr_[turn], vstr_[turn].size(), vstr_[turn].capacity());
    PoolBlocks(vstrBeg_[turn], vstrBeg_[turn].size(), vstrBeg_[turn].capacity());
    ClearStack();
  }


  void RePair::ReplaceNonRepeating
  (
   const uint8_t turn,
   const std::pair<uint32_t, uint32_t> & bigram,
   const uint32_t idx,
   const uint32_t left,
   const uint32_t right
   ) {
    // std::cout << __func__ << ": turn = " << (int)turn << ", idx = " << idx << ", left = " << left << ", right = " << right << std::endl;
    auto & vstr1 = vstr_[turn];
    auto & vbeg1 = vstrBeg_[turn];
    auto & vstr2 = vstr_[1 - turn];

    uint32_t i = vbeg1[idx];
    const uint32_t end = vbeg1[idx+1];
    uint32_t llen = (left < UINT32_MAX - 1) && (bigram.first == RightMC(left));
    uint32_t rlen = 0;
    uint32_t nlen = 0;
    bool isPopoutToLeft = (left == UINT32_MAX);

    while (i < end) {
      const auto temp = ReadBlock(vstr1, i);
      i += 1 + (temp.second > 1);
      if (temp.first == bigram.first) {
        llen += temp.second;
        if (rlen) {
          PushBlock(vstr2, std::make_pair(bigram.second, rlen));
          rlen = 0;
        } else if (nlen && llen > 1) {
          PushBlock(vstr2, std::make_pair(repairVarId_, nlen));
          nlen = 0;
        }
        isPopoutToLeft = false;
      } else if (temp.first == bigram.second) {
        if (llen) {
          rlen = temp.second - 1;
          ++nlen;
          if (--llen) {
            PushBlock(vstr2, std::make_pair(bigram.first, llen));
            llen = 0;
          }
          if (nlen && rlen) {
            PushBlock(vstr2, std::make_pair(repairVarId_, nlen));
            nlen = 0;
          }
        } else {
          if (nlen) {
            PushBlock(vstr2, std::make_pair(repairVarId_, nlen));
            nlen = 0;
          }
          rlen += temp.second - isPopoutToLeft;
        }
        isPopoutToLeft = false;
      } else {
        if (nlen) {
          PushBlock(vstr2, std::make_pair(repairVarId_, nlen));
          nlen = 0;
        }
        if (llen) {
          PushBlock(vstr2, std::make_pair(bigram.first, llen));
          llen = 0;
        }
        if (rlen) {
          PushBlock(vstr2, std::make_pair(bigram.second, rlen));
          rlen = 0;
        }
        PushBlock(vstr2, temp);
        isPopoutToLeft = false;
      }
    }

    rlen += (right < UINT32_MAX - 1) && (bigram.second == LeftMC(right));
    rlen -= (isPopoutToLeft && rlen);
    if (rlen == 1) {
      if (llen) {
        if (--llen) {
          PushBlock(vstr2, std::make_pair(bigram.first, llen));
        }
        PushBlock(vstr2, std::make_pair(repairVarId_, nlen + 1));
      } else {
        if (nlen) {
          PushBlock(vstr2, std::make_pair(repairVarId_, nlen));
        }
        PushVal(vstr2, GetCharVal(bigram.second));
      }
    } else {
      if (nlen) {
        PushBlock(vstr2, std::make_pair(repairVarId_, nlen));
      }
      if (llen) {
        llen -= (right == UINT32_MAX); // popout to right
        if (llen) {
          PushBlock(vstr2, std::make_pair(bigram.first, llen));
        }
      }
      if (rlen) {
        PushBlock(vstr2, std::make_pair(bigram.second, rlen));
      }
    }

    PoolBlocks(vstr1, vbeg1[idx], end);
    PoolBlocks(vbeg1, idx, idx+1);
  }


  void RePair::ReplaceRepeating
  (
   const uint8_t turn,
   const uint32_t ch
   ) {
    const uint32_t numNonnull = GetNumNonnull();
    auto & vbeg2 = vstrBeg_[1 - turn];
    auto & vstr2 = vstr_[1 - turn];
    PushVal(vbeg2, 0);

    if (numNonnull) {
      nonnullRank_.clear();
      uint32_t inner_i = 0;
      uint32_t leafBit_i = 0;
      uint32_t leaf_i = 0;
      uint32_t leafBit_wi = 0;
      uint32_t leaf_wi = 0;
      uint32_t skelton_wi = 0;
      uint32_t numNonnull = 0;
      const uint32_t length = skelton_.size();
      for (uint32_t i = 0; i < length; ++i) {
        skelton_.writeBit(skelton_.readBit(i), skelton_wi++);
        if (IsOpenBit(i)) { // leaf
          const bool lb = leafBit_.readBit(leafBit_i++);
          if (lb) {
            const uint32_t var = leafVar_[leaf_i++];
            if (nonnullRank_.readBit(var)) {
              leafVar_[leaf_wi++] = nonnullRank_.rank_1(var) - 1;
              leafBit_.writeBit(true, leafBit_wi++);
            } else {
              leafBit_.writeBit(false, leafBit_wi++);
            }
            stack_.push(var);
          } else {
            leafBit_.writeBit(false, leafBit_wi++);
            stack_.push(UINT32_MAX);
          }
        } else {
          const auto right = stack_.top();
          stack_.pop();
          const auto left = stack_.top();
          stack_.pop();
          stack_.push(inner_i);
          const auto beg = vstr2.size();
          ReplaceRepeating(turn, ch, inner_i, left, right);
          if ((left != UINT32_MAX && nonnullRank_.readBit(left)) ||
              (right != UINT32_MAX && nonnullRank_.readBit(right)) ||
              beg < vstr2.size()) { // nonnull
            PushVal(vbeg2, vstr2.size());
            vocc_[numNonnull++] = vocc_[inner_i];
            nonnullRank_.appendBit(true);
          } else { // variable becomes null
            nonnullRank_.appendBit(false);
            skelton_wi -= 2;
            leafBit_wi -= 2;
            leafBit_.writeBit(false, leafBit_wi++);
          }
          ++inner_i;
        }
      }
      skelton_.resize(skelton_wi);
      leafVar_.resize(leaf_wi);
      leafBit_.resize(leafBit_wi);
      // top string
      const auto start = inner_i - 1;
      if (numNonnull) {
        ReplaceRepeating(turn, ch, inner_i, UINT32_MAX - 1, start);
        PushVal(vbeg2, vstr2.size());
        ReplaceRepeating(turn, ch, inner_i+1, start, UINT32_MAX - 1);
        PushVal(vbeg2, vstr2.size());
      } else {
        auto & vbeg1 = vstrBeg_[turn];
        auto & vstr1 = vstr_[turn];
        uint32_t chlen = 0;
        for (uint32_t j = 0; j < 2; ++j) {
          for (uint32_t i = vbeg1[inner_i + j]; i < vbeg1[inner_i + j + 1]; ++i) {
            const auto temp = ReadBlock(vstr1, i);
            i += 1 + (temp.second > 1);
            if (temp.first != ch) {
              if (chlen) {
                PushBlock(vstr2, std::make_pair(repairVarId_, chlen / 2));
                if (chlen % 2) {
                  PushVal(vstr2, GetCharVal(ch));
                }
              }
              PushBlock(vstr2, temp);
              chlen = 0;
            } else {
              chlen += temp.second;
            }
          }
          if (j == 0) {
            const auto temp = LeftMBlock(start);
            chlen += (temp.second >> 1);
          }
        }
        if (chlen) {
          PushBlock(vstr2, std::make_pair(repairVarId_, chlen / 2));
          if (chlen % 2) {
            PushVal(vstr2, GetCharVal(ch));
          }
        }
        PushVal(vbeg2, vstr2.size());
        PoolBlocks(vstr1, vbeg1[inner_i], vbeg1[inner_i + 2]);
        PoolBlocks(vbeg1, inner_i, vbeg1.size());
      }
    } else {
      ReplaceRepeating(turn, ch, 0, UINT32_MAX - 1, UINT32_MAX - 1);
      PushVal(vbeg2, vstr2.size());
    }
    PoolBlocks(vstr_[turn], vstr_[turn].size(), vstr_[turn].capacity());
    PoolBlocks(vstrBeg_[turn], vstrBeg_[turn].size(), vstrBeg_[turn].capacity());
    ClearStack();
  }


  void RePair::ReplaceRepeating
  (
   const uint8_t turn,
   const uint32_t ch,
   const uint32_t idx,
   const uint32_t left,
   const uint32_t right
   ) {
    // std::cout << __func__ << ": turn = " << (int)turn << ", idx = " << idx << ", left = " << left << ", right = " << right << std::endl;
    auto & vstr1 = vstr_[turn];
    auto & vbeg1 = vstrBeg_[turn];
    auto & vstr2 = vstr_[1 - turn];

    uint32_t i = vbeg1[idx];
    const uint32_t end = vbeg1[idx+1];
    uint32_t chlen = UINT32_MAX; // UNT32_MAX indicates that it is left-open
    if (left < UINT32_MAX - 1) {
      const auto temp = RightMBlock(left);
      if (temp.first != ch) {
        chlen = 0; // left-close
      } else if (temp.second & 1) { // stopping
        chlen = (temp.second >> 1);
      }
    } else if (left == UINT32_MAX - 1) {
      chlen = 0; // left-close
    }

    if (chlen) { // left-open or "chlen" ch are popping in from left var
      if (i < end) {
        const auto temp = ReadBlock(vstr1, i);
        i += 1 + (temp.second > 1);
        if (temp.first == ch) {
          if (chlen != UINT32_MAX) {
            chlen += temp.second;
          } else if (i < end) { // left-open and vstr contains another char other than "ch"
            chlen = 0; // ignore the leftmost block of "ch"
          }
        } else {
          if (chlen != UINT32_MAX) {
            if (chlen >= 2) {
              PushBlock(vstr2, std::make_pair(repairVarId_, chlen / 2));
            }
            if (chlen % 2) {
              PushVal(vstr2, GetCharVal(ch));
            }
          }
          PushBlock(vstr2, temp);
          chlen = 0;
        }
      }
    }

    while (i < end) {
      const auto temp = ReadBlock(vstr1, i);
      i += 1 + (temp.second > 1);
      if (temp.first != ch) {
        if (chlen) {
          if (chlen >= 2) {
            PushBlock(vstr2, std::make_pair(repairVarId_, chlen / 2));
          }
          if (chlen % 2) {
            PushVal(vstr2, GetCharVal(ch));
          }
        }
        PushBlock(vstr2, temp);
        chlen = 0;
      } else {
        chlen = temp.second;
      }
    }

    if (chlen != UINT32_MAX && right != UINT32_MAX) {
      if (right < UINT32_MAX - 1) {
        const auto temp = LeftMBlock(right);
        if (temp.first == ch) {
          if (temp.second & 1) { // right-close
            chlen += (temp.second >> 1);
          } else { //
            chlen = 0;
          }
        }
      }
      if (chlen) {
        if (chlen >= 2) {
          PushBlock(vstr2, std::make_pair(repairVarId_, chlen / 2));
        }
        if (chlen % 2) {
          PushVal(vstr2, GetCharVal(ch));
        }
      }
    }

    PoolBlocks(vstr1, vbeg1[idx], end);
    PoolBlocks(vbeg1, idx, idx+1);
  }


  void RePair::InitDictRule()
  {
    dict_rule_ = (RULE*)malloc(sizeof(RULE)*INIT_DICTIONARY_SIZE);
    for (uint32_t i = 0; i < INIT_DICTIONARY_SIZE; i++) {
      dict_rule_[i].left = DUMMY_CODE;
      dict_rule_[i].right = DUMMY_CODE;
    }
    for (uint32_t i = 0; i < CHAR_SIZE+1; i++) {
      dict_rule_[i].left  = (CODE)i;
      dict_rule_[i].right = DUMMY_CODE;
    }
  }


  void RePair::RePairRecompression(size_t input_size, uint32_t turn_point)
  {
    std::cout << "RePair 1st step..." << std::endl;
    auto start = std::chrono::system_clock::now();

    std::cout << "turning textsize : " << input_size/turn_point << "(bytes)" << std::endl;

    InitDictRule();

    // PrintStatus(0); // debug
    // leafVar_.printStatistics(true);
    freqTable_.Init();
    KickOffRecomp();

    uint8_t turn = 0; // switching between 0 and 1 by turns
    // debug
    std::vector<uint32_t> naive_str;

    const size_t step = 100; // print status every step characters
    size_t next_print_step = step;
    size_t cur_step = 0;
    size_t naive_text_size = input_size;

    while (freqTable_.GetMaxOcc() > 1 and naive_text_size > input_size/turn_point) {
      const auto bigram = freqTable_.GetMaxPair();
      dict_rule_[cur_step+257].left = bigram.first;
      dict_rule_[cur_step+257].right = bigram.second;
      naive_text_size -= freqTable_.GetMaxOcc();
      cur_step++;

      if (cur_step == next_print_step) { // debug
        // PrintStatus(turn);
        auto end = std::chrono::system_clock::now();      
        auto dur = end - start;        
        auto sec = std::chrono::duration_cast<std::chrono::seconds>(dur).count();
        cout << cur_step << ": " << sec << " sec, numNonnull = " << GetNumNonnull() 
             << ", seqsize = " << vstr_[turn].size() << ", textsize = " << naive_text_size << ": ";
        std::cout << repairVarId_ << " -> (" << bigram.first << ", "
                  << bigram.second << "), freq = " << freqTable_.GetMaxOcc() << std::endl;
        next_print_step = cur_step + step;
      }
      freqTable_.Clear();

      if (bigram.first != bigram.second) {
        ReplaceNonRepeating(turn, bigram);
      } else {
        ReplaceRepeating(turn, bigram.first);
      }

      { // clear old blocks
        // const auto size1 = vstrBeg_[turn].size();
        // if (size1 % BlockVecT::kBlockSize) {
        //   PoolBlocks(vstrBeg_[turn], size1 - 1, size1 / BlockVecT::kBlockSize * BlockVecT::kBlockSize + BlockVecT::kBlockSize);
        // }
        // const auto size2 = vstr_[turn].size();
        // if (size2 % BlockVecT::kBlockSize) {
        //   PoolBlocks(vstr_[turn], size2 - 1, size2 / BlockVecT::kBlockSize * BlockVecT::kBlockSize + BlockVecT::kBlockSize);
        // }
        vstrBeg_[turn].freeBlock();
        vstr_[turn].freeBlock();
      }

      ++repairVarId_;
      turn = 1 - turn;
      BuildLRTBL(vstrBeg_[turn], vstr_[turn]);

      // { // debug correctness
      //   // if (repairVarId_ == 262) {
      //   //   const uint32_t numNonnull = GetNumNonnull();
      //   //   std::vector<uint32_t> lrIds(2 * (numNonnull + 2));
      //   //   lrIds.resize(2 * (numNonnull + 2));
      //   //   CalcLrIds(lrIds);
      //   //   std::vector<uint32_t> len(numNonnull + 2);
      //   //   len.resize(numNonnull + 2);
      //   //   CalcLen(turn, lrIds, len);
      //   //   for (uint32_t i = 571200; i < 571200 + 10; ++i) {
      //   //     const uint32_t left = lrIds[2 * i];
      //   //     const uint32_t right = lrIds[2 * i + 1];
      //   //     PrintStatus(turn, i, left, right, 0, 0);
      //   //     if (left != UINT32_MAX) {
      //   //       std::cout << "left  ";
      //   //       PrintStatus(turn, left, lrIds[2 * left], lrIds[2 * left + 1], 0, 0);
      //   //     }
      //   //     if (right != UINT32_MAX) {
      //   //       std::cout << "right ";
      //   //       PrintStatus(turn, right, lrIds[2 * right], lrIds[2 * right + 1], 0, 0);
      //   //     }
      //   //   }
      //   // }

      //   // {// print naive before
      //   //   PrintStatus(turn);
      //   //   for (uint32_t i = 0; i < naive_str.size(); ++i) {
      //   //     std::cout << naive_str[i] << ", ";
      //   //   }
      //   //   std::cout << "naive before" << std::endl;
      //   // }
      //   {// replace naively
      //     NaiveReplace(naive_str, bigram, repairVarId_-1);
      //   }
      //   // {// print naive after
      //   //   for (uint32_t i = 0; i < naive_str.size(); ++i) {
      //   //     std::cout << naive_str[i] << ", ";
      //   //   }
      //   //  std::cout << "naive after" << std::endl;
      //   // }
      //   std::vector<uint32_t> gra_str;
      //   CalcSequence_Debug(turn, gra_str);
      //   // {// print gra_str
      //   //   for (uint32_t i = 0; i < gra_str.size(); ++i) {
      //   //     std::cout << gra_str[i] << ", ";
      //   //   }
      //   //   std::cout << "grammar" << std::endl;
      //   // }
      //   for (uint32_t i = 0; i < naive_str.size(); ++i) {
      //     if (naive_str[i] != gra_str[i]) {
      //       std::cout << "error: i = " << i << ", naive_str = " << naive_str[i] << ", gra_str = " << gra_str[i] << std::endl;
      //       for (uint32_t j = (i < 100) ? 0 : i - 100; j < i; ++j) {
      //         std::cout << naive_str[j] << ", ";
      //       }
      //       std::cout << "[" << naive_str[i] << "], ";
      //       uint32_t end = (naive_str.size() < i + 100) ? naive_str.size() : i + 100;
      //       for (uint32_t j = i + 1; j < end; ++j) {
      //         std::cout << naive_str[j] << ", ";
      //       }
      //       std::cout << "naive" << std::endl;
      //       for (uint32_t j = (i < 100) ? 0 : i - 100; j < i; ++j) {
      //         std::cout << gra_str[j] << ", ";
      //       }
      //       std::cout << "[" << gra_str[i] << "], ";
      //       end = (gra_str.size() < i + 100) ? gra_str.size() : i + 100;
      //       for (uint32_t j = i + 1; j < end; ++j) {
      //         std::cout << gra_str[j] << ", ";
      //       }
      //       std::cout << "grammar" << std::endl;
      //       {
      //         const uint32_t numNonnull = GetNumNonnull();
      //         std::vector<uint32_t> lrIds(2 * (numNonnull + 2));
      //         lrIds.resize(2 * (numNonnull + 2));
      //         CalcLrIds(lrIds);
      //         std::vector<uint32_t> len(numNonnull + 2);
      //         len.resize(numNonnull + 2);
      //         CalcLen(turn, lrIds, len);
      //         {
      //           const auto varpos = CalcVar(i, lrIds, len);
      //           const uint32_t left = lrIds[2 * varpos.first];
      //           const uint32_t right = lrIds[2 * varpos.first + 1];
      //           std::cout << "pos in str = " << varpos.second << std::endl;
      //           PrintStatus(turn, varpos.first, left, right, 0, 0);
      //           if (left != UINT32_MAX) {
      //             PrintStatus(turn, left, lrIds[2 * left], lrIds[2 * left + 1], 0, 0);
      //           }
      //           if (right != UINT32_MAX) {
      //             PrintStatus(turn, right, lrIds[2 * right], lrIds[2 * right + 1], 0, 0);
      //           }
      //         }
      //         {
      //           const auto varpos = CalcVar(i - 1, lrIds, len);
      //           const uint32_t left = lrIds[2 * varpos.first];
      //           const uint32_t right = lrIds[2 * varpos.first + 1];
      //           std::cout << "pos in str = " << varpos.second << std::endl;
      //           PrintStatus(turn, varpos.first, left, right, 0, 0);
      //         }
      //       }
      //       // PrintStatus(turn);
      //       exit(1);
      //     }
      //   }
      // }

      ComputeFreq(vstrBeg_[turn], vstr_[turn]);
    }
    freqTable_.Clear();

    naive_seq_ = (SEQ*)malloc(sizeof(SEQ)*naive_text_size);
    CalcSequence(turn);

    auto end = std::chrono::system_clock::now();      
    auto dur = end - start;      
    auto sec = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
    std::cout << "SLP to RePair 1st step done in " << sec/1000.0 << " sec" << std::endl << std::endl;
  }   

  void RePair::PrintStatus
  (
   const uint8_t turn,
   const uint32_t nodeId,
   const uint32_t left,
   const uint32_t right,
   const uint32_t lparenth,
   const uint32_t rparenth
   ) const {
    static char oc[2];
    oc[kOpen] = 'l';
    oc[kClose] = 'i';

    std::cout << nodeId << " -> ["
              << left << "(" << oc[lparenth] << "), "
              << right << "(" << oc[rparenth] << ")]"
              << ", vocc = " << vocc_[nodeId];
    if (lrTableBeg_.size()) {
      const auto lmb = LeftMBlock(nodeId);
      if (!(lmb.second & 1)) {
        std::cout << ", SB = " << lmb.first << "^" << (lmb.second >> 1);
      } else {
        std::cout << ", LMB = " << lmb.first << "^" << (lmb.second >> 1) << "; ";
        for (uint32_t j = vstrBeg_[turn][nodeId]; j < vstrBeg_[turn][nodeId+1]; ) {
          const auto temp = ReadBlock(vstr_[turn], j);
          j += 1 + (temp.second > 1);
          std::cout << temp.first << "^" << temp.second << ", ";
        }
        const auto rmb = RightMBlock(nodeId);
        std::cout << "RMB = " << rmb.first << "^" << (rmb.second >> 1);
      }
    }
    std::cout << std::endl;
  }


  void RePair::PrintStatus(const uint8_t turn) const {
    const auto numNonnull = GetNumNonnull();
    static char oc[2];
    oc[kOpen] = 'l';
    oc[kClose] = 'i';
    std::cout << __func__  << ", numNonnull = " << numNonnull << std::endl;

    std::stack<uint32_t, std::vector<uint32_t>> sv;
    std::stack<uint8_t, std::vector<uint8_t>> soc; // stack of open close
    uint32_t inner_i = 0;
    uint32_t leafBit_i = 0;
    uint32_t leaf_i = 0;
    uint32_t length = skelton_.size();
    
    for (uint32_t i = 0; i < length; ++i) {
      if (IsOpenBit(i)) { // leaf
        soc.push(kOpen);
        if (leafBit_.size() == 0 || leafBit_.readBit(leafBit_i++)) {
          sv.push(leafVar_[leaf_i++]);
        } else {
          sv.push(UINT32_MAX);
        }
      } else {
        auto right = sv.top();
        auto rightP = soc.top();
        sv.pop();
        soc.pop();
        auto left = sv.top();
        auto leftP = soc.top();
        sv.pop();
        soc.pop();
        sv.push(inner_i);
        soc.push(kClose);
        PrintStatus(turn, inner_i, left, right, leftP, rightP);
        ++inner_i;
      }
    }
    std::cout << "topleft string: vocc = " << vocc_[numNonnull] << ": ";
    for (uint32_t j = vstrBeg_[turn][numNonnull]; j < vstrBeg_[turn][numNonnull+1]; ) {
      const auto temp = ReadBlock(vstr_[turn], j);
      j += 1 + (temp.second > 1);
      std::cout << temp.first << "^" << temp.second << ", ";
    }
    std::cout << std::endl << "topright string: vocc = " << vocc_[numNonnull+1] << ": ";
    for (uint32_t j = vstrBeg_[turn][numNonnull+1]; j < vstrBeg_[turn][numNonnull+2]; ) {
      const auto temp = ReadBlock(vstr_[turn], j);
      j += 1 + (temp.second > 1);
      std::cout << temp.first << "^" << temp.second << ", ";
    }
    std::cout << std::endl;

    std::cout << "print skelton_" << std::endl;
    skelton_.printStatistics(true);
    std::cout << "print leafBit_" << std::endl;
    leafBit_.printStatistics(true);
    std::cout << "print leafVar_" << std::endl;
    leafVar_.printStatistics(true);
    std::cout << "print vstrBeg_" << std::endl;
    vstrBeg_[turn].printStatistics(std::cout, true); 
    std::cout << "print vstr_" << std::endl;
    vstr_[turn].printStatistics(std::cout, true);
    std::cout << "print lrTableBeg_" << std::endl;
    lrTableBeg_.printStatistics(std::cout, true);
    std::cout << "print lrTable_" << std::endl;
    lrTable_.printStatistics(std::cout, true);
  }


  void RePair::CalcSequence
  (
   const uint8_t turn,
   std::vector<uint32_t> & aux,
   uint32_t & processedNum,
   const uint32_t nodeId
   ) {
    const BlockVecT & vbeg = vstrBeg_[turn];
    const BlockVecT & vstr = vstr_[turn];

    const auto left = aux[2 * nodeId];
    aux[2 * nodeId] = idx_seq_;
    if (left != UINT32_MAX) {
      if (left+1 <= processedNum) {
        // Copy the text from previous occurrence.
        for (uint32_t j = aux[2 * left]; j < aux[2 * left + 1]; ++j) {
          naive_seq_[idx_seq_].code = naive_seq_[j].code;
          naive_seq_[idx_seq_].next = DUMMY_POS;
          naive_seq_[idx_seq_].prev = DUMMY_POS;
          idx_seq_++;
        }
      } else {
        CalcSequence(turn, aux, processedNum, left);
      }
    }
    {
      uint32_t k = vbeg[nodeId];
      while (k < vbeg[nodeId+1]) {
        const auto temp = ReadBlock(vstr, k);
        k += 1 + (temp.second > 1);
        for (uint32_t j = 0; j < temp.second; ++j) {
          naive_seq_[idx_seq_].code = temp.first;
          naive_seq_[idx_seq_].next = DUMMY_POS;
          naive_seq_[idx_seq_].prev = DUMMY_POS;
          idx_seq_++;
        }
      }
    }
    const uint32_t right = aux[2 * nodeId + 1];
    if (right != UINT32_MAX) {
      if (right+1 <= processedNum) {
        // Copy the text from previous occurrence.
        for (uint32_t j = aux[2 * right]; j < aux[2 * right + 1]; ++j) {
          naive_seq_[idx_seq_].code = naive_seq_[j].code;
          naive_seq_[idx_seq_].next = DUMMY_POS;
          naive_seq_[idx_seq_].prev = DUMMY_POS;
          idx_seq_++;
        }
      } else {
        CalcSequence(turn, aux, processedNum, right);
      }
    }
    aux[2 * nodeId + 1] = idx_seq_;
    ++processedNum;
  }


  void RePair::CalcSequence
  (
   const uint8_t turn
   ) {
    const BlockVecT & vbeg = vstrBeg_[turn];
    const BlockVecT & vstr = vstr_[turn];

    const uint32_t numNonnull = GetNumNonnull();
    {
      uint32_t k = vbeg[numNonnull];
      while (k < vbeg[numNonnull+1]) {
        const auto temp = ReadBlock(vstr, k);
        k += 1 + (temp.second > 1);
        for (uint32_t j = 0; j < temp.second; ++j) {
          naive_seq_[idx_seq_].code = temp.first;
          naive_seq_[idx_seq_].next = DUMMY_POS;
          naive_seq_[idx_seq_].prev = DUMMY_POS;
          idx_seq_++;
        }
      }
    }

    if (numNonnull) {
      // std::cout << __func__ << ": numNonnull = " << numNonnull << std::endl;
      uint32_t inner_i = 0;
      uint32_t leafBit_i = 0;
      uint32_t leaf_i = 0;
      std::vector<uint32_t> aux(2 * numNonnull);
      aux.resize(2 * numNonnull);
      for (uint32_t i = 0; i < skelton_.size(); ++i) {
        if (IsOpenBit(i)) { // leaf
          if (leafBit_.readBit(leafBit_i++)) {
            stack_.push(leafVar_[leaf_i++]);
          } else {
            stack_.push(UINT32_MAX);
          }
        } else {
          auto right = stack_.top();
          stack_.pop();
          auto left = stack_.top();
          stack_.pop();
          aux[2 * inner_i] = left;
          aux[2 * inner_i + 1] = right;
          stack_.push(inner_i);
          ++inner_i;
        }
      }
      ClearStack();

      uint32_t processedNum = 0;
      CalcSequence(turn, aux, processedNum, numNonnull-1);

      uint32_t k = vbeg[numNonnull+1];
      while (k < vbeg[numNonnull+2]) {
        const auto temp = ReadBlock(vstr, k);
        k += 1 + (temp.second > 1);
        for (uint32_t j = 0; j < temp.second; ++j) {
          naive_seq_[idx_seq_].code = temp.first;
          naive_seq_[idx_seq_].next = DUMMY_POS;
          naive_seq_[idx_seq_].prev = DUMMY_POS;
          idx_seq_++;
        }
      }
    }
  }


  // void RePair::CalcSequence
  // (
  //  const uint8_t turn,
  //  std::vector<uint32_t> & v
  //  ) {
  //   const BlockVecT & vbeg = vstrBeg_[turn];
  //   const BlockVecT & vstr = vstr_[turn];

  //   const uint32_t numNonnull = GetNumNonnull();
  //   {
  //     uint32_t k = vbeg[numNonnull];
  //     while (k < vbeg[numNonnull+1]) {
  //       const auto temp = ReadBlock(vstr, k);
  //       k += 1 + (temp.second > 1);
  //       for (uint32_t j = 0; j < temp.second; ++j) {
  //         v.push_back(temp.first);
  //       }
  //     }
  //   }

  //   if (numNonnull) {
  //     // std::cout << __func__ << ": numNonnull = " << numNonnull << std::endl;
  //     uint32_t inner_i = 0;
  //     uint32_t leafBit_i = 0;
  //     uint32_t leaf_i = 0;
  //     std::vector<uint32_t> aux(2 * numNonnull);
  //     aux.resize(2 * numNonnull);
  //     for (uint32_t i = 0; i < skelton_.size(); ++i) {
  //       if (IsOpenBit(i)) { // leaf
  //         if (leafBit_.readBit(leafBit_i++)) {
  //           stack_.push(leafVar_[leaf_i++]);
  //         } else {
  //           stack_.push(UINT32_MAX);
  //         }
  //       } else {
  //         auto right = stack_.top();
  //         stack_.pop();
  //         auto left = stack_.top();
  //         stack_.pop();
  //         aux[2 * inner_i] = left;
  //         aux[2 * inner_i + 1] = right;
  //         stack_.push(inner_i);
  //         ++inner_i;
  //       }
  //     }
  //     ClearStack();

  //     // for (uint32_t i = 0; i < aux.size(); ++i) {
  //     //   std::cout << aux[i] << ", ";
  //     // }
  //     // std::cout << std::endl;

  //     uint32_t id = 0;
  //     for (uint32_t lmpath = numNonnull - 1; lmpath != UINT32_MAX; lmpath = aux[2 * lmpath]) {
  //       stack_.push(lmpath);
  //     }
  //     while (!stack_.empty()) {
  //       const auto parent = stack_.top();
  //       if (!(parent & (UINT32_C(1) << 31))) {
  //         // std::cout << "f: parent = " << parent << ", id = " << id << std::endl;
  //         const auto left = aux[2 * parent];
  //         stack_.pop();
  //         stack_.push(parent | (UINT32_C(1) << 31)); // mark that next time right child to be processed
  //         if (left < id) { // NOTE: it is not a mistake of "left <= id"
  //           // Copy the text from previous occurrence.
  //           for (uint32_t j = aux[2 * left]; j < aux[2 * left + 1]; ++j) {
  //             v.push_back(v[j]);
  //           }
  //         }
  //         aux[2 * parent] = v.size(); // overwrite to starting position
  //         if (left != UINT32_MAX) {
  //           aux[2 * parent] -= (aux[2 * left + 1] - aux[2 * left]);
  //         }
  //         {
  //           uint32_t k = vbeg[parent];
  //           while (k < vbeg[parent+1]) {
  //             const auto temp = ReadBlock(vstr, k);
  //             k += 1 + (temp.second > 1);
  //             for (uint32_t j = 0; j < temp.second; ++j) {
  //               v.push_back(temp.first);
  //             }
  //           }
  //         }
  //         const uint32_t right = aux[2 * parent + 1];
  //         if (right != UINT32_MAX && right <= id) {
  //           // Copy the text from previous occurrence.
  //           for (uint32_t j = aux[2 * right]; j < aux[2 * right + 1]; ++j) {
  //             v.push_back(v[j]);
  //           }
  //         } else {
  //           for (uint32_t path = right; path != UINT32_MAX && path > id; path = aux[2 * path]) {
  //             stack_.push(path);
  //           }
  //         }
  //       } else {
  //         stack_.pop();
  //         id = parent & itmmti::ctcbits::UINTW_MAX(31); // remove the MSB bit that is used for flag
  //         aux[2 * id + 1] = v.size();
  //         std::cout << "done: id = " << id << ", v.size() = " << v.size() << ", start = " << aux[2 * id] << ", len = " << aux[2 * id + 1] - aux[2 * id] << std::endl;
  //       }
  //     }
  //     ClearStack();

  //     uint32_t k = vbeg[numNonnull+1];
  //     while (k < vbeg[numNonnull+2]) {
  //       const auto temp = ReadBlock(vstr, k);
  //       k += 1 + (temp.second > 1);
  //       for (uint32_t j = 0; j < temp.second; ++j) {
  //         v.push_back(temp.first);
  //       }
  //     }
  //   }
  // }


  void RePair::CalcSequence_Debug
  (
   const uint8_t turn,
   std::vector<uint32_t> & v,
   std::vector<uint32_t> & aux,
   uint32_t & processedNum,
   const uint32_t nodeId
   ) {
    const BlockVecT & vbeg = vstrBeg_[turn];
    const BlockVecT & vstr = vstr_[turn];

    const auto left = aux[2 * nodeId];
    aux[2 * nodeId] = v.size();
    if (left != UINT32_MAX) {
      if (left+1 <= processedNum) {
        // Copy the text from previous occurrence.
        for (uint32_t j = aux[2 * left]; j < aux[2 * left + 1]; ++j) {
          v.push_back(v[j]);
        }
      } else {
        CalcSequence_Debug(turn, v, aux, processedNum, left);
      }
    }
    {
      uint32_t k = vbeg[nodeId];
      while (k < vbeg[nodeId+1]) {
        const auto temp = ReadBlock(vstr, k);
        k += 1 + (temp.second > 1);
        for (uint32_t j = 0; j < temp.second; ++j) {
          v.push_back(temp.first);
        }
      }
    }
    const uint32_t right = aux[2 * nodeId + 1];
    if (right != UINT32_MAX) {
      if (right+1 <= processedNum) {
        // Copy the text from previous occurrence.
        for (uint32_t j = aux[2 * right]; j < aux[2 * right + 1]; ++j) {
          v.push_back(v[j]);
        }
      } else {
        CalcSequence_Debug(turn, v, aux, processedNum, right);
      }
    }
    aux[2 * nodeId + 1] = v.size();
    ++processedNum;

    {// check correctness of data structures
      {
        uint32_t k = vbeg[nodeId];
        if (k > vbeg[nodeId + 1]) {
          std::cout << "error: nodeId = " << nodeId << ", k = " << k << ", vbeg[nodeId + 1] = " << vbeg[nodeId + 1] << std::endl;
          exit(1);
        } else if (k == vbeg[nodeId + 1]) {
          {
            const auto temp = LeftMBlockInStr(vbeg, vstr, nodeId);
            if (temp.first != UINT32_MAX || temp.second != 0) {
              std::cout << "error: lmbstr nodeId = " << nodeId << ", temp.first = " << temp.first << ", temp.second = " << temp.second << std::endl;
              exit(1);
            }
          }
          {
            const auto temp = RightMBlockInStr(vbeg, vstr, nodeId);
            if (temp.first != UINT32_MAX || temp.second != 0) {
              std::cout << "error: rmbstr nodeId = " << nodeId << ", temp.first = " << temp.first << ", temp.second = " << temp.second << std::endl;
              exit(1);
            }
          }
        } else {
          uint32_t ch;
          {
            const auto temp = ReadBlock(vstr, k);
            k += 1 + (temp.second > 1);
            ch = temp.first;
          }
          while (k < vbeg[nodeId + 1]) {
            const auto temp = ReadBlock(vstr, k);
            k += 1 + (temp.second > 1);
            if (ch == temp.first) {
              std::cout << "error: consecutive nodeId = " << nodeId << ", k = " << k << ", temp.first = " << temp.first << std::endl;
              exit(1);
            }
            ch = temp.first;
          }
          if (k > vbeg[nodeId + 1]) {
            std::cout << "error: nodeId = " << nodeId << ", k = " << k << ", vbeg[nodeId + 1] = " << vbeg[nodeId + 1] << std::endl;
            exit(1);
          }
        }
      }
      {// check leftmost block
        const auto temp = LeftMBlock(nodeId);
        const uint32_t ch = temp.first;
        uint32_t len = temp.second >> 1;
        for (uint32_t j = aux[2 * nodeId]; j < aux[2 * nodeId + 1]; ++j) {
          if (ch != v[j]) {
            break;
          }
          --len;
        }
        if (len) {
          std::cout << "error: nodeId = " << nodeId << ", lmb block = (" << temp.first << ", " << temp.second << "), len = " << len
                    << ", left = " << left << ", right = " << right << std::endl;
          if (left != UINT32_MAX) {
            std::cout << "leftlen = " << aux[2 * left + 1] - aux[2 * left] << std::endl;
          }
          for (uint32_t j = vbeg[nodeId]; j < vbeg[nodeId + 1]; ++j) {
            std::cout << vstr[j] << " ";
          }
          std::cout << std::endl;
          for (uint32_t j = aux[2 * nodeId]; j < aux[2 * nodeId + 1]; ++j) {
            std::cout << v[j] << " ";
          }
          std::cout << std::endl;
          exit(1);
        }
      }
      {// check rightmost block
        const auto temp = RightMBlock(nodeId);
        const uint32_t ch = temp.first;
        uint32_t len = temp.second >> 1;
        for (uint32_t j = aux[2 * nodeId + 1]; j > aux[2 * nodeId]; --j) {
          if (ch != v[j - 1]) {
            break;
          }
          --len;
        }
        if (len) {
          std::cout << "error: nodeId = " << nodeId << ", rmb block = (" << temp.first << ", " << temp.second << "), len = " << len << std::endl;
          exit(1);
        }
      }
      // {// check if LeftMBlockInStr works
      //   const auto temp = LeftMBlockInStr(vstr, vbeg, nodeId);
      //   const uint32_t ch = temp.first;
      //   uint32_t len = temp.second >> 1;
      //   uint32_t j = aux[2 * nodeId];
      //   if (left != UINT32_MAX) {
      //     j += aux[2 * left + 1] - aux[2 * left];
      //   }
      //   const uint32_t end = j;
      //   while (j < end) {
      //     if (ch != v[j]) {
      //       break;
      //     }
      //     --len;
      //     ++j;
      //   }
      //   if (len) {
      //     std::cout << "error: nodeId = " << nodeId << ", lmb block in str = (" << temp.first << ", " << temp.second << "), len = " << len << std::endl;
      //     exit(1);
      //   }
      // }
      // {// check if RightMBlockInStr works
      //   const auto temp = LeftMBlockInStr(vstr, vbeg, nodeId);
      //   const uint32_t ch = temp.first;
      //   uint32_t len = temp.second >> 1;
      //   uint32_t j = aux[2 * nodeId + 1];
      //   if (right != UINT32_MAX) {
      //     j -= aux[2 * right + 1] - aux[2 * right];
      //   }
      //   const uint32_t beg = j;
      //   while (j > beg) {
      //     if (ch != v[j - 1]) {
      //       break;
      //     }
      //     --len;
      //     --j;
      //   }
      //   if (len) {
      //     std::cout << "error: nodeId = " << nodeId << ", rmb block in str = (" << temp.first << ", " << temp.second << "), len = " << len << std::endl;
      //     exit(1);
      //   }
      // }
    }
  }


  void RePair::CalcSequence_Debug
  (
   const uint8_t turn,
   std::vector<uint32_t> & v
   ) {
    const BlockVecT & vbeg = vstrBeg_[turn];
    const BlockVecT & vstr = vstr_[turn];

    const uint32_t numNonnull = GetNumNonnull();
    {
      uint32_t k = vbeg[numNonnull];
      while (k < vbeg[numNonnull+1]) {
        const auto temp = ReadBlock(vstr, k);
        k += 1 + (temp.second > 1);
        for (uint32_t j = 0; j < temp.second; ++j) {
          v.push_back(temp.first);
        }
      }
    }

    if (numNonnull) {
      // std::cout << __func__ << ": numNonnull = " << numNonnull << std::endl;
      uint32_t inner_i = 0;
      uint32_t leafBit_i = 0;
      uint32_t leaf_i = 0;
      std::vector<uint32_t> aux(2 * numNonnull);
      aux.resize(2 * numNonnull);
      for (uint32_t i = 0; i < skelton_.size(); ++i) {
        if (IsOpenBit(i)) { // leaf
          if (leafBit_.readBit(leafBit_i++)) {
            stack_.push(leafVar_[leaf_i++]);
          } else {
            stack_.push(UINT32_MAX);
          }
        } else {
          auto right = stack_.top();
          stack_.pop();
          auto left = stack_.top();
          stack_.pop();
          aux[2 * inner_i] = left;
          aux[2 * inner_i + 1] = right;
          stack_.push(inner_i);
          ++inner_i;
        }
      }
      ClearStack();

      uint32_t processedNum = 0;
      CalcSequence_Debug(turn, v, aux, processedNum, numNonnull-1);

      uint32_t k = vbeg[numNonnull+1];
      while (k < vbeg[numNonnull+2]) {
        const auto temp = ReadBlock(vstr, k);
        k += 1 + (temp.second > 1);
        for (uint32_t j = 0; j < temp.second; ++j) {
          v.push_back(temp.first);
        }
      }
    }
  }


  void RePair::CalcLrIds
  (
   std::vector<uint32_t> & lrIds //!< capacity should be >= 2 * (numNonnull + 2)
   ) {
    const uint32_t numNonnull = GetNumNonnull();

    lrIds[2 * numNonnull] = UINT32_MAX;
    lrIds[2 * numNonnull + 1] = numNonnull - 1; // UINT32_MAX if numNonnull == 0
    if (numNonnull) {
      // std::cout << __func__ << ": numNonnull = " << numNonnull << std::endl;
      uint32_t inner_i = 0;
      uint32_t leafBit_i = 0;
      uint32_t leaf_i = 0;
      for (uint32_t i = 0; i < skelton_.size(); ++i) {
        if (IsOpenBit(i)) { // leaf
          if (leafBit_.readBit(leafBit_i++)) {
            stack_.push(leafVar_[leaf_i++]);
          } else {
            stack_.push(UINT32_MAX);
          }
        } else {
          auto right = stack_.top();
          stack_.pop();
          auto left = stack_.top();
          stack_.pop();
          lrIds[2 * inner_i] = left;
          lrIds[2 * inner_i + 1] = right;
          stack_.push(inner_i);
          ++inner_i;
        }
      }
      ClearStack();
      lrIds[2 * numNonnull + 2] = numNonnull;
      lrIds[2 * numNonnull + 3] = UINT32_MAX;
      lrIds.resize(2 * (numNonnull + 2));
    } else {
      lrIds.resize(2);
    }
  }


  void RePair::CalcLen
  (
   const uint8_t turn,
   const std::vector<uint32_t> & lrIds,
   std::vector<uint32_t> & len //!< capacity should be >= (numNonnull + 2)
   ) {
    const BlockVecT & vbeg = vstrBeg_[turn];
    const BlockVecT & vstr = vstr_[turn];

    for (uint32_t nodeId = 0; nodeId < lrIds.size() / 2; ++nodeId) {
      len[nodeId] = 0;
      const auto left = lrIds[2 * nodeId];
      if (left != UINT32_MAX) {
        len[nodeId] += len[left];
      }
      {
        uint32_t k = vbeg[nodeId];
        while (k < vbeg[nodeId+1]) {
          const auto temp = ReadBlock(vstr, k);
          k += 1 + (temp.second > 1);
          len[nodeId] += temp.second;
        }
      }
      const uint32_t right = lrIds[2 * nodeId + 1];
      if (right != UINT32_MAX) {
        len[nodeId] += len[right];
      }
    }
  }


  std::pair<uint32_t, uint32_t> RePair::CalcVar
  (
   uint32_t pos,
   const std::vector<uint32_t> & lrIds,
   const std::vector<uint32_t> & len
   ) {
    uint32_t nodeId = lrIds.size() / 2 - 1;

    while (true) {
      const auto left = lrIds[2 * nodeId];
      const auto right = lrIds[2 * nodeId + 1];
      const auto llen = (left != UINT32_MAX) ? len[left] : 0;
      if (pos < llen) {
        nodeId = left;
      } else if (right != UINT32_MAX && pos >= len[nodeId] - len[right]) {
        pos -= (len[nodeId] - len[right]);
        nodeId = right;
      } else {
        return {nodeId, pos - llen};
      }
    }
  }


  // void RePair::SaveSequence(const string outputFileName, const BlockVecT & vbeg, const BlockVecT & vstr) {
  //   ofstream outputfile(outputFileName, ios::app);

  //   const uint32_t numNonnull = GetNumNonnull();
  //   {
  //     uint32_t k = vbeg[numNonnull];
  //     while (k < vbeg[numNonnull+1]) {
  //       const auto temp = ReadBlock(vstr, k);
  //       k += 1 + (temp.second > 1);
  //       for (uint32_t j = 0; j < temp.second; ++j) {
  //         outputfile << temp.first << " ";
  //       }
  //     }
  //   }
  //   if (numNonnull) {
  //     uint32_t inner_i = 0;
  //     uint32_t leafBit_i = 0;
  //     const uint32_t length = skelton_.size();
  //     std::vector<uint32_t> order(numNonnull + 1);
  //     order.resize(numNonnull + 1);
  //     order[numNonnull] = numNonnull+1;
  //     for (uint32_t i = 0; i < length; ++i) {
  //       if (IsOpenBit(i)) { // leaf
  //         stack_.push(leafBit_i++);
  //       } else {
  //         auto right = stack_.top();
  //         stack_.pop();
  //         auto left = stack_.top();
  //         stack_.pop();
  //         order[left] = inner_i++;
  //         stack_.push(right);
  //       }
  //     }
  //     ClearStack();

  //     uint32_t leaf_i = 0;
  //     leafBit_i = 0;
  //     for (uint32_t i = 0; i < numNonnull+1; ++i) {
  //       if (leafBit_.readBit(leafBit_i++)) {
  //         outputfile << LeftMC(leafVar_[leaf_i++]) << " ";
  //       }
  //       const uint32_t id = order[i];
  //       uint32_t k = vbeg[id];
  //       while (k < vbeg[id+1]) {
  //         const auto temp = ReadBlock(vstr, k);
  //         k += 1 + (temp.second > 1);
  //         for (uint32_t j = 0; j < temp.second; ++j) {
  //           outputfile << temp.first << " ";
  //         }
  //       }
  //     }
  //   }

  //   outputfile.close();
  // }

  uint RePair::GetSeqSize(){
    return idx_seq_;
  }

  uint RePair::GetRepairVarId(){
    return repairVarId_;
  }

  SEQ *RePair::GetNaiveSeq(){
    return std::move(naive_seq_);
  }

  RULE *RePair::GetDictRule(){
    return std::move(dict_rule_);
  }
}
