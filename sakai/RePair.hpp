#ifndef REPAIR_HPP_
#define REPAIR_HPP_

#define BIT_FLAG_C 1
#define BIT_FLAG_N 0

#include <regex>
#include <unordered_map>
#include <chrono>
#include <list>

#include "solca.hpp"
#include "hash_table_rp.hpp"
#include "BlockVec.hpp"
#include "BitVec.hpp"
#include "WBitsVec.hpp"
#include "RankVec.hpp"
#include "rp/repair.h"

namespace slp_repair
{

  class RePair{
  public:
    static constexpr uint8_t kPoolMax{8};
    uint8_t kOpen{1};
    uint8_t kClose{0};

  private:
    using BlockVecT = itmmti::BlockVec<uint32_t, 4096>;

  private:
    uint32_t repairVarId_;
    solca_comp::HashTableRP freqTable_;

    std::vector<uint32_t> vocc_;
    itmmti::RankVec<> nonnullRank_;
    itmmti::BitVec<> skelton_;
    itmmti::BitVec<> leafBit_;
    itmmti::WBitsVec leafVar_;
    BlockVecT vstrBeg_[2]; 
    BlockVecT vstr_[2];
    BlockVecT lrTableBeg_;
    BlockVecT lrTable_;

    std::stack<uint32_t, std::vector<uint32_t>> stack_;
    uint32_t * pool_[kPoolMax];
    uint8_t poolTop_;

    SEQ * naive_seq_;
    uint32_t idx_seq_;
    RULE * dict_rule_;

  public: 
    RePair(uint32_t numRules, itmmti::BitVec<> && skelton, itmmti::WBitsVec && leaf)
      : repairVarId_(257),
        freqTable_(),
        vocc_(numRules + 2, 0),
        nonnullRank_(numRules),
        skelton_(std::move(skelton)),
        leafBit_(numRules + 1),
        leafVar_(std::move(leaf)),
        poolTop_(0),
        naive_seq_(),
        idx_seq_(0),
        dict_rule_()
    {}
    ~RePair(){
      for(uint32_t i=0; i < poolTop_; ++i) {
        delete[] pool_[i];
      }
      for(uint32_t i=0; i < vstrBeg_[0].getNumBlocks(); ++i) {
        delete[] vstrBeg_[0].getBlockPtr(i);
      }
      for(uint32_t i=0; i < vstrBeg_[1].getNumBlocks(); ++i) {
        delete[] vstrBeg_[1].getBlockPtr(i);
      }
      for(uint32_t i=0; i < vstr_[0].getNumBlocks(); ++i) {
        delete[] vstr_[0].getBlockPtr(i);
      }
      for(uint32_t i=0; i < vstr_[1].getNumBlocks(); ++i) {
        delete[] vstr_[1].getBlockPtr(i);
      }
      for(uint32_t i=0; i < lrTableBeg_.getNumBlocks(); ++i) {
        delete[] lrTableBeg_.getBlockPtr(i);
      }
      for(uint32_t i=0; i < lrTable_.getNumBlocks(); ++i) {
        delete[] lrTable_.getBlockPtr(i);
      }
    }
    
    void RePairRecompression(size_t input_size, uint32_t turn_point);
    uint32_t GetNumNonnull() const;

    uint GetSeqSize();
    uint GetRepairVarId();
    SEQ *GetNaiveSeq();
    RULE *GetDictRule();


  private:
    void ClearStack();

    bool IsSB(uint32_t val) const;
    uint32_t SetIsSB(uint32_t val) const;
    uint32_t UnsetIsSB(uint32_t val) const;
    uint32_t IsChar(uint32_t val) const;
    uint32_t GetCharVal(uint32_t ch) const;
    uint32_t GetExponentVal(uint32_t exp) const;
    uint32_t GetActualVal(uint32_t val) const;
    bool IsOpenBit(uint32_t bit_idx) const;
    bool IsCloseBit(uint32_t bit_idx) const;

    void ShrinkBV(BlockVecT & bv);
    void PoolBlockPtr(uint32_t * ptr);
    void PoolBlocks(BlockVecT & bv, const uint32_t beg, const uint32_t end);
    void PushVal(BlockVecT & bv, const uint32_t val);
    void PushBlock(BlockVecT & bv, const std::pair<uint32_t, uint32_t> & val);
    std::pair<uint32_t, uint32_t> ReadBlock(const BlockVecT & bv, uint32_t idx) const;

    uint32_t LeftMC(const uint32_t idx) const;
    uint32_t RightMC(const uint32_t idx) const;
    std::pair<uint32_t, uint32_t> LeftMBlock(const uint32_t idx) const;
    std::pair<uint32_t, uint32_t> RightMBlock(const uint32_t idx) const;
    std::pair<uint32_t, uint32_t> LeftMBlockInStr
    (
     const BlockVecT & begbv,
     const BlockVecT & bv,
     const uint32_t idx
     ) const;
    std::pair<uint32_t, uint32_t> RightMBlockInStr
    (
     const BlockVecT & begbv,
     const BlockVecT & bv,
     const uint32_t idx
     ) const;

    void CalcVocc();
    void KickOffRecomp();

    void BuildLRTBL
    (
     const BlockVecT & vbeg,
     const BlockVecT & vstr
     );
    void BuildLRTBL
    (
     const BlockVecT & vbeg,
     const BlockVecT & vstr,
     const uint32_t idx,
     const uint32_t left,
     const uint32_t right
     );

    void ComputeFreq
    (
     const BlockVecT & vbeg,
     const BlockVecT & vstr
     );
    void ComputeFreq
    (
     const BlockVecT & vbeg,
     const BlockVecT & vstr,
     const uint32_t idx,
     const uint32_t left,
     const uint32_t right,
     bool lCountFlag, // when it is 0, skip counting repeating bigram
     bool rCountFlag // when it is 0, skip counting repeating bigram
     );

    void ReplaceNonRepeating
    (
     const uint8_t turn,
     const std::pair<uint32_t, uint32_t> & bigram
     );
    void ReplaceNonRepeating
    (
     const uint8_t turn,
     const std::pair<uint32_t, uint32_t> & bigram,
     const uint32_t idx,
     const uint32_t left,
     const uint32_t right
     );

    void ReplaceRepeating
    (
     const uint8_t turn,
     const uint32_t ch
     );
    void ReplaceRepeating
    (
     const uint8_t turn,
     const uint32_t ch,
     const uint32_t idx,
     const uint32_t left,
     const uint32_t right
     );

    void CalcSequence
    (
     const uint8_t turn,
     std::vector<uint32_t> & aux,
     uint32_t & processedNum,
     const uint32_t nodeId
     );

    void CalcSequence
    (
     const uint8_t turn
     );

    void CalcSequence_Debug
    (
     const uint8_t turn,
     std::vector<uint32_t> & v,
     std::vector<uint32_t> & aux,
     uint32_t & processedNum,
     const uint32_t nodeId
     );

    void CalcSequence_Debug
    (
     const uint8_t turn,
     std::vector<uint32_t> & v
     );

    void CalcLrIds
    (
     std::vector<uint32_t> & lrIds
     );

    void CalcLen
    (
     const uint8_t turn,
     const std::vector<uint32_t> & lrIds,
     std::vector<uint32_t> & len
     );

    std::pair<uint32_t, uint32_t> CalcVar
    (
     uint32_t pos,
     const std::vector<uint32_t> & lrIds,
     const std::vector<uint32_t> & len
     );

    void PrintStatus(const uint8_t turn) const;

    void PrintStatus
    (
     const uint8_t turn,
     const uint32_t nodeId,
     const uint32_t left,
     const uint32_t right,
     const uint32_t lparenth,
     const uint32_t rparenth
     ) const;

    // void SaveSequence(const std::string outputFileName, const BlockVecT & vbeg, const BlockVecT & vstr);

    void InitDictRule();
  };

  inline void RePair::ClearStack() {
    while (!stack_.empty()) {
      stack_.pop();
    }
  }

  inline uint32_t RePair::GetNumNonnull() const {
    return (skelton_.size() - 1) / 2;
  }

  inline bool RePair::IsSB(uint32_t val) const {
    return (val & UINT32_C(1) << 31);
  }

  inline uint32_t RePair::SetIsSB(uint32_t val) const {
    return (val | UINT32_C(1) << 31);
  }

  inline uint32_t RePair::UnsetIsSB(uint32_t val) const {
    return (val & itmmti::ctcbits::UINTW_MAX(31));
  }

  inline uint32_t RePair::IsChar(uint32_t val) const {
    return val & 1;
  }

  inline uint32_t RePair::GetCharVal(uint32_t ch) const {
    return (ch << 1) + 1;
  }

  inline uint32_t RePair::GetExponentVal(uint32_t exp) const {
    return (exp << 1);
  }

  inline uint32_t RePair::GetActualVal(uint32_t val) const {
    return (val >> 1);
  }

  inline bool RePair::IsOpenBit(uint32_t bit_idx) const {
    return (skelton_.readBit(bit_idx) == kOpen);
  }

  inline bool RePair::IsCloseBit(uint32_t bit_idx) const {
    return (skelton_.readBit(bit_idx) == kClose);
  }

  inline void RePair::ShrinkBV(BlockVecT & bv) {
    const uint32_t size = bv.size();
    const uint32_t numBlocks = (size) ? (size - 1) / BlockVecT::kBlockSize + 1 : 0;
    // std::cout << __func__ << ": size = " << size << ", numBlocks = " << numBlocks << std::endl;
    // bv.printStatistics(std::cout, false);
    for (uint32_t i = numBlocks; i < bv.getNumBlocks(); ++i) {
      PoolBlockPtr(bv.getBlockPtr(i));
    }
    bv.reduceNumBlocks(numBlocks);
  }

  inline void RePair::PoolBlockPtr(uint32_t * ptr) {
    if (poolTop_ < kPoolMax) {
      pool_[poolTop_++] = ptr;
    } else {
      delete[] ptr;
    }
  }

  inline void RePair::PoolBlocks(BlockVecT & bv, const uint32_t beg, const uint32_t end) {
    for (uint32_t i = beg / BlockVecT::kBlockSize; i < end / BlockVecT::kBlockSize; ++i) {
      PoolBlockPtr(bv.getBlockPtr(i));
    }
  }

  inline void RePair::PushVal(BlockVecT & bv, const uint32_t val) {
    const auto size = bv.size();
    if (size >= bv.capacity()) {
      if (poolTop_) {
        bv.appendBlock(pool_[--poolTop_]);
      } else {
        bv.appendBlock();
      }
    }
    bv.resize(size + 1);
    bv[size] = val;
  }

  inline void RePair::PushBlock(BlockVecT & bv, const std::pair<uint32_t, uint32_t> & val) {
    PushVal(bv, GetCharVal(val.first));
    if (val.second > 1) {
      PushVal(bv, GetExponentVal(val.second));
    }
  }

  inline std::pair<uint32_t, uint32_t> RePair::ReadBlock(const BlockVecT & bv, uint32_t idx) const {
    const auto ch = GetActualVal(bv[idx]);
    if (idx+1 == bv.size()) {
      return std::make_pair(ch, 1);
    }
    const auto val = bv[idx + 1];
    if (IsChar(val)) {
      return std::make_pair(ch, 1);
    } else {
      return std::make_pair(ch, GetActualVal(val));
    }
  }

  inline uint32_t RePair::LeftMC(const uint32_t idx) const {
    return GetActualVal(lrTable_[UnsetIsSB(lrTableBeg_[idx])]);
  }

  inline uint32_t RePair::RightMC(const uint32_t idx) const {
    auto i = UnsetIsSB(lrTableBeg_[idx + 1]) - 1;
    auto val = lrTable_[i];
    if (IsChar(val)) {
      return GetActualVal(val);
    } else {
      return GetActualVal(lrTable_[i - 1]);
    }
  }

  inline std::pair<uint32_t, uint32_t> RePair::LeftMBlock(const uint32_t idx) const {
    const auto tblBeg = lrTableBeg_[idx];
    const auto pair = ReadBlock(lrTable_, UnsetIsSB(tblBeg));
    return {pair.first, (pair.second << 1) | !IsSB(tblBeg)};
  }

  inline std::pair<uint32_t, uint32_t> RePair::RightMBlock(const uint32_t idx) const {
    auto tblBeg = lrTableBeg_[idx];
    if (IsSB(tblBeg)) {
      return LeftMBlock(idx);
    }
    tblBeg = UnsetIsSB(tblBeg);
    const auto pair = ReadBlock(lrTable_, tblBeg + 1 + (!IsChar(lrTable_[tblBeg+1])));
    return {pair.first, (pair.second << 1) + 1};
  }

  inline std::pair<uint32_t, uint32_t> RePair::LeftMBlockInStr
  (
   const BlockVecT & begbv,
   const BlockVecT & bv,
   const uint32_t idx
   ) const {
    const auto beg = begbv[idx];
    const auto size = begbv[idx+1] - beg;
    if (size == 0) {
      return {UINT32_MAX, 0};
    }
    const auto val = bv[beg];
    const auto ch = GetActualVal(val);
    if (size == 1) {
      return {ch, 1 << 1};
    }
    const auto val2 = bv[beg + 1];
    if (IsChar(val2)) {
      return {ch, (1 << 1) + 1};
    } else {
      return {ch, (GetActualVal(val2) << 1) | (size > 2)}; // If LSB is raised it means repeating is terminated.
    }
  }

  inline std::pair<uint32_t, uint32_t> RePair::RightMBlockInStr
  (
   const BlockVecT & begbv,
   const BlockVecT & bv,
   const uint32_t idx
   ) const {
    const auto end = begbv[idx+1];
    const auto size = end - begbv[idx];
    if (size == 0) {
      return {UINT32_MAX, 0};
    }
    const auto val2 = bv[end - 1];
    const auto aval2 = GetActualVal(val2);
    if (IsChar(val2)) {
      return {aval2, (1 << 1) | (size > 1)};
    } else {
      return {GetActualVal(bv[end - 2]), (aval2 << 1) | (size > 2)};
    }
  }
}

#endif
