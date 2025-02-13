//
// Created by Enno Adler on 25.11.24.
//

#ifndef GRAMMAR_CUTTER_CONNECTOR_H
#define GRAMMAR_CUTTER_CONNECTOR_H
#include "RePair.hpp"

extern "C" {
void rekursion_sequence(DICT *dict, uint64_t from, uint64_t to, itmmti::BitVec<>* skelton, uint64_t* post_order_index, itmmti::WBitsVec* leaf, uint64_t* leaf_index, itmmti::BitVec<>* outputted);
void rekursion_grammar(DICT *dict, CODE symbol, itmmti::BitVec<>* skelton, uint64_t* post_order_index, itmmti::WBitsVec* leaf, uint64_t* leaf_index, itmmti::BitVec<>* outputted);

void rekursion_sequence(DICT *dict, uint64_t from, uint64_t to, itmmti::BitVec<>* skelton, uint64_t* post_order_index, itmmti::WBitsVec* leaf, uint64_t* leaf_index, itmmti::BitVec<>* outputted)
{
    if (from == to - 1)
    {
        rekursion_grammar(dict, dict->comp_seq[from], skelton, post_order_index, leaf, leaf_index, outputted);
    }
    else {
        uint64_t mid = ((to - from) / 2) + from;
        // left
        if (from < mid) {
            rekursion_sequence(dict, from, mid, skelton, post_order_index, leaf, leaf_index, outputted);
        }
        // right
        if (mid < to) {
            rekursion_sequence(dict, mid, to, skelton, post_order_index, leaf, leaf_index, outputted);
        }
        // parent
        skelton->writeBit(true, (*post_order_index)++);
    }
}

void rekursion_grammar(DICT *dict, CODE symbol, itmmti::BitVec<>* skelton, uint64_t* post_order_index, itmmti::WBitsVec* leaf, uint64_t* leaf_index, itmmti::BitVec<>* outputted)
{
    // First Case: Symbol is variable, but the subtree was already written (POPPT)
    // Second Case: Symbol is terminal
    if (symbol < 256 || outputted->readBit(symbol-256) == 1)
    {
        skelton->writeBit(false, (*post_order_index)++);
        leaf->write(symbol, (*leaf_index)++);
        // TODO: There might be gaps in the naming of the rules --> Is this a problem?
    }
    else
    {
        // Otherwise: We need to write the tree in Post-order
        rekursion_grammar(dict, dict->rule[symbol].left, skelton, post_order_index, leaf, leaf_index, outputted);
        rekursion_grammar(dict, dict->rule[symbol].right, skelton, post_order_index, leaf, leaf_index, outputted);
        skelton->writeBit(true, (*post_order_index)++);
        outputted->writeBit(true, symbol-256); // Set a symbol we output their children to 1.
    }
}

__attribute__((visibility("default")))
void recompress(DICT *dict, uint32_t turn_point) {
    printf("Converting to POPPT Bit representation...\n");
    uint64_t num_rules = dict->num_rules - 256;
    uint64_t seq_len = dict->seq_len;

    itmmti::BitVec<> skelton;
    itmmti::WBitsVec leaf;
    itmmti::BitVec<> printed_rule;
    skelton.resize(seq_len + 2 * num_rules + 1);
    leaf.convert(itmmti::bits::bitSize(num_rules + 256), num_rules + 1);
    // leaf.resize(num_rules + 1);
    printed_rule.resize(2 * num_rules);
    // Step 1: Build grammar fully binary
    // Step 2: Do post-order traversal of grammar
    uint64_t skelton_pos = 0, leaf_pos = 0;
    rekursion_sequence(dict, 0, seq_len, &skelton, &skelton_pos, &leaf, &leaf_pos, &printed_rule);
    skelton.writeBit(true, skelton_pos++); // Super roots 1, need to be added always.
//
//    // copy solca grammar
//    uint32_t inner_i = 256;
//    uint32_t leaf_i = 0;
//    const uint32_t length = solca.Length() - 1; // -1 to remove bit for super root
//    for (uint32_t i = 0; i < length; ++i) {
//        const uint8_t bit = solca.GetBit(i);
//        skelton.writeBit(bit, i);
//        if (bit) { // leaf
//            leaf[leaf_i] = solca.GetLeaf(leaf_i, inner_i);
//            ++leaf_i;
//        } else { // internal node
//            ++inner_i;
//        }
//    }


// Truncate skeleton to remove unused space
    skelton.resize(skelton_pos);
    leaf.resize(leaf_pos);

    size_t size = num_rules;  // Assuming `size()` returns the number of bits
    std::cout << "Used Rules: ";
    for (size_t i = 0; i < size; ++i) {
        if (printed_rule.readBit(i))
        {
            printf("-%zu- ", i + 256);
        }
    }
    std::cout << std::endl;

    size = skelton.size();  // Assuming `size()` returns the number of bits
    std::cout << "BitVec B: ";
    for (size_t i = 0; i < size; ++i) {
        std::cout << skelton.readBit(i);
    }
    std::cout << std::endl;

    size = leaf.size();
    std::cout << "Leaf L:   ";
    for (size_t i = 0; i < size; ++i) {
        uint64_t c = leaf.read(i);
        if (c < 256) {
            printf("\'%c\' ", (unsigned char) c);
        } else {
            printf("-%llu- ", c);
        }
    }
    std::cout << std::endl;

    printf("Doing Recompression...\n");
    slp_repair::RePair compressor(num_rules, std::move(skelton), std::move(leaf));
    compressor.RePairRecompression(dict->seq_len, turn_point);

    // TODO: Translate result back.
}
}
#endif //GRAMMAR_CUTTER_CONNECTOR_H
