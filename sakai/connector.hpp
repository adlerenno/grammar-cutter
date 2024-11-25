//
// Created by Enno Adler on 25.11.24.
//

#ifndef GRAMMAR_CUTTER_CONNECTOR_H
#define GRAMMAR_CUTTER_CONNECTOR_H
#include "RePair.hpp"

extern "C" {
void recompress(DICT *dict, uint32_t turn_point) {
    uint64_t num_rules = dict->num_rules - 256; // TODO: correct?
    itmmti::BitVec<> skelton;
    itmmti::WBitsVec leaf;

    skelton.resize(2 * num_rules + 1);
    leaf.convert(itmmti::bits::bitSize(num_rules + 256), num_rules + 1);
    leaf.resize(num_rules + 1);
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

    // TODO: This is ChatGPT generated code. Doesn't work yet, but I'm missing a good description of the representation of Sakai+ structure to represent input data. The above is the conversation from Solca to their format.
    uint32_t inner_i = 0;  // Counter for internal nodes
    uint32_t leaf_i = 0;   // Counter for leaf nodes
    uint32_t skelton_pos = 0;

    for (uint i = 256; i < num_rules; ++i) {
        const RULE &rule = dict->rule[i];

// Check if the rule is a leaf (both left and right point to terminal symbols or are invalid)
        bool is_left_leaf = (rule.left < 256);  // Assuming terminal symbols are in the range [0, 255]
        bool is_right_leaf = (rule.right < 256);

        if (is_left_leaf) {
            skelton.writeBit(true, skelton_pos++); // Leaf
            leaf[leaf_i] = rule.left;
            ++leaf_i;
        } else {
            skelton.writeBit(false, skelton_pos++); // Internal node
            ++inner_i;
        }

        if (is_right_leaf) {
            skelton.writeBit(true, skelton_pos++); // Leaf
            leaf[leaf_i] = rule.right;
            ++leaf_i;
        } else {
            skelton.writeBit(false, skelton_pos++); // Internal node
            ++inner_i;
        }
    }

// Truncate skeleton to remove unused space
    skelton.resize(skelton_pos);

    slp_repair::RePair compressor(num_rules, std::move(skelton), std::move(leaf));
    compressor.RePairRecompression(dict->seq_len, turn_point);
}
}
#endif //GRAMMAR_CUTTER_CONNECTOR_H
