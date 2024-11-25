//
// Created by Enno Adler on 20.11.24.
//

#include "recompression.h"

uint replacePairs(RDS *rds, PAIR *max_pair, CODE new_code)
{
    uint i, j;
    uint num_replaced = 0;
    SEQ *seq = rds->seq;

    i = max_pair->f_pos;
    while (i != DUMMY_POS) {
        j = seq[i].next;
        if (j == rightPos_SQ(rds, i)) {
            j = seq[j].next;
        }
        updateBlock_SQ(rds, new_code, i);
        i = j;
        num_replaced++;
    }

    if (max_pair->freq != 1) {
        destructPair(rds, max_pair);
    }
    resetPQ(rds, 1);
    return num_replaced;
}
