//
// Created by Enno Adler on 14.11.24.
//

#ifndef GRAMMAR_CUTTER_CUTTER_H
#define GRAMMAR_CUTTER_CUTTER_H

#include "repair.h"

extern void recompress(DICT* input_file, uint32_t turn_point); // Needed to compile against the library.

DICT *get_excerpt_from_grammar(DICT *d, uint from, uint to);

#endif //GRAMMAR_CUTTER_CUTTER_H
