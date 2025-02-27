/* 
 *  Copyright (c) 2011 Shirou Maruyama
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


#ifndef ENCODERINCLUDED
#define ENCODERINCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bits.h"
#include "basics.h"
#include "repair.h"

typedef struct EncodeDictionary
{
  uint txt_len;
  uint seq_len;
  uint num_rules;
  CODE *comp_seq;
  RULE *rule;
  CODE *tcode;
} EDICT;

EDICT *ReadCFG(FILE *input);
void EncodeCFG(EDICT *dict, FILE *output);
void encodeCFG_rec(CODE code, EDICT *dict, BITOUT *bitout);
void putLeaf(uint numcode, CODE lcode, BITOUT *bitout);
void putParen(uchar b, BITOUT *bitout);
EDICT *convertDict(DICT *dict);
DICT *convertEDict(EDICT *dict);
void DestructEDict(EDICT *dict);
#endif
