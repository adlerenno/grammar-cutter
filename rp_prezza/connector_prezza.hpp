//
// Created by Enno Adler on 16.12.24.
//

#ifndef GRAMMAR_CUTTER_CONNECTOR_PREZZA_H
#define GRAMMAR_CUTTER_CONNECTOR_PREZZA_H

#include <cstdint>
#include <vector>

typedef uint32_t itype;
vector<itype> A; //alphabet (mapping int->ascii)
vector<pair<itype, itype> > G; //grammar
vector<itype> T_vec;// compressed text

extern "C" {
void repair_prezza() { // TODO: parameters.
    compute_repair(in);

    packed_gamma_file3<> out_file(out);
    //compress the grammar with Elias' gamma-encoding and store it to file
    out_file.compress_and_store(A,G,T_vec);
}
}
#endif //GRAMMAR_CUTTER_CONNECTOR_PREZZA_H
