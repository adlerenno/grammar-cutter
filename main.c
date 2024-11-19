//
// Created by Enno Adler on 14.11.24.
//

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include "main.h"
#include "repair.h"
#include "encoder.h"
#include "cutter.h"

#define check_mode(mode_compress, mode_read, compress_expected, option_name) \
do { \
	if(compress_expected) { \
		if(mode_read) { \
			fprintf(stderr, "option '-%s' not allowed when compressing\n", option_name); \
			return -1; \
		} \
		mode_compress = true; \
	} \
	else { \
		if(mode_compress) { \
			fprintf(stderr, "option '-%s' not allowed when reading the compressed graph\n", option_name); \
			return -1; \
		} \
		mode_read = true; \
	} \
} while(0)

void do_compress(char *target_filename, char *output_filename)
{
    FILE *input, *output;
    DICT *dict;
    EDICT *edict;

    input  = fopen(target_filename, "r");
    if (input == NULL) {
        puts("File open error at the beginning.");
        exit(1);
    }

    output = fopen(output_filename, "wb");
    if (output == NULL) {
        puts("File open error at the beginning.");
        exit(1);
    }

    dict = RunRepair(input);
    OutputGeneratedCFG(dict, output);
    DestructDict(dict);

    fclose(input);
    fclose(output);
    exit(0);
}

void do_excerpt(char *target_filename, uint from, uint to)
{
    FILE *input;
//    DICT *d = malloc(sizeof(DICT));
//    uint i;

    input  = fopen(target_filename, "rb");
    if (input == NULL) {
        puts("File open error at the beginning.");
        exit(1);
    }

//    fread(&d->txt_len, sizeof(uint), 1, input);
//    fread(&d->num_rules, sizeof(uint), 1, input);
//    fread(&d->seq_len, sizeof(uint), 1, input);
//
//    d->rule = (RULE*)malloc(sizeof(RULE)*d->num_rules);
//    if (d->rule == NULL) {
//        puts("Memory allocate error at InitRuleionary.");
//        exit(1);
//    }
//
//    for (i = 0; i <= CHAR_SIZE; i++) {
//        d->rule[i].left = i;
//        d->rule[i].right = DUMMY_CODE;
//    }
//
//    fread(d->rule+CHAR_SIZE+1, sizeof(RULE), d->num_rules-(CHAR_SIZE+1), input);
//
//    d->comp_seq = (CODE*)malloc(sizeof(CODE)*d->seq_len);
//    fread(d->comp_seq, sizeof(CODE), d->seq_len, input);
    EDICT *ed = ReadCFG(input);

    printf("Doing Excerpt...");
    DICT *d = convertEDict(ed);
    DICT *r = get_excerpt_from_grammar(d, from, to);

    printf("\rResult start rule: ");
    for (uint i = 0; i < r->seq_len; i++)
    {
        if (r->comp_seq[i] < 256)
        {
            printf("\'%c\' ", (unsigned char) r->comp_seq[i]);
        }
        else
        {
            printf("-%d- ", r->comp_seq[i]);
        }
    }
    printf("\n");

    fclose(input);
    DestructDict(d);
    DestructDict(r);
}

int main(int argc, char** argv) {
    int opt;
    char * filename = "";
    char *output_filename = "";
    uint from = 0, to = 0;
    bool mode_compress = false;
    bool mode_read = false;
    while ((opt = getopt(argc, argv, "hi:o:f:t:")) != -1) {
        switch (opt) {
            case 'i':
                filename = optarg;
                if (access(filename, F_OK) != 0) {
                    printf("Invalid input file.");
                    return EXIT_FAILURE;
                }
                break;
            case 'o':
                check_mode(mode_compress, mode_read, true, "o");
                output_filename = optarg;
                if (strcmp(output_filename, "") == 0)
                {
                    printf("Output file not specified correctly.");
                    return -1;
                }
                break;
            case 'f':
                check_mode(mode_compress, mode_read, false, "f");
                from = atoi(optarg);
                if (from < 0)
                {
                    printf("Invalid position.");
                    return EXIT_FAILURE;
                }
                break;
            case 't':
                check_mode(mode_compress, mode_read, false, "t");
                to = atoi(optarg);
                if (to < 0)
                {
                    printf("Invalid position.");
                    return EXIT_FAILURE;
                }
                break;
            case 'h':
            default:
                printf("Usage: \n\t%s [-i <input_file>] [-o <output_file>] [-f <from_position>] [-t <to_position>] \nor\n\t%s -h\n\n", argv[0], argv[0]);
                printf("Default parameters: \nrun_length: %d\n\n", -1);
                return 0;
        }
    }

    if (mode_compress)
    {
        do_compress(filename, output_filename);
    }
    if (mode_read)
    {
        do_excerpt(filename, from , to);
    }
}
