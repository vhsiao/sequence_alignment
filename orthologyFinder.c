//
//  main.c
//  orthologyFinder
//
//  Created by Vivian Hsiao on 10/7/12.
//  Copyright (c) 2012 Vivian Hsiao. All rights reserved.


// Uses Fitting Alignment to determine whether any
// genes in one collection of sequences has orthologs in another. 
//

#include <stdio.h>
#include "fitting.h"

//structs
typedef struct orthology {
    char *seq1_header;
    char *seq2_header;
} orthology_t;

//constants
#define ORTHO_MATCH_SCORE 1
#define ORTHO_MISMATCH_SCORE -1
#define ORTHO_GAP_PENALTY -1
#define SIMILARITY_THRESHOLD 85

//function declarations
int getFilesOrtho(int argc, int correctArgc, const char *argv[], FILE **fasta1, FILE **fasta2);
int orthologyFinder(FILE *fasta1, FILE *fasta2);

//function definitions
int getFilesOrtho(int argc, int correctArgc, const char *argv[], FILE **fasta1, FILE **fasta2) {
    if (argc != correctArgc) {
        printf("Usage: orthologyFinder <shorter-sequence-fasta> <longer-sequence-fasta> (e.g., ./orthologyFinder mouse.fa human.fa)");
        exit(1);
    } else {
        *fasta1 = fopen(argv[1], "r");
        if (!*fasta1) {
            printf("Trouble opening input file %s", argv[1]);
            exit(1);
        }
        *fasta2 = fopen(argv[2], "r");
        if (!*fasta2) {
            printf("Trouble opening input file %s", argv[2]);
            exit(2);
        }
    }
    return 0;
}
int orthologyFinder(FILE *fasta1, FILE *fasta2) {
    sequence_t *sequenceArray1 = (sequence_t*) malloc(sizeof(sequence_t) * MAX_SEQUENCES_PER_FILE);
    sequence_t *sequenceArray2 = (sequence_t*) malloc(sizeof(sequence_t) * MAX_SEQUENCES_PER_FILE);
    int numSeqs1 = getSequences(fasta1, sequenceArray1);
    int numSeqs2 = getSequences(fasta2, sequenceArray2);
    int score[3] = {ORTHO_MATCH_SCORE, ORTHO_MISMATCH_SCORE, ORTHO_GAP_PENALTY};
    float threshold = ((float)SIMILARITY_THRESHOLD)/100;
    for (int i=0; i<numSeqs1; i++) {
        for (int j=0; j<numSeqs2; j++) {
            fittingAlignment(&sequenceArray1[i], &sequenceArray2[j], score, threshold);
        }
    }
    return 0;
}
int main(int argc, const char * argv[]) {
    FILE *fasta1 = NULL;
    FILE *fasta2 = NULL;
    getFilesOrtho(argc, 3, argv, &fasta1, &fasta2);
    printf("Now finding orthologs between files %s and %s", argv[1], argv[2]);
    orthologyFinder(fasta1, fasta2);
    return 0;
}