//
//  fitting.c
//  fitting
//
//  Created by Vivian Hsiao on 10/4/12.
//  Copyright (c) 2012 Vivian Hsiao. All rights reserved.
//

#include "fitting.h"
int main(int argc, const char * argv[])
{
    FILE *fasta1 = NULL;
    FILE *fasta2 = NULL;
    int scoring[3];
    getFiles(argc, 6, argv, &fasta1, &fasta2, scoring);
    //int scoring[3] = {MATCH_SCORE, MISMATCH_SCORE, GAP_SCORE};
    alignTwoSeqsFitting(fasta1, fasta2, scoring);
    return 0;
}