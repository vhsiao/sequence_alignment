//
//  global.c
//  global
//
//  Created by Vivian Hsiao on 10/4/12.
//  Copyright (c) 2012 Vivian Hsiao. All rights reserved.
//
#include "alignment.h"
// Function declarations
void globalAlignment(FILE *fasta1, FILE *fasta2, int score[3]);
void fillMatrixGlobal(int dp_rows, int dp_cols, intMatrix dist, directionMatrix backtrace, char *seq1, char *seq2, int score[3]);
// Function definitions
// Takes 2 fasta files and scoring matrix. Prints optimal alignment score and alignment.
void globalAlignment(FILE *fasta1, FILE *fasta2, int score[3]) {
    int numSeqs1;
    int numSeqs2;
    sequence_t *seq1 = readSequence(fasta1, &numSeqs1);
    sequence_t *seq2 = readSequence(fasta2, &numSeqs2);
    int dp_rows = seq1->sequence_length+1; //length of sequence + null character + 1. Need extra space for base case.
    int dp_cols = seq2->sequence_length+1;
    direction **bt = allocateMemDirection(dp_rows, dp_cols);
    int **distance = allocateMemInt(dp_rows, dp_cols);
    fillMatrixGlobal(dp_rows, dp_cols, distance, bt, seq1->sequence, seq2->sequence, score);
    //printMatrix(dp_rows, dp_cols, distance, bt); //print DP table for debugging
    //start backtrace at the lower right of the dp table.
    int row = seq1->sequence_length;
    int col = seq2->sequence_length;
    backtraceAlignment(row, col, seq1, seq2, bt, -3); // recover and print alignment.
}//Populate the distances matrix and the backtrace pointer matrix according to seq1 and seq1. Score holds the scoring scheme.
void fillMatrixGlobal(int dp_rows, int dp_cols, intMatrix dist, directionMatrix backtrace, char *seq1, char *seq2, int score[3]) {
    //Initialization
    int gap = score[2];
    for (int row=0; row<dp_rows; row++) {
        dist[row][0] = row*gap;
        backtrace[row][0] = prestart;
    }
    for (int col=0; col<dp_cols; col++) {
        dist[0][col] = col*gap;
        backtrace[0][col] = prestart;
    }
    // Recurrence:
    // dist[row][col] = max(
    //      dist[row-1][col] + gap <backtrace: north>,
    //      dist[row][col-1] + gap <backtrace: west>,
    //      dist[row-1][col-1]+ match <backtrace: nw>,
    //      dist[row-1][col-1] + mismatch <backtrace: nw>)
    // )
    // Filling in the DP tables (distance and backtrace)
    int best[2] = {0, 0}; //{max score, traceback pointer}
    for (int row=1; row<dp_rows; row++) {
        for (int col=1; col<dp_cols; col++) {
            get_best_score_global(row, col, dp_cols, dist, best, score, (seq1[row-1] == seq2[col-1]));
            dist[row][col] = best[0];
            backtrace[row][col] = best[1];
        }
    }
    printScore(dist[dp_rows-1][dp_cols-1], -3);
}
int main(int argc, const char * argv[])
{
    FILE *fasta1 = NULL;
    FILE *fasta2 = NULL;
    int scoring[3];
    getFiles(argc, 6, argv, &fasta1, &fasta2, scoring);
    globalAlignment(fasta1, fasta2, scoring);
}