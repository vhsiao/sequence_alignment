//
//  local.c
//  local
//
//  Created by Vivian Hsiao on 10/4/12.
//  Copyright (c) 2012 Vivian Hsiao. All rights reserved.
//

#include <stdio.h>
#include "alignment.h"
//function declarations
int localAlignment(FILE *fasta1, FILE *fasta2, int scoring[3]);
void fillMatrixLocal(int dp_rows, int dp_cols, intMatrix distance, directionMatrix bt, char *seq1, char *seq2, int score[]);
void best_score_for_cell_local(int row, int col, int dp_cols, intMatrix distance, int best[2], int score[3], int match);
//function definitions
// Takes 2 fasta files and scoring matrix. Prints optimal alignment score and local alignment.
int localAlignment(FILE *fasta1, FILE *fasta2, int score[3]) {
    int numSeqs1;
    int numSeqs2;
    sequence_t *seq1 = readSequence(fasta1, &numSeqs1);
    sequence_t *seq2 = readSequence(fasta2, &numSeqs2);
    int dp_rows = seq1->sequence_length+1; //length of sequence + null character + 1. Need extra space for base case.
    int dp_cols = seq2->sequence_length+1;
    direction **bt = allocateMemDirection(dp_rows, dp_cols); //backtrace table
    int **distance = allocateMemInt(dp_rows, dp_cols); //table of scores
    fillMatrixLocal(dp_rows, dp_cols, distance, bt, seq1->sequence, seq2->sequence, score);
    //printMatrix(dp_rows, dp_cols, distance, bt); //print DP table for debugging
    //print out the highest score in the table.
    int highest_score_coordinates[2] = {0, 0};
    int maxScore = max_score_in_dp(dp_rows, dp_cols, distance, highest_score_coordinates);
    printScore(maxScore, -2);
    //void backtrace(int srow, int scol, sequence_t *seq1, sequence_t *seq2, int dp_cols, direction bt[][dp_cols])
    backtraceAlignment(highest_score_coordinates[0], highest_score_coordinates[1], seq1, seq2, bt, -2);
    return 0;
}
//Populate the distances matrix and the backtrace pointer matrix according to seq1 and seq2. Score holds the scoring scheme.
void fillMatrixLocal(int dp_rows, int dp_cols, intMatrix distance, directionMatrix bt, char *seq1, char *seq2, int score[]) {
    //Initializiation
    for (int row=0; row<dp_rows; row++) {
        distance[row][0] = 0;
        bt[row][0] = prestart;
    }
    for (int col=0; col<dp_cols; col++) {
        distance[0][col] = 0;
        bt[0][col] = prestart;
    }
    //Recurrence:
    //
    //dist[row][col] = max (
    //  0 <backtrace: pre-start>
    //  dist[row-1][col] + gap <backtrace: north>
    //  dist[row][col-1] + gap <backtrace: west>
    //  dist[row-1][col-1] + mismatch <backtrace: nw>
    //  dist[row-1][col-1] + match <backtrace: nw>
    //)
    int best[2] = {0, 0};
    for (int row=1; row<dp_rows; row++) {
        for (int col=1; col<dp_cols; col++) {
            best_score_for_cell_local(row, col, dp_cols, distance, best, score,  (seq1[row-1] == seq2[col-1]));
            distance[row][col] = best[0];
            bt[row][col] = best[1];
        }
    }
}
//Finds the max score possible for a single cell of the DP table
//as the cell is being filled in.
void best_score_for_cell_local(int row, int col, int dp_cols, intMatrix dist, int best[2], int score[3], int match) {
    get_best_score_global(row, col, dp_cols, dist, best, score, match);
    if (0 > best[0]) { //0 is always an option with local alignment <backtrace: pre-start>
        best[0] = 0;
        best[1] = prestart;
    }
}
int main(int argc, const char * argv[])
{
    FILE *fasta1 = NULL;
    FILE *fasta2 = NULL;
    int scoring[3];
    getFiles(argc, 6, argv, &fasta1, &fasta2, scoring);
    //int scoring[3] = {MATCH_SCORE, MISMATCH_SCORE, GAP_SCORE};
    localAlignment(fasta1, fasta2, scoring);
    return 0;
}