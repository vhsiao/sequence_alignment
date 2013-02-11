//
//  fittingFunctions.c
//  Alignment
//
//  Created by Vivian Hsiao on 10/7/12.
//  Copyright (c) 2012 Vivian Hsiao. All rights reserved.
//

#include "fitting.h"
//function definitions
void alignTwoSeqsFitting(FILE *fasta1, FILE *fasta2, int score[3]) {
    sequence_t *sequenceArray1 = (sequence_t*) malloc(sizeof(sequence_t) * MAX_SEQUENCES_PER_FILE);
    sequence_t *sequenceArray2 = (sequence_t*) malloc(sizeof(sequence_t) * MAX_SEQUENCES_PER_FILE);
    getSequences(fasta1, sequenceArray1);
    getSequences(fasta2, sequenceArray2);
    sequence_t *seq1 = &sequenceArray1[0];
    sequence_t *seq2 = &sequenceArray2[0];
    fittingAlignment(seq1, seq2, score, -1.0);
}
int fittingAlignment(sequence_t *seq1, sequence_t *seq2, int score[3], float printThreshold) {
    sequence_t *shorter;
    sequence_t *longer;
    shorter = seq1;
    longer = seq2;
    if(seq1->sequence_length > seq2->sequence_length) {
        longer = seq1;
        shorter = seq2;
    }
    //arbitrarily, let the shorter sequence by indexed by the rows and the longer sequence by the columns.
    int dp_rows = shorter->sequence_length+1;
    int dp_cols = longer->sequence_length+1;
    direction **bt = allocateMemDirection(dp_rows, dp_cols); //backtrace table
    int **distance = allocateMemInt(dp_rows, dp_cols); //table of scores
    fillMatrixFitting(dp_rows, dp_cols, distance, bt, shorter->sequence, longer->sequence, score);
    //print out the highest score in the table.
    int highest_score_coordinates[2] = {0, 0};
    int maxScore = max_score_in_row(dp_rows-1, dp_rows, dp_cols, distance, highest_score_coordinates);
    printScore(maxScore, printThreshold);
    backtraceAlignment(highest_score_coordinates[0], highest_score_coordinates[1], shorter, longer, bt, printThreshold);
    freeMemDirection(dp_rows, dp_cols, bt);
    freeMemInt(dp_rows, dp_cols, distance);
    return 0;
}
void fillMatrixFitting(int dp_rows, int dp_cols, intMatrix distance, directionMatrix bt, char *shorter, char *longer, int score[]) {
    //Initializiation
    for (int row=0; row<dp_rows; row++) {
        distance[row][0] = score[2]*row;
        bt[row][0] = prestart;
    }
    for (int col=0; col<dp_cols; col++) {
        distance[0][col] = 0;
        bt[0][col] = prestart;
    }
    //Recurrence:
    //
    //dist[row][col] = max (
    //  0 <backtrace: pre-start> if row is 0
    //  dist[row-1][col] + gap <backtrace: north>
    //  dist[row][col-1] + gap <backtrace: west>
    //  dist[row-1][col-1] + mismatch <backtrace: nw>
    //  dist[row-1][col-1] + match <backtrace: nw>
    //)
    int best[2] = {0, 0};
    for (int row=1; row<dp_rows; row++) {
        for (int col=1; col<dp_cols; col++) {
            char a = shorter[row-1];
            char b = longer[col-1];
            best_score_for_cell_fitting(row, col, dp_cols, distance, best, score,  (a == b));
            distance[row][col] = best[0];
            bt[row][col] = best[1];
        }
    }
}
void best_score_for_cell_fitting(int row, int col, int dp_cols, intMatrix dist, int best[2], int score[3], int match) {
    get_best_score_global(row, col, dp_cols, dist, best, score, match);
    if ((0 > best[0]) && (row==1)) { //0 is an option with fitting alignment when we are at the start of the shorter string.<backtrace: pre-start>
        best[0] = 0;
        best[1] = prestart;
    }
}
int max_score_in_row(int row, int dp_rows, int dp_cols, intMatrix distance, int highest_score_coordinates[2]) {
    int max = -1;
    dp_rows+=0; //ignore this line of code.
    for (int col=0; col<dp_cols; col++) {
        if (distance[row][col] > max) {
            max = distance[row][col];
            highest_score_coordinates[0] = row;
            highest_score_coordinates[1] = col;
        }
    }
    return max;
}

