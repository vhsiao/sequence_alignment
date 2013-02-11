//
//  main.c
//  globalAffine
//
//  Created by Vivian Hsiao on 10/6/12.
//  Copyright (c) 2012 Vivian Hsiao. All rights reserved.
//
#include <stdio.h>
#include <assert.h>
#include "alignment.h"
//function declarations
void globalAffineAlignment(FILE *fasta1, FILE *fasta2, int score[]);
void fillMatricesGlobalAffine(int dp_rows, int dp_cols, intMatrix scoreMain, intMatrix scoreVGap, intMatrix scoreHGap, directionMatrix backtrace, char *seq1, char *seq2, int score[4]);
void fillEntryGlobalAffine(int row, int col, intMatrix scoreMain, intMatrix scoreVGap, intMatrix scoreHGap, directionMatrix backtrace,int score[4], int match);
int max (int a, int b);
//function definitions
void globalAffineAlignment(FILE *fasta1, FILE *fasta2, int score[]){
    int numSeqs1;
    int numSeqs2;
    sequence_t *seq1 = readSequence(fasta1, &numSeqs1);
    sequence_t *seq2 = readSequence(fasta2, &numSeqs2);
    int dp_rows = seq1->sequence_length+1; //length of sequence + null character + 1. Need extra space for base case.
    int dp_cols = seq2->sequence_length+1;
    direction **backtrace = allocateMemDirection(dp_rows, dp_cols);
    int **scoreMain = allocateMemInt(dp_rows, dp_cols); //best score of prefixes
    int **scoreHGap = allocateMemInt(dp_rows, dp_cols);//best score of prefixes ending in horizontal gap (in sequence 2)
    int **scoreVGap = allocateMemInt(dp_rows, dp_cols);//best score of prefixes ending in vertical gap (in sequence 2)
    fillMatricesGlobalAffine(dp_rows, dp_cols, scoreMain, scoreHGap, scoreVGap, backtrace, seq1->sequence, seq2->sequence, score);
    //printMatrix(dp_rows, dp_cols, scoreMain, backtrace); //print DP table for debugging
    //start backtrace at the lower right of the main dp table.
    int row = seq1->sequence_length;
    int col = seq2->sequence_length;
    backtraceAlignment(row, col, seq1, seq2, backtrace, -3); // recover and print alignment.
}
void fillMatricesGlobalAffine(int dp_rows, int dp_cols, intMatrix scoreMain, intMatrix scoreVGap, intMatrix scoreHGap, directionMatrix backtrace, char *seq1, char *seq2, int score[4]){
        //Initialization
    
    int gapOpen = score[2];
    int gapExtend = score[3];
        for (int row=1; row<dp_rows; row++) {
            int x = gapOpen + row*gapExtend;
            scoreMain[row][0] = x;
            scoreHGap[row][0] = x;
            scoreVGap[row][0] = x;
            backtrace[row][0] = prestart;
        }
        for (int col=1; col<dp_cols; col++) {
            int x = gapOpen + col*gapExtend;
            scoreMain[0][col] = x;
            scoreHGap[0][col] = x;
            scoreVGap[0][col] = x;
            backtrace[0][col] = prestart;
        }
        // Recurrence: see comments of get_best_score_global_affine
        for (int row=1; row<dp_rows; row++) {
            for (int col=1; col<dp_cols; col++) {
                fillEntryGlobalAffine(row, col, scoreMain, scoreHGap, scoreVGap, backtrace, score, (seq1[row-1] == seq2[col-1]));
            }
        }
        printScore(scoreMain[dp_rows-1][dp_cols-1], -3);
}
// RECURRENCE:
// A three level edit graph keeps the runtime quadratic.
// scoreMain[row][col] = max(
//      dist[row-1][col-1]+ match <backtrace: nw>,
//      dist[row-1][col-1] + mismatch <backtrace: nw>)
//      scoreVGap[row][col]
//      scoreHGap[row][col]
//      )
// scoreVGap = max (
//      scoreVGap[row-1][col] + gapextend
//      scoreMain[row-1][col] + gapextend + gapopen
// scoreHGap = max (
//      scoreHGap[row][col-1] + gapextend
//      scoreMain[row][col-1] + gapextend + gapopen
//
// score: 0:match 1:mismatch 2:gapopen 3:gapextend
// At the end, the entries at [row][col] in the three edit graphs are set correctly, as well as the corresponding entry in the
// backtrace table.
void fillEntryGlobalAffine(int row, int col, intMatrix scoreMain, intMatrix scoreVGap, intMatrix scoreHGap, directionMatrix backtrace,int score[4], int match) {
    int matchScore = score[0];
    int mismatch = score[1];
    int gapopen = score[2];
    int gapextend = score[3];
    scoreVGap[row][col] = max(scoreVGap[row-1][col] + gapextend, scoreMain[row-1][col] + gapextend + gapopen);
    scoreHGap[row][col] = max(scoreHGap[row][col-1] + gapextend, scoreMain[row][col-1] + gapextend + gapopen);
    // scoreMain[row][col] = max(
    //      dist[row-1][col-1]+  matchScore <backtrace: nw>,
    //      dist[row-1][col-1] + mismatch <backtrace: nw>)
    //      scoreVGap[row][col]
    //      scoreHGap[row][col]
    //      )
    int x = match?matchScore:mismatch;
    int possibleScores[3] = {scoreVGap[row][col], scoreHGap[row][col], scoreMain[row-1][col-1] + x};
    int bestScore = findMaxInArray(possibleScores, 3);
    scoreMain[row][col] =  bestScore;
    if (bestScore == scoreVGap[row][col]) {
        backtrace[row][col] = north;
    } else if (bestScore == scoreHGap[row][col]) {
        backtrace[row][col] = west;
    } else { //bestScore == scoreMain[row][col] + x
        assert(bestScore == scoreMain[row][col]);
        backtrace[row][col] = nw;
    }
}
int max (int a, int b) {
    return (a>b)?a:b;
}
int main(int argc, const char * argv[]){
    FILE *fasta1 = NULL;
    FILE *fasta2 = NULL;
    int scoring[4];
    //<affineGlobal> <fasta-1> <fasta-2> <match-score> <mismatch-score> <gap-open-penalty> <gap-extend-penalty>
    //scoring: match, mismatch, gapopen, gapextend
    getFiles(argc, 7, argv, &fasta1, &fasta2, scoring);
    globalAffineAlignment(fasta1, fasta2, scoring);
}