//
//  fitting.h
//  Alignment
//
//  Created by Vivian Hsiao on 10/7/12.
//  Copyright (c) 2012 Vivian Hsiao. All rights reserved.
//

#ifndef Alignment_fitting_h
#define Alignment_fitting_h

#include <stdio.h>
#include "alignment.h"
//function declarations
void alignTwoSeqsFitting(FILE *fasta1, FILE *fasta2, int scoring[3]);
void fillMatrixFitting(int dp_rows, int dp_cols, intMatrix distance, directionMatrix bt, char *seq1, char *seq2, int score[]);
void best_score_for_cell_fitting(int row, int col, int dp_cols, intMatrix distance, int best[2], int score[3], int match);
int max_score_in_row(int col, int dp_rows, int dp_cols, intMatrix distance, int highest_score_coordinates[2]);

#endif
