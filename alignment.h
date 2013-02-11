//
//  alignment.h
//  Alignment
//
//  Created by Vivian Hsiao on 10/5/12.
//  Copyright (c) 2012 Vivian Hsiao. All rights reserved.
//

#ifndef Alignment_alignment_h
#define Alignment_alignment_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

//Structs, enums
typedef enum dir {
    prestart, west, north, nw
} direction;

typedef struct sequence {
	int sequence_length;
	int header_length;
	char *sequence;
	char *header;
} sequence_t ;

typedef direction **directionMatrix;
typedef int **intMatrix;

//constants
#define MAX_SEQUENCES_PER_FILE 200
#define BUFFER_SIZE 1000
//REDO THIS PART LATER
#define MATCH_SCORE 1
#define MISMATCH_SCORE -1
#define GAP_SCORE -1
#define FASTA_MAXLINE 512	/* Requires FASTA file lines to be <512 characters */

//Function declarations
int fittingAlignment(sequence_t *seq1, sequence_t *seq2, int score[3], float printThreshold);
int getFiles(int argc, int correctArgc, const char *argv[], FILE **fasta1, FILE **fasta2, int score[]);
void getScoringMatrix(const char *argv[], int numConstants, int score[]);
int getSequences(FILE* fafile, sequence_t sequences[]);
sequence_t* readSequence(FILE *fasta, int *numSeqs);
int read(FILE *fp, char *buffer, char **ret_seq, char **ret_name, int *ret_L);
double backtraceAlignment(int srow, int scol, sequence_t *seq1, sequence_t *seq2, directionMatrix backtrace, float printThreshold);
void printScore(int score, float printThreshold);
void printAlignment(char *align, int largest_index);
void printMatrix (int dp_rows, int dp_cols, intMatrix distance, directionMatrix backtrace);
int get_best_score_global(int row, int col, int dp_cols, intMatrix dist, int best[2], int score[], int match);
int findMaxInArray(int scores[], int num_scores);
int max_score_in_dp(int dp_rows, int dp_cols, intMatrix distance, int highest_score_coordinates[2]);
direction** allocateMemDirection(int nrows, int ncolumns);
int** allocateMemInt(int nrows, int ncolumns);
void freeMemInt(int rows, int columns, intMatrix array);
void freeMemDirection(int rows, int columns, directionMatrix array);

#endif
