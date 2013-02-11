//
//  alignment.c
//  Alignment
//
//  Created by Vivian Hsiao on 10/5/12.
//  Copyright (c) 2012 Vivian Hsiao. All rights reserved.
//
//  Contains function and struct definitions for methods shared between alignment algorithms.
//

#include <stdio.h>
#include "alignment.h"// Function definitions.
// Checks validity of arguments.
// Reads filenames into FILE *
int getFiles(int argc, int correctArgc, const char *argv[], FILE **fasta1, FILE **fasta2, int score[]) {
    if (argc != correctArgc) {
        printf("Usage: (global | local | fitting) <fasta-1> <fasta-2> <match-score> <mismatch-score> <gap-penalty>\n or <affineGlobal> <fasta-1> <fasta-2> <match-score> <mismatch-score> <gap-open-penalty> <gap-extend-penalty> Note that mismatch and gap  penalties should be entered as nonnegative integers");
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
    getScoringMatrix(argv, correctArgc-3, score);
    return 0;
}

void getScoringMatrix(const char *argv[], int numConstants, int score[]) {
    score[0] = atoi(argv[3]);
    for (int i=1; i<numConstants; i++) {
        score[i] = -atoi(argv[i+3]);
    }
}
// BACKTRACE
// int dp_rows = seq1->sequence_length+1;
// int dp_cols = seq2->sequence_length+1;
// Follow the traceback pointers stored in the backtrace table in order to regenerate an alignment.
/*
 * ALGORITHM: given 2 sequences: seq1 (length dp_rows-1) and seq2 (length dp_cols-1)
 * __________
 * Start from the bottom right.
 * direction = backtrace[row, col]
 * while (row, col > 0) {
 *       if direction=northwest-> take seq1[row] and seq2[col]. Decrement row, col.
 *       if direction=west------> take seq1[row] and gap seq2. Decrement row.
 *       if direction=north-----> gap seq1 and take seq2[col]. Decrement col.
 * }
 * while (row>0) {take seq1[row] and gap seq2. Decrement row.}
 * while (col>0) {gap seq1 and take seq2[col]. Decrement col.}
 * 
 * return value is sequence identity.
 */

//printThreshold: -3 global/globalAffine, -2 fitting/local
double backtraceAlignment(int srow, int scol, sequence_t *seq1, sequence_t *seq2, directionMatrix bt, float printThreshold) {
    int i = 0; //alignment index
    int numMatches = 0;
    int matchLength = 0;
    int max_length = seq1->sequence_length + seq2->sequence_length;
    char *seq1_buffer = (char *)(malloc(max_length * sizeof(char)));
    char *seq2_buffer = (char *)(malloc(max_length * sizeof(char)));
    memset(seq1_buffer, '\0', max_length); //[0 0 0 ... 0]
    memset(seq2_buffer, '\0', max_length); //[0 0 0 ... 0]
    //seq1->sequence[row-1] sitting:g
    //seq2->sequence[col-1] kitten:n
    //printf("row is %d: %s\n", row, &(seq1->sequence[row-1])); //sitting->g
    //printf("col is %d: %s\n", col, &(seq2->sequence[col-1])); //kitten->n
    int row = srow;
    int col = scol;
    direction dir = bt[row][col];
    char seq1_char;
    char seq2_char;
    while (dir != prestart) {
        matchLength += 1;
        if(seq1_char == seq2_char) {
            numMatches += 1;
        }
        seq1_char = seq1->sequence[row-1];
        seq2_char = seq2->sequence[col-1];
        if (dir == nw) {
            seq1_buffer[i] = seq1_char;
            seq2_buffer[i] = seq2_char;
            row--;
            col--;
        } else if (dir == north) {
            seq1_buffer[i] = seq1_char;
            seq2_buffer[i] = '-';
            row--;
        } else { // (dir == west)
            seq1_buffer[i] = '-';
            seq2_buffer[i] = seq2_char;
            col--;
        }
        dir = bt[row][col];
        i++;
    }
    if (printThreshold<-2) {
        //while (row>0) {take seq1[row] and gap seq2. Decrement row.}
        while(row>0) {
            seq1_buffer[i] = seq1->sequence[row-1];
            seq2_buffer[i] = '-';
            i++;
            row--;
            matchLength += 1;
        }
        //while (col>0) {gap seq1 and take seq2[col]. Decrement col.}
        while(col>0) {
            seq1_buffer[i] = '-';
            seq2_buffer[i] = seq2->sequence[col-1];
            i++;
            col--;
            matchLength += 1;
        }
    }
    //i now holds the largest index of the alignment +1.
    double similarity = ((double)numMatches)/matchLength;
    if (similarity >= printThreshold) {
        printf("\n**********************************\n");
        printf("Sequence 1\n: %s", seq1->header);
        printf("Sequence 2: %s\n", seq2->header);
        printf("Similarity: %.2f\n\n", similarity);
        printf("Alignment:\n");
        printAlignment(seq1_buffer, i-1);
        printAlignment(seq2_buffer, i-1);
        printf("**********************************\n");
    }
    return similarity;
}
void printMatrix (int dp_rows, int dp_cols, intMatrix distance, directionMatrix bt) {
    printf("Score:\n");
    for (int i=0; i<dp_rows; i++) {
        for (int j=0; j<dp_cols; j++) {
            printf("%3d", distance[i][j]);
        }
        printf("\n");
    }
    printf("Backtrace:\n");
    for (int i=0; i<dp_rows; i++) {
        for (int j=0; j<dp_cols; j++) {
            printf("%3d", bt[i][j]);
        }
        printf("\n");
    }
}
void printScore(int score, float printThreshold) {
    if (printThreshold < 0){
        printf("SCORE: %d\n", score);
    }
    //printf("SCORE: %d\n", score);
}
void printAlignment(char *align, int largest_index) {
    for (int i=largest_index; i>-1; i--) {
        printf("%2c", align[i]);
    }
    printf ("\n");
}

// Modifies best. best[0] is max score, best[1] is best direction.
int get_best_score_global(int row, int col, int dp_cols, intMatrix dist, int best[2], int score[], int match) {
    int north_score = dist[row-1][col] + score[2];        //dist[row-1][col] + gap
    int west_score = dist[row][col-1] + score[2];         //dist[row][col-1] + gap
    int nw_score;
    if (match) {
        nw_score = dist[row-1][col-1] + score[0];   //dist[row-1][col-1] + match
    } else {
        nw_score = dist[row-1][col-1] + score[1];   //dist[row-1][col-1] + mismatch
    }
    int scores[3] = {north_score, west_score, nw_score};
    int max = findMaxInArray(scores, 3);
    if (max == north_score) {
        best[0] = north_score;
        best[1] = north;
    } else if (max == west_score) {
        best[0] = west_score;
        best[1] = west;
    } else { //(max == nw_score)
        best[0] = nw_score;
        best[1] = nw;
    }
    return max;
}

int findMaxInArray(int scores[], int num_scores) {
    int max = scores[0];
    for (int i=0; i<num_scores; i++) {
        if (scores[i] > max) {
            max = scores[i];
        }
    }
    return max;
}
//Finds the max score over every entry in the complete DP table
int max_score_in_dp(int dp_rows, int dp_cols, intMatrix distance, int highest_score_coordinates[2]) {
    int max = distance[0][0];
    for (int row=0; row<dp_rows; row++) {
        for(int col=0; col<dp_cols; col++) {
            if (distance[row][col] > max) {
                max = distance[row][col];
                highest_score_coordinates[0] = row;
                highest_score_coordinates[1] = col;
            }
        }
    }
    return max;
}

direction** allocateMemDirection(int nrows, int ncolumns) {
	direction **array;
	array = malloc(nrows * sizeof(direction *));
	if(array == NULL) {
		fprintf(stderr, "out of memory\n");
		exit(1);
	}
	for(int i = 0; i < nrows; i++) {
		array[i] = malloc(ncolumns * sizeof(direction));
		if(array[i] == NULL) {
			fprintf(stderr, "out of memory\n");
			exit(1);
		}
	}
	return array;
}

int** allocateMemInt(int nrows, int ncolumns) {
	int **array;
	array = malloc(nrows * sizeof(int *));
	if(array == NULL) {
		fprintf(stderr, "out of memory\n");
		exit(1);
	}
	for(int i = 0; i < nrows; i++) {
		array[i] = malloc(ncolumns * sizeof(int));
		if(array[i] == NULL) {
			fprintf(stderr, "out of memory\n");
			exit(1);
		}
	}
	return array;
}
/*
 * Functions for reading fasta files
 */
int read(FILE *fp, char *buffer, char **ret_seq, char **ret_name, int *ret_L) {
	char *s;
	char *name;
	char *seq;
	int   n;
	int   nalloc;
	
	// Peek at the lookahead buffer; see if it starts with '>'.
	if (buffer[0] != '>') return 0;
	
	// Parse out the name: the first non-whitespace token after the >
    
	name = malloc(sizeof(char) * (strlen(buffer)+1));
	strcpy(name, buffer);
	
	/* Everything until next '>' is the sequence.
	 * Seq is dynamically allocated.
	 */
	seq = malloc(sizeof(char) * 128);     // allocate seq in blocks of 128
	nalloc = 128;
	n = 0;
	while (fgets(buffer, FASTA_MAXLINE, fp))
    {
		if (buffer[0] == '>') break;        // Another '>' has been reached
		for (s = buffer; *s != '\0'; s++)
		{
			if (! isalpha(*s)) continue;    // accept any alphabetic character
			seq[n] = *s;                    // store the character, bump length n
			n++;
			if (nalloc == n)                // if out of room in seq, if so, expand
			{                               // need space for the final '\0'
				nalloc += 128;
				seq = realloc(seq, sizeof(char) * nalloc);
			}
		}
    }
	seq[n] = '\0';
	*ret_name = name;
	*ret_seq  = seq;
	*ret_L    = n;
	return 1;
}
// Reads the contents of fafile into an array of sequence_t.
int getSequences(FILE* fafile, sequence_t sequences[]) {
	int numSeqs = 0;
	char *seq;
	char *name;
	int   L;
    char  buffer[512];
    if (fafile == NULL) {
        printf("Invalid file. Exiting.");
        exit(1);
    }
	fgets(buffer, FASTA_MAXLINE, fafile);
	while (read(fafile, buffer, &seq, &name, &L)) {
		sequence_t* sequence = &(sequences[numSeqs]);
		sequence->sequence_length = (int) strlen(seq);
		sequence->header_length = (int) strlen(name);
		sequence->sequence = seq;
		sequence->header = name;
		numSeqs++;
        if (numSeqs >= MAX_SEQUENCES_PER_FILE) {
            printf("Exceeded max number of sequences per file (%d). Exiting. (The max number can be increased by editing the value of the macro MAX_SEQUENCES_PER_FILE in alignment.h)", MAX_SEQUENCES_PER_FILE);
            exit(1);
        }
    }
	return numSeqs;
}
sequence_t* readSequence(FILE *fasta, int *numSeqs) {
    sequence_t *seq = (sequence_t*) malloc(sizeof(sequence_t) * MAX_SEQUENCES_PER_FILE); //Assume no more than 200 sequences per file.
    *numSeqs = getSequences(fasta, seq);
    return seq;
}
void freeMemInt(int rows, int columns, intMatrix array) {
    for (int i=0; i<rows; i++) {
        free(array[i]);
    }
}
void freeMemDirection(int rows, int columns, directionMatrix array) {
    for (int i=0; i<rows; i++) {
        free(array[i]);
    }
}
