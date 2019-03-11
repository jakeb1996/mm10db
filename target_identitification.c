/**
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
* *                                                     * *
* *        Port of: target_identitification_viaC.py     * *
* *                                                     * *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
*
*
* Purpose: identify all potential targets in exons
*
* Dependencies:
*			- PCRE2 (Perl-compatible Regular Expressions) 
*				Installation: 	http://www.linuxfromscratch.org/blfs/view/svn/general/pcre2.html
*				Sample:			https://github.com/luvit/pcre2/blob/master/src/pcre2demo.c
*				Documentation:	https://www.pcre.org/current/doc/html/pcre2.html
*
* Inputs:
*           - one file with all the exon sequences, one exon per line
*           - the list of exons of interest
*
* Output:
*           - selected targets
*
* Compile:
*			clear && gcc -std=c99 target_identitification.c -o target_identitification -g -lpcre2-8 && ./target_identitification 1 1 all
*
* Porting notes:
*			- RNAfold sequence not correct
*
**/

// Flags
#define PCRE2_CODE_UNIT_WIDTH 8
#define _GNU_SOURCE

// `
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <time.h>
#include <pcre2.h> 
#include <stdbool.h>

// Function headers
char * rc(char * strIn);
char * transToDNA(char * rna);
double AT_percentage(char * seq);
int getSimilarityScore(char first, char second);
double NeedlemanWunsch(char * seq1, char * seq2);
char * itoa(int intChar);
char * prestrcat(char * strToPrepend, char * strToAppend);
char * getCurrentDt();
char * safeStrCat(char * s1, char * s2);



// Defining the patterns used to detect sequences
char * pattern_forward = "([ACG][ACGT]{19}[ACGT]GG)"; //"(?=([ACG][ACGT]{19}[ACGT]GG))";
char * pattern_forward_offsite = "(?=([ACG][ACGT]{19}[ACGT][AG]G))";
char * pattern_reverse = "(CC[ACGT][ACGT]{19}[TGC])"; //"(?=(CC[ACGT][ACGT]{19}[TGC]))";
char * pattern_reverse_offsite = "(?=(C[CT][ACGT][ACGT]{19}[TGC]))";

// Defining the patterns used for secondary structures
char * pattern_RNAstructure = "(.{28}\\({4}\\.{4}\\){4}\\.{3}\\){4}.{21}\\({4}\\.{4}\\){4}\\({7}\\.{3}\\){7}\\.{3})\\s\\((.+)\\)";
char * pattern_RNAenergy = "(-[\\d]*.[\\d]*)"; //"\\s\\((.+)\\)";

// Thresholds used when processing secondary structures
int low_energy_threshold = -30;
int high_energy_threshold = -18;

// Threshold when looking at the off-target sites
int offtarget_threshold = 75;

// Guide RNA
char * guide = "GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUU";

// Directory settings
char * dir_seq = "./mm10_input/";
char * dir_trg = "./mm10_input/";
char * dir_list = "./mm10_input/";

// IO file settings (these are modified in main() if a gene list is supplied)
char * default_inputFile_exon = "exon_sequences.txt";
char * default_inputFile_chr = "all_sequences.txt";
char * default_outputFile = "potentialTargets.txt";
char * default_out_RNAfold = "RNAfold_output.txt";
char * default_out_targetsToScore = "targetsToScore.txt";
char * default_in_targetScores = "targets_scored.txt";
char * default_exon_list = "exonList.txt";
// (and the following need dir_trg prepended)
char * default_accepted_targets = "accepted_targets.txt";
char * default_accepted_targets_sorted = "accepted_targets_sortedByGene.txt";
char * default_accepted_targets_Excel = "accepted_targets_ExcelFriendly.tsv";
char * default_rejected_targets = "rejected_targets.txt";
char * default_offTargetSites = "offtargetSites.txt";

// the actual IO file settings
char * inputFile_exon;
char * inputFile_chr;
char * outputFile;
char * out_RNAfold;
char * out_targetsToScore;
char * in_targetScores;
char * exon_list;
char * accepted_targets;
char * accepted_targets_sorted;
char * accepted_targets_Excel;
char * rejected_targets;
char * offTargetSites;

// Temp IO file settings
char * tempTargetFile = "reads.txt";
char * alignmentFile = "alignedReads.txt";

char * tempRnaFoldInputFile = "RNAfoldInput.txt";

// Name of mismatch detection binary
char * C_program = "./findMismatches_threads";

// Threading settings (defaults)
int nb_threads_C = 16;
int nb_threads_Bowtie = 16;

// List of mm10 chromosomes
char * chromosomes[22] = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "M", "X", "Y"};

// Exon size padding
// Gigestion site must be on exon, but the whole target does not have to.
// Note: if this value is different than that used to create exon sequences, we might have an error with trying to locate targets.
int padding = 17; 

typedef struct casTarget {
	char * seq_23;
	char * seq_20;
	int removed;
	int isReverseComplement;
	struct casTarget * next;
		
	// store the exon ID's here
	int * exons; 
	// keep track of how many exons it was found in. Used for memory (re)alloc.
	int exonCount; 
} casTarget;

void freeCasTarget(casTarget ** freeCasTarget);
bool removeNode(casTarget ** head, casTarget ** tail, casTarget * itemToRemove);

void freeCasTarget(casTarget ** itemToFree) {
	if (*itemToFree != NULL) {
		if ((*itemToFree)->seq_23 != NULL)
			free((*itemToFree)->seq_23);
		if ((*itemToFree)->seq_20 != NULL)
			free((*itemToFree)->seq_20);
		if ((*itemToFree)->exons != NULL)
			free((*itemToFree)->exons);
		free((*itemToFree));
	}
}
	
// START custom for C port
int EXON_SEQ_LINE_LEN = 50000;
int INT_REMOVE_MULTIPLE_MATCH_EXON = 1;
int INT_REMOVE_MULTIPLE_MATCH_GENOME = 2;
int INT_REMOVE_GC_CONTENT = 3;
int INT_REMOVE_TTTT = 4;
int INT_REMOVE_TTT = 5;
int INT_REMOVE_REVERSE_PRIMER = 6;
int INT_REMOVE_SECONDARY_STRUCTURE_THRESHOLD = 7;
int INT_REMOVE_SCORE = 8;
// END custom for C port

// Because itoa is not implemented natively, here it is
char * itoa(int intChar) {
	char * strOut;
	strOut = malloc(sizeof(char) * 2);
	snprintf(strOut, 2, "%c", intChar);
	return strOut;
}

char * prestrcat(char * strToPrepend, char * strToAppend) {
	char * strOut;
	int intTotalStrLength = (strlen(strToPrepend) + strlen(strToAppend) + 1);
	strOut = malloc(sizeof(char) * intTotalStrLength);
	snprintf(strOut, intTotalStrLength, "%s%s", strToPrepend, strToAppend);
	return strOut;
}

// returns nicely formatted date time
char * getCurrentDt() {
	char * strOut; 
	strOut = malloc(sizeof(char) * 12); // "23:59:59 PM\0"
	time_t curtime;
	struct tm * loc_time;
	curtime = time(NULL);
	loc_time = localtime(&curtime);
	strftime(strOut, 12, "%I:%M:%S %p", loc_time);
	return strOut;
}

casTarget * doRegex(char * strPattern, char * strSearchString, int offset) {
	// Declare some variables
	pcre2_code * re; 		// for the compiled re
	PCRE2_SPTR pattern;		// the regex pattern to compile
	PCRE2_SPTR subject;		// the search string	
	PCRE2_SIZE erroroffset;	// error offset
	PCRE2_SIZE * ovector = NULL;	// pointer to the output vector, where string offsets are stored.
	int errornumber;
	int rc;
	size_t subject_length;
	pcre2_match_data * match_data;
	
	// cast the pattern and searchString to PCRE2 native types
	pattern = (PCRE2_SPTR)strPattern;
	subject = (PCRE2_SPTR)strSearchString;
	
	subject_length = strlen((char *)subject);
	
	// compile the pattern
	re = pcre2_compile(pattern, PCRE2_ZERO_TERMINATED, 0, &errornumber, &erroroffset, NULL);
	
	// make sure the RE compile did not fail
	if (re == NULL) {
		PCRE2_UCHAR buffer[256];
		pcre2_get_error_message(errornumber, buffer, sizeof(buffer));
		//printf("PCRE2 compilation failed at offset %d: %s\n", (int)erroroffset,	buffer);
		return NULL;
	}
	
	match_data = pcre2_match_data_create_from_pattern(re, NULL);

	// execute the RE
	PCRE2_SIZE start_offset;
	casTarget * matchesHead = NULL;
	casTarget * matchesTail = NULL;
	casTarget * tempMatch;
	while (1) {
		if (ovector == NULL) {
			start_offset = 0;
		} else {
			start_offset = ovector[1] - offset;
		}
		rc = pcre2_match(
			re, 
			subject, 
			subject_length, 
			start_offset, 
			0, 
			match_data, 
			NULL
		);
		
		// something went wrong
		if (rc <= 0) {
			pcre2_match_data_free(match_data);   /* Release memory used for the match */
			pcre2_code_free(re);                 /* data and the compiled pattern. */
			return matchesHead;
		}

		ovector = pcre2_get_ovector_pointer(match_data);

		for (int i = 0; i < (rc - 1); i++) {		
			tempMatch = malloc(sizeof(casTarget));
			tempMatch->seq_23 = NULL;
			tempMatch->seq_20 = NULL;
			tempMatch->exons = NULL;
			tempMatch->exonCount = 0;
			
			// how long is the match?
			size_t substring_length = ovector[2*i+1] - ovector[2*i];
			tempMatch->seq_23 = malloc(sizeof(char) * (substring_length + 1));
			PCRE2_SPTR substring_start = subject + ovector[2*i];
			snprintf(tempMatch->seq_23, substring_length+1, "%.*s", (int)substring_length+1, substring_start);

			tempMatch->next = NULL;
			// for each result, allocate space to store
			if (matchesHead == NULL) {
				matchesHead = tempMatch;
				matchesTail = tempMatch;
			} else {
				matchesTail->next = tempMatch;
				matchesTail = matchesTail->next;
				matchesTail->next = NULL;
			}

		}
		
		if (offset == 0) {
			break;
		}
	}
	pcre2_code_free(re);
	pcre2_match_data_free(match_data);
	return matchesHead;	
}

/*
* * * * * * * * * * * * * * *
* * * * * * * * * * * * * * *
* *                       * *   
* *       Main method     * *
* *                       * *   
* * * * * * * * * * * * * * *
* * * * * * * * * * * * * * *
*/


int main(int argc, char * argv[])  {
	
	if (argc == 4) {
		
		nb_threads_C = atoi(argv[1]);
		nb_threads_Bowtie = atoi(argv[2]);

		char * gene;
		gene = malloc(strlen(argv[3]) + 1);
		strcpy(gene, argv[3]);
		
		inputFile_exon = malloc(strlen(default_inputFile_exon) + 1);
		strcpy(inputFile_exon, default_inputFile_exon);
		
		inputFile_chr = malloc(strlen(default_inputFile_chr) + 1);
		strcpy(inputFile_chr, default_inputFile_chr);
		
		outputFile = malloc(strlen(default_outputFile) + 1);
		strcpy(outputFile, default_outputFile);
		
		out_RNAfold = malloc(strlen(default_out_RNAfold) + 1);
		strcpy(out_RNAfold, default_out_RNAfold);
		
		out_targetsToScore = malloc(strlen(default_out_targetsToScore) + 1);
		strcpy(out_targetsToScore, default_out_targetsToScore);
		
		in_targetScores = malloc(strlen(default_in_targetScores) + 1);
		strcpy(in_targetScores, default_in_targetScores);
				
		offTargetSites = malloc(strlen(default_offTargetSites) + 1);
		strcpy(offTargetSites, default_offTargetSites);
		
		if (strcmp(gene, "all") != 0) {
			// add the gene list name to these filenames
			char * fileNameWithoutExtension;
			
			int extraMemReq = 4;
			
			// inputFile_exon
			fileNameWithoutExtension = strtok(inputFile_exon, ".");
			free(inputFile_exon);
			inputFile_exon = malloc(strlen(dir_trg) + strlen(fileNameWithoutExtension) + strlen(gene) + extraMemReq);
			sprintf(inputFile_exon, "%s%s_%s.txt", dir_trg, fileNameWithoutExtension, gene);
				
			// exon_list
			fileNameWithoutExtension = strtok(exon_list, ".");
			exon_list = malloc(strlen(dir_trg) + strlen(fileNameWithoutExtension) + strlen(gene) + extraMemReq);
			sprintf(exon_list, "%s%s_%s.txt", dir_trg, fileNameWithoutExtension, gene);

			// outputFile
			fileNameWithoutExtension = strtok(outputFile, ".");
			free(outputFile);
			outputFile = malloc(strlen(dir_trg) + strlen(fileNameWithoutExtension) + strlen(gene) + extraMemReq);
			sprintf(outputFile, "%s%s_%s.txt", dir_trg, fileNameWithoutExtension, gene);

			// accepted_targets
			fileNameWithoutExtension = strtok(accepted_targets, ".");
			accepted_targets = malloc(strlen(dir_trg) + strlen(fileNameWithoutExtension) + strlen(gene) + extraMemReq);
			sprintf(accepted_targets, "%s%s_%s.txt", dir_trg, fileNameWithoutExtension, gene);
			
			// accepted_targets_sorted
			fileNameWithoutExtension = strtok(accepted_targets_sorted, ".");
			accepted_targets_sorted = malloc(strlen(dir_trg) + strlen(fileNameWithoutExtension) + strlen(gene) + extraMemReq);
			sprintf(accepted_targets_sorted, "%s%s_%s.txt", dir_trg, fileNameWithoutExtension, gene);
			
			// accepted_targets_Excel
			fileNameWithoutExtension = strtok(accepted_targets_Excel, ".");
			accepted_targets_Excel = malloc(strlen(dir_trg) + strlen(fileNameWithoutExtension) + strlen(gene) + extraMemReq);
			sprintf(accepted_targets_Excel, "%s%s_%s.txt", dir_trg, fileNameWithoutExtension, gene);
			
			// rejected_targets
			fileNameWithoutExtension = strtok(rejected_targets, ".");
			rejected_targets = malloc(strlen(dir_trg) + strlen(fileNameWithoutExtension) + strlen(gene) + extraMemReq);
			sprintf(rejected_targets, "%s%s_%s.txt", dir_trg, fileNameWithoutExtension, gene);

		} 
		else {

			// (the following need dir_trg prepended)
			accepted_targets = malloc(strlen(dir_trg) + strlen(default_accepted_targets) + 1);
			sprintf(accepted_targets, "%s%s", dir_trg, default_accepted_targets);
			
			accepted_targets_sorted = malloc(strlen(dir_trg) + strlen(default_accepted_targets_sorted) + 1);
			sprintf(accepted_targets_sorted, "%s%s", dir_trg, default_accepted_targets_sorted);
			
			accepted_targets_Excel = malloc(strlen(dir_trg) + strlen(default_accepted_targets_Excel) + 1);
			sprintf(accepted_targets_Excel, "%s%s", dir_trg, default_accepted_targets_Excel);
			
			rejected_targets = malloc(strlen(dir_trg) + strlen(default_rejected_targets) + 1);
			sprintf(rejected_targets, "%s%s", dir_trg, default_rejected_targets);
			
			free(offTargetSites);
			offTargetSites = malloc(strlen(dir_trg) + strlen(default_offTargetSites) + 1);
			sprintf(offTargetSites, "%s%s", dir_trg, default_offTargetSites);
			
			printf("Running the method on the whole genome.\nWARNING: the program might crash if the number of potential targets exceeds the memory available on this computer.\n");
		}
		
		free(gene);
		
	} else {
		fprintf(stderr, "Wrong number of arguments: 3 expected, %d given.\nUsage: ./target_identitification <nb_threads_C:int> <nb_threads_Bowtie:int> <genes:all or gene name or set name>", argc);
		return -1;
	}
	
	
// ###################################
// ##   Processing the input file   ##
// ###################################
// 

	char * strCurrentDt = NULL;
	strCurrentDt = getCurrentDt();
	printf("%s:\tGetting ready to process %s\n", strCurrentDt, inputFile_exon);
	free(strCurrentDt);

	int POTENTIAL_GUIDE_COUNT = 0;
	
	FILE * fp = NULL;
	char * line = NULL;
	line = malloc(sizeof(char) * EXON_SEQ_LINE_LEN);
	size_t len = 0;
	size_t read;
	
	char * fileName = prestrcat(dir_seq, inputFile_exon);
	fp = fopen(fileName, "r");
	if (fp == NULL) {
		printf("Error opening file: %s\n\n", prestrcat(dir_seq, inputFile_exon));
		exit(1);
	}
	free(fileName);
	
	int exonLineNumber = 0;
	
	casTarget * possibleTargetsHead = NULL;
	casTarget * possibleTargetsTail = NULL;
	casTarget * possibleTargetsCurrent = NULL;
	
	casTarget * matchesInExon = NULL;
	casTarget * matchesInExonCurrent = NULL;
	casTarget * tempMatch = NULL;
	casTarget * nextTemp = NULL;
	bool doesMatchExistAlready = false;
	
	while ((read = getline(&line, &len, fp)) != -1) {
		exonLineNumber++;
		
		///// FORWARD 
		matchesInExon = doRegex(pattern_forward, line, 22);
		
		// Iterate over the matches for the current exon
		matchesInExonCurrent = matchesInExon;
		while (matchesInExonCurrent != NULL) {
			
			// Make sure we haven't found this target already
			possibleTargetsCurrent = possibleTargetsHead;
			doesMatchExistAlready = false;
			while (possibleTargetsCurrent != NULL) {
				if (strcmp(matchesInExonCurrent->seq_23, possibleTargetsCurrent->seq_23) == 0) {
					doesMatchExistAlready = true;
					break;
				}
				possibleTargetsCurrent = possibleTargetsCurrent->next;
			}
			
			if (doesMatchExistAlready) {
				// We have already found this target: instead of creating a new 
				// obj to handle it, just add the exon number to the list of exons
				// which it exists within.
				possibleTargetsCurrent->exons = realloc(possibleTargetsCurrent->exons, ((possibleTargetsCurrent->exonCount + 1) * sizeof(int)));
				possibleTargetsCurrent->exonCount++;
				possibleTargetsCurrent->exons[possibleTargetsCurrent->exonCount - 1] = exonLineNumber;
			
			} else {
				POTENTIAL_GUIDE_COUNT++;
				// We have never seen this target: setup a new obj for it
				tempMatch = malloc(sizeof(casTarget));
				tempMatch->seq_23 = NULL;
				tempMatch->seq_20 = NULL;
				tempMatch->exons = NULL;
				
				tempMatch->seq_23 = malloc(sizeof(char) * (strlen(matchesInExonCurrent->seq_23) + 1));
				strcpy(tempMatch->seq_23, matchesInExonCurrent->seq_23);
				//sprintf(tempMatch->seq_23, "%s", matchesInExonCurrent->seq_23);
				
				tempMatch->seq_20 = malloc(sizeof(char) * 21);
				strncpy(tempMatch->seq_20, tempMatch->seq_23, 20);
				tempMatch->seq_20[20] = '\0';
				
				tempMatch->exons = malloc(sizeof(int) * 1);
				tempMatch->exons[0] = exonLineNumber;
				tempMatch->exonCount = 1;
				tempMatch->removed = 0;
				tempMatch->next = NULL;
				
				if (possibleTargetsHead == NULL) {
					possibleTargetsHead = tempMatch;
					possibleTargetsTail = tempMatch;
				} else {		
					possibleTargetsTail->next = tempMatch;
					possibleTargetsTail = possibleTargetsTail->next;
				}
			}

			matchesInExonCurrent = matchesInExonCurrent->next;
		}		
		
		// clear matches
		matchesInExonCurrent = matchesInExon;
		while (matchesInExonCurrent != NULL) {
			nextTemp = matchesInExonCurrent->next;
			freeCasTarget(&matchesInExonCurrent);
			matchesInExonCurrent = nextTemp;
		}

		///// REVERSE		
		doesMatchExistAlready = false;
		
		matchesInExon = NULL;
		matchesInExon = doRegex(pattern_reverse, line, 22);
		// Iterate over the matches for the current exon
		matchesInExonCurrent = matchesInExon;

		char * rcTemp;
		while (matchesInExonCurrent != NULL) {
			
			// Make sure we haven't found this target already
			possibleTargetsCurrent = possibleTargetsHead;
			doesMatchExistAlready = false;
			rcTemp = rc(matchesInExonCurrent->seq_23);
			while (possibleTargetsCurrent != NULL) {
				if (strcmp(rcTemp, possibleTargetsCurrent->seq_23) == 0) {
					doesMatchExistAlready = true;
					break;
				}
				possibleTargetsCurrent = possibleTargetsCurrent->next;
			}
			
			if (doesMatchExistAlready) {
				// We have already found this target: instead of creating a new 
				// obj to handle it, just add the exon number to the list of exons
				// which it exists within.
				possibleTargetsCurrent->exons = realloc(possibleTargetsCurrent->exons, ((possibleTargetsCurrent->exonCount + 1) * sizeof(int)));
				possibleTargetsCurrent->exonCount++;
				possibleTargetsCurrent->exons[possibleTargetsCurrent->exonCount - 1] = exonLineNumber;
			
			} else {
				POTENTIAL_GUIDE_COUNT++;
				// We have never seen this target: setup a new obj for it
				tempMatch = malloc(sizeof(casTarget));
				tempMatch->seq_23 = NULL;
				tempMatch->seq_20 = NULL;
				tempMatch->exons = NULL;
				
				tempMatch->seq_23 = malloc(sizeof(char) * strlen(rcTemp) + 1);
				//sprintf(tempMatch->seq_23, "%s", rcTemp);
				strcpy(tempMatch->seq_23, rcTemp);
				
				tempMatch->seq_20 = malloc(sizeof(char) * 21);
				strncpy(tempMatch->seq_20, tempMatch->seq_23, 20);
				tempMatch->seq_20[20] = '\0';
				
				tempMatch->exons = malloc(sizeof(int) * 1);
				tempMatch->exons[0] = exonLineNumber;
				tempMatch->exonCount = 1;
				tempMatch->removed = 0;
				tempMatch->next = NULL;
				
				if (possibleTargetsHead == NULL) {
					possibleTargetsHead = tempMatch;
					possibleTargetsTail = tempMatch;
				} else {		
					possibleTargetsTail->next = tempMatch;
					possibleTargetsTail = possibleTargetsTail->next;
				}
			}

			free(rcTemp);
			matchesInExonCurrent = matchesInExonCurrent->next;
		}
		
		// clear matches
		matchesInExonCurrent = matchesInExon;
		while (matchesInExonCurrent != NULL) {
			//printf("free...\n");
			nextTemp = matchesInExonCurrent->next;
			freeCasTarget(&matchesInExonCurrent);
			matchesInExonCurrent = nextTemp;
		}

	}
	free(line);
	fclose(fp);
	//

	strCurrentDt = getCurrentDt();
	printf("\n%s:\t%i potential targets have been identified.\n", strCurrentDt, POTENTIAL_GUIDE_COUNT);
	free(strCurrentDt);

//  ##############################################################
//  ##   Removing targets that have multiple matches in exons   ##
//  ##############################################################

	strCurrentDt = getCurrentDt();
	printf("\n%s:\tRemoving all targets that have been observed more than once.\n", strCurrentDt);
	free(strCurrentDt);
	
	casTarget * possibleTargetsCurrentRC = NULL;
	char * possibleTargetsCurrentRcSeq = NULL;
	int total_occurrences = 0;
	bool reverse_also_exists = false;
	
	// For every guide, check if there were multiple perfect matches identified
	possibleTargetsCurrent = possibleTargetsHead;
	while (possibleTargetsCurrent != NULL) { 
		total_occurrences = possibleTargetsCurrent->exonCount;
		
		// Generate the reverse complement string too
		possibleTargetsCurrentRcSeq = rc(possibleTargetsCurrent->seq_23);
		
		reverse_also_exists = false;
		
		// Check if the reverse complement was found in the previous step
		possibleTargetsCurrentRC = possibleTargetsHead; 
		while (possibleTargetsCurrentRC != NULL) {
			if (strcmp(possibleTargetsCurrentRcSeq, possibleTargetsCurrentRC->seq_23) == 0) {
				total_occurrences += possibleTargetsCurrentRC->exonCount;
				reverse_also_exists = true;
				break;
			}
			possibleTargetsCurrentRC = possibleTargetsCurrentRC->next;
		}

		// Either the guide or its RC had multiple perfect matches
		if (total_occurrences > 1) {
			POTENTIAL_GUIDE_COUNT--;
			possibleTargetsCurrent->removed = INT_REMOVE_MULTIPLE_MATCH_EXON;
			if (reverse_also_exists) {
				possibleTargetsCurrentRC->removed = INT_REMOVE_MULTIPLE_MATCH_EXON;
			}
		}
		
		free(possibleTargetsCurrentRcSeq);
		possibleTargetsCurrent = possibleTargetsCurrent->next;
	}

	printf("\t\t%i potential targets are selected for the next step.\n", POTENTIAL_GUIDE_COUNT);
	
//  ###############################################
//  ##   Using Bowtie to find multiple matches   ##
//  ###############################################

	strCurrentDt = getCurrentDt();
	printf("\n%s:\tPreparing file for Bowtie analysis.\n", strCurrentDt);
	free(strCurrentDt);

	fp = fopen(tempTargetFile, "w");
	if (fp == NULL) {
		printf("Error opening file: %s\n\n", tempTargetFile);
		exit(1);
	}

	casTarget * tempTargetDict_offset = NULL;
	casTarget * tempTargetDict_offset_tail = NULL;
	
	char threeBasePamSequences[31] = "AGGCGGGGGTGGAAGCAGGAGTAG";
	char tempPamSeq[4];
	char * tempFullSeq;

	// For every spacer sequence (first 20 nucleotides), append all PAM variants.
	// Construct a file to pass to Bowtie.
	// Bowtie will find how many times each of these guides aligns onto the genome (with imperfect matches)
	possibleTargetsCurrent = possibleTargetsHead; 
	while (possibleTargetsCurrent != NULL) {
		if (possibleTargetsCurrent->removed == 0) {
			for (int i = 0; i < strlen(threeBasePamSequences); i += 3) {
				
				strncpy(tempPamSeq, threeBasePamSequences + i, 3);
				tempPamSeq[3] = '\0';
				
				//tempFullSeq = strcat(possibleTargetsCurrent->seq_20, tempPamSeq);
				tempFullSeq = malloc(sizeof(char) * 24);
				sprintf(tempFullSeq, "%s%s", possibleTargetsCurrent->seq_20, tempPamSeq);
				//printf("tempFullSeq: %s\n", tempFullSeq);
				fprintf(fp, "%s\n", tempFullSeq);
				
				if (tempTargetDict_offset == NULL) {
					tempTargetDict_offset = malloc(sizeof(casTarget));
					tempTargetDict_offset->seq_23 = NULL;
					tempTargetDict_offset->seq_20 = NULL;
					tempTargetDict_offset->exons = NULL;
					tempTargetDict_offset->next = NULL;
						
					
					tempTargetDict_offset->seq_23 = tempFullSeq;
					tempTargetDict_offset_tail = tempTargetDict_offset;
				} else {
					tempTargetDict_offset_tail->next = malloc(sizeof(casTarget));
					tempTargetDict_offset_tail = tempTargetDict_offset_tail->next;
					
					tempTargetDict_offset_tail->seq_23 = NULL;
					tempTargetDict_offset_tail->seq_20 = NULL;
					tempTargetDict_offset_tail->exons = NULL;
					tempTargetDict_offset_tail->next = NULL;
					
					tempTargetDict_offset_tail->seq_23 = tempFullSeq;
				}
			}
		}
		possibleTargetsCurrent = possibleTargetsCurrent->next;
	}

	fclose(fp);

	strCurrentDt = getCurrentDt();
	printf("\n%s:\tFile ready. Calling Bowtie.\n\n", strCurrentDt);
	free(strCurrentDt);

	// replace this with a safe method to prevent buffer overflows. for now, just allocate a ridiculous amount of space.
	char cmd[500];
	sprintf(cmd, "bowtie2 -x mm10_input/chr19.fa -p %i --reorder --no-hd -t -r -U %s -S %s", nb_threads_Bowtie, tempTargetFile, alignmentFile);
	printf("%s\n", cmd);

	system(cmd);

	// START PROCESSING BOWTIE OUTPUT
	strCurrentDt = getCurrentDt();
	printf("\n%s:\tStarting to process the Bowtie results.\n", strCurrentDt);
	free(strCurrentDt);
	
	fp = fopen(alignmentFile, "r");
	if (fp == NULL) {
		printf("Error opening file: %s\n\n", alignmentFile);
		exit(1);
	}

	int nb_occurences = 0;

	line = malloc(sizeof(char) * 3000); // how many characters is a Bowtie line????
	char * lineCopy = malloc(sizeof(char) * 3000); // how many characters is a Bowtie line????
	char * token; // does not need free'ing because strtok just passes a pointer. no mem allocated to this ptr
	int bowtieRowCount = 0;
	bool notAddedYet = 1;
	char * seq = NULL;
	
	// read each line of the bowtie output
	casTarget * curr = NULL;
	char * readRC = NULL;
	char readSeq[23];
	while ((read = getline(&line, &len, fp)) != -1) {
		strcpy(lineCopy, line);
		
		// first obtain the read from the bowtie results
		// pretty sure we only do this for the first line out of the set
		if (bowtieRowCount == 0) {
			
			
			token = strtok(line, "\t");
			int tokenNumber = 1;
			while (token != NULL) {
				token = strtok(NULL, "\t");
				if (tokenNumber == 9) {
					break;
				}
				tokenNumber++;
			}
			strcpy(readSeq, token);
			
			//seq = read;
			// get the RC of the read
			readRC = rc(readSeq);
		}
		
		// check if the read exists in the tempTargetDict
		// edit: we're going to go ahead and assume that the alignment exists
		//  curr = tempTargetDict_offset;
		//  while (curr != NULL) { 
		//  	if (strcmp(curr->seq_23, read) == 0 || strcmp(curr->seq_23, readRC) == 0) {
		//  		seq = malloc(sizeof(char) * 25);
		//  		strcpy(seq, curr->seq_23);
		//  		break;
		//  	}
		//  	curr = curr->next;
		//  }

		// we count how many of the eight reads for this target have a perfect alignment
		if (strstr(lineCopy, "XM:i:0")) {
			nb_occurences++;
			if (strstr(lineCopy, "XS:i:0")) {
				nb_occurences++;
			}
		}
		
		// make sure the current target is only added once
		if (notAddedYet) {
			if (nb_occurences > 1) {
				notAddedYet = 0;
				
				// find the guide and mark it to be deleted
				casTarget * currGuide = possibleTargetsHead;
				while (currGuide != NULL) {
					if (currGuide->removed == 0) {
						if (strcmp(currGuide->seq_23, readSeq) == 0) {
							printf("removed guide\n");
							POTENTIAL_GUIDE_COUNT--;
							currGuide->removed = INT_REMOVE_MULTIPLE_MATCH_GENOME;
						}
					}
					currGuide = currGuide->next;
				}
				
				
			}
		}

		
		// reset the counts if we're on to the next target
		if (bowtieRowCount == strlen(threeBasePamSequences) / 3 || nb_occurences > 1) {
			//free(readRC);
			bowtieRowCount = 0;
			nb_occurences = 0;
			notAddedYet = 1;
		} else {
			bowtieRowCount++;
		}
	}
	free(line);
	free(lineCopy);
	
	curr = tempTargetDict_offset;
	while (curr != NULL) { 
		nextTemp = curr->next;
		freeCasTarget(&curr);
		curr = nextTemp;
	}
	
	printf("\t\t%i potential targets are selected for the next step.\n", POTENTIAL_GUIDE_COUNT);

//  ##############################################
//  ##   Removing targets that have AT% < 45%   ##
//  ##############################################

	strCurrentDt = getCurrentDt();
	printf("\n%s:\tRemoving all targets that have AT strictly below 45.\n", strCurrentDt);
	free(strCurrentDt);

	possibleTargetsCurrent = possibleTargetsHead;

	while (possibleTargetsCurrent != NULL) {
		if (possibleTargetsCurrent->removed == 0) {
			if (AT_percentage(possibleTargetsCurrent->seq_20) < 45.0) {
				POTENTIAL_GUIDE_COUNT--;
				possibleTargetsCurrent->removed = INT_REMOVE_GC_CONTENT;
			}
		}
		possibleTargetsCurrent = possibleTargetsCurrent->next;
	}

	printf("\t\t%i potential targets are selected for the next step.\n", POTENTIAL_GUIDE_COUNT);

//  ############################################
//  ##   Removing targets that contain TTTT   ##
//  ############################################

	strCurrentDt = getCurrentDt();
	printf("\n%s:\tRemoving all targets that contain TTTT.\n", strCurrentDt);
	free(strCurrentDt);

	possibleTargetsCurrent = possibleTargetsHead;
	while (possibleTargetsCurrent != NULL) {
		if (possibleTargetsCurrent->removed == 0) {
			if (strstr(possibleTargetsCurrent->seq_23, "TTTT")) {
				POTENTIAL_GUIDE_COUNT--;
				possibleTargetsCurrent->removed = INT_REMOVE_TTTT;
			}
		}
		possibleTargetsCurrent = possibleTargetsCurrent->next;
	}

	printf("\t\t%i potential targets are selected for the next step.\n", POTENTIAL_GUIDE_COUNT);

//  ############################################
//  ##   Removing targets that contain TTT   ##
//  ############################################

	strCurrentDt = getCurrentDt();
	printf("\n%s:\tRemoving all targets that contain TTT.\n", strCurrentDt);
	free(strCurrentDt);
	
	printf("\t\tSKIPPED!\n");
	possibleTargetsCurrent = possibleTargetsHead;
	
	//while (possibleTargetsCurrent != NULL) {
	//	if (possibleTargetsCurrent->removed == 0) {
	//		if (strstr(possibleTargetsCurrent->seq_23, "TTT")) {
	//			POTENTIAL_GUIDE_COUNT--;
	//			possibleTargetsCurrent->removed = INT_REMOVE_TTT;
	//		}
	//	}
	//	possibleTargetsCurrent = possibleTargetsCurrent->next;
	//}

	printf("\t\t%i potential targets are selected for the next step.\n", POTENTIAL_GUIDE_COUNT);

//  ################################################################
//  ##   Removing targets that are too close to reverser primer   ##
//  ################################################################

	strCurrentDt = getCurrentDt();
	printf("\n%s:\tRemoving all targets that are too close to reverse primer.\n", strCurrentDt);
	free(strCurrentDt);

	possibleTargetsCurrent = possibleTargetsHead;
	while (possibleTargetsCurrent != NULL) {
		if (possibleTargetsCurrent->removed == 0) {
			double nwScore = NeedlemanWunsch(possibleTargetsCurrent->seq_20, "AAAAGCACCGACTCGGTGCC");
			//printf("%s %i\n", possibleTargetsCurrent->seq_20, (int)nwScore);
			if ((int)nwScore  > 60) {
				POTENTIAL_GUIDE_COUNT--;
				possibleTargetsCurrent->removed = INT_REMOVE_REVERSE_PRIMER;
			}
		}
		possibleTargetsCurrent = possibleTargetsCurrent->next;
	}

	printf("\t\t%i potential targets are selected for the next step.\n", POTENTIAL_GUIDE_COUNT);

//  ###############################################
//  ##   Prepare file for RNAfold                ##
//  ###############################################

	strCurrentDt = getCurrentDt();
	printf("\n%s:\tPreparing file for RNAfold analysis.\n", strCurrentDt);
	free(strCurrentDt);

	fp = fopen(tempRnaFoldInputFile, "w");
	if (fp == NULL) {
		printf("Error opening file: %s\n\n", tempRnaFoldInputFile);
		exit(1);
	}
	
	char * structure = malloc(sizeof(char) * (1 + 20 + strlen(guide) + 1)); // (G + spacer + guide) + nullTerminator
	char targetOneToTwenty[20];
	possibleTargetsCurrent = possibleTargetsHead; 
	while (possibleTargetsCurrent != NULL) {
		if (possibleTargetsCurrent->removed == 0) {
			
			memcpy(targetOneToTwenty, &possibleTargetsCurrent->seq_20[1], 20);
			targetOneToTwenty[20] = '\0';
			
			sprintf(structure, "G%s%s", targetOneToTwenty, guide);
			fprintf(fp, "%s\n", structure);

		}
		possibleTargetsCurrent = possibleTargetsCurrent->next;
	}

	fclose(fp);

	strCurrentDt = getCurrentDt();
	printf("\n%s:\tFile ready. Next: call RNAfold.\n\n", strCurrentDt);
	free(strCurrentDt);

	
//  ##########################################
//  ##   Calculating secondary structures   ##
//  ##########################################

	strCurrentDt = getCurrentDt();
	printf("\n%s:\tCalculating secondary structures.\n\n", strCurrentDt);
	free(strCurrentDt);
	
	// WARNING: removing any existing version of the RNAfold output file.
	sprintf(cmd, "rm -f %s", out_RNAfold);
	system(cmd);

	//int temp_counter = 0;
    //
	//char targetOneToTwenty[20];
    //
	//possibleTargetsCurrent = possibleTargetsHead;
	//while (possibleTargetsCurrent != NULL) {
	//	if (possibleTargetsCurrent->removed == 0) {
	//		
	//		memcpy(targetOneToTwenty, &possibleTargetsCurrent->seq_20[1], 20);
	//		targetOneToTwenty[20] = '\0';
	//		
	//		sprintf(structure, "G%s%s", targetOneToTwenty, guide);
	//		
	//		// nb: cmd was defined earlier in the script
	//		sprintf(cmd, "echo %s | RNAfold --noPS >> %s", structure, out_RNAfold);
	//		//printf("RNAfold: %s\n", cmd);
	//		system(cmd);
    //
	//		temp_counter++;
	//		if (temp_counter % 100000 == 0) {
	//			strCurrentDt = getCurrentDt();
	//			printf("%s:\t\t%i targets processed.\n", strCurrentDt, temp_counter);
	//			free(strCurrentDt);
	//		}
	//	}
	//	possibleTargetsCurrent = possibleTargetsCurrent->next;
	//}
	
	// nb: cmd was defined earlier in the script
	sprintf(cmd, "RNAfold -i %s --noPS >> %s", tempRnaFoldInputFile, out_RNAfold);
	system(cmd);
	
	free(structure);
	int total_number_structures = POTENTIAL_GUIDE_COUNT;

//  #########################################
//  ##   Processing secondary structures   ##
//  #########################################

	strCurrentDt = getCurrentDt();
	printf("\n%s:\tProcessing secondary structures.\n", strCurrentDt);
	free(strCurrentDt);

	fp = fopen(out_RNAfold, "r");
	if (fp == NULL) {
		printf("Error opening file: %s\n\n", alignmentFile);
		exit(1);
	}

	line = malloc(sizeof(char) * 1000); // how many characters is a Bowtie line????
	
	// read each line of the RNAfold output
	// nb: we need to read lines in pairs
	
	char L1[1000];
	char L2[1000];
	char target[21];

	casTarget * matchesHead = NULL;
	casTarget * matchesHeadCurr = NULL;
	possibleTargetsCurrent = possibleTargetsHead;
	double energy = 0.0;
	while (possibleTargetsCurrent != NULL) {
		if (possibleTargetsCurrent->removed == 0) {
			
			// first line of pair
			read = getline(&line, &len, fp);
			strcpy(L1, line);
			sprintf(target, "%*s", 20, L1 + (strlen(L1) - 20));
			target[20] = '\0';

			// second line of pair
			read = getline(&line, &len, fp);
			strcpy(L2, line);
			
			// TODO: wrap things with transToDNA: see python source code
			// port note: the massive if statement with quit() has been skipped. we'll forget about error checking here for now :')
			if (possibleTargetsCurrent->removed == 0) {
				matchesHead = doRegex(pattern_RNAstructure, L2, 0);
				if (matchesHead != NULL) {
					
					// clear the matches
					matchesHeadCurr = matchesHead;
					while (matchesHeadCurr != NULL) {
						nextTemp = matchesHeadCurr->next;
						freeCasTarget(&matchesHeadCurr);
						matchesHeadCurr = nextTemp;
					}
					
					// performing the regular expression again isn't in the python code
					// but for some reason, the second group isn't matching in pattern_RNAstructure
					// whereas pattern_RNAenergy performs the match successfully. so we'll just 
					// do that.
					matchesHead = doRegex(pattern_RNAenergy, L2, 0);
					
					// The structure is correct, we only reject if the energy is too low
					energy = atof(matchesHead->seq_23);
					if (energy < low_energy_threshold) {
						POTENTIAL_GUIDE_COUNT--;
						possibleTargetsCurrent->removed = INT_REMOVE_SECONDARY_STRUCTURE_THRESHOLD;
					}
					
					// clear the matches
					matchesHeadCurr = matchesHead;
					while (matchesHeadCurr != NULL) {
						nextTemp = matchesHeadCurr->next;
						freeCasTarget(&matchesHeadCurr);
						matchesHeadCurr = nextTemp;
					}
					
					
				} else {
					matchesHead = doRegex(pattern_RNAenergy, L2, 0);
					
					if (matchesHead != NULL) {
						
						
						// The structure is not correct, we only reject if the energy is not high enough
						energy = atof(matchesHead->seq_23);
						if (energy <= high_energy_threshold) {
							POTENTIAL_GUIDE_COUNT--;
							possibleTargetsCurrent->removed = INT_REMOVE_SECONDARY_STRUCTURE_THRESHOLD;
						}
						
						
						// clear the matches
						matchesHeadCurr = matchesHead;
						while (matchesHeadCurr != NULL) {
							nextTemp = matchesHeadCurr->next;
							freeCasTarget(&matchesHeadCurr);
							matchesHeadCurr = nextTemp;
						}
					}
				}
			}
		}
		possibleTargetsCurrent = possibleTargetsCurrent->next;
	}	

	free(line);
	
	printf("\t\t%i potential targets are selected for the next step.\n", POTENTIAL_GUIDE_COUNT);

//  #########################
//  ##   Scoring targets   ##
//  #########################

	strCurrentDt = getCurrentDt();
	printf("\n%s:\tScoring targets.\n\n", strCurrentDt);
	free(strCurrentDt);
	
	strCurrentDt = getCurrentDt();
	printf("\n%s:\tWriting targets to file.\n\n", strCurrentDt);
	free(strCurrentDt);
	
	int i = 0;

	fp = fopen(out_targetsToScore, "w");
	if (fp == NULL) {
		printf("Error opening file: %s\n\n", out_targetsToScore);
		exit(1);
	}

	int temp_counter = 0;
	possibleTargetsCurrent = possibleTargetsHead;
	while (possibleTargetsCurrent != NULL) {
		if (possibleTargetsCurrent->removed == 0) {
			fprintf(fp, "%s\n", possibleTargetsCurrent->seq_20);
			temp_counter++;
		}
		possibleTargetsCurrent = possibleTargetsCurrent->next;
	}
	fclose(fp);

	int total_number_scores = temp_counter;

	strCurrentDt = getCurrentDt();
	printf("\n%s:\tCalling C program.\n\n", strCurrentDt);
	free(strCurrentDt);

	sprintf(cmd, "%s %i %s %s %s %i", C_program, nb_threads_C, out_targetsToScore, offTargetSites, in_targetScores, offtarget_threshold);
	printf("%s\n", cmd);
	system(cmd);

	strCurrentDt = getCurrentDt();
	printf("\n%s:\tReading the results from file, and processing.\n", strCurrentDt);
	free(strCurrentDt);

	FILE * inFile;
	inFile = fopen(in_targetScores, "r");
	if (inFile == NULL) {
		printf("Error opening file: %s\n\n", in_targetScores);
		exit(1);
	}
	
	FILE * outFile;
	outFile = fopen(tempTargetFile, "w");
	if (outFile == NULL) {
		printf("Error opening file: %s\n\n", tempTargetFile);
		exit(1);
	}
	
	line = malloc(sizeof(char) * 1000); // how many characters is a C-program line????

	possibleTargetsCurrent = possibleTargetsHead;
	char * tempGuide; // = malloc(sizeof(char) * 24);
	char * score;// = malloc(sizeof(char) * 8);
	while (possibleTargetsCurrent != NULL) {
		if (possibleTargetsCurrent->removed == 0) {
			if (read = getline(&line, &len, inFile) != -1) {
				
				tempGuide = strtok(line, "\t");
				score = strtok(NULL, "\t");
				
				// skipping the error check that is in the python edition
				double dblScore = atof(score);
				if (dblScore < offtarget_threshold) {
					possibleTargetsCurrent->removed = INT_REMOVE_SCORE;
				} else {
					fprintf(outFile, "%s\n", possibleTargetsCurrent->seq_23);
				}
			}
		}

		possibleTargetsCurrent = possibleTargetsCurrent->next;
	}

	fclose(inFile);
	fclose(outFile);

	
	printf("\t\t%i potential targets are selected for the next step.\n", POTENTIAL_GUIDE_COUNT);
	
	
	
	strCurrentDt = getCurrentDt();
	printf("\n%s:\tCalculating their exact position using Bowtie.\n\n", strCurrentDt);
	free(strCurrentDt);
	
	sprintf(cmd, "bowtie2 -x mm10_input/chr19.fa -p %i --reorder --no-hd -t -r -U %s -S %s", nb_threads_Bowtie, tempTargetFile, alignmentFile);
	printf("%s\n", cmd);
	system(cmd);
	
	strCurrentDt = getCurrentDt();
	printf("\n%s:\t(TODO) Saving the results.\n\n", strCurrentDt);
	free(strCurrentDt);

	
	
	
	/*********************************************************/
	/************      BEGIN FREEING MEMORY     **************/
	/*********************************************************/
	/**/	
	/**/	// Free possible targets
	/**/	possibleTargetsCurrent = possibleTargetsHead;
	/**/	while (possibleTargetsCurrent != NULL) {
	/**/		nextTemp = possibleTargetsCurrent->next;
	/**/		freeCasTarget(&possibleTargetsCurrent);
	/**/		possibleTargetsCurrent = nextTemp;
	/**/	}
	/**/	
	/**/	free(accepted_targets);
	/**/	free(accepted_targets_sorted);
	/**/	free(accepted_targets_Excel);
	/**/	free(rejected_targets);
	/**/	free(offTargetSites);
	/**/	free(in_targetScores);
	/**/	free(out_targetsToScore);
	/**/	free(out_RNAfold);
	/**/	free(outputFile);
	/**/	free(inputFile_chr);
	/**/	free(inputFile_exon);
	/**/	
	/**/	return 1;
	/**/	
	/*********************************************************/
	/*********************************************************/
	

	
	FILE * inFile_score;
	inFile_score = fopen(in_targetScores, "r");
	if (inFile_score == NULL) {
		printf("Error opening file: %s\n\n", in_targetScores);
		exit(1);
	}
	
	FILE * inFile_Bowtie;
	inFile_Bowtie = fopen(alignmentFile, "r");
	if (inFile_Bowtie == NULL) {
		printf("Error opening file: %s\n\n", alignmentFile);
		exit(1);
	}
	
	FILE * outFileAccepted;
	outFileAccepted = fopen(accepted_targets, "w");
	if (outFileAccepted == NULL) {
		printf("Error opening file: %s\n\n", accepted_targets);
		exit(1);
	}
	
	FILE * inFileExonList;
	inFileExonList = fopen(prestrcat(dir_list, exon_list), "r");
	if (inFileExonList == NULL) {
		printf("Error opening file: %s\n\n", prestrcat(dir_list, exon_list));
		exit(1);
	}
	fclose(inFileExonList);
	
	
	i = 0;
	possibleTargetsCurrent = possibleTargetsHead;
	char * output_line = malloc(sizeof(char) * 5000); // surely thats enough for now
	int j = 0;
	while (possibleTargetsCurrent != NULL) {
		if (possibleTargetsCurrent->removed == 0) {
			memset(output_line, 0, sizeof(output_line));

			strcat(output_line, possibleTargetsCurrent->seq_20);
			strcat(output_line, "\t");
			
			while (i < total_number_structures) {
				
				i++;
			}
			
			if (i == total_number_scores) {
				printf("Error? target not found in RNAFold output\n");
				return 0;
			}
			
			j = 0;
			while (j < total_number_scores) {
				
				j++;
			}
			
			if (j == total_number_scores) {
				printf("Error? target not found in targetScores file\n");
				return 0;
			}
			
			fprintf(outFile, "%s", output_line);
		}

		possibleTargetsCurrent = possibleTargetsCurrent->next;
	}
	
	fclose(inFile_score);
	fclose(inFile_Bowtie);
	fclose(outFile);
}
// https://stackoverflow.com/a/8465083
char * safeStrCat(char * s1, char * s2) {
	int length = strlen(s1) + strlen(s2) + 1;
	char * result = malloc(length); // +1 for the null-terminator
    // in real code you would check for errors in malloc here
    strcpy(result, s1);
    strcat(result, s2);
	result[length] = '\0';
    return result;
}

bool removeNode(casTarget ** head, casTarget ** tail, casTarget * itemToRemove) {
	// case: where itemToRemove is the head
	//	- replace the head with the next item
	if (itemToRemove == *head) {
		if ((*head)->next != NULL) {
			*head = (*head)->next;
		} else {	
			*head = NULL;
		}
		
		return true;
	}
	
	// case: where the itemToRemove is the tail
	//	- replace the tail with the second last item
	if (itemToRemove == *tail) {	
		casTarget * curr = *head;
		while (curr->next != *tail) {
			curr = curr->next;
		}
		printf("tail: %p\n", *tail);
		*tail = curr;
		(*tail)->next = NULL;
		
		return true;
	}
	
	// case: where the itemToRemove is in the middle somewhere
	// 	- find where the itemToRemove is, then link the next item to the previous item
	casTarget * current = (*head)->next; // we already know it's not the head
	casTarget * previous = *head;
	while (current != NULL) {
		if (itemToRemove == current) {
			previous->next = current->next;
			
			return true;
		}
		
		previous = current;
		current = current->next;
	}
	
	// case: the node was not found in the list, uh oh!
	return false;
}


char * rc(char * dna) {
	char rcHash[122];
	rcHash[97]  = 't'; 	// a	ie: 97 is 'a', reverse complement is 't'
	rcHash[99]  = 'g'; 	// c
	rcHash[103] = 'c'; 	// g
	rcHash[116] = 'a'; 	// t
	rcHash[114] = 'y'; 	// r
	rcHash[121] = 'r'; 	// y
	rcHash[109] = 'k'; 	// m
	rcHash[107] = 'm'; 	// k
	rcHash[98]  = 'v'; 	// b
	rcHash[100] = 'h'; 	// d
	rcHash[104] = 'd'; 	// h
	rcHash[118] = 'b'; 	// v
	rcHash[65]  = 'T'; 	// A
	rcHash[67]  = 'G'; 	// C
	rcHash[71]  = 'C'; 	// G
	rcHash[84]  = 'A'; 	// T
	rcHash[82]  = 'Y'; 	// R
	rcHash[89]  = 'R'; 	// Y
	rcHash[77]  = 'K'; 	// M
	rcHash[75]  = 'M'; 	// K
	rcHash[66]  = 'V'; 	// B
	rcHash[68]  = 'H'; 	// D
	rcHash[72]  = 'D'; 	// H
	rcHash[86]  = 'B'; 	// V

	char * strOut = NULL;
	int strOutLen = strlen(dna) + 1;
	strOut = malloc(strOutLen);

	for (int i = 0; i < strOutLen; i++) {
		strOut[i] = rcHash[dna[strOutLen - 2 - i]];
	}
	strOut[strOutLen - 1] = 0;
	
	return strOut;
}

char * transToDNA(char * rna) {
	char * strOut = malloc(strlen(rna));
	
	for (int i = 0; i < strlen(rna); i++) {
		if (rna[i] == 85) { // 85 = 'U'
			sprintf(strOut, "%s%c", strOut, 'T');
		} else {
			sprintf(strOut, "%s%c", strOut, rna[i]);
		}
	}
	
	return strOut;
}

double AT_percentage(char * seq) {
    int total = 0;
    int length = strlen(seq);
    for (int i = 0; i < length; i++) {
		// (65 = A, 97 = a, 84 = T, 116 = t)
        if (seq[i] == 65 || seq[i] == 97 || seq[i] == 84 || seq[i] == 116) {
            total++;
		}
	}
    return (100.0 * total / length);
}

int getSimilarityScore(char first, char second) {
	// what is a tidier way to do this?
	// how can I make a hashtable effectively for this?
	// hmm
	switch (first) {
		case 65: // A
			switch (second) {
				case 65: // A
					return 10;
				case 71: // G
					return -1;
				case 67: // C
					return -3;
				case 84: // T
					return -4;
			}
		case 71: // G
			switch (second) {
				case 65: // A
					return -1;
				case 71: // G
					return 7;
				case 67: // C
					return -5;
				case 84: // T
					return -3;
			}
		case 67: // C
			switch (second) {
				case 65: // A
					return -3;
				case 71: // G
					return -5;
				case 67: // C
					return 9;
				case 84: // T
					return 0;
			}
		case 84: // T
			switch (second) {
				case 65: // A
					return -4;
				case 71: // G
					return -3;
				case 67: // C
					return 0;
				case 84: // T
					return 8;
			}
	}
	return 0;
}


// Function that aligns two given sequences and returns their similarity
double NeedlemanWunsch(char * seq1, char * seq2) {
	// TODO: cleanup the string concat methods
	int d = -5;

	char * A = malloc((strlen(seq1) + 1) * sizeof(char));
	strcpy(A, seq1);
	
	char * B = malloc((strlen(seq2) + 1) * sizeof(char));
	strcpy(B, seq2);

	int I = strlen(seq1);
	int J = strlen(seq2);

	int F[I][J];
	for (int counterI = 0; counterI < I; counterI++) {
		for (int counterJ = 0; counterJ < I; counterJ++) {
			F[counterI][counterJ] = 0;
		}
	}

	// Initialization
	for (int counterI = 0; counterI < I; counterI++) {
		F[counterI][0] = d * counterI;
	}
	
	for (int counterJ = 0; counterJ < J; counterJ++) {
		F[0][counterJ] = d * counterJ;
	}

	// Scoring
	for (int counterI = 1; counterI < I; counterI++) {

		for (int counterJ = 1; counterJ < J; counterJ++) {

			int Match = F[counterI - 1][counterJ - 1] + getSimilarityScore(A[counterI], B[counterJ]);
			int Delete = F[counterI - 1][counterJ] + d;
			int Insert = F[counterI][counterJ - 1] + d;

			if (Match >= Delete && Match >= Insert) {
				F[counterI][counterJ] = Match;
			}
			else if (Delete >= Match && Delete >= Insert) {
				F[counterI][counterJ] = Delete;
			}
			else if (Insert >= Match && Insert >= Delete) {
				F[counterI][counterJ] = Insert;
			}
		}
	}
	
	// Traceback
	int i = strlen(seq1) - 1;
	int j = strlen(seq2) - 1;

	char * AlignmentA;
	int AlignmentASize = i + 1;
	AlignmentA = malloc(sizeof(char) * (AlignmentASize));
	strcpy(AlignmentA, "");
	
	char * AlignmentB;
	int AlignmentBSize = j + 1;
	AlignmentB = malloc(sizeof(char) * (AlignmentBSize));
	strcpy(AlignmentB, "");
	
	char * itoaTemp;
	char * alignmentTemp;
	
	while (i > 0 && j > 0) {
		int Score = F[i][j];
		int ScoreDiag = F[i - 1][j - 1];
		int ScoreUp = F[i][j - 1];
		int ScoreLeft = F[i - 1][j];
		
		if (Score == ScoreDiag + getSimilarityScore(A[i], B[j])) {
			itoaTemp = itoa(A[i]);
			alignmentTemp = prestrcat(itoaTemp, AlignmentA);
			free(AlignmentA);
			free(itoaTemp);
			AlignmentA = alignmentTemp;
			
			itoaTemp = itoa(B[j]);
			alignmentTemp = prestrcat(itoaTemp, AlignmentB);
			free(AlignmentB);
			free(itoaTemp);
			AlignmentB = alignmentTemp;
			
			i--;
			j--;
			
		} else if (Score == ScoreLeft + d) {
			
			itoaTemp = itoa(A[i]);
			alignmentTemp = prestrcat(itoaTemp, AlignmentA);
			free(AlignmentA);
			free(itoaTemp);
			AlignmentA = alignmentTemp;			
			
			alignmentTemp = AlignmentB;
			AlignmentB = prestrcat("-", AlignmentB);
			free(alignmentTemp);
			
			i--;
		} else if (Score == ScoreUp + d) {			
			alignmentTemp = AlignmentA;
			AlignmentA = prestrcat("-", AlignmentA);
			free(alignmentTemp);			
			
			itoaTemp = itoa(B[j]);
			alignmentTemp = prestrcat(itoaTemp, AlignmentB);
			free(AlignmentB);
			free(itoaTemp);
			AlignmentB = alignmentTemp;
			
			j--;
		} else {
			printf("algorithm error?\n");
		}
	}

	while (i > 0) {
		itoaTemp = itoa(A[i]);
		alignmentTemp = prestrcat(itoaTemp, AlignmentA);
		free(AlignmentA);
		free(itoaTemp);
		AlignmentA = alignmentTemp;
		
		alignmentTemp = AlignmentB;
		AlignmentB = prestrcat("-", AlignmentB);
		free(alignmentTemp);
		i--;
	}
	
	while (j > 0) {
		alignmentTemp = AlignmentA;
		AlignmentA = prestrcat("-", AlignmentA);
		free(alignmentTemp);
		
		itoaTemp = itoa(B[j]);
		alignmentTemp = prestrcat(itoaTemp, AlignmentB);
		free(AlignmentB);
		free(itoaTemp);
		AlignmentB = alignmentTemp;
		j--;
	}
	
	int lenA;
	lenA = strlen(AlignmentA);
	int lenB;
	lenB = strlen(AlignmentB);

	char * sim1;
	sim1 = malloc(sizeof(char) * (lenA + 1));

	char * sim2;
	sim2 = malloc(sizeof(char) * (lenB + 1));
	
	int len0 = 0;
	int k = 0;
	float total = 0.0;
	float similarity = 0.0;

	len0 = (lenA > lenB) ? lenA : lenB;
	
	strncpy(sim1, AlignmentA, lenA+1);
	
	strncpy(sim2, AlignmentB, lenB+1);

	while (k < len0) {
		if (sim1[k] == sim2[k]) {
			total++;
		}
		k++;
	}
	
	similarity = total / len0 * 100.0;
	
	free(A);
	free(B);
	free(AlignmentA);
	free(AlignmentB);
	free(sim1);
	free(sim2);

	return similarity;
}













/*
// ###################################
// ##   Processing the input file   ##
// ###################################
// 

	char * strCurrentDt = NULL;
	strCurrentDt = getCurrentDt();
	printf("%s:\tGetting ready to process %s\n", strCurrentDt, inputFile_exon);
	free(strCurrentDt);

	int POTENTIAL_GUIDE_COUNT = 0;
	
	FILE * fp = NULL;
	char * line = NULL;
	line = malloc(sizeof(char) * EXON_SEQ_LINE_LEN);
	size_t len = 0;
	size_t read;
	
	char * fileName = prestrcat(dir_seq, inputFile_exon);
	fp = fopen(fileName, "r");
	if (fp == NULL) {
		printf("Error opening file: %s\n\n", prestrcat(dir_seq, inputFile_exon));
		exit(1);
	}
	free(fileName);
	
	int exonLineNumber = 0;
	
	casTarget * possibleTargetsHead = NULL;
	casTarget * possibleTargetsTail = NULL;
	casTarget * possibleTargetsCurrent = NULL;
	
	casTarget * matchesInExon = NULL;
	casTarget * matchesInExonCurrent = NULL;
	casTarget * tempMatch = NULL;
	casTarget * nextTemp = NULL;
	bool doesMatchExistAlready = false;
	
	while ((read = getline(&line, &len, fp)) != -1) {
		exonLineNumber++;
		
		///// FORWARD 
		matchesInExon = doRegex(pattern_forward, line, 22);
		
		// Iterate over the matches for the current exon
		matchesInExonCurrent = matchesInExon;
		while (matchesInExonCurrent != NULL) {
			//printf("%p -> %p\n", matchesInExonCurrent, matchesInExonCurrent->next);
			// Make sure we haven't found this target already
			possibleTargetsCurrent = possibleTargetsHead;
			doesMatchExistAlready = false;
			while (possibleTargetsCurrent != NULL) {
				if (strcmp(matchesInExonCurrent->seq_23, possibleTargetsCurrent->seq_23) == 0 && possibleTargetsCurrent->exonCount > 0) {
					doesMatchExistAlready = true;
					break;
				}
				possibleTargetsCurrent = possibleTargetsCurrent->next;
			}

			if (doesMatchExistAlready) {
				//printf("%s %s\n", matchesInExonCurrent->seq_23, possibleTargetsCurrent->seq_23);
				//printf("old\n");
				// We have already found this target: instead of creating a new 
				// obj to handle it, just add the exon number to the list of exons
				// which it exists within.
				possibleTargetsCurrent->exons = realloc(possibleTargetsCurrent->exons, ((possibleTargetsCurrent->exonCount + 1) * sizeof(int)));
				possibleTargetsCurrent->exonCount++;
				possibleTargetsCurrent->exons[possibleTargetsCurrent->exonCount - 1] = exonLineNumber;
			
			} else {
				//printf("new\n");
				POTENTIAL_GUIDE_COUNT++;
				
				matchesInExonCurrent->seq_20 = malloc(sizeof(char) * 21);
				strncpy(matchesInExonCurrent->seq_20, matchesInExonCurrent->seq_23, 20);
				matchesInExonCurrent->seq_20[20] = '\0';
				
				matchesInExonCurrent->exons = malloc(sizeof(int) * 1);
				matchesInExonCurrent->exons[0] = exonLineNumber;
				matchesInExonCurrent->exonCount = 1;
				matchesInExonCurrent->removed = 0;
				
				if (possibleTargetsHead == NULL) {
					possibleTargetsHead = matchesInExonCurrent;
					possibleTargetsTail = matchesInExonCurrent;
				} else {		
					possibleTargetsTail->next = matchesInExonCurrent;
					possibleTargetsTail = possibleTargetsTail->next;
				}
			}

			matchesInExonCurrent = matchesInExonCurrent->next;
		}		
		possibleTargetsTail->next = NULL;
		// clear matches
		//  matchesInExonCurrent = matchesInExon;
		//  while (matchesInExonCurrent != NULL) {
		//  	nextTemp = matchesInExonCurrent->next;
		//  	freeCasTarget(&matchesInExonCurrent);
		//  	matchesInExonCurrent = nextTemp;
		//  }

		
		
		strCurrentDt = getCurrentDt();
		printf("\n%s:\t%i potential targets have been identified (FORWARD).\n", strCurrentDt, POTENTIAL_GUIDE_COUNT);
		free(strCurrentDt);
		
		
		///// REVERSE		
		doesMatchExistAlready = false;
		
		matchesInExon = NULL;
		matchesInExon = doRegex(pattern_reverse, line, 22);
		// Iterate over the matches for the current exon
		matchesInExonCurrent = matchesInExon;

		char * rcTemp;
		int reverseRegexMatches = 0;
		while (matchesInExonCurrent != NULL) {
			reverseRegexMatches++;
			// Make sure we haven't found this target already
			possibleTargetsCurrent = possibleTargetsHead;
			doesMatchExistAlready = false;
			rcTemp = rc(matchesInExonCurrent->seq_23);
			int tempCounter = 0;
			while (possibleTargetsCurrent != NULL) {
				tempCounter++;
				if (strcmp(rcTemp, possibleTargetsCurrent->seq_23) == 0) {
					doesMatchExistAlready = true;
					break;
				}
				possibleTargetsCurrent = possibleTargetsCurrent->next;
			}
			printf("tempCounter: %i\n", tempCounter);
			if (doesMatchExistAlready) {
				// We have already found this target: instead of creating a new 
				// obj to handle it, just add the exon number to the list of exons
				// which it exists within.
				possibleTargetsCurrent->exons = realloc(possibleTargetsCurrent->exons, ((possibleTargetsCurrent->exonCount + 1) * sizeof(int)));
				possibleTargetsCurrent->exonCount++;
				possibleTargetsCurrent->exons[possibleTargetsCurrent->exonCount - 1] = exonLineNumber;
			
			} else {
				POTENTIAL_GUIDE_COUNT++;

				matchesInExonCurrent->seq_20 = malloc(sizeof(char) * 21);
				strncpy(matchesInExonCurrent->seq_20, matchesInExonCurrent->seq_23, 20);
				matchesInExonCurrent->seq_20[20] = '\0';
				
				matchesInExonCurrent->exons = malloc(sizeof(int) * 1);
				matchesInExonCurrent->exons[0] = exonLineNumber;
				matchesInExonCurrent->exonCount = 1;
				matchesInExonCurrent->removed = 0;
				
				if (possibleTargetsHead == NULL) {
					possibleTargetsHead = matchesInExonCurrent;
					possibleTargetsTail = matchesInExonCurrent;
				} else {		
					possibleTargetsTail->next = matchesInExonCurrent;
					possibleTargetsTail = possibleTargetsTail->next;
				}
			}

			free(rcTemp);
			matchesInExonCurrent = matchesInExonCurrent->next;
		}
		printf("reverseRegexMatches: %i\n", reverseRegexMatches);
		possibleTargetsTail->next = NULL;
		// clear matches
		// matchesInExonCurrent = matchesInExon;
		// while (matchesInExonCurrent != NULL) {
		// 	//printf("free...\n");
		// 	nextTemp = matchesInExonCurrent->next;
		// 	freeCasTarget(&matchesInExonCurrent);
		// 	matchesInExonCurrent = nextTemp;
		// }

	}
	free(line);
	fclose(fp);
	//

	strCurrentDt = getCurrentDt();
	printf("\n%s:\t%i potential targets have been identified.\n", strCurrentDt, POTENTIAL_GUIDE_COUNT);
	free(strCurrentDt);
	*/