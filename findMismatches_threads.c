
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h>




#define MAX_LINE_LENGTH 30
#define SEQ_LENGTH 20
#define LINES_PREALLOCATED_TARGETS 64
#define LINES_PREALLOCATED_OFFSITE 134217728
#define MAX_NUMBER_MISMATCHES 4




/**************************************/
/**            arg_struct            **/
/**************************************/

/*
 Structure used to pass arguments to threads
 */

struct arg_struct {
    double* global_scores;
    char** targets;
    char** sites;
    int nb_targets;
    int nb_sites;
    int targetOffset;
    int nbThreads;
    int threshold;
};




/**************************************/
/**         printTimestamp()         **/
/**************************************/

/*
 Prints a custom string with the current time and a given message
 Input: the message we want to include
 Output: none
 */

int printTimestamp(char* infoString) {
    time_t rawtime;
    struct tm * timeinfo;
    char buffer1 [15];
    char buffer2 [80];

    time (&rawtime);
    timeinfo = localtime (&rawtime);
    strftime (buffer1,80,"%H:%M:%S ->\t",timeinfo);
    
    strcpy(buffer2, buffer1);
    strcat(buffer2, infoString);
    puts (buffer2);

    return 0;
}




/**************************************/
/**          single_score()          **/
/**************************************/

/*
 Computes the score of a single off-target site
 Input: array containing all the mismatches between this site and the target, length of the array
 Output: score
 */

double single_score(int* mismatch_array, int length) {
    int i;
    double T1=1.0, T2, T3, d=0.0, score;
    double M[] = {0.0, 0.0, 0.014, 0.0, 0.0, 0.395, 0.317, 0.0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583};

    /* 1st term */
    for(i=0; i<length; ++i)
        T1 = T1*(1.0-M[mismatch_array[i]]);

    /* 2nd term */
    if(length==1)
        d = 19.0;
    else {
        for(i=0; i<length-1; ++i)
            d += mismatch_array[i+1]-mismatch_array[i];
        d = d/(length-1);
    }
    T2 = 1.0 / ((19.0-d)/19.0 * 4.0 + 1);

    /* 3rd term */
    T3 = 1.0 / (length*length);

    /* Total score */
    score = T1*T2*T3*100;
    return score;
}




/**************************************/
/**         thread_scoring()         **/
/**************************************/

/*
 Computes the global score of targets (using threads)
 Input: data from the main function (targets, off-target sites, etc.)
 Output: updated data
 */

void *thread_scoring(void *arguments) {
    struct arg_struct *args = arguments;
    double* global_scores = args->global_scores;
    char** targets = args->targets;
    char** sites = args->sites;
    int nb_targets = args->nb_targets;
    int nb_sites = args->nb_sites;
    int targetOffset = args->targetOffset;
    int nb_threads = args->nbThreads;
    int threshold = args->threshold;
    int* mismatches;
    int nb_mismatches;
    double sum_local_scores;
    int i,j,k;
    int nb_sites_that_are_close_enough = 0;
    double total_score = 1.0;
    
    double maximum_sum = (10000.0 - threshold*100) / threshold;

    mismatches = (int*)malloc(SEQ_LENGTH*sizeof(int));
    
    for(i=targetOffset; i<nb_targets; i+=nb_threads) {
    
        sum_local_scores = 0.0;
        nb_sites_that_are_close_enough = 0;
        total_score = 1.0;
        
        for(j=0; j<nb_sites; ++j) {
            
            /* Find mismatches */
            nb_mismatches = 0;
            k = 0;
            while(k<SEQ_LENGTH && nb_mismatches<=MAX_NUMBER_MISMATCHES) {
                if(targets[i][k]!=sites[j][k]) {
                    mismatches[nb_mismatches] = k;
                    ++nb_mismatches;
                }
                ++k;
            }
            
            /* Score the mismatches (if they are similar enough to target) */
            if(nb_mismatches>0 && nb_mismatches<=MAX_NUMBER_MISMATCHES) {
                sum_local_scores += single_score(mismatches,nb_mismatches);
                
                /* if the current sum is already higher than the limit,  */
                /*   then the score will always be below the threshold,  */
                /*      so we can stop early.                            */
                if(sum_local_scores > maximum_sum)
                    break;
            }
        }
        
        global_scores[i] = 10000.0/(100.0+sum_local_scores);

    }
    free(mismatches);
    
    pthread_exit(NULL);
}




/**************************************/
/**              main()              **/
/**************************************/

int main(int argc, char* argv[]) {
    int lines_allocated_targets = LINES_PREALLOCATED_TARGETS;
    int lines_allocated_offsite = LINES_PREALLOCATED_OFFSITE;
    int max_line_len = MAX_LINE_LENGTH;
    int i;
    FILE *fp;
    int nb_targets=0, nb_sites=0, global_score_threshold;
    char *file_targets, *file_offsite, *file_output, **targets, **offsite;
    struct arg_struct *args;
    pthread_attr_t attr;
    pthread_t *thread;
    int rc;
    int t;
    void *status;
    double *global_scores;
    int nb_threads;

    if(argc!=6) {
        fprintf(stderr,"*** Error. Exactly five arguments are needed.\n");
        exit(6);
    }
    
    nb_threads = atoi(argv[1]);
    file_targets = argv[2];
    file_offsite = argv[3];
    file_output = argv[4];
    global_score_threshold = atoi(argv[5]);
    printTimestamp("Ready to start.");
    printf("\t\tThe scoring will run on %d threads.\n", nb_threads);
    
    
    

    
    /**************************************/
    /**       READING THE TARGETS        **/
    /**************************************/


    printTimestamp("Starting to read targets to score.");
    
    /* Allocate lines of text */
    targets = (char **)malloc(sizeof(char*)*lines_allocated_targets);
    if (targets==NULL) {
        fprintf(stderr,"*** Out of memory (1).\n");
        exit(1);
    }
    
    fp = fopen(file_targets, "r");
    if (fp == NULL) {
        fprintf(stderr,"*** Error opening file with sequences.\n");
        exit(2);
    }
    
    for (i=0;1;i++) {
        
        /* Have we gone over our line allocation? */
        if (i >= lines_allocated_targets) {
            int new_size;
            
            /* Double our allocation and re-allocate */
            new_size = lines_allocated_targets*2;
            fprintf(stdout,"\t\tReallocating space for target sequences.\n");
            targets = (char **)realloc(targets,sizeof(char*)*new_size);
            if (targets==NULL) {
                fprintf(stderr,"*** Out of memory.\n");
                exit(3);
            }
            lines_allocated_targets = new_size;
        }
        
        /* Allocate space for the next line */
        targets[i] = (char*)malloc(max_line_len*sizeof(char));
        if (targets[i]==NULL) {
            fprintf(stderr,"*** Out of memory (3).\n");
            exit(4);
        }
        if (fgets(targets[i],max_line_len-1,fp)==NULL)
            break;
        
        /* Get rid of CR or LF at end of line */
        targets[i][SEQ_LENGTH]='\0';
        ++nb_targets;
    }
    
    fclose(fp);
    

    
    /**************************************/
    /**   READING THE OFF-TARGET SITES   **/
    /**************************************/
    

    printTimestamp("Starting to read off-target sites.");

    /* Allocate lines of text */
    offsite = (char **)malloc(sizeof(char*)*lines_allocated_offsite);
    if (offsite==NULL) {
        fprintf(stderr,"*** Out of memory (1).\n");
        exit(1);
    }
    
    fp = fopen(file_offsite, "r");
    if (fp == NULL) {
        fprintf(stderr,"*** Error opening file with off-target sites.\n");
        exit(2);
    }
    
    for (i=0;1;i++) {
        
        /* Have we gone over our line allocation? */
        if (i >= lines_allocated_offsite) {
            int new_size;
            
            /* Double our allocation and re-allocate */
            new_size = lines_allocated_offsite*2;
            fprintf(stdout,"\t\tReallocating space for off-target sites.\n");
            offsite = (char **)realloc(offsite,sizeof(char*)*new_size);
            if (offsite==NULL) {
                fprintf(stderr,"*** Out of memory.\n");
                exit(3);
            }
            lines_allocated_offsite = new_size;
        }
        
        /* Allocate space for the next line */
        offsite[i] = (char*)malloc(max_line_len*sizeof(char));
        if (offsite[i]==NULL) {
            fprintf(stderr,"*** Out of memory (3).\n");
            exit(4);
        }
        if (fgets(offsite[i],max_line_len-1,fp)==NULL)
            break;
        
        /* Get rid of CR or LF at end of line */
        offsite[i][SEQ_LENGTH]='\0';
        ++nb_sites;
    }
    
    fclose(fp);
    

    
    /**************************************/
    /**            PROCESSING            **/
    /**************************************/


    printTimestamp("Ready to proceed.");
    printf("\t\tData: %d targets, %d off-target sites.\n", nb_targets, nb_sites);

    printTimestamp("Starting to score.");
    
    thread = (pthread_t*)malloc(nb_threads*sizeof(pthread_t));
    global_scores = (double*)malloc(nb_targets*sizeof(double));
    
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    args=(struct arg_struct*)malloc(nb_threads*sizeof(struct arg_struct));
    for(t=0; t<nb_threads; t++) {
        args[t].global_scores = global_scores;
        args[t].targets = targets;
        args[t].sites = offsite;
        args[t].nb_targets = nb_targets;
        args[t].nb_sites = nb_sites;
        args[t].targetOffset = t;
        args[t].nbThreads = nb_threads;
        args[t].threshold = global_score_threshold;
    }
    
    for(t=0; t<nb_threads; t++) {
        rc = pthread_create(&thread[t], &attr, thread_scoring, (void *)&args[t]);
        if (rc) {
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
    }
    
    /* Free attribute and wait for the other threads */
    pthread_attr_destroy(&attr);
    for(t=0; t<nb_threads; t++) {
        rc = pthread_join(thread[t], &status);
        if (rc) {
            printf("ERROR; return code from pthread_join() is %d\n", rc);
            exit(-1);
        }
    }
    

    fp = fopen(file_output, "w");
    if (fp == NULL) {
        fprintf(stderr,"*** Error opening output file.\n");
        exit(2);
    }
    
    for(i=0; i<nb_targets; ++i) {
        printf("%s\t%f\n",targets[i],global_scores[i]);
        fprintf(fp,"%s\t%f\n",targets[i],global_scores[i]);
    }
    
    fclose(fp);
    

    /**************************************/
    /**           FINAL STEPS            **/
    /**************************************/
    
    printTimestamp("Freeing the memory.");
    
    /* Good practice to free memory */
    free(args);
    free(thread);
    free(global_scores);
    for (i=nb_targets-1;i>=0;i--)
        free(targets[i]);
    free(targets);
    for (i=nb_sites-1;i>=0;i--)
        free(offsite[i]);
    free(offsite);

    printTimestamp("Done.");
    return 0;
}
