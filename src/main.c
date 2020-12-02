//
//  main.c
//  NETEXPO-CORE
//
//  Created by Tuan Amith
//  
//  Copyright (c) 2020 School of Biomedical Informatics,
//  The University of Texas Health Science Center at Houston. All Rights Reserved.
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//
//  **************
//  main.c demonstrates how to wield networkexposure.c and affiliationexposure.c and also
//  perform an evaluation of its performance to handle large network files. The coding 
//  below uses either multi-threading of exposure computation to handle multiple files or
//  just a single thread for one network or a handful of network files. See the main function
//  for further details on executing either the multi-threaded or single threaded execution.
//  ************** 

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <pthread.h>
#include "networkexposure.h"
#include "affiliationexposure.h"


pthread_mutex_t lock;

typedef enum{
    NETWORK,
    AFFILIATION,
    NETWORK_THREAD,
    AFFILIATION_THREAD
} test_type;

struct thread_config{
    bool *header_config;
    char netFile[1024];
    char yFile[1024];
    char first_mode_flag[1024];
    char file_output[1024];
    char target_file[1024];
};

void *network_exposure_worker(void *args){
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    struct thread_config *config_args = args;
    bool HAS_HEADER = config_args->header_config;
    
    pthread_mutex_lock(&lock);
    int matrix_dim = getDimension(config_args->yFile, HAS_HEADER);
    pthread_mutex_unlock(&lock);
    
    double **yArray = malloc(sizeof(double*) * matrix_dim);
    for(int i=0; i < matrix_dim ; ++i){
        yArray[i] = malloc(sizeof(double) * 2);
    }
    
    double **wMatrix = malloc(sizeof(double*) * matrix_dim);
    for(int i=0; i < matrix_dim; ++i){
        wMatrix[i] = malloc(sizeof(double) * matrix_dim);
    }
    
    init_w_matrix(wMatrix, &matrix_dim);
    
    double **eArray = malloc(sizeof(double*) * matrix_dim);
    for(int i=0; i<matrix_dim; ++i){
        eArray[i] = malloc(sizeof(double) * 2);
    }
    
    //read the files
    
    //... read and store the y file
    pthread_mutex_lock(&lock);
    read_y_file(config_args->yFile, yArray, HAS_HEADER);
    pthread_mutex_unlock(&lock);
    //... read and store the net file
    pthread_mutex_lock(&lock);
    read_net_file(config_args->netFile, wMatrix, HAS_HEADER);
    pthread_mutex_unlock(&lock);
    //compute
    calc_network_exposure(wMatrix, yArray, eArray, matrix_dim);
    
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    
    //uncomment to see output
    print_net_exposure_values(yArray, eArray, matrix_dim);
    
    FILE* outputFile = fopen(config_args->file_output, "a");
    fprintf(outputFile, "%s,%ld,%ld,%ld,%ld,%ld\n",config_args->target_file,start.tv_nsec, end.tv_nsec, start.tv_nsec/1000000, end.tv_nsec/1000000,
            (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000);
    fclose(outputFile);
    
    //deallocate
    for(int i=0; i< matrix_dim; ++i){
        free(yArray[i]);
    }
    free(yArray);
    
    for(int i=0; i< matrix_dim; ++i){
        free(wMatrix[i]);
    }
    free(wMatrix);
    
    for(int i=0; i<matrix_dim; ++i){
        free(eArray[i]);
    }
    free(eArray);

    return NULL;
}

void *affiliation_exposure_worker(void *args){
    struct timespec start, end;

    struct thread_config *config_args = args;
    
    int i;
    bool HAS_HEADER = config_args->header_config;
    struct ModeInfo *mode_data_thread = malloc(sizeof *mode_data_thread);
    mode_data_thread->first_mode_flag = config_args->first_mode_flag;
    
    pthread_mutex_lock(&lock);
    setup_mode_data(config_args->yFile, HAS_HEADER, mode_data_thread);
    pthread_mutex_unlock(&lock);
    
    int *first_mode_array = (int*)calloc(mode_data_thread->first_mode_indices, sizeof(int));
    int *second_mode_array = (int*)calloc(mode_data_thread->second_mode_indices, sizeof(int));
    
    double** yArray = malloc(sizeof(double*) * mode_data_thread->first_mode_indices);
    for(i=0; i < mode_data_thread->first_mode_indices; ++i){
        yArray[i] = malloc(sizeof(double) * 2);
    }
    
    double** eArray = malloc(sizeof(double*) * mode_data_thread->first_mode_indices);
    for (i=0; i< mode_data_thread->first_mode_indices; ++i) {
        eArray[i] = malloc(sizeof(double) * 2);
    }
    
    pthread_mutex_lock(&lock);
    read_y_mode_file(config_args->yFile, yArray, first_mode_array, second_mode_array, HAS_HEADER, mode_data_thread->first_mode_flag);
    pthread_mutex_unlock(&lock);
    
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    
    pthread_mutex_lock(&lock);
    double **aMatrix = generateAMatrix(config_args->netFile, first_mode_array, second_mode_array, HAS_HEADER, mode_data_thread);
    pthread_mutex_unlock(&lock);
    
    double **aPrimeMatrix = generateAPrimeMatrix(aMatrix, mode_data_thread);
    
    double **cMatrix = transposeAMatrix(aMatrix, aPrimeMatrix, mode_data_thread);
    
    calc_affiliation_exposure(cMatrix, eArray, yArray, mode_data_thread);
    
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    
    //save_time(config_args->target_file,&start, &end, config_args->file_output);
    
    FILE* outputFile = fopen(config_args->file_output, "a");
    fprintf(outputFile, "%s,%ld,%ld,%ld,%ld,%ld\n",config_args->target_file,start.tv_nsec, end.tv_nsec, start.tv_nsec/1000000, end.tv_nsec/1000000,
            (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000);
    fclose(outputFile);
    
    //uncomment to see output
    print_aff_exposure_values(yArray, eArray, mode_data_thread);

    //deallocate
    for(i=0; i< mode_data_thread->first_mode_indices; ++i){
        free(yArray[i]);
    }
    free(yArray);
    
    
    free(first_mode_array);
    free(second_mode_array);
    
    for(i=0; i< mode_data_thread->first_mode_indices; ++i){
        free(aMatrix[i]);
    }
    free(aMatrix);
    
    
    for(i=0; i < mode_data_thread->second_mode_indices; ++i){
        free(aPrimeMatrix[i]);
    }
    free(aPrimeMatrix);
    
    
    for(i=0; i< mode_data_thread->first_mode_indices; ++i){
        free(cMatrix[i]);
    }
    free(cMatrix);
    
    
    for(i=0; i< mode_data_thread->first_mode_indices; ++i){
        free(eArray[i]);
    }
    free(eArray);
    
    
    free(mode_data_thread);
    config_args = NULL;

    return NULL;
    
}

void execute_affiliation_exposure(bool FILE_HAS_HEADER, char *netFile, char *yFile, char *first_mode_flag){
    int i;
    bool HAS_HEADER = FILE_HAS_HEADER;
    struct ModeInfo mode_data;
    
    mode_data.first_mode_flag = first_mode_flag;
    
    setup_mode_data(yFile, HAS_HEADER, &mode_data);
    
    int *first_mode_array = (int*)calloc(mode_data.first_mode_indices, sizeof(int));
    int *second_mode_array = (int*)calloc(mode_data.second_mode_indices, sizeof(int));
    
    double** yArray = malloc(sizeof(double*) * mode_data.first_mode_indices);
    for(i=0; i < mode_data.first_mode_indices; ++i){
        yArray[i] = malloc(sizeof(double) * 2);
    }
    
    double** eArray = malloc(sizeof(double*) * mode_data.first_mode_indices);
    for (i=0; i< mode_data.first_mode_indices; ++i) {
        eArray[i] = malloc(sizeof(double) * 2);
    }
    
    //... read and store y file
    read_y_mode_file(yFile, yArray, first_mode_array, second_mode_array, HAS_HEADER, mode_data.first_mode_flag);
    
    double **aMatrix = generateAMatrix(netFile, first_mode_array, second_mode_array, HAS_HEADER, &mode_data);
    
    double **aPrimeMatrix = generateAPrimeMatrix(aMatrix, &mode_data);
    
    double **cMatrix = transposeAMatrix(aMatrix, aPrimeMatrix, &mode_data);
    
    
    calc_affiliation_exposure(cMatrix, eArray, yArray, &mode_data);
    
    print_aff_exposure_values(yArray, eArray, &mode_data);
    
    //deallocate
    for(i=0; i< mode_data.first_mode_indices; ++i){
        free(yArray[i]);
    }
    free(yArray);
    
    
    free(first_mode_array);
    free(second_mode_array);
    
    for(i=0; i< mode_data.first_mode_indices; ++i){
        free(aMatrix[i]);
    }
    free(aMatrix);
    
    
    for(i=0; i < mode_data.second_mode_indices; ++i){
        free(aPrimeMatrix[i]);
    }
    free(aPrimeMatrix);
    
    
    for(i=0; i<mode_data.first_mode_indices; ++i){
        free(cMatrix[i]);
    }
    free(cMatrix);
    
    
    for(i=0; i<mode_data.first_mode_indices; ++i){
        free(eArray[i]);
    }
    free(eArray);
    
    
}

void execute_network_exposure(bool FILE_HAS_HEADER, char *netFile, char *yFile){
    bool HAS_HEADER = FILE_HAS_HEADER;
    int matrix_dim = getDimension(yFile, HAS_HEADER);
    
    //intialize the arrays and matrices - dynamic allocation
    
    double **yArray = malloc(sizeof(double*) * matrix_dim);
    for(int i=0; i < matrix_dim ; ++i){
        yArray[i] = malloc(sizeof(double) * 2);
    }
    
    double **wMatrix = malloc(sizeof(double*) * matrix_dim);
    for(int i=0; i < matrix_dim; ++i){
        wMatrix[i] = malloc(sizeof(double) * matrix_dim);
    }
    
    init_w_matrix(wMatrix, &matrix_dim);
    
    double **eArray = malloc(sizeof(double*) * matrix_dim);
    for(int i=0; i<matrix_dim; ++i){
        eArray[i] = malloc(sizeof(double) * 2);
    }
    
    //read the files
    
    //... read and store the y file
    read_y_file(yFile, yArray, HAS_HEADER);
    
    //... read and store the net file
    read_net_file(netFile, wMatrix, HAS_HEADER);
    
    //compute
    calc_network_exposure(wMatrix, yArray, eArray, matrix_dim);
    
    //uncomment to see output
    print_net_exposure_values(yArray, eArray, matrix_dim);
    
    //deallocate
    for(int i=0; i< matrix_dim; ++i){
        free(yArray[i]);
    }
    free(yArray);
    
    for(int i=0; i< matrix_dim; ++i){
        free(wMatrix[i]);
    }
    free(wMatrix);
    
    for(int i=0; i<matrix_dim; ++i){
        free(eArray[i]);
    }
    free(eArray);
    
}

void save_time(char *targetFile,struct timespec *start,struct timespec *end, char* fileName){
    
    FILE* outputFile = fopen(fileName, "a");
    fprintf(outputFile, "%s,%ld,%ld,%ld,%ld,%ld\n",targetFile,start->tv_nsec, end->tv_nsec, start->tv_nsec/1000000, end->tv_nsec/1000000,
            (end->tv_sec - start->tv_sec) * 1000000 + (end->tv_nsec - start->tv_nsec) / 1000);
    fclose(outputFile);
    
}

void init_save_file(char* fileName){
    FILE* outputFile = fopen(fileName, "a");
    fprintf(outputFile, "file_name, start(nano), end(nano), start(milisec), end(milisec), completed(microsec)\n");
    fclose(outputFile);
}

void runBatchTest(char *target_directory, char *fileName, test_type test, char *first_mode_flag, bool files_has_headers){
    
    
    int file_count =0;
    struct thread_config *config_array = NULL;
    if(test == AFFILIATION_THREAD || test == NETWORK_THREAD){
        config_array = (struct thread_config*)calloc(1,sizeof(struct thread_config));
    }
    
    DIR *data_folder;
    struct dirent *file;
    struct stat filestat;
    data_folder = opendir(target_directory);
    
    struct timespec start, end;
    
    char* yfile_suffix = NULL;
    if(test == AFFILIATION || test == AFFILIATION_THREAD){
        yfile_suffix = "_y_mode.csv";
    }
    else if(test == NETWORK || test == NETWORK_THREAD){
        yfile_suffix = "_y.csv";
    }
    
    init_save_file(fileName);
    
    if(data_folder == NULL)
    {
        perror("Unable to read directory");
        exit(0);

    }
    
    char netPath[1024];
    char yPath[1024];
    while( (file=readdir(data_folder)) )
    {
        
        lstat(file->d_name,&filestat);
        if( S_ISDIR(filestat.st_mode) )
        {
            if(strcmp(".",file->d_name) == 0 ||
               strcmp("..",file->d_name) == 0 ||
               strcmp(".DS_Store",file->d_name) == 0)
                continue;
            
            
            if(strstr(file->d_name, yfile_suffix) == NULL){
                printf("File : %s\n",file->d_name);

                char y_file_name[file->d_namlen -4];
                memcpy(y_file_name, file->d_name, file->d_namlen-4);
                y_file_name[file->d_namlen-4] = '\0';
                
                strcat(y_file_name, yfile_suffix);

                strcpy(netPath, target_directory);
                strcpy(yPath, target_directory);
                
                strcat(netPath, file->d_name);
                strcat(yPath,y_file_name);
                
                
                if(test == NETWORK){
                    
                    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
                    
                    execute_network_exposure(files_has_headers, netPath, yPath);
                    
                    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
                    save_time(file->d_name,&start, &end, fileName);
                }
                else if(test == NETWORK_THREAD){
                    config_array[file_count].header_config = &files_has_headers;
                    snprintf(config_array[file_count].netFile, sizeof(config_array[file_count].netFile), "%s", netPath);
                    snprintf(config_array[file_count].yFile, sizeof(config_array[file_count].yFile), "%s", yPath);
                    snprintf(config_array[file_count].file_output, sizeof(config_array[file_count].file_output), "%s", fileName);
                    snprintf(config_array[file_count].target_file, sizeof(config_array[file_count].target_file),"%s", file->d_name);
                    file_count++;
                    config_array = (struct thread_config*)realloc(config_array, (file_count+1)*sizeof(struct thread_config));
                }
                else if (test == AFFILIATION){
                   
                    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
                    
                    execute_affiliation_exposure(files_has_headers, netPath, yPath, first_mode_flag);
                    
                    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
                    save_time(file->d_name,&start, &end, fileName);

                }
                else if(test == AFFILIATION_THREAD){

                    config_array[file_count].header_config = &files_has_headers;
                    snprintf(config_array[file_count].netFile, sizeof(config_array[file_count].netFile), "%s", netPath);
                    snprintf(config_array[file_count].yFile, sizeof(config_array[file_count].yFile), "%s", yPath);
                    snprintf(config_array[file_count].first_mode_flag, sizeof(config_array[file_count].first_mode_flag), "%s", first_mode_flag);
                    snprintf(config_array[file_count].file_output, sizeof(config_array[file_count].file_output), "%s", fileName);
                    snprintf(config_array[file_count].target_file, sizeof(config_array[file_count].target_file),"%s", file->d_name);
                    file_count++;
                    config_array = (struct thread_config*)realloc(config_array, (file_count+1)*sizeof(struct thread_config));
                }
                
                
            }
            

        }
        

    }
    
    closedir(data_folder);
    
    if(test == AFFILIATION_THREAD){
        pthread_t worker[file_count];
        for(int x =0; x < file_count; ++x){
            pthread_create(&worker[x], NULL, affiliation_exposure_worker, &config_array[x]);
        }
    }
    if(test == NETWORK_THREAD){
        pthread_t worker[file_count];
        for(int x=0; x < file_count; ++x){
            pthread_create(&worker[x], NULL, network_exposure_worker, &config_array[x]);
        }
    }
    /*
    for(int x = 0; x < file_count; ++x){
        pthread_join(worker[x], NULL);
    }
    */
}



int main(void){

    //Uncomment the following lines 5 lines to utilize multi-threading
    /*
    if (pthread_mutex_init(&lock, NULL) != 0)
    {
        printf("\n thread initializing failed\n");
        return 1;
    }
    */

    //Define the path for output of the performance results
    char fileName [] = "path_to_performance_output.csv";
    
    //Define the path location of the network files
    char target_directory [] = "path_to_folder_network_file_and_attribute_file";
    
    //runBatchTest(target_directory, fileName, NETWORK_THREAD, NULL, true);
    //runBatchTest(target_directory, fileName, AFFILIATION_THREAD, "1",true);
    
    runBatchTest(target_directory, fileName, NETWORK, NULL, true);
    //runBatchTest(target_directory, fileName, AFFILIATION, "1",true);
    
    //Uncomment the following lines 2 lines to utilize multi-threading
    /*
    pthread_mutex_destroy(&lock);
    pthread_exit(NULL);
    */

    return(0);
    
}

