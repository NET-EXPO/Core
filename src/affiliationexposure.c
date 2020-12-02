//
//  affiliationexposure.c
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
//


#include "affiliationexposure.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>


void read_y_mode_file(char* path, double** yArray, int* first_mode_array, int* second_mode_array, bool HAS_HEADER, char *first_mode_flag){
    bool hasHeader = HAS_HEADER;
    int first_index = -1;
    int second_index = -1;
    FILE * yFile = fopen(path, "r");
    char line[1024];

    while(fgets(line, 1024, yFile)){
        if(hasHeader == false)
        {
            char * line_temp = strdup(line);
            //remove new line
            if(line_temp[strlen(line_temp)-1] == '\n'){
                line_temp[strlen(line_temp)-1] = '\0';
            }
            
            
            if(line_temp != NULL){
                
                char * tokens = strtok(line_temp, ",");
                char * token_array[3];
                int i=0;
                while(tokens != NULL){
                    if(i<3) token_array[i++] = tokens;
                    tokens = strtok(NULL, ",");
                }
                double id = atof(token_array[0]);
                double y_value = atof(token_array[1]);
                //if row belongs to first mode, add it
                char* mode_token = token_array[2];
                if(mode_token[strlen(mode_token)-1] == '\n'){
                    mode_token[strlen(mode_token)-1] = '\0';
                }
                
                if(mode_token[strlen(mode_token)-1] == '\r'){
                    mode_token[strlen(mode_token)-1] = '\0';
                }
                
                if(strncmp(first_mode_flag, mode_token, strlen(mode_token)) == 0){
                    
                    first_index++;
                    yArray[first_index][0] = id;
                    yArray[first_index][1] = y_value;
                    
                    first_mode_array[first_index] = id;

                }
                else{

                    second_index++;
                    second_mode_array[second_index] = id;
                }
                
                free(tokens);
                
            }

            free(line_temp);
            
        }
        else
        {
            hasHeader = false;
        }
    }
    fclose(yFile);
    
}

void init_a_matrix(double** wMatrix, struct ModeInfo *mode_data){
    
    for(int row=0; row < mode_data->first_mode_indices; ++row){
        for(int column=0; column< mode_data->second_mode_indices; ++column){
            wMatrix[row][column] = 0;

        }
    }
}

void setup_mode_data(char* yPathFile, bool HAS_HEADER, struct ModeInfo *mode_data){
    //get number of first mode and number of second mode
    int num_modes = 0;;
    bool hasHeader = HAS_HEADER;
    FILE * yFile = fopen(yPathFile, "r");
    char line[1024];
    int num_rows = 0;
    
    while(fgets(line, 1024,yFile)){
        if(hasHeader == false){
            num_rows++;
            char * line_temp = strdup(line);
            //remove new line
            if(line_temp[strlen(line_temp)-1] == '\n')
                line_temp[strlen(line_temp)-1] = '\0';
            
            //printf("line: %s\n", line_temp);
            char* token = strtok(line_temp, ",");
            char* last_token = NULL;
            while(token != NULL){
                
                token = strtok(NULL, ",");
                if(token != NULL){
                    
                    last_token = token;
                    
                }
                else{
                    //using strncmp because last token has mysterious added character
                    
                    if(last_token[strlen(last_token)-1] == '\n'){
                        last_token[strlen(last_token)-1] = '\0';
                    }
                    
                    if(last_token[strlen(last_token)-1] == '\r'){
                        last_token[strlen(last_token)-1] = '\0';
                    }
                    
                    if(strncmp(mode_data->first_mode_flag, last_token, strlen(last_token)) == 0){
                        num_modes++;
                    }
                    
                }
            }
            
            free(token);
            free(line_temp);
            
        }
        else{
            hasHeader = false;
        }
    }
    
    mode_data->second_mode_indices = num_rows - num_modes;
    mode_data->first_mode_indices = num_modes;

}



void populate_a_matrix(double** wMatrix, char* netPathFile, bool HAS_HEADER, struct ModeInfo *mode_data, int* f_mode_array, int* s_mode_array){
   
    int idx = 0;
    int i; int j;
    int *row_ = NULL; int *column_;
    int row; int column;
    
    bool hasHeader = HAS_HEADER;
    FILE *nFile = fopen(netPathFile, "r");
    char line[1024];
    while(fgets(line,1024,nFile)){
        if(hasHeader == false){
            char *line_temp = strdup(line);
            //remove new line
            if(line_temp[strlen(line_temp)-1] == '\n')
                line_temp[strlen(line_temp)-1] = '\0';
            
            if(line_temp != NULL){
                char *tokens = strtok(line_temp, ",");
                char *token_array[2];
                idx=0;
                while(tokens != NULL){
                    token_array[idx++] = tokens;
                    tokens = strtok(NULL, ",");
                }
                
                row = atoi(token_array[0]);
                column = atoi(token_array[1]);
                //puts("iteriating through...\n");
                for(i=0; i < mode_data->first_mode_indices; ++i){
                    row_ = &f_mode_array[i];
                    if(row == *row_){
                      for(j=0; j< mode_data->second_mode_indices; ++j){
                        column_ = &s_mode_array[j];
                          if(column == *column_){
                              wMatrix[i][j]++;
                          }
                      }
                    }
                    else if(column == *row_){
                       for(j=0; j< mode_data->second_mode_indices; ++j){
                           column_ = &s_mode_array[j];
                           if(row == *column_){
                               wMatrix[i][j]++;
                           }
                       }
                    }
                    
                    
                }

                free(tokens);
                
            }
            
            free(line_temp);
            
        }
        else{
            hasHeader = false;
        }
    }
    
    fclose(nFile);
    
}

double** generateAMatrix(char* netPathFile, int* f_mode_array, int* s_mode_array, bool HAS_HEADER, struct ModeInfo *mode_data){

    double** aMatrix = malloc(sizeof(double*) * mode_data->first_mode_indices);
    for(int i=0; i < mode_data->first_mode_indices; ++i){
        aMatrix[i] = malloc(sizeof(double) * mode_data->second_mode_indices);
    }

    init_a_matrix(aMatrix, mode_data);
    
    
    populate_a_matrix(aMatrix, netPathFile, HAS_HEADER, mode_data, f_mode_array, s_mode_array);
    
    
    return aMatrix;
}

double** generateAPrimeMatrix(double** aMatrix, struct ModeInfo *mode_data){
    
    double** aPrimeMatrix = malloc(sizeof(double*) * mode_data->second_mode_indices);
    for(int i=0; i < mode_data->second_mode_indices; ++i){
        aPrimeMatrix[i] = malloc(sizeof(double) * mode_data->first_mode_indices);
    }
    
    for(int i=0; i < mode_data->second_mode_indices; ++i){
        for(int j=0; j< mode_data->first_mode_indices; ++j){
            
            aPrimeMatrix[i][j] = aMatrix[j][i];
            
        }
    }
    
    return aPrimeMatrix;
    
}

double** transposeAMatrix(double** aMatrix, double** aPrimeMatrix, struct ModeInfo *mode_data){
    
    
    double** cMatrix = malloc(sizeof(double*) * mode_data->first_mode_indices);
    for(int i=0; i < mode_data->first_mode_indices; ++i){
        cMatrix[i] = malloc(sizeof(double) * mode_data->first_mode_indices);
    }
    
    for(int i=0; i < mode_data->first_mode_indices; ++i){
        for(int j=0; j< mode_data->first_mode_indices; ++j){
            cMatrix[i][j] = 0;
            for(int k=0; k < mode_data->second_mode_indices; ++k){
                cMatrix[i][j] += aMatrix[i][k] * aPrimeMatrix[k][j];
            }
        }
    }
    
    
    return cMatrix;
}

void calc_affiliation_exposure(double** cMatrix, double** eArray, double** yArray, struct ModeInfo *mode_data){

    int i; int j;
    //zero out the diag
    for(i=0; i < mode_data->first_mode_indices; ++i){
        for(j=0; j < mode_data->first_mode_indices; ++j){
            if(i==j){
                cMatrix[i][j] = 0;
            }
        }
    }
    
    double s_row = 0;
    double s_mul = 0;
    
    for(i=0; i < mode_data->first_mode_indices; ++i){
        s_row = 0;
        s_mul = 0;
        for (j = 0; j < mode_data->first_mode_indices; j++) {
            s_row += cMatrix[i][j];
            s_mul += cMatrix[i][j] * yArray[j][1];
            
        }
        if(s_row != 0){
            eArray[i][1] = s_mul / s_row;
        }
        else{
            eArray[i][1] = 0;
        }
    }
    
}

void print_aff_exposure_values(double** y_array, double** e_array, struct ModeInfo *mode_data){
    for(int i=0; i < mode_data->first_mode_indices; ++i)
    {
        printf("Exposure of %d:\t %f \n", (int)y_array[i][0],e_array[i][1]);
    }
}

