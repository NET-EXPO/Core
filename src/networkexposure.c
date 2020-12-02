//
//  networkexposure.c
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

#include "networkexposure.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>




void calc_network_exposure(double** wMatrix, double ** yArray, double ** eArray, int matrix_dim){
    
    for(int row=0; row < matrix_dim; ++row){
        double s_row = 0;
        double s_mul = 0;
        for(int column = 0 ; column < matrix_dim; ++column){
            s_row += wMatrix[row][column];
            s_mul += wMatrix[row][column] * yArray[column][1];
        }
        
        eArray[row][1] = s_mul / s_row;
    }
    
}

void read_net_file(char * path, double** wMatrix, bool HAS_HEADER){
    bool hasHeader = HAS_HEADER;
    
    FILE* nFile = fopen(path, "r");
    char line[1024];
    int row; int column; int i=0;
    while(fgets(line,1024,nFile))
    {
        if(hasHeader == false){
            char * line_temp = strdup(line);
            
            //remove new line
            if(line_temp[strlen(line_temp)-1] == '\n')
                line_temp[strlen(line_temp)-1] = '\0';
            
            if(line_temp != NULL){
                char * tokens = strtok(line_temp, ",");
                char * token_array[2];
                i=0;
                while(tokens != NULL){
                    token_array[i++] = tokens;
                    tokens = strtok(NULL, ",");
                }
                
                row = atoi(token_array[0]);
                column = atoi(token_array[1]);
                wMatrix[row][column] += 1;
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

int getDimension(char *pathYFile, bool HAS_HEADER){
    bool hasHeader = HAS_HEADER;
    int line_count = 0;
    FILE * yFile = fopen(pathYFile, "r");
    char line[1024];
    
    while(fgets(line, 1024, yFile)){
        if(hasHeader == false){
            char * line_temp = strdup(line);
            if(line_temp != NULL){
                line_count ++;
            }
            free(line_temp);
        }
        else{
            hasHeader = false;
        }
    }
    
    fclose(yFile);
    return line_count;
}

void read_y_file(char * path, double** yArray, bool HAS_HEADER){
    bool hasHeader = HAS_HEADER;
    int index = -1;
    FILE * yFile = fopen(path, "r");
    char line[1024];
    double id; double y_value;
    while(fgets(line, 1024, yFile)){
        if(hasHeader == false)
        {
            char * line_temp = strdup(line);
            //remove new line
            if(line_temp[strlen(line_temp)-1] == '\n')
                line_temp[strlen(line_temp)-1] = '\0';
            
            if(line_temp != NULL){
                index++;
                char * tokens = strtok(line_temp, ",");
                char * token_array[2];
                int i=0;
                while(tokens != NULL){
                    if(i<2) token_array[i++] = tokens;
                    tokens = strtok(NULL, ",");
                }
                id = atof(token_array[0]);
                y_value = atof(token_array[1]);
                
                yArray[index][0] = id;
                yArray[index][1] = y_value;
                
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

void init_w_matrix(double** wMatrix, int *matrix_dim){
    for(int row=0; row < *matrix_dim; ++row){
        for(int column=0; column<*matrix_dim; ++column){
            wMatrix[row][column] = 0;
        }
    }
}

void print_net_exposure_values(double** y_array, double** e_array, int matrix_dim){
    for(int i=0; i < matrix_dim; ++i){
        printf("Exposure of %d: %f \n", (int)y_array[i][0],e_array[i][1]);
    }
}

void get_net_exposure_values(double** y_array, double** e_array, int matrix_dim){
    
}
