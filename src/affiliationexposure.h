//
//  affiliationexposure.h
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

#ifndef affiliationexposure_h
#define affiliationexposure_h

#include <stdio.h>
#include <stdbool.h>

struct ModeInfo{
    int first_mode_indices;
    int second_mode_indices;
    char *first_mode_flag;
};

void read_y_mode_file(char* path, double** yArray, int* first_mode_array, int* second_mode_array, bool HAS_HEADER, char *first_mode_flag);
void init_a_matrix(double** wMatrix, struct ModeInfo *mode_data);
void setup_mode_data(char* yPathFile, bool HAS_HEADER, struct ModeInfo *mode_data);
int findCountPair(int *source, int *target, char* netPathFile, bool HAS_HEADER);
double** generateAMatrix(char* netPathFile, int* f_mode_array, int* s_mode_array, bool HAS_HEADER, struct ModeInfo *mode_data);
double** generateAPrimeMatrix(double** aMatrix, struct ModeInfo *mode_data);
double** transposeAMatrix(double** aMatrix, double** aPrimeMatrix, struct ModeInfo *mode_data);
void calc_affiliation_exposure(double** cMatrix, double** eArray, double** yArray, struct ModeInfo *mode_data);
void print_aff_exposure_values(double** y_array, double** e_array, struct ModeInfo *mode_data);



#endif /* affiliationexposure_h */
