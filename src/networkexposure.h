//
//  networkexposure.h
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

#ifndef networkexposure_h
#define networkexposure_h
#include <stdbool.h>

void calc_network_exposure(double** wMatrix, double ** yArray, double ** eArray, int matrix_dim);
void read_net_file(char * path, double** wMatrix, bool HAS_HEADER);
int getDimension(char *pathYFile, bool HAS_HEADER);
void read_y_file(char * path, double** yArray, bool HAS_HEADER);
void init_w_matrix(double** wMatrix, int *matrix_dim);
void print_net_exposure_values(double** y_array, double** e_array, int matrix_dim);

#endif /* networkexposure_h */
