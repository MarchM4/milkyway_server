#ifndef NBODY_PARAMETERS_H
#define NBODY_PARAMETERS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct nbody_parameters {
	double parameters_version;
	int number_nbody_parameters;
	int number_nbodies;
    int number_extra;
    
    double* min_orbit_time;
    double* max_orbit_time;
    double min_simulation_time;
    double max_simulation_time;
    
	double* initial_x_min;
    double* initial_x_max;
	double* initial_y_min;
    double* initial_y_max;
	double* initial_z_min;
    double* initial_z_max;
    
	double* initial_dx_min;
    double* initial_dx_max;
	double* initial_dy_min;
    double* initial_dy_max;
	double* initial_dz_min;
    double* initial_dz_max;
    
	double* radius_1_min;
    double* radius_1_max;
	double* radius_2_min;
    double* radius_2_max;
    
	double* mass_1_min;
    double* mass_1_max;
	double* mass_2_min;
    double* mass_2_max;
    
    double* other_min;
    double* other_max;
} NBODY_PARAMETERS;

int read_nbody_parameters(const char* file, NBODY_PARAMETERS *np);
void fread_nbody_parameters(FILE* file, NBODY_PARAMETERS *np);
void free_nbody_parameters(NBODY_PARAMETERS *np);
int get_min_parameters(NBODY_PARAMETERS *np, double ** result);
int get_max_parameters(NBODY_PARAMETERS *np, double ** result);
#endif
