#ifdef BOINC_APPLICATION
#ifdef _WIN32
#include "boinc_win.h"
#else
#include "config.h"
#endif

#ifndef _WIN32
#include <cstdio>
#include <cctype>
#include <ctime>
#include <cstring>
#include <cstdlib>
#include <csignal>
#include <unistd.h>
#endif

#include "diagnostics.h"
#include "util.h"
#include "filesys.h"
#include "boinc_api.h"
#include "mfile.h"
#endif

#include "nbody_parameters.hxx"

int fread_double_array(FILE *file, const char *array_name, double** array_t) {
    int i, size;
    fscanf(file, array_name);
    fscanf(file, "[%d]: ", &size);
    
    (*array_t) = (double*)malloc(sizeof(double) * size);
    
    for (i = 0; i < size; i++) {
        if (fscanf(file, "%lf", &(*array_t)[i]) != 1) {
            fprintf(stderr,"Error reading into %s\n",array_name);
            exit(-1);
        }
        if (i < size-1) fscanf(file, ", ");
            }
    fscanf(file, "\n");
    return size;
}


void free_nbody_parameters(NBODY_PARAMETERS* np){
    free(np->min_orbit_time);
    free(np->max_orbit_time);
    free(np->min_simulation_time);
    free(np->max_simulation_time);
    free(np->initial_x_min);
    free(np->initial_x_max);
    free(np->initial_y_min);
    free(np->initial_y_max);
    free(np->initial_z_min);
    free(np->initial_z_max);
    free(np->initial_dx_min);
    free(np->initial_dx_max);
    free(np->initial_dy_min);
    free(np->initial_dy_max);
    free(np->initial_dz_min);
    free(np->initial_dz_max);
    free(np->radius_1_min);
    free(np->radius_1_max);
    free(np->radius_2_min);
    free(np->radius_2_max);
    free(np->mass_1_min);
    free(np->mass_1_max);
    free(np->mass_2_min);
    free(np->mass_2_max);
}

int read_nbody_parameters(const char* filename, NBODY_PARAMETERS *np) {
#ifdef BOINC_APPLICATION
    char input_path[512];
    int retval = boinc_resolve_filename(filename, input_path, sizeof(input_path));
    
    if (retval) {
        fprintf(stderr, "APP: error resolving parameters file [%s], %d\n", filename, retval);
        return retval;
    }
    
    FILE* data_file = boinc_fopen(input_path, "r");
#else
    FILE* data_file = fopen(filename, "r");
#endif
    if (!data_file) {
        fprintf(stderr, "Couldn't find input file [%s] to read nbody parameters.\n", filename);
        return 1;
    }
    
    fread_nbody_parameters(data_file, np);
    if (np->parameters_version < 0) {
        fprintf(stderr, "Input file [%s] did not specify parameter file version.\n", filename);
        return 1;
    }
    fclose(data_file);
    return 0;
}

void fread_nbody_parameters(FILE* file, NBODY_PARAMETERS *np) {
    int i, retval;
    
    //    retval = fscanf(file, "parameters_version: %lf\n", &ap->parameters_version);
    if (retval < 1) {
        np->parameters_version = 0.01;
        //      fprintf(stderr, "Error reading nbody parameters file. Parameters version not specified\n");
        //      return;
    }
    
    fscanf(file, "number_nbodies: %d\n", &np->number_nbodies);
    fscanf(file, "multi_stage: %d\n", &np->multi_stage);
    int number_orbit_parameters;
    if (np->multi_stage){
        number_orbit_parameters=np->number_nbodies;
    }
    else{
        number_orbit_parameters=1;
    }
    /*
    fread_double_array(file, "background_parameters", &np->background_parameters);
    fread_double_array(file, "background_step", &np->background_step);
    fread_double_array(file, "background_min", &np->background_min);
    fread_double_array(file, "background_max", &np->background_max);
    fread_int_array(file, "optimize_parameter", &np->background_optimize);
     */
    np->min_orbit_time              = (double*)malloc(sizeof(double) * number_orbit_parameters);
    np->max_orbit_time              = (double*)malloc(sizeof(double) * number_orbit_parameters);
    np->min_simulation_time         = (double*)malloc(sizeof(double) * number_orbit_parameters);
    np->max_simulation_time         = (double*)malloc(sizeof(double) * number_orbit_parameters);
    
    np->initial_x_min              = (double*)malloc(sizeof(double) * np->number_nbodies);
    np->initial_x_max              = (double*)malloc(sizeof(double) * np->number_nbodies);
    np->initial_y_min              = (double*)malloc(sizeof(double) * np->number_nbodies);
    np->initial_y_max              = (double*)malloc(sizeof(double) * np->number_nbodies);
    np->initial_z_min              = (double*)malloc(sizeof(double) * np->number_nbodies);
    np->initial_z_max              = (double*)malloc(sizeof(double) * np->number_nbodies);
    np->initial_dx_min              = (double*)malloc(sizeof(double) * np->number_nbodies);
    np->initial_dx_max              = (double*)malloc(sizeof(double) * np->number_nbodies);
    np->initial_dy_min              = (double*)malloc(sizeof(double) * np->number_nbodies);
    np->initial_dy_max              = (double*)malloc(sizeof(double) * np->number_nbodies);
    np->initial_dz_min              = (double*)malloc(sizeof(double) * np->number_nbodies);
    np->initial_dz_max              = (double*)malloc(sizeof(double) * np->number_nbodies);
    np->radius_1_min               = (double*)malloc(sizeof(double) * np->number_nbodies);
    np->radius_1_max               = (double*)malloc(sizeof(double) * np->number_nbodies);
    np->radius_2_min          =   (double*)malloc(sizeof(double)   * np->number_nbodies);
    
    np->radius_2_max              = (double*)malloc(sizeof(double) * np->number_nbodies);
    np->mass_1_min                 = (double*)malloc(sizeof(double) * np->number_nbodies);
    np->mass_1_max                  = (double*)malloc(sizeof(double) * np->number_nbodies);
    np->mass_2_min                  = (double*)malloc(sizeof(double) * np->number_nbodies);
    np->mass_2_max             =   (double*)malloc(sizeof(int)   * np->number_nbodies);
    
    fread_double_array(file, "min_orbit_time", &np->min_orbit_time);
    fread_double_array(file, "max_orbit_time", &np->max_orbit_time);
    fread_double_array(file, "min_simulation_time", &np->min_simulation_time);
    fread_double_array(file, "max_simulation_time", &np->max_simulation_time);
    fread_double_array(file, "min_x", &np->initial_x_min);
    fread_double_array(file, "max_x", &np->initial_x_max);
    fread_double_array(file, "min_y", &np->initial_y_min);
    fread_double_array(file, "max_y", &np->initial_y_max);
    fread_double_array(file, "min_z", &np->initial_z_min);
    fread_double_array(file, "max_z", &np->initial_z_max);
    fread_double_array(file, "min_dx", &np->initial_dx_min);
    fread_double_array(file, "max_dx", &np->initial_dx_max);
    fread_double_array(file, "min_dy", &np->initial_dy_min);
    fread_double_array(file, "max_dy", &np->initial_dy_max);
    fread_double_array(file, "min_dz", &np->initial_dz_min);
    fread_double_array(file, "max_dz", &np->initial_dz_max);
    fread_double_array(file, "min_radius_1", &np->radius_1_min);
    fread_double_array(file, "max_radius_1", &np->radius_1_max);
    fread_double_array(file, "min_radius_2", &np->radius_2_min);
    fread_double_array(file, "max_radius_2", &np->radius_2_max);
    fread_double_array(file, "min_mass_1", &np->mass_1_min);
    fread_double_array(file, "max_mass_1", &np->mass_1_max);
    fread_double_array(file, "min_mass_2", &np->mass_2_min);
    fread_double_array(file, "max_mass_2", &np->mass_2_max);
    
    /*
    fscanf(file, "convolve: %d\n", &ap->convolve);
    fscanf(file, "sgr_coordinates: %d\n", &ap->sgr_coordinates);
     */
    //    if (ap->parameters_version > 0.01) {
    //       fscanf(file, "aux_bg_profile: %d\n", &ap->aux_bg_profile); //vickej2_bg
    //       } else {
    //ap->aux_bg_profile = 0;
    //    }
    //fscanf(file, "wedge: %d\n", &ap->wedge);
    /*
    ap->integral = (INTEGRAL**)malloc(sizeof(INTEGRAL*));
    ap->integral[0] = (INTEGRAL*)malloc(sizeof(INTEGRAL));
    
    fscanf(file, "r[min,max,steps]: %lf, %lf, %d\n", &ap->integral[0]->r_min, &ap->integral[0]->r_max, &ap->integral[0]->r_steps);
    fscanf(file, "mu[min,max,steps]: %lf, %lf, %d\n", &ap->integral[0]->mu_min, &ap->integral[0]->mu_max, &ap->integral[0]->mu_steps);
    fscanf(file, "nu[min,max,steps]: %lf, %lf, %d\n", &ap->integral[0]->nu_min, &ap->integral[0]->nu_max, &ap->integral[0]->nu_steps);
    ap->integral[0]->r_step_size = (ap->integral[0]->r_max - ap->integral[0]->r_min)/(double)ap->integral[0]->r_steps;
    ap->integral[0]->mu_step_size = (ap->integral[0]->mu_max - ap->integral[0]->mu_min)/(double)ap->integral[0]->mu_steps;
    ap->integral[0]->nu_step_size = (ap->integral[0]->nu_max - ap->integral[0]->nu_min)/(double)ap->integral[0]->nu_steps;
    ap->integral[0]->min_calculation = 0;
    ap->integral[0]->max_calculation = ap->integral[0]->r_steps * ap->integral[0]->mu_steps * ap->integral[0]->nu_steps;
    
    fscanf(file, "number_cuts: %d\n", &ap->number_integrals);
    ap->number_integrals++;
    if (ap->number_integrals > 1) {
        ap->integral = (INTEGRAL**)realloc(ap->integral, sizeof(INTEGRAL*) * ap->number_integrals);
        for (i = 1; i < ap->number_integrals; i++) {
            int temp;
            ap->integral[i] = (INTEGRAL*)malloc(sizeof(INTEGRAL));
            fscanf(file, "r_cut[min,max,steps][%d]: %lf, %lf, %d\n", &temp, &ap->integral[i]->r_min, &ap->integral[i]->r_max, &ap->integral[i]->r_steps);
            fscanf(file, "mu_cut[min,max,steps][%d]: %lf, %lf, %d\n", &temp, &ap->integral[i]->mu_min, &ap->integral[i]->mu_max, &ap->integral[i]->mu_steps);
            fscanf(file, "nu_cut[min,max,steps][%d]: %lf, %lf, %d\n", &temp, &ap->integral[i]->nu_min, &ap->integral[i]->nu_max, &ap->integral[i]->nu_steps);
            ap->integral[i]->r_step_size = (ap->integral[i]->r_max - ap->integral[i]->r_min)/(double)ap->integral[i]->r_steps;
            ap->integral[i]->mu_step_size = (ap->integral[i]->mu_max - ap->integral[i]->mu_min)/(double)ap->integral[i]->mu_steps;
            ap->integral[i]->nu_step_size = (ap->integral[i]->nu_max - ap->integral[i]->nu_min)/(double)ap->integral[i]->nu_steps;
            ap->integral[i]->min_calculation = 0;
            ap->integral[i]->max_calculation = ap->integral[i]->r_steps * ap->integral[i]->mu_steps * ap->integral[i]->nu_steps;
        }
    }*/
}

int get_min_parameters(NBODY_PARAMETERS *np, double ** result){
    int size=np->number_nbodies*10;
    size += 2;
    if (np->multi_stage){
        size += 2;
    }
    (*result) = (double*)malloc(sizeof(double) * size);
    if (np->multi_stage){
        for (int i=0; i<np->number_nbodies; ++i){
            (*result)[12*i]=np->min_orbit_time[i];
            (*result)[12*i+1]=np->min_simulation_time[i];
            (*result)[12*i+2]=np->initial_x_min[i];
            (*result)[12*i+3]=np->initial_y_min[i];
            (*result)[12*i+4]=np->initial_z_min[i];
            (*result)[12*i+5]=np->initial_dx_min[i];
            (*result)[12*i+6]=np->initial_dy_min[i];
            (*result)[12*i+7]=np->initial_dz_min[i];
            (*result)[12*i+8]=np->radius_1_min[i];
            (*result)[12*i+9]=np->radius_2_min[i];
            (*result)[12*i+10]=np->mass_1_min[i];
            (*result)[12*i+11]=np->mass_2_min[i];
        }
    }
    else{
        (*result)[0]=np->min_orbit_time[0];
        (*result)[1]=np->min_simulation_time[0];
        for (int i=0; i<np->number_nbodies; ++i){
            (*result)[10*i+2]=np->initial_x_min[i];
            (*result)[10*i+3]=np->initial_y_min[i];
            (*result)[10*i+4]=np->initial_z_min[i];
            (*result)[10*i+5]=np->initial_dx_min[i];
            (*result)[10*i+6]=np->initial_dy_min[i];
            (*result)[10*i+7]=np->initial_dz_min[i];
            (*result)[10*i+8]=np->radius_1_min[i];
            (*result)[10*i+9]=np->radius_2_min[i];
            (*result)[10*i+10]=np->mass_1_min[i];
            (*result)[10*i+11]=np->mass_2_min[i];

        }
    }
    return size;
}

int get_max_parameters(NBODY_PARAMETERS *np, double **result){
    int size=np->number_nbodies*10;
    size += 2;
    if (np->multi_stage){
        size += 2;
    }
    (*result) = (double*)malloc(sizeof(double) * size);
    if (np->multi_stage){
        for (int i=0; i<np->number_nbodies; ++i){
            (*result)[12*i]=np->max_orbit_time[i];
            (*result)[12*i+1]=np->max_simulation_time[i];
            (*result)[12*i+2]=np->initial_x_max[i];
            (*result)[12*i+3]=np->initial_y_max[i];
            (*result)[12*i+4]=np->initial_z_max[i];
            (*result)[12*i+5]=np->initial_dx_max[i];
            (*result)[12*i+6]=np->initial_dy_max[i];
            (*result)[12*i+7]=np->initial_dz_max[i];
            (*result)[12*i+8]=np->radius_1_max[i];
            (*result)[12*i+9]=np->radius_2_max[i];
            (*result)[12*i+10]=np->mass_1_max[i];
            (*result)[12*i+11]=np->mass_2_max[i];
        }
    }
    else{
        (*result)[0]=np->max_orbit_time[0];
        (*result)[1]=np->max_simulation_time[0];
        for (int i=0; i<np->number_nbodies; ++i){
            (*result)[10*i+2]=np->initial_x_max[i];
            (*result)[10*i+3]=np->initial_y_max[i];
            (*result)[10*i+4]=np->initial_z_max[i];
            (*result)[10*i+5]=np->initial_dx_max[i];
            (*result)[10*i+6]=np->initial_dy_max[i];
            (*result)[10*i+7]=np->initial_dz_max[i];
            (*result)[10*i+8]=np->radius_1_max[i];
            (*result)[10*i+9]=np->radius_2_max[i];
            (*result)[10*i+10]=np->mass_1_max[i];
            (*result)[10*i+11]=np->mass_2_max[i];
        }
    }
    return size;

}
