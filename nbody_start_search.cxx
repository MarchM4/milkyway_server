/*
 * Copyright 2008, 2009 Travis Desell, Dave Przybylo, Nathan Cole,
 * Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail
 * and Rensselaer Polytechnic Institute.
 *
 * This file is part of Milkway@Home.
 *
 * Milkyway@Home is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Milkyway@Home is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
 * */

#include <vector>
#include <string>
#include <fstream>

#include <boost/filesystem.hpp>

#include "nbody_parameters.hxx"

//From BOINC
#include "boinc_db.h"
#include "error_numbers.h"
#include "backend_lib.h"
#include "parse.h"
#include "util.h"

//From TAO
#include "tao/boinc/workunit_information.hxx"
#include "tao/asynchronous_algorithms/particle_swarm_db.hxx"
#include "tao/asynchronous_algorithms/differential_evolution_db.hxx"
#include "tao/undvc_common/arguments.hxx"
#include "tao/undvc_common/file_io.hxx"


#define WORKUNIT_XML "milkyway_nbody_wu_new.xml"
#define ALTERNATE_WORKUNIT_XML "milkyway_nbody_multistage_wu_new.xml"
#define RESULT_XML "milkyway_nbody_result.xml"



using std::vector;
using std::string;

int main(int argc, char **argv) {
    vector<string> arguments(argv, argv + argc);

    //need to initialize the BOINC database
    int retval = config.parse_file();
    if (retval) {
        log_messages.printf(MSG_CRITICAL, "Can't parse config.xml: %s\n", boincerror(retval));
        exit(1);
    }

    retval = boinc_db.open(config.db_name, config.db_host, config.db_user, config.db_passwd);
    if (retval) {
        log_messages.printf(MSG_CRITICAL, "can't open db\n");
        exit(1);
    }

    DB_APP app;
    const char* app_name = "milkyway_nbody";
    char buf[256];
    sprintf(buf, "where name='%s'", app_name);
    if (app.lookup(buf)) {
        log_messages.printf(MSG_CRITICAL, "can't find app %s\n", app_name);
        exit(1);
    }

    if (argument_exists(arguments, "--create_tables")) {
        try {
            WorkunitInformation::create_table(boinc_db.mysql);
            ParticleSwarmDB::create_tables(boinc_db.mysql);
            DifferentialEvolutionDB::create_tables(boinc_db.mysql);
        } catch (string err_msg) {
            cerr << "Error creating the workunit information tables in the database." << endl;
            cerr << "threw message: '" << err_msg << "'" << endl;
        }
        
        exit(0);
    }


    //1. get search id
    //      a. create search from name and input parameters
    double min_simulation_time, max_simulation_time;
    double min_orbit_time, max_orbit_time;
    double min_radius_1, max_radius_1;
    double min_radius_2, max_radius_2;
    double min_mass_1, max_mass_1;
    double min_mass_2, max_mass_2;
    std::vector<double> other_min;
    std::vector<double> other_max;
            NBODY_PARAMETERS *np = (NBODY_PARAMETERS*)malloc(sizeof(NBODY_PARAMETERS));
    string extra_parameters;
    string first_input;
    uint64_t n_bodies;
    uint64_t n_dwarves=0;
    uint64_t multi_stage=0;
    get_argument(arguments, "--n_dwarves", false, n_dwarves);
    get_argument(arguments, "--multi_stage", false, multi_stage);
    if (n_dwarves){
        get_argument(arguments, "--extra_parameters", true, extra_parameters);
        read_nbody_parameters( extra_parameters.c_str(), np );

    }
    else if(multi_stage){
        get_argument(arguments, "--first_input", true, first_input);
        get_argument(arguments, "--extra_parameters", true, extra_parameters);
        read_nbody_parameters( extra_parameters.c_str(), np );

    }
    else{
        get_argument(arguments, "--min_simulation_time", true, min_simulation_time);
        get_argument(arguments, "--max_simulation_time", true, max_simulation_time);
        get_argument(arguments, "--min_orbit_time", true, min_orbit_time);
        get_argument(arguments, "--max_orbit_time", true, max_orbit_time);
        get_argument(arguments, "--min_radius_1", true, min_radius_1);
        get_argument(arguments, "--max_radius_1", true, max_radius_1);
        get_argument(arguments, "--min_radius_2", true, min_radius_2);
        get_argument(arguments, "--max_radius_2", true, max_radius_2);
        get_argument(arguments, "--min_mass_1", true, min_mass_1);
        get_argument(arguments, "--max_mass_1", true, max_mass_1);
        get_argument(arguments, "--min_mass_2", true, min_mass_2);
        get_argument(arguments, "--max_mass_2", true, max_mass_2);
    }
    get_argument(arguments, "--n_bodies", true, n_bodies);
    //2. get app id
    int app_id = app.id;         /*  The NBody simulation application's id is 7. */
    //5. get input filenames (from command line)
    //      a. copy to download directory if needed
    string parameters_filename, histogram_filename;

    get_argument(arguments, "--parameters", true, parameters_filename);
    get_argument(arguments, "--histogram", true, histogram_filename);

    // Copy the input files to the right place in the download directory hierarchy
    copy_file_to_download_dir(parameters_filename);
    copy_file_to_download_dir(histogram_filename);
    vector<string> input_filenames;
    input_filenames.push_back( parameters_filename.substr( parameters_filename.find_last_of('/') + 1) );
    input_filenames.push_back( histogram_filename.substr( histogram_filename.find_last_of('/') + 1) );
    if(multi_stage){
        copy_file_to_download_dir(first_input);
        input_filenames.push_back( first_input.substr( first_input.find_last_of('/') + 1) );
    }

    cout << "input_filenames[0]: " << input_filenames[0] << endl;
    cout << "input_filenames[1]: " << input_filenames[1] << endl;
    if (multi_stage){
        cout << "input_filenames[2]: " << input_filenames[2] << endl;
    }

    //6. get any command line options
    //      -f nbody_parameters.lua -h histogram.txt
    string command_line_options = "-f nbody_parameters.lua -h histogram.txt";
    if(multi_stage){
        command_line_options += " -F -I first_input.lua";
    }


    //7. get any extra xml:
    //      <credit>13297.805963949866</credit>, <rsc_fpops_est>4.925113319981432E15</rsc_fpops_est>, <rsc_fpops_bound>4.925113319981432E19</rsc_fpops_bound>, <rsc_disk_bound>5.24288E7</rsc_disk_bound>
    //      a. calculate rsc_fpops_est
    //      b. calculate rsc_fpops_bound
    //      c. calculate rsc_disk_bound

    vector<double> min_bound;
    vector<double> max_bound;
   /* if (multi_stage){
        min_bound=other_min;
        max_bound=other_max;
    }
    else if (n_dwarves){
        min_bound.push_back(min_simulation_time);
        min_bound.push_back(min_orbit_time);
        max_bound.push_back(max_simulation_time);
        max_bound.push_back(max_orbit_time);
        for (int i=0; i<other_min.size(); ++i){
            min_bound.push_back(other_min[i]);
            max_bound.push_back(other_max[i]);
        }
    }
    */
    if(n_dwarves || multi_stage){
        double *min_b, *max_b;
        int size_min=get_min_parameters(np, &min_b);
        int size_max=get_max_parameters(np, &max_b);
        assert(size_min==size_max);
        vector<double> n_bound(min_b, min_b + size_min);
        vector<double> x_bound(max_b, max_b + size_max);
        min_bound=n_bound;
        max_bound=x_bound;
        free_nbody_parameters(np);
        free(np);

    }
    else{
        double min_b[] = { min_simulation_time, min_orbit_time, min_radius_1, min_radius_2, min_mass_1, min_mass_2};
        double max_b[] = { max_simulation_time, max_orbit_time, max_radius_1, max_radius_2, max_mass_1, max_mass_2};
        vector<double> n_bound(min_b, min_b + 6);
        vector<double> x_bound(max_b, max_b + 6);
        min_bound=n_bound;
        max_bound=x_bound;



    }
    /*

    cout << "min_bounds: ";
    for (int i=0; i<min_bound.size(); ++i){
        cout << min_bound[i] << " ";
    }
    cout << "\nmax_bounds: ";
    for (int j=0; j<max_bound.size();++j){
        cout << max_bound[j] << " ";
    }
    cout << "\n";
     */
    double rsc_disk_bound = 50 * 1024 * 1024; // 50MB

    cout.precision(15);

    ostringstream oss;
    oss << "<rsc_disk_bound>"   << rsc_disk_bound   << "</rsc_disk_bound>" << endl;
    oss << "<n_bodies>"         << n_bodies         << "</n_bodies>" << endl;       //this is needed to calculate the fpops based on the input parameters
    if (multi_stage){
        oss << "<multi_stage>" << 1 << "</multi_stage>" << endl;
    } //this is needed to calculate fpops

    string extra_xml = oss.str();

    //Put the search in the database
    EvolutionaryAlgorithmDB *ea = NULL;
    string search_name;
    try {
        get_argument(arguments, "--search_name", true, search_name);

        if (search_name.substr(0,3).compare("ps_") == 0) {
            ea = new ParticleSwarmDB(boinc_db.mysql, app.id, min_bound, max_bound, arguments);
        } else if (search_name.substr(0,3).compare("de_") == 0) {
            ea = new DifferentialEvolutionDB(boinc_db.mysql, app.id, min_bound, max_bound, arguments);
        } else {
            cerr << "Improperly specified search name: '" << search_name <<"'" << endl;
            cerr << "Search name must begin with either:" << endl;
            cerr << "    'de_'     -       for differential evolution" << endl;
            cerr << "    'ps_'     -       for particle swarm optimization" << endl;
            exit(1);
        }
    } catch (string err_msg) {
            cerr << "Error putting the search into the database." << endl;
            cerr << "threw message: '" << err_msg << "'" << endl;
    }
    if(!multi_stage){

        try{
            WorkunitInformation(boinc_db.mysql,
                                ea->get_name(),
                                app_id,
                                WORKUNIT_XML,
                                RESULT_XML,
                                input_filenames,
                                command_line_options,
                                extra_xml);
            } catch (string err_msg) {
                cerr << "Error putting the workunit information into the database." << endl;
                cerr << "threw message: '" << err_msg << "'" << endl;
            }
    }
    else{
        try{
            WorkunitInformation(boinc_db.mysql,
                                ea->get_name(),
                                app_id,
                                ALTERNATE_WORKUNIT_XML,
                                RESULT_XML,
                                input_filenames,
                                command_line_options,
                                extra_xml);
        } catch (string err_msg) {
            cerr << "Error putting the workunit information into the database." << endl;
            cerr << "threw message: '" << err_msg << "'" << endl;
        }
    }
    return 0;
}
