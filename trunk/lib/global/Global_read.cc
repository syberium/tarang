/* Tarang-2
 *
 * Copyright (C) 2008, 2009  Mahendra K. Verma
 *
 * Mahendra K. Verma
 * Indian Institute of Technology, Kanpur-208016
 * UP, India
 *
 * mkv@iitk.ac.in
 *
 * This file is part of Tarang-2 .
 *
 * Tarang-2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * Tarang-2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Tarang-2; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, U
 */

/*! \file  Global.h
 * 
 * @brief  Class constructor of Global.
 *
 * @author  M. K. Verma, A. G. Chatterjee
 * @version 4.0  MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */

#include "Global.h"
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>

	//*********************************************************************************************
template<typename T>
void operator >> (const YAML::Node& node, vector<complex<T> >& vec) {
	string str;
	string str_tmp;
	int n=node.size();  
	node[0] >> str;
	for (int i=1; i<n; i++) {
		node[i] >> str_tmp;
		str+=", "+str_tmp;
	}
	
	stringstream ss;
	complex<T> c;
	
	string Complex_str="";
	bool real=false;
	bool imag=false;
	bool separator=false;
	
	for (int i=0; i< str.length(); i++){
		switch (str[i]){
			case '(':
				if (imag)
					cerr<<"para.yaml: "<<str<<" is not a complex number";
				real=true;
				imag=false;
				break;
			case ',':
				if (real){
					separator=true;
					real=false;
					imag=true;
				}
				else{
					separator=false;
					real=false;
					imag=false;
				}
				break;
			case ')':
				if (imag)
					Complex_str+=")";
				else
					cerr<<"para.yaml: "<<str<<" is not a complex number";
				imag=false;
				real=false;
				separator=false;
				ss<<Complex_str;
				ss>>c;
				vec.push_back(c);
				Complex_str="";
				ss.flush();
				break;
		}
		if ((real || imag) && (str[i] != ' '))
			Complex_str+=str[i];
	}
}

template<typename T>
void operator >> (const YAML::Node& node, Array<T,1 >& bl_arr) {
	vector<T> vec;
	node >> vec;
	bl_arr.resize(vec.size());
	for (int i=0; i<vec.size(); i++)
		bl_arr(i)=vec[i];
}


template<typename T>
void operator >> (const YAML::Node& node, Array<complex<T>,1 >& bl_arr) {
	string str;
	string str_tmp;
	vector<complex<T> > vec;
	int n=node.size();
	node[0] >> str;
	for (int i=1; i<n; i++) {
		node[i] >> str_tmp;
		str+=", "+str_tmp;
	}
	
	stringstream ss;
	complex<T> c;
	
	string Complex_str="";
	bool real=false;
	bool imag=false;
	bool separator=false;
	
	for (int i=0; i< str.length(); i++){
		switch (str[i]){
			case '(':
				if (imag)
					cerr<<"para.yaml: "<<str<<" is not a complex number";
				real=true;
				imag=false;
				break;
			case ',':
				if (real){
					separator=true;
					real=false;
					imag=true;
				}
				else{
					separator=false;
					real=false;
					imag=false;
				}
				break;
			case ')':
				if (imag)
					Complex_str+=")";
				else
					cerr<<"para.yaml: "<<str<<" is not a complex number";
				imag=false;
				real=false;
				separator=false;
				ss<<Complex_str;
				ss>>c;
				vec.push_back(c);
				Complex_str="";
				ss.flush();
				break;
		}
		if ((real || imag) && (str[i] != ' '))
			Complex_str+=str[i];
	}
	
	bl_arr.resize(vec.size());
	//cout<<"No of C numbers = "<<vec.size()<<endl;
	for (int i=0; i< vec.size(); i++)
		{
		bl_arr(i)=vec[i];
		//cout<<i<<": "<<bl_arr(i)<<endl;
		}
}


template < typename T >
void Global::Assign_if_input_provided(const YAML::Node& node, const  string parameter, T &var, T default_value) {
	if (Input_provided(node, parameter))
		node[parameter] >> var;
	else
		var = default_value;
}

//*********************************************************************************************


void Global::Parse(int argc, char** argv, bool parse_para)
{
	mpi.master=(mpi.my_id==mpi.master_id);
	mpi.num_p_rows=1;
	mpi.num_p_cols=mpi.numprocs;

	int opt = 0;
	int longIndex = 0;

	const char *optString = "vhr:";

	const struct option longOpts[] = {
		{ "version", no_argument, NULL, 'v' },
		{ "help", no_argument, NULL, 'h' },
		{ "num_rows", required_argument, NULL, 'r' },
		{ NULL, no_argument, NULL, 0 }
	};

	bool stop=false;
	bool help_message_printed=false;

	stringstream help_message;
	help_message << "Usage: mpirun -np num_procs " << argv[0] << "  [-v,-version] [-h,-help] [-r,-num_rows <num_row_procs>] [/path/to/data]";

	while( (opt = getopt_long_only(argc, argv, optString, longOpts, &longIndex )) != -1 ) {
		switch( opt ) {
			case 'v':
				if (mpi.master)
					cout << "Compiled with:\n" << CONFIGURE_OPTIONS << "\ntarang " << program.version << endl;
				stop=true;
				help_message_printed=true;
				break;
				
			case 'h':
				if (mpi.master)
					cout << help_message.str() << endl;
				stop=true;
				break;
				
			case 'r':
				mpi.num_p_rows=atoi(optarg);
				if (mpi.num_p_rows>0){
					if (mpi.numprocs%mpi.num_p_rows==0){
						mpi.num_p_cols = mpi.numprocs/mpi.num_p_rows;
					
					}
					else{
						if (mpi.master)
							cerr << "mpi.numprocs(="<<mpi.numprocs<<") must be divisible by mpi.num_p_rows(="<<mpi.num_p_rows<<")." << endl;
						stop=true;
					}

				}
				else{
					if (mpi.master)
						cerr << "num_p_rows must be greater than 0." << endl;
					stop=true;
				}
				break;
				
			default:
				if (mpi.master)
					cerr << "Global_Parse getopt_long default" << endl;
				stop=true;
				break;
		}
	}

	if (optind>=argc && parse_para && !help_message_printed) {
		if (mpi.master)
			cout << help_message.str() << endl;
		stop=true;
	}

	if (stop){
		MPI_Finalize();
		exit(0);
	}   
	
	//Parse the parameter file
	if (optind<argc){
		struct stat st;

		ifstream para_yaml;
		string para_path="";

		if (stat(argv[optind], &st) != -1) {
			if (S_ISDIR(st.st_mode)) {
				io.data_dir = argv[optind];
				para_path = io.data_dir+"/in/para.yaml";
			}
			else if (S_ISREG(st.st_mode) || S_ISLNK(st.st_mode) ) {
				para_path = argv[optind];
			}
			else{
				cerr << "Global::Parse: '" + string(argv[optind]) + "' is neither a file nor a folder." << endl;
				exit(1);
			}
		}
		else{
				cerr << "Global::Parse: '" + string(argv[optind]) + "' is neither a file nor a folder." << endl;
				exit(1);
		}
		if (parse_para) {
			para_yaml.open(para_path.c_str());
			if (para_yaml.is_open()){
				try {
					YAML::Parser parser(para_yaml);
					parser.GetNextDocument(para);
				}
				catch(YAML::ParserException& e) {
					cerr << "Global::Parse: Error reading parameter file: \n" << e.what() << endl;
				}

			}
			else {
				if (mpi.master)
					cerr << "Global::Parse: Unable to open '" + para_path + "'." << endl;
				exit(1);
			}
		}
	} 
}

//*********************************************************************************************

void Global::Read()
{
	// program
	para["program"]["kind"] >> program.kind;
	para["program"]["basis_type"] >> program.basis_type;
	para["program"]["iter_or_diag"] >> program.iter_or_diag;
	para["program"]["alias_option"] >> program.alias_option;
	para["program"]["integration_scheme"] >> program.integration_scheme;
	para["program"]["LES_switch"] >> program.LES_switch;
	para["program"]["apply_strong_realitycond_alltime_switch"] >> program.apply_strong_realitycond_alltime_switch;
	para["program"]["apply_weak_realitycond_alltime_switch"] >> program.apply_weak_realitycond_alltime_switch;
	para["program"]["low_dimensional_switch"] >> program.low_dimensional_switch;
	para["program"]["two_and_half_dimension"] >> program.two_and_half_dimension;
	para["program"]["two_dimension"] >> program.two_dimension;
	para["program"]["dt_option"] >> program.dt_option;
	para["program"]["helicity_switch"] >> program.helicity_switch;
	para["program"]["sincostr_switch"] >> program.sincostr_switch;
	if (Input_provided(para, "PHYSICS")){
		Assign_if_input_provided<string>(para["PHYSICS"], "Pr_option", PHYSICS.Pr_option, string("PRLARGE"));
		Assign_if_input_provided<string>(para["PHYSICS"], "Uscaling", PHYSICS.Uscaling, string("USMALL"));
		Assign_if_input_provided<Real>(para["PHYSICS"], "Rayleigh", PHYSICS.Rayleigh, 2000.0);
		Assign_if_input_provided<Real>(para["PHYSICS"], "Prandtl", PHYSICS.Prandtl, 1.0);
		Assign_if_input_provided<int>(para["PHYSICS"], "temperature_grad", PHYSICS.temperature_grad, 1);
		Assign_if_input_provided<Real>(para["PHYSICS"], "Chandrasekhar", PHYSICS.Chandrasekhar, 0.0);
		Assign_if_input_provided<Real>(para["PHYSICS"], "Prandtl_mag", PHYSICS.Prandtl_mag, 1.0);
		Assign_if_input_provided<Real>(para["PHYSICS"], "Prandtl_c", PHYSICS.Prandtl_c, 1.0);
		Assign_if_input_provided<Real>(para["PHYSICS"], "Reynolds", PHYSICS.Reynolds, 1.0);
		Assign_if_input_provided<Real>(para["PHYSICS"], "Reynolds_mag", PHYSICS.Reynolds_mag, 1.0);
		Assign_if_input_provided<Real>(para["PHYSICS"], "Peclet", PHYSICS.Peclet, 1.0);
		Assign_if_input_provided<Real>(para["PHYSICS"], "Peclet_c", PHYSICS.Peclet_c, 1.0);
		
	//  if ( (PHYSICS.temperature_grad != 1) && (PHYSICS.temperature_grad != -1) )
	//      Show_error("PHYSICS.temperature_grad can be either +1 or -1");
	}
	if (program.kind == "MRBC"){
		para["MRBC"]["Pr_option"] >> MRBC.Pr_option;
		para["MRBC"]["Uscaling"] >> MRBC.Uscaling;
		
		para["MRBC"]["Pr"] >> MRBC.Pr;

		if (Input_provided(para["MRBC"] ,"RaD"))
			para["MRBC"]["RaD"] >> MRBC.RaD;
		
		if (Input_provided(para["MRBC"] ,"RaM"))
			para["MRBC"]["RaM"] >> MRBC.RaM;
		
		if (Input_provided(para["MRBC"] ,"SSD"))
			para["MRBC"]["SSD"] >> MRBC.SSD;
		
		if (Input_provided(para["MRBC"] ,"CSA"))
			para["MRBC"]["CSA"] >> MRBC.CSA;
	}
	
	// field
	para["field"]["incompressible"] >> field.incompressible;
	para["field"]["waveno_switch"] >> field.waveno_switch;
	para["field"]["anisotropy_dirn"] >> field.anisotropy_dirn;
	
	field.N[0]=0;
	para["field"]["N"][0] >> field.N[1];
	para["field"]["N"][1] >> field.N[2];
	para["field"]["N"][2] >> field.N[3];

	if ( (Input_provided(para["field"], "kfactor")) && (para["field"]["kfactor"].size()==3) ) {
		field.kfactor.resize(4);
		field.kfactor[0]=0;
		para["field"]["kfactor"][0] >> field.kfactor[1];
		para["field"]["kfactor"][1] >> field.kfactor[2];
		para["field"]["kfactor"][2] >> field.kfactor[3];
	}
	else if ((Input_provided(para["field"], "L")) && (para["field"]["L"].size()==3)){
		field.L.resize(4);
		field.L[0]=0;
		para["field"]["L"][0] >> field.L[1];
		para["field"]["L"][1] >> field.L[2];
		para["field"]["L"][2] >> field.L[3];

	}
	else {
		if (mpi.master) 
			cout << "WARNING: Improper field.kfactor or field.L provided, default values will be taken." << endl;
	}
	
	int no_diss_coeff = para["field"]["diss_coefficients"].size();
	field.diss_coefficients.resize(no_diss_coeff);
	for (int i=0; i<no_diss_coeff; i++) 
		para["field"]["diss_coefficients"][i] >> field.diss_coefficients[i];
	
	field.hyper_diss_coefficients.resize(no_diss_coeff);
	field.hyper_diss_exponents.resize(no_diss_coeff);
	
	if (Input_provided(para["field"],"hyper_diss_coefficients")) 
		for (int i=0; i<no_diss_coeff; i++)
			para["field"]["hyper_diss_coefficients"][i] >> field.hyper_diss_coefficients[i];
	
	if (Input_provided(para["field"],"hyper_diss_exponents")) 
		for (int i=0; i<no_diss_coeff; i++)
			para["field"]["hyper_diss_exponents"][i] >> field.hyper_diss_exponents[i];
	
	// time
	para["time"]["init"] >> time.init;
	para["time"]["final"] >> time.final;
	para["time"]["dt_fixed"] >> time.dt_fixed;
	para["time"]["Courant_no"] >> time.Courant_no;
	Assign_if_input_provided(para["time"], "job_time", time.job_time, string(""));
	time.dt = time.dt_fixed;
	
	// force
	para["force"]["U_switch"] >> force.U_switch;
	Assign_if_input_provided(para["force"], "W_switch", force.W_switch, false);
	Assign_if_input_provided(para["force"], "T_switch", force.T_switch, false);
	Assign_if_input_provided(para["force"], "C_switch", force.C_switch, false);
	
	para["force"]["field_procedure"] >> force.field_procedure;
	
	para["force"]["int_para"] 		>> force.int_para;
	para["force"]["double_para"] 	>> force.double_para;
	para["force"]["string_para"] 	>> force.string_para;
//the following input is used by case No 24 and 25
	if (Input_provided(para["force"],"force_I0")) { 
		para["force"]["force_I0"] 				>> force.force_I0;
		force.force_alpha_is_adaptive = true;
	} else {
		force.force_alpha_is_adaptive = false;
	}
	if (Input_provided(para["force"],"force_I0_helical")) { 
		para["force"]["force_I0_helical"] 		>> force.force_I0_helical;
		force.force_beta_is_adaptive = true;
	} else {
		force.force_beta_is_adaptive = false;
	}
	if (Input_provided(para["force"],"force_I0_crosshelical")) { 
		para["force"]["force_I0_crosshelical"] 	>> force.force_I0_crosshelical;
		force.force_gamma_is_adaptive = true;
	} else {
		force.force_gamma_is_adaptive = false;
	}

	if (Input_provided(para["force"],"mag_force_I0")) { 
		para["force"]["mag_force_I0"] 				>> force.mag_force_I0;
		force.mag_force_alpha_is_adaptive = true;
	} else {
		force.mag_force_alpha_is_adaptive = false;
	}
	if (Input_provided(para["force"],"mag_force_I0_helical")) { 
		para["force"]["mag_force_I0_helical"] 		>> force.mag_force_I0_helical;
		force.mag_force_beta_is_adaptive = true;
	} else {
		force.mag_force_beta_is_adaptive = false;
	}
	if (Input_provided(para["force"],"mag_force_I0_crosshelical")) { 
		para["force"]["mag_force_I0_crosshelical"]	>> force.mag_force_I0_crosshelical;	
		force.mag_force_gamma_is_adaptive = true; 
	} else { 
		force.mag_force_gamma_is_adaptive = false; 
	}
// it is used by forces 24 and 25	
    force.force_alpha = force.double_para(2);
    force.force_beta  = force.double_para(3);
    force.force_gamma = force.double_para(4);

    force.mag_force_alpha = force.double_para(7);
    force.mag_force_beta  = force.double_para(8);
    force.mag_force_gamma = force.double_para(9);

    force.force_T 				= 100.0;
    force.force_T_helical 		= 100.0;
    force.force_T_crosshelical 	= 100.0;
    
    force.mag_force_T 				= 100.0;
    force.mag_force_T_helical 		= 100.0;
    force.mag_force_T_crosshelical 	= 100.0;
// it is used by force 25 & 26
    force.update_time = time.init-1;
    force.delta_t = 0.0005;

    force.total_injected_kin_energy = force.total_injected_kin_helicity = force.total_injected_kin_crosshel = 0;
    force.total_injected_mag_energy = force.total_injected_mag_helicity = force.total_injected_mag_crosshel = 0;


	force.modes.number=para["force"]["modes"].size();
	force.modes.number_components = no_components_table[program.kind];
	force.modes.coords.resize(force.modes.number,4);

	// force coord and modes
	if (program.basis_type == "SSS") {
		Array<int, 1> blitz_int_temp_array;
		Array<Real, 1> blitz_Real_temp_array;
		force.modes.field_array_real.resize(force.modes.number,force.modes.number_components);
		for (int i=0; i<force.modes.number; i++) {
			para["force"]["modes"][i]["coord"]>>blitz_int_temp_array;
			para["force"]["modes"][i]["mode"]>>blitz_Real_temp_array;
			
			force.modes.coords(i, Range::all()) = blitz_int_temp_array;
			force.modes.field_array_real(i, Range::all()) = blitz_Real_temp_array;
		}
	}
	else {
		Array<int, 1> blitz_int_temp_array;
		Array<Complex, 1> blitz_Complex_temp_array;
		force.modes.field_array_complex.resize(force.modes.number,force.modes.number_components);
		for (int i=0; i<force.modes.number; i++) {
			para["force"]["modes"][i]["coord"]>>blitz_int_temp_array;
			para["force"]["modes"][i]["mode"]>>blitz_Complex_temp_array;
			
			force.modes.coords(i, Range::all()) = blitz_int_temp_array;
			force.modes.field_array_complex(i, Range::all()) = blitz_Complex_temp_array;
			}
		}
	
	// IO
	if (Input_provided(para["io"], "data_dir"))
		para["io"]["data_dir"] >> io.data_dir;
	
	para["io"]["input_field_procedure"] >> io.input_field_procedure;
	para["io"]["input_vx_vy_switch"] >> io.input_vx_vy_switch;
	para["io"]["output_vx_vy_switch"] >> io.output_vx_vy_switch;
	
	if (Input_provided(para["io"],"diagnostic_procedures"))
		para["io"]["diagnostic_procedures"] >> io.diagnostic_procedures;
	
	if ( (io.input_field_procedure==2) && (Input_provided(para["io"],"N_in_reduced") ) ) {
		io.N_in_reduced.resize(3);
		para["io"]["N_in_reduced"] >> io.N_in_reduced;
	}
	
	if (Input_provided(para["io"],"N_out_reduced") ) {
		io.N_out_reduced.resize(3);
		para["io"]["N_out_reduced"] >> io.N_out_reduced;
	}

	para["io"]["int_para"] >> io.int_para;
	para["io"]["double_para"] >> io.double_para;
	para["io"]["string_para"] >> io.string_para;


	if (Input_provided(para["io"], "slice_save")) {
		io.slice_save.resize(para["io"]["slice_save"].size());

		string expression_str;
		for (int i=0; i<para["io"]["slice_save"].size(); i++) {
			para["io"]["slice_save"][i] >> expression_str;
			io.slice_save[i] = h5::Expression(expression_str);
		}
	}
	
		// init cond modes
	if (io.input_field_procedure==4) {
		io.init_cond_modes.number=para["io"]["init_cond_modes"].size();
		io.init_cond_modes.number_components = no_components_table[program.kind];
		io.init_cond_modes.coords.resize(io.init_cond_modes.number, 4);
		
		if  (program.basis_type == "SSS") {
			Array<int, 1> blitz_int_temp_array;
			Array<Real, 1> blitz_Real_temp_array;
			
			io.init_cond_modes.field_array_real.resize(io.init_cond_modes.number,io.init_cond_modes.number_components);
			for (int i=0; i<io.init_cond_modes.number; i++) {
				para["io"]["init_cond_modes"][i]["coord"]>>blitz_int_temp_array;
				para["io"]["init_cond_modes"][i]["mode"]>>blitz_Real_temp_array;
				
				io.init_cond_modes.coords(i,0) = 0;
				io.init_cond_modes.coords(i, Range(1,3)) = blitz_int_temp_array;
				io.init_cond_modes.field_array_real(i, Range::all()) = blitz_Real_temp_array;
			}
		}
		else {
			Array<int, 1> blitz_int_temp_array;
			Array<Complex, 1> blitz_Complex_temp_array;
			
			io.init_cond_modes.field_array_complex.resize(io.init_cond_modes.number,io.init_cond_modes.number_components);
			for (int i=0; i<io.init_cond_modes.number; i++)
			{
				para["io"]["init_cond_modes"][i]["coord"]>>blitz_int_temp_array;
				para["io"]["init_cond_modes"][i]["mode"]>>blitz_Complex_temp_array;
				
				io.init_cond_modes.coords(i,0) = 0;
				io.init_cond_modes.coords(i, Range(1,3)) = blitz_int_temp_array;
				io.init_cond_modes.field_array_complex(i, Range::all()) = blitz_Complex_temp_array;
			}
		}
	}

		// Probes
	io.probes.spectral_space.number = para["io"]["probes"]["spectral_space"].size();
	io.probes.spectral_space.coords.resize(io.probes.spectral_space.number,4);
	for (int i=0; i< io.probes.spectral_space.number; i++) {
		io.probes.spectral_space.coords(i,0) = 0;
		para["io"]["probes"]["spectral_space"][i]["coord"][0] >> io.probes.spectral_space.coords(i,1);
		para["io"]["probes"]["spectral_space"][i]["coord"][1] >> io.probes.spectral_space.coords(i,2);
		para["io"]["probes"]["spectral_space"][i]["coord"][2] >> io.probes.spectral_space.coords(i,3);
	}

	io.probes.real_space.number = para["io"]["probes"]["real_space"].size();
	io.probes.real_space.coords.resize(io.probes.real_space.number,4);
	for (int i=0; i< io.probes.real_space.number; i++) {
		io.probes.real_space.coords(i,0) = 0;
		para["io"]["probes"]["real_space"][i]["coord"][0] >> io.probes.real_space.coords(i,1);
		para["io"]["probes"]["real_space"][i]["coord"][1] >> io.probes.real_space.coords(i,2);
		para["io"]["probes"]["real_space"][i]["coord"][2] >> io.probes.real_space.coords(i,3);
	}
	
	Assign_if_input_provided(para["io"]["time"], "global_save_first"        , io.time.global_save_next, myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "complex_field_save_first" , io.time.complex_field_save_next , myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "field_frequent_save_first", io.time.field_frequent_save_next, myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "field_reduced_save_first" , io.time.field_reduced_save_next , myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "real_field_save_first"    , io.time.real_field_save_next, myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "slice_save_first"         , io.time.slice_save_next , myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "field_k_save_first"       , io.time.field_k_save_next       , myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "field_r_save_first"       , io.time.field_r_save_next       , myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "spectrum_save_first"      , io.time.spectrum_save_next      , myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "pressure_save_first"      , io.time.pressure_save_next      , myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "pressure_spectrum_save_first", io.time.pressure_spectrum_save_next, myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "flux_save_first"          , io.time.flux_save_next          , myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "shell_to_shell_save_first", io.time.shell_to_shell_save_next, myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "ring_spectrum_save_first" , io.time.ring_spectrum_save_next , myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "ring_to_ring_save_first"  , io.time.ring_to_ring_save_next  , myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "cylindrical_ring_spectrum_save_first", io.time.cylindrical_ring_spectrum_save_next, myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "cylindrical_ring_to_ring_save_first", io.time.cylindrical_ring_to_ring_save_next, myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "structure_fn_save_first"  , io.time.structure_fn_save_next  , myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "Tk_shell_spectrum_save_first"  , io.time.Tk_shell_spectrum_save_next  , myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "Tk_ring_spectrum_save_first"  , io.time.Tk_ring_spectrum_save_next  , myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "Tk_cylindrical_ring_spectrum_save_first"  , io.time.Tk_cylindrical_ring_spectrum_save_next  , myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "cout_save_first"          , io.time.cout_save_next          , myconstant.INF_TIME);

	//make save_first relative to time.init
	io.time.global_save_next                        +=time.init;
	io.time.complex_field_save_next                 +=time.init;
	io.time.field_frequent_save_next                +=time.init;
	io.time.field_reduced_save_next                 +=time.init;
	io.time.real_field_save_next                    +=time.init;
	io.time.slice_save_next                         +=time.init;
	io.time.field_k_save_next                       +=time.init;
	io.time.field_r_save_next                       +=time.init;
	io.time.spectrum_save_next                      +=time.init;
	io.time.pressure_save_next                      +=time.init;
	io.time.pressure_spectrum_save_next             +=time.init;
	io.time.flux_save_next                          +=time.init;
	io.time.shell_to_shell_save_next                +=time.init;
	io.time.ring_spectrum_save_next                 +=time.init;
	io.time.ring_to_ring_save_next                  +=time.init;
	io.time.cylindrical_ring_spectrum_save_next     +=time.init;
	io.time.cylindrical_ring_to_ring_save_next      +=time.init;
	io.time.structure_fn_save_next                  +=time.init;
	io.time.Tk_shell_spectrum_save_next             +=time.init;
	io.time.Tk_ring_spectrum_save_next              +=time.init;
	io.time.Tk_cylindrical_ring_spectrum_save_next  +=time.init;
	io.time.cout_save_next                          +=time.init;

	
	Assign_if_input_provided(para["io"]["time"], "global_save_interval"        , io.time.global_save_interval        , myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "complex_field_save_interval" , io.time.complex_field_save_interval , myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "field_frequent_save_interval", io.time.field_frequent_save_interval, myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "field_reduced_save_interval" , io.time.field_reduced_save_interval , myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "real_field_save_interval"    , io.time.real_field_save_interval    , myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "slice_save_interval"         , io.time.slice_save_interval         , myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "field_k_save_interval"       , io.time.field_k_save_interval       , myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "field_r_save_interval"       , io.time.field_r_save_interval       , myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "spectrum_save_interval"      , io.time.spectrum_save_interval      , myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "pressure_save_interval"      , io.time.pressure_save_interval      , myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "pressure_spectrum_save_interval", io.time.pressure_spectrum_save_interval, myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "flux_save_interval"          , io.time.flux_save_interval          , myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "shell_to_shell_save_interval", io.time.shell_to_shell_save_interval, myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "ring_spectrum_save_interval" , io.time.ring_spectrum_save_interval , myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "ring_to_ring_save_interval"  , io.time.ring_to_ring_save_interval  , myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "cylindrical_ring_spectrum_save_interval", io.time.cylindrical_ring_spectrum_save_interval, myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "cylindrical_ring_to_ring_save_interval", io.time.cylindrical_ring_to_ring_save_interval, myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "structure_fn_save_interval"       , io.time.structure_fn_save_interval       , myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "Tk_shell_spectrum_save_interval"  , io.time.Tk_shell_spectrum_save_interval  , myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "Tk_ring_spectrum_save_interval"   , io.time.Tk_ring_spectrum_save_interval   , myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "Tk_cylindrical_ring_spectrum_save_interval"  , io.time.Tk_cylindrical_ring_spectrum_save_interval  , myconstant.INF_TIME);
	Assign_if_input_provided(para["io"]["time"], "cout_save_interval"           , io.time.cout_save_interval         , myconstant.INF_TIME);
	
	//last
	Assign_if_input_provided(para["io"]["time"], "global_save_last"        , io.time.global_save_last        , false);
	Assign_if_input_provided(para["io"]["time"], "complex_field_save_last" , io.time.complex_field_save_last , false);
	Assign_if_input_provided(para["io"]["time"], "field_frequent_save_last", io.time.field_frequent_save_last, false);
	Assign_if_input_provided(para["io"]["time"], "field_reduced_save_last" , io.time.field_reduced_save_last , false);
	Assign_if_input_provided(para["io"]["time"], "real_field_save_last"    , io.time.real_field_save_last    , false);
	Assign_if_input_provided(para["io"]["time"], "field_k_save_last"       , io.time.field_k_save_last       , false);
	Assign_if_input_provided(para["io"]["time"], "field_r_save_last"       , io.time.field_r_save_last       , false);
	Assign_if_input_provided(para["io"]["time"], "spectrum_save_last"      , io.time.spectrum_save_last      , false);
	Assign_if_input_provided(para["io"]["time"], "pressure_save_last"      , io.time.pressure_save_last      , false);
	Assign_if_input_provided(para["io"]["time"], "pressure_spectrum_save_last", io.time.pressure_spectrum_save_last, false);
	Assign_if_input_provided(para["io"]["time"], "flux_save_last"          , io.time.flux_save_last          , false);
	Assign_if_input_provided(para["io"]["time"], "shell_to_shell_save_last", io.time.shell_to_shell_save_last, false);
	Assign_if_input_provided(para["io"]["time"], "ring_spectrum_save_last" , io.time.ring_spectrum_save_last , false);
	Assign_if_input_provided(para["io"]["time"], "ring_to_ring_save_last"  , io.time.ring_to_ring_save_last  , false);
	Assign_if_input_provided(para["io"]["time"], "cylindrical_ring_spectrum_save_last", io.time.cylindrical_ring_spectrum_save_last, false);
	Assign_if_input_provided(para["io"]["time"], "cylindrical_ring_to_ring_save_last", io.time.cylindrical_ring_to_ring_save_last, false);
	Assign_if_input_provided(para["io"]["time"], "structure_fn_save_last"  , io.time.structure_fn_save_last  , false);
	Assign_if_input_provided(para["io"]["time"], "Tk_shell_spectrum_save_last"  , io.time.Tk_shell_spectrum_save_last  , false);
	Assign_if_input_provided(para["io"]["time"], "Tk_ring_spectrum_save_last"  , io.time.Tk_ring_spectrum_save_last  , false);
	Assign_if_input_provided(para["io"]["time"], "Tk_cylindrical_ring_spectrum_save_last"  , io.time.Tk_cylindrical_ring_spectrum_save_last  , false);
	Assign_if_input_provided(para["io"]["time"], "cout_save_last"          , io.time.cout_save_last          , false);

	// SPECTRUM
	Assign_if_input_provided(para["spectrum"]["shell"], "turnon", spectrum.shell.turnon, true);
	
	Assign_if_input_provided(para["spectrum"]["ring"], "turnon", spectrum.ring.turnon, false);
	
	if (spectrum.ring.turnon) {
		para["spectrum"]["ring"]["no_sectors"] >> spectrum.ring.no_sectors;
		Assign_if_input_provided<string>(para["spectrum"]["ring"],"sector_option",spectrum.ring.sector_option,"EQUISPACED");
		
		if (spectrum.ring.sector_option == "USER_DEFINED")  {
			if (para["spectrum"]["ring"]["sector_angles"].size() == (spectrum.ring.no_sectors+1)){
				para["spectrum"]["ring"]["sector_angles"] >> spectrum.ring.sector_angles;
				// provide the boundaries (no+1)
			}
			else {
				if (mpi.master) 
					cout << "WARNING: size(spectrum.ring.sector_angles) != number supplied.  Computer proceeds to make EQUISPACED sector_angles " << endl;
				spectrum.ring.sector_option = "EQUISPACED";
			}
		}
	}
	
	Assign_if_input_provided(para["spectrum"]["cylindrical_ring"], "turnon", spectrum.cylindrical_ring.turnon, false);
	
	if (spectrum.cylindrical_ring.turnon) {
		para["spectrum"]["cylindrical_ring"]["no_slabs"] >> spectrum.cylindrical_ring.no_slabs;
		Assign_if_input_provided<string>(para["spectrum"]["cylindrical_ring"], "kpll_option", spectrum.cylindrical_ring.kpll_option, "EQUISPACED");
		
		if (spectrum.cylindrical_ring.kpll_option == "USER_DEFINED") {
			if (para["spectrum"]["cylindrical_ring"]["kpll_array"].size() == (spectrum.cylindrical_ring.no_slabs+1)) {
				para["spectrum"]["cylindrical_ring"]["kpll_array"] >> spectrum.cylindrical_ring.kpll_array;
					// provide the boundaries (no+1)
			}
			else {
				if (mpi.master) 
					cout << "WARNING: size(spectrum.cylindrical_ring.kpll_array) != number supplied.  Computer proceeds to make EQUISPACED kpll_array " << endl;
				spectrum.cylindrical_ring.kpll_option = "EQUISPACED";
			}
		}
	}
	
	
	// EnergyTr
	Assign_if_input_provided(para["energy_transfer"], "turnon", energy_transfer.turnon, false);
	
	if (energy_transfer.turnon) {
		Assign_if_input_provided(para["energy_transfer"], "helicity_flux_switch", energy_transfer.helicity_flux_switch, false);
		
		Assign_if_input_provided(para["energy_transfer"], "helicity_shell_to_shell_switch", energy_transfer.helicity_shell_to_shell_switch, false);
      
        Assign_if_input_provided(para["energy_transfer"], "helicity_ring_to_ring_switch", energy_transfer.helicity_ring_to_ring_switch, false);
    
		
		Assign_if_input_provided(para["energy_transfer"], "Elsasser", energy_transfer.Elsasser, true);
		
		Assign_if_input_provided(para["energy_transfer"], "Vpll_switch", energy_transfer.Vpll_switch, false);
		
			// flux radii assign
		Assign_if_input_provided(para["energy_transfer"]["flux"], "turnon", energy_transfer.flux.turnon, false);
		
		if (energy_transfer.flux.turnon) {
			Assign_if_input_provided(para["energy_transfer"]["flux"], "no_spheres", energy_transfer.flux.no_spheres, 0);
			
			if ( Input_provided(para["energy_transfer"]["flux"],"radii") ) {
				if (para["energy_transfer"]["flux"]["radii"].size() == energy_transfer.flux.no_spheres + 1){
						para["energy_transfer"]["flux"]["radii"] >> energy_transfer.flux.radii;
				}
				else {
					if (mpi.master) 
						cout << "WARNING: size(energy_transfer.flux.radii) != number supplied + 1.  Computer builds the array " << endl;
				}
			}
		}
		
			// SHELL-to-SHELL
		Assign_if_input_provided(para["energy_transfer"]["shell_to_shell"], "turnon", energy_transfer.shell_to_shell.turnon, false);
		
		if (energy_transfer.shell_to_shell.turnon) {
			Assign_if_input_provided(para["energy_transfer"]["shell_to_shell"], "no_shells", energy_transfer.shell_to_shell.no_shells, 0);
			
			if (Input_provided(para["energy_transfer"]["shell_to_shell"],"radii") ) {
				if (para["energy_transfer"]["shell_to_shell"]["radii"].size() == energy_transfer.shell_to_shell.no_shells) {
					para["energy_transfer"]["shell_to_shell"]["radii"]>> energy_transfer.shell_to_shell.radii;
				}
				
				else {
					if (mpi.master) 
						cout << "WARNING: size(energy_transfer.shell_to_shell.radii) != number supplied.  Computer builds the array " << endl;
				}
			}
		}

		
		// ring-to-ring
		Assign_if_input_provided(para["energy_transfer"]["ring_to_ring"], "turnon", energy_transfer.ring_to_ring.turnon, false);
		
		if (energy_transfer.ring_to_ring.turnon) {
				// Read ring radii
			Assign_if_input_provided(para["energy_transfer"]["ring_to_ring"], "no_shells", energy_transfer.ring_to_ring.no_shells, 0);
			
			if (Input_provided(para["energy_transfer"]["ring_to_ring"],"radii") ) {
				if (para["energy_transfer"]["ring_to_ring"]["radii"].size() == energy_transfer.ring_to_ring.no_shells) {
						para["energy_transfer"]["ring_to_ring"]["radii"] >> energy_transfer.ring_to_ring.radii;
				}
				
				else {
					if (mpi.master) 
						cout << "WARNING: size(energy_transfer.ring_to_ring.radii) != number supplied.  Computer builds the array " << endl;
				}
			}
			
			para["energy_transfer"]["ring_to_ring"]["no_sectors"] >> energy_transfer.ring_to_ring.no_sectors;
			
			Assign_if_input_provided<string>(para["energy_transfer"]["ring_to_ring"],"sector_option", energy_transfer.ring_to_ring.sector_option, "EQUISPACED");

			if (energy_transfer.ring_to_ring.sector_option == "USER_DEFINED") {
				if (para["energy_transfer"]["ring_to_ring"]["sector_angles"].size() == (energy_transfer.ring_to_ring.no_sectors+1)) {
					para["energy_transfer"]["ring_to_ring"]["sector_angles"] >> energy_transfer.ring_to_ring.sector_angles;
					// provide the boundaries (no+1)
				}           
				else {
					if (mpi.master) 
						cout << "WARNING: size(energy_transfer.ring_to_ring.sector_angles) != number supplied.  Computer proceeds to make EQUISPACED sector_angles " << endl;
					energy_transfer.ring_to_ring.sector_option = "EQUISPACED";
				}
			}
		}
		
		// Cylindrical ring to ring
		Assign_if_input_provided(para["energy_transfer"]["cylindrical_ring_to_ring"], "turnon", energy_transfer.cylindrical_ring_to_ring.turnon, false);
		
		if (energy_transfer.cylindrical_ring_to_ring.turnon) {
				// Read ring radii
			Assign_if_input_provided(para["energy_transfer"]["cylindrical_ring_to_ring"], "no_shells", energy_transfer.cylindrical_ring_to_ring.no_shells, 0);
			
			if (Input_provided(para["energy_transfer"]["cylindrical_ring_to_ring"],"radii") ) {
				if (para["energy_transfer"]["cylindrical_ring_to_ring"]["radii"].size() == energy_transfer.cylindrical_ring_to_ring.no_shells) {
					para["energy_transfer"]["cylindrical_ring_to_ring"]["radii"] >> energy_transfer.cylindrical_ring_to_ring.radii;
				}           
				else {
					if (mpi.master) 
						cout << "WARNING: size(energy_transfer.cylindrical_ring_to_ring.radii) != number supplied.  Computer builds the array " << endl;
				}
			}
			
			para["energy_transfer"]["cylindrical_ring_to_ring"]["no_slabs"] >> energy_transfer.cylindrical_ring_to_ring.no_slabs;
			
			Assign_if_input_provided<string>(para["energy_transfer"]["cylindrical_ring_to_ring"],"kpll_option",energy_transfer.cylindrical_ring_to_ring.kpll_option, "EQUISPACED");
			
			
			if (energy_transfer.cylindrical_ring_to_ring.kpll_option == "USER_DEFINED") {
				if (para["energy_transfer"]["cylindrical_ring_to_ring"]["kpll_array"].size() == (energy_transfer.cylindrical_ring_to_ring.no_slabs+1)) {
					para["energy_transfer"]["cylindrical_ring_to_ring"]["kpll_array"] >> energy_transfer.cylindrical_ring_to_ring.kpll_array;
				}
				else {
					if (mpi.master) 
						cout << "WARNING: size(energy_transfer.cylindrical_ring_to_ring.kpll_array) != number supplied.  Computer proceeds to make EQUISPACED kpll_array" << endl;
					energy_transfer.cylindrical_ring_to_ring.kpll_option = "EQUISPACED";
				}
			}
		
		}
	}
}

