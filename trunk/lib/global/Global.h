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
			  

#ifndef _H_Global
#define _H_Global

#include "def_vars.h"

#include <yaml-cpp/yaml.h>
#include <vector>
#include <map>

//*********************************************************************************************

class Global
{
public:

	YAML::Node para;            //stores all the para in the para.yaml file

	map<string, map<string, map<string, string> > > basis_table;
	map<string, int> no_components_table;

	map<string, int> global_data_packet_size_table;
	map<string, map<string, int> > spectral_probe_packet_size_table;
	
	map<string, int> real_probe_packet_size_table;
	

	struct program {
		const string version;
		
		string kind;
		string iter_or_diag;
		string alias_option;
		string integration_scheme;
		string basis_type;
		bool LES_switch;
		bool T_exists;
		bool W_exists;
		bool apply_strong_realitycond_alltime_switch;
		bool apply_weak_realitycond_alltime_switch;
		bool low_dimensional_switch;
		bool two_and_half_dimension;
		bool two_dimension;
		bool helicity_switch;
		bool adaptive_viscosity;

		int dt_option;  // 0=fixed; 1=dx/umax; 
		string sincostr_switch;             // SINX, COSX, SCC, CSC, CCS, SSC, CSS, SCS
		string sincostr_switch_Vx, sincostr_switch_Vy, sincostr_switch_Vz;
		string sincostr_switch_F;
		string sincostr_switch_Visqr;
		string sincostr_switch_VxVy, sincostr_switch_VxVz, sincostr_switch_VyVz;
		string sincostr_switch_FVx, sincostr_switch_FVy, sincostr_switch_FVz;
		string sincostr_switch_divergence;	
		
		

		program();
		void Print(int my_level);
	} program;

	struct PHYSICS {
		string	Pr_option;					// Prandtl number switch (PRLARGNE..) for RB
		string	Uscaling;					// UBscaling (ULARGE... ) for RB
		Real Rayleigh;							// Rayleigh number
		Real Prandtl;							// Prandtl number	
		int temperature_grad;			// +1 for convection; -x for stratification;
		
		Real Chandrasekhar; // Chandrashekhar number
		Real Prandtl_mag; //Magnetic Prandlt number
		Real Prandtl_c;
		Real Reynolds;
		Real Reynolds_mag;
		Real Peclet;
		Real Peclet_c;
		// factor for u3 in temperature eqn
		void Print(int my_level);
	} PHYSICS;

	struct MRBC {
		string	Pr_option;					// Prandtl number switch (PRLARGNE..) for RB
		string	Uscaling;					// UBscaling (ULARGE... ) for RB
		Real		Pr;							// Prandtl number
		
		Real		RaD;
		Real		RaM;
		Real		SSD;						// Surface saturation deficit
		Real		CSA;						// Condensate in saturated ascent
		void Print(int my_level);
	} MRBC;
	
	struct field {
		//string SPACE;	// REAL or SPECTRAL
		bool incompressible;
		bool waveno_switch;
		int anisotropy_dirn;  //1,2,3
		
		int N[4];
		int Nx, Ny, Nz;
		int howmany;  // how many fields, e.g. RBC has 4.

		vector<Real> kfactor;
		Real xfactor[4];
		vector<Real> L;		// Box size
		Real Delta_x[4];	// Delta x along each direction
		vector<Real> diss_coefficients;
		vector<Real> hyper_diss_coefficients;
		vector<int> hyper_diss_exponents;
		
		TinyVector<int,3> shape_complex_array;
		TinyVector<int,3> shape_real_array;
		
		Real twobyL1;			// for Chebyshev basis
		
		int maxlx, maxly, maxlz;

		//
		int lx_start;
		int ly_start;
		int lz_start;

		int maxrx;
		int maxry;
		int maxrz;

		int rx_start;
		int ry_start;
		int rz_start;
		//
		
		void Print(int my_level);
	} field;

	struct mpi {
		int my_id;
		int numprocs;
		int master_id;
		bool master;
		
		// for pencil

		int num_p_rows;             // processors along the vertical direction
		int num_p_cols;             // processors along the horizontal direction  
		
		int num_x_procs;
		int num_y_procs;
		int num_z_procs;

		int num_x_procs_real;
		int num_y_procs_real;
		int num_z_procs_real;

		int my_x_pcoord;
		int my_y_pcoord;
		int my_z_pcoord;

		int my_x_pcoord_real;
		int my_y_pcoord_real;
		int my_z_pcoord_real;
		
		MPI_Comm MPI_COMM_ROW;
		MPI_Comm MPI_COMM_COL;
		
		mpi();
		void Print(int my_level);
	} mpi;


	struct time{
		Real init; 
		Real final;
		Real dt_fixed;
		Real Courant_no;
		string job_time;
		
		Real dt;				// variable time including CFL
		Real now;
		Real previous;
		Real keplerian;
		bool dt_computation_done;
		
		clock_t job_time_final;
		void Print(int my_level);
	} time;

	struct force{
		bool U_switch;
		bool W_switch;
		bool T_switch;
		bool C_switch;
		bool configuration_done;

		int field_procedure;

		Array<int,1> int_para;
		Array<Real,1> double_para;
		Array<string,1> string_para;

		bool force_alpha_is_adaptive;
		Real force_alpha;
		Real force_injection;
		Real force_T;
		Real force_I0;

		bool force_beta_is_adaptive;
		Real force_beta;
		Real force_injection_helical;
		Real force_T_helical;
		Real force_I0_helical;

		bool force_gamma_is_adaptive;
		Real force_gamma;
		Real force_injection_crosshelical;
		Real force_T_crosshelical;
		Real force_I0_crosshelical;

		bool mag_force_alpha_is_adaptive;
		Real mag_force_alpha;
		Real mag_force_injection;
		Real mag_force_T;
		Real mag_force_I0;

		bool mag_force_beta_is_adaptive;
		Real mag_force_beta;
		Real mag_force_injection_helical;
		Real mag_force_T_helical;
		Real mag_force_I0_helical;

		bool mag_force_gamma_is_adaptive;
		Real mag_force_gamma;
		Real mag_force_injection_crosshelical;
		Real mag_force_T_crosshelical;
		Real mag_force_I0_crosshelical;

		Real A1, A2, A3;
		Real B1, B2, B3;
		Real C1, C2, C3;
//For RANDOM FORCE by Titov V
		Real update_time, delta_t;
		bool force_U_lock,force_B_lock;
		bool empty_force;
		TinyVector<Real, 3> Ku, Kb;
//end RANDOM FORCE by Titov V
		Real total_injected_kin_energy, total_injected_kin_helicity, total_injected_kin_crosshel;
		Real total_injected_mag_energy, total_injected_mag_helicity, total_injected_mag_crosshel;

		struct modes{
			bool read_done;
			int number;
			int number_components;  // 2 for fluid, 4 for U,W (incompressiblity)
			Array<int,2> coords;
			Array<Real,2> field_array_real;
			Array<Complex,2> field_array_complex;
			void Print(int my_level);
		} modes;
		void Print(int my_level);
	} force;

	struct io {
		string data_dir;
		int input_field_procedure;
		bool input_vx_vy_switch;
		bool output_vx_vy_switch;
		const int output_precision;  // 6 for float and 12 for double in ascii format

		bool output_real_field_done;
		bool output_field_k_done;
		bool output_pressure_done;
		bool output_pressure_spectrum_done;
		bool real_space_field_available;
		
		bool output_nlin_magnitude_done;
		
		vector<int>	diagnostic_procedures;

		vector<int> N_in_reduced;
		vector<int> N_out_reduced;

		Array<int,1> int_para;
		Array<Real,1> double_para;
		Array<string,1> string_para;

		vector<h5::Expression> slice_save;
		
		struct init_cond_modes {
			int number;
			int number_components;		// if incompressible: subtract 1 for V, W
			Array<int,2> coords;
			Array<Real,2> field_array_real;
			Array<Complex,2> field_array_complex;
			void Print(int my_level);
			int In_table(int kx, int ky, int kz);
		} init_cond_modes;
		
		struct global_data {
			unsigned long buffer_size;
			unsigned long buffer_index;		// initialize to zero
			unsigned long packet_size;

			Array<Real,1> buffer;
		} global_data;

		struct probes{
			struct spectral_space{
				int number;
				unsigned long buffer_size;
				Array<int,2> coords;
				unsigned long buffer_index;		// initialize to zero
				int packet_size;		// for each time frame (define)
				// (global.io.probes.spectral_space.number * 6) +1 + Tk for each field; 

				Array<Real,1> field_buffer;  // in all nodes
				void Print(int my_level);
			} spectral_space;

			struct real_space {
				int number;
				unsigned long buffer_size;
				Array<int,2> coords;
				unsigned long buffer_index;		// initialize to zero
				int packet_size;

				Array<Real,1> field_buffer;
				void Print(int my_level);
			} real_space;

			void Print(int my_level);
		} probes;

		struct time{
			Real global_save_next;
			Real complex_field_save_next;
			Real field_frequent_save_next;
			Real field_reduced_save_next;
			Real real_field_save_next;
			Real slice_save_next;
			Real field_k_save_next;
			Real field_r_save_next;
			Real spectrum_save_next;
			Real pressure_save_next;
			Real pressure_spectrum_save_next;
			Real flux_save_next;
			Real shell_to_shell_save_next;
			Real ring_spectrum_save_next;
			Real ring_to_ring_save_next;
			Real cylindrical_ring_spectrum_save_next;
			Real cylindrical_ring_to_ring_save_next;
			Real structure_fn_save_next;
			Real Tk_shell_spectrum_save_next;
			Real Tk_ring_spectrum_save_next;
			Real Tk_cylindrical_ring_spectrum_save_next;
			Real cout_save_next;


			Real global_save_interval;	
			Real complex_field_save_interval; 
			Real field_frequent_save_interval;
			Real field_reduced_save_interval;
			Real real_field_save_interval;
			Real slice_save_interval; 
			Real field_k_save_interval;				
			Real field_r_save_interval;	
			Real pressure_save_interval;
			Real spectrum_save_interval;
			Real pressure_spectrum_save_interval;
			Real flux_save_interval;
			Real shell_to_shell_save_interval;
			Real ring_spectrum_save_interval;
			Real ring_to_ring_save_interval;
			Real cylindrical_ring_spectrum_save_interval;
			Real cylindrical_ring_to_ring_save_interval;
			Real structure_fn_save_interval;
			Real Tk_shell_spectrum_save_interval;
			Real Tk_ring_spectrum_save_interval;
			Real Tk_cylindrical_ring_spectrum_save_interval;
			Real cout_save_interval;
			
			bool global_save_last;	
			bool complex_field_save_last; 
			bool field_frequent_save_last;
			bool field_reduced_save_last;
			bool real_field_save_last;
			bool field_k_save_last;				
			bool field_r_save_last;	
			bool pressure_save_last;
			bool spectrum_save_last;
			bool pressure_spectrum_save_last;
			bool flux_save_last;
			bool shell_to_shell_save_last;
			bool ring_spectrum_save_last;
			bool ring_to_ring_save_last;
			bool cylindrical_ring_spectrum_save_last;
			bool cylindrical_ring_to_ring_save_last;
			bool structure_fn_save_last;
			bool Tk_shell_spectrum_save_last;
			bool Tk_ring_spectrum_save_last;
			bool Tk_cylindrical_ring_spectrum_save_last;
			bool cout_save_last;
			
			void Print(int my_level);
		} time;

		io();
		void Print(int my_level);
		void Dump_buffer(ofstream& file, Array<Real,1> buffer, unsigned long& buffer_index, unsigned long packet_size);
	} io;


	struct  spectrum {
		struct shell {
			bool turnon;
			int no_shells;
			void Print(int my_level);
		} shell;

		struct ring {
			bool turnon;
			int no_shells;     // prog defined
			int no_sectors;   // to be read
			string sector_option;
			Array<Real,1> sector_angles;
			void Print(int my_level);
		} ring;

		struct cylindrical_ring {
			bool turnon;
			int no_shells;		// prog define
			int no_slabs;		 // to be read
			string kpll_option;
			Array<Real,1> kpll_array;  // z coords for the sections
			void Print(int my_level);
		} cylindrical_ring;

		void Print(int my_level);
	} spectrum;

	struct energy_transfer {
		bool turnon;
		bool helicity_flux_switch;
		bool helicity_shell_to_shell_switch;
        bool helicity_ring_to_ring_switch;
		bool Elsasser;
		bool Vpll_switch;  // for energy transfer of parallel diagnostics

		struct flux {
			bool turnon;
			int no_spheres;
			Array<Real,1> radii;
			void Print(int my_level);
		} flux;

		struct shell_to_shell {
			bool turnon;
			int no_shells;
			Array<Real,1> radii;
			void Print(int my_level);
		} shell_to_shell;

		struct ring_to_ring {
			bool turnon;
			int no_shells;
			int no_sectors;
			string sector_option;
			Array<Real,1> radii;
			Array<Real,1> sector_angles;
			void Print(int my_level);
		} ring_to_ring;

		struct cylindrical_ring_to_ring {
			bool turnon;
			int no_shells;
			int no_slabs;
			Array<Real,1> radii;
			string kpll_option;
			Array<Real,1> kpll_array;  // kpll coords for the sections
			void Print(int my_level);
		} cylindrical_ring_to_ring;
		void Print(int my_level);
	} energy_transfer;

	struct structure_fn {
		bool turnon;
		bool box_switch;
		bool planar_switch;
		bool approx_switch;

		int qmin;
		int qmax;
		int array_size;
		int rmax_div_st_fn_r_max;
		int nearest_neighbour_index;  // points_to_skip+1
		int dist_farthest_proc;			// for i1 index...
	} structure_fn;


	struct temp_array {
		Array<Complex,3> X_transform;
		
		Array<Complex,3> X;	

		// Used in div calc & in MHD (Compute_nlin_offdiag(W))
		Array<Complex,3> X2;					

		//!  \f$ (local_{N2}, N_1, N_3/2+1) \f$. 	
		Array<Real,3> Xr;
		Array<Real,3> Xr_slab; 	//Used for reading and writing HDF5

		// Ar not reqd NOW.. CHANGE to Xr2
		// BOTH Ar and Xr are reqd. Xr is the temp array..Ar is the real space array.
		//!  temp array \f$ (local_{Ny}, N_x, N_z/2+1) \f$.
		// For transforms.. MUST
		// Also in Compute_nlin_offdiag(W)
		Array<Real,3> Xr2;


		// for the function satisfy reality condition
		Array<Complex,2> plane_xy;  // Ny,Nx 	
		Array<Complex,2> plane_xy_inproc;  // local_Nx,Ny orlocal_Nx_hor,Ny 
		
		// for Helmholtz solver: Chebyshev basis
		Array<Complex,1> in_helm_complex;    // Nx
		Array<Complex,1> out_helm_complex;    // Nx
   
		Array<Complex,1> V1_x;
		
		// pl remove to
		Array<Real,1> in_helm_real;    // Nx
		Array<Real,1> out_helm_real;    // Nx

		
		Array<Real,3> pressure_plus,pressure_minus;   // local_Nz,Ny,Nx
		Array<Real,1> pressure_derivative_real_1d;
		
		Array<Real,3> vx_plus, vx_minus;      // local_Nz,Ny, Nx
		Array<Real,3> Xreal;
		
		Array<Real,4> influence_matrix;       // local_Nz,Ny,2,2
		//ChFF end
		
		/////////////////////////////////////////////////////////////////////////////
		/////////////// ONLY FOR ADAPTIVE FORCE!!!///////////////////////////////////

		// Array<Complex,3> ForceA1;	
		// Array<Complex,3> ForceA2;	
		// Array<Complex,3> ForceA3;	
		
		// Array<Complex,3> ForceB1;	
		// Array<Complex,3> ForceB2;	
		// Array<Complex,3> ForceB3;	
		
		// Array<Complex,3> ForceC1;	
		// Array<Complex,3> ForceC2;	
		// Array<Complex,3> ForceC3;	

		/////////////// END OF ADAPTIVE FORCE!!!///////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////

	} temp_array;


	struct myconstant{
		const Complex I;

		/// minusI = -sqrt(-1).			
		const Complex minusI;

		/// minus2I = -2*sqrt(-1).			
		const Complex minus2I;	

		// Number of digits for output files
		// NOTE:  for double only.. For float put MY_PRECISION = 6
		//const int MY_PRECISION = 12;


		const Real MYEPS;
		const Real MYEPS2;

		const int MY_MAX_INT;
		// cut off while reading diagnostic_procedure() array from input file and similar ops

		/// Infinite radius.. All the modes outside -- for flux and shelltr calc
		const Real INF_RADIUS;

		const Real INF_TIME; 
		
		const int MAX_NO_GLOB_BUFFER_PACKETS;
		const int MAX_NO_PROBE_PACKETS;

		myconstant();
	} myconstant;

	Global();

	void Parse(int argc, char** argv, bool parse_para=true);
	void Read();
	void Init_defaults();
	void Process_basic_vars();
	void Process_advanced_vars();
	void Assign_parameters();
	void Alias_some_global_vars();
	bool Input_provided(const YAML::Node& node, const string parameter);
	
	template <typename T> 
	void Assign_if_input_provided(const YAML::Node& node, const string parameter, T &var, T default_value);
	
	void Show_error(string message);

	void Print();
};

#endif

//******************************* End of Global.h ********************************
