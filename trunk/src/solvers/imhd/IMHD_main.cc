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

/*! \file Ifluid_main.cc 
 * 
 * @brief Main program executing fluid turbulence
 *
 * @author  M. K. Verma, A. G. Chatterjee
 * @date 2 October 2008
 * 
 * @bugs  No known bug
 */ 

#include "IMHD_main.h"


//****************************************************************************************					

int IMHD_main()
{
    // ITERATION...
    if (global.program.iter_or_diag == "ITERATION") {
        
        fluidIO_incompress.Open_files();
        fluidIO_incompress.Init_energy_transfer();
        
        FluidVF  U(global.field.diss_coefficients[0], global.field.hyper_diss_coefficients[0], global.field.hyper_diss_exponents[0], global.force.U_switch, "U");
        
        FluidVF  B(global.field.diss_coefficients[1], global.field.hyper_diss_coefficients[1], global.field.hyper_diss_exponents[1], global.force.W_switch, "B");
   
        FluidVF helicalU("helicalU");
        FluidVF helicalW("helicalW");

        Pressure P;

        FORCE  Force;

        Time_advance_incompress  time_advance_incompress;

        fluidIO_incompress.Read_init_cond(U, B);
        
        Real total_abs_div;
        U.Compute_divergence_field(global.temp_array.X2, total_abs_div, true);
        // true mean print nonzero div modes
            if (total_abs_div > MYEPS2) {
                cout << "max(abs(Divergence)) of U = " << total_abs_div << "is large. " << '\n' << "Therefore exiting the program." << endl;
                MPI_Abort(MPI_COMM_WORLD, 1); 
            }
        
        B.Compute_divergence_field(global.temp_array.X2, total_abs_div, true);
            if (total_abs_div > MYEPS2) {
                cout << "max(abs(Divergence)) of B = " << total_abs_div << "is large. " << '\n' << "Therefore exiting the program." << endl;
                MPI_Abort(MPI_COMM_WORLD, 1); 
            }

        // Real A1, B1, C1, A2, B2, C2, A3, B3, C3;
        // A1 = B1 = C1 = A2 = B2 = C2 = A3 = B3 = C3 = 0;
        // Real temp_value;

        int  iter=0;  // iterations 

        global.time.now = global.time.init;

        Real viscosity_update_time=global.time.now;
        
        int stop_flag=0;

        global.force.force_U_lock = global.force.force_B_lock = true;

        Force.Compute_force(U, B);

        fluidIO_incompress.Output_all_inloop(U, B, P, helicalU, helicalW);  // for initial cond
        
        if (my_id == master_id) { 
            cout << endl << "STARTING THE SIMULATION NOW" << endl;
        }

        do 	{

            //if (global.mpi.my_id == global.mpi.master_id) cout<<"started time step"<<". now = "<<global.time.now<< ", dt = " << global.time.dt <<endl;
            
            global.time.dt_computation_done = false;
            //global.io.output_real_field_done = false;
            global.io.output_field_k_done = false;
            global.io.output_pressure_spectrum_done = false;
            global.io.output_pressure_done = false;

            iter++;

            time_advance_incompress.Time_advance_step(U, B, P, Force);

            if (global.program.kind == "KEPLERIAN")  {
                time_advance_incompress.Make_field_incompressible(U);
                time_advance_incompress.Make_field_incompressible(B);
            }

            U.Compute_divergence_field(global.temp_array.X2, total_abs_div, true);
            // true mean print nonzero div modes
            if ((total_abs_div > MYEPS2) && (my_id == master_id)) {
                cout << "max(abs(Divergence)) of U = " << total_abs_div << " is large. " << '\n' << "Therefore exiting the program." << endl;
                MPI_Abort(MPI_COMM_WORLD, 1); 
            }
            
            B.Compute_divergence_field(global.temp_array.X2, total_abs_div, true);
            // true mean print nonzero div modes
            if ((total_abs_div > MYEPS2) && (my_id == master_id)) {
                cout << "max(abs(Divergence)) of B = " << total_abs_div << " is large. " << '\n' << "Therefore exiting the program." << endl;
                MPI_Abort(MPI_COMM_WORLD, 1); 
            }

            U.Make_incomplete_shells_zero();
            B.Make_incomplete_shells_zero();


            global.force.force_U_lock = global.force.force_B_lock = false;

            Force.Compute_force(U, B);

            if (U.force_switch) {
                universal->Compute_vorticity(U.cvf.V1, U.cvf.V2, U.cvf.V3, helicalU.cvf.V1, helicalU.cvf.V2, helicalU.cvf.V3, 0, universal->Max_radius_inside());
                global.force.force_injection = Correlation::Get_energy_injection_rate(U);
                global.force.force_injection_helical = Correlation::Get_injection_rate(U, helicalU);
                global.force.force_injection_crosshelical = Correlation::Get_injection_rate(U, B);

                global.force.total_injected_kin_energy += global.time.dt * global.force.force_injection;
                global.force.total_injected_kin_helicity += global.time.dt * global.force.force_injection_helical;
                global.force.total_injected_kin_crosshel += global.time.dt * global.force.force_injection_crosshelical;

                if (global.force.force_alpha_is_adaptive){
                    global.force.force_alpha = global.force.force_alpha +//* exp(
                        (global.time.dt / global.force.force_T) * (global.force.force_I0 - global.force.force_injection) // / global.force.force_I0
                        ;//);
                }
                if (global.force.force_beta_is_adaptive) {
                    global.force.force_beta = global.force.force_beta +
                        (global.time.dt / global.force.force_T_helical) * (global.force.force_I0_helical - global.force.force_injection_helical) // / global.force.force_I0_helical
                        ;
                }
                if (global.force.force_gamma_is_adaptive) {
                    global.force.force_gamma = global.force.force_gamma +
                        (global.time.dt / global.force.force_T_crosshelical) * (global.force.force_I0_crosshelical - global.force.force_injection_crosshelical) // / global.force.force_I0_crosshelical
                        ;
                }

                // if (my_id == master_id) cout<< "U.F = "<< global.force.force_injection<<endl;
            }

            if (B.force_switch) {
                universal->Compute_vorticity(B.cvf.V1, B.cvf.V2, B.cvf.V3, helicalW.cvf.V1, helicalW.cvf.V2, helicalW.cvf.V3, 0, universal->Max_radius_inside());
                helicalW.cvf.Divide_ksqr();
                global.force.mag_force_injection = Correlation::Get_energy_injection_rate(B);
                global.force.mag_force_injection_helical = Correlation::Get_injection_rate(B, helicalW);
                global.force.mag_force_injection_crosshelical = Correlation::Get_injection_rate(B, U);

                global.force.total_injected_mag_energy += global.time.dt * global.force.mag_force_injection;
                global.force.total_injected_mag_helicity += global.time.dt * global.force.mag_force_injection_helical;
                global.force.total_injected_mag_crosshel += global.time.dt * global.force.mag_force_injection_crosshelical;

                if (global.force.mag_force_alpha_is_adaptive) {
                    global.force.mag_force_alpha = global.force.mag_force_alpha +//* exp(
                        (global.time.dt / global.force.mag_force_T) * (global.force.mag_force_I0 - global.force.mag_force_injection) // / global.force.mag_force_I0
                        ;//);
                }

                if (global.force.mag_force_beta_is_adaptive) {
                    global.force.mag_force_beta = global.force.mag_force_beta +
                        (global.time.dt / global.force.mag_force_T_helical) * (global.force.mag_force_I0_helical - global.force.mag_force_injection_helical) // / global.force.mag_force_I0_helical
                        ;
                }
                if (global.force.mag_force_gamma_is_adaptive) {
                    global.force.mag_force_gamma = global.force.mag_force_gamma +
                        (global.time.dt / global.force.mag_force_T_crosshelical) * (global.force.mag_force_I0_crosshelical - global.force.mag_force_injection_crosshelical) // / global.force.mag_force_I0_crosshelical
                        ;
                }

            }

            fluidIO_incompress.Output_all_inloop(U, B, P, helicalU, helicalW);
            if (global.force.empty_force) {
                U.Force1 = 0;  U.Force2 = 0;  U.Force3 = 0;
                B.Force1 = 0;  B.Force2 = 0;  B.Force3 = 0;
                global.force.empty_force = false;
            }

            if ( (isnan(U.cvf.total_energy) || (isnan(B.cvf.total_energy))))  {
                cout << "ERROR: Numerical Overflow " << endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }

            if ((my_id==master_id) && (clock() > global.time.job_time_final)) {
                stop_flag = 1;
            }

            MPI_Bcast(&stop_flag,1, MPI_INT, master_id, MPI_COMM_WORLD);
        } 
        // while ( (global.time.now < global.time.final) && (clock() < global.time.job_time_final) );
        while ( (global.time.now < global.time.final) && (!stop_flag));
        
        fluidIO_incompress.Output_last(U, B, P, helicalU, helicalW);
        
        fluidIO_incompress.Close_files();
    }

    else if (global.program.iter_or_diag == "FILTER") {

        fluidIO_incompress.Open_files();
        fluidIO_incompress.Init_energy_transfer();
        
        FluidVF  U(global.field.diss_coefficients[0], global.field.hyper_diss_coefficients[0], global.field.hyper_diss_exponents[0], global.force.U_switch, "U");
        
        FluidVF  B(global.field.diss_coefficients[1], global.field.hyper_diss_coefficients[1], global.field.hyper_diss_exponents[1], global.force.W_switch, "B");
   
        FluidVF helicalU("helicalU");
        FluidVF helicalW("helicalW");

        Pressure P;

        FORCE  Force;

        Time_advance_incompress  time_advance_incompress;

        fluidIO_incompress.Read_init_cond(U, B);

        global.time.now = 0;

        if (my_id == master_id) cout<<"STARTED FILTERING MODE..."<<endl;

        Real kr = global.io.int_para(0); kr=(kr==0 ? 67.0 : kr);
        Real fsize=global.io.int_para(1); fsize=(fsize==0 ? 3.0 : fsize);
        int filter_num=global.io.int_para(2); filter_num=(filter_num==0 ? 3 : filter_num);
        Real kf;

        if (my_id == master_id) cout<<"Cube will be filtered with kr = "<<kr<<", "<<filter_num<<" times."<<endl;

        U.Write_real_field("U");
        B.Write_real_field("B");

        helicalU.cvf.Set_field_name("helicalU");
        helicalW.cvf.Set_field_name("helicalW");
    
        universal->Compute_vorticity(B.cvf.V1, B.cvf.V2, B.cvf.V3, helicalW.cvf.V1, helicalW.cvf.V2, helicalW.cvf.V3, 0, universal->Max_radius_inside());
        helicalW.Write_real_field("J");
        universal->Compute_vorticity(U.cvf.V1, U.cvf.V2, U.cvf.V3, helicalU.cvf.V1, helicalU.cvf.V2, helicalU.cvf.V3, 0, universal->Max_radius_inside());
        helicalU.Write_real_field("W");


        for(int iii=1; iii<=filter_num; iii++) {

            kf = fsize*iii;
            if (my_id == master_id) cout<<"Stated filter with kf = "<<kf<<"."<<endl;
            //I use helical fields as a buffer for filtering

            helicalU.cvf.Set_field_name("uFaround");
            helicalW.cvf.Set_field_name("bFaround");
    
            helicalU.dissipation_coefficient = U.dissipation_coefficient;
            helicalW.dissipation_coefficient = B.dissipation_coefficient;
    
            U.Filter_field_and_write_real(kr, 1, 2, kf, "uFaround",helicalU.cvf);
            B.Filter_field_and_write_real(kr, 1, 2, kf, "bFaround",helicalW.cvf);
            
            // I need some spectra for this filtered fields
    
            if (global.mpi.master) {
                fluidIO_incompress.spectrum_file << "\t!!! Filtered spectra !!!" << "\n";
                fluidIO_incompress.spectrum_file << "%% Time = " << global.time.now << "\n \n";
            }
            
            Correlation::Compute_shell_spectrum(helicalU);
            fluidIO_incompress.Print_array(fluidIO_incompress.spectrum_file, "Eu(k)", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
            fluidIO_incompress.Print_array(fluidIO_incompress.spectrum_file, "Du(k)", Correlation::shell_dissk1, Correlation::shell_dissk2, Correlation::shell_dissk3);
            
            Correlation::Compute_shell_spectrum(helicalW);
            fluidIO_incompress.Print_array(fluidIO_incompress.spectrum_file, "Eb(k)", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
            fluidIO_incompress.Print_array(fluidIO_incompress.spectrum_file, "Db(k)", Correlation::shell_dissk1, Correlation::shell_dissk2, Correlation::shell_dissk3);
            
            Correlation::Compute_shell_spectrum(helicalU,helicalW);
            fluidIO_incompress.Print_array(fluidIO_incompress.spectrum_file, "Hc(k)", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
    
            Correlation::Compute_shell_spectrum_dissk(helicalU,helicalW);
            fluidIO_incompress.Print_array(fluidIO_incompress.spectrum_file, "Dhc(k)", Correlation::shell_ek1, Correlation::shell_ek2, Correlation::shell_ek3);
    
            if (master)
                fluidIO_incompress.spectrum_file.flush();
    
            helicalU.cvf.Set_field_name("helicalU");
            helicalW.cvf.Set_field_name("helicalW");
    
            universal->Compute_vorticity(B.cvf.V1, B.cvf.V2, B.cvf.V3, helicalW.cvf.V1, helicalW.cvf.V2, helicalW.cvf.V3, 0, universal->Max_radius_inside());
            // helicalW.Write_real_field("J");
            universal->Compute_vorticity(U.cvf.V1, U.cvf.V2, U.cvf.V3, helicalU.cvf.V1, helicalU.cvf.V2, helicalU.cvf.V3, 0, universal->Max_radius_inside());
            // helicalU.Write_real_field("W");
    
            helicalU.Filter_field_and_write_real(kr, 1, 2, kf, "WFaround");
            helicalW.Filter_field_and_write_real(kr, 1, 2, kf, "JFaround");

            global.time.now += 1;
        }

        fluidIO_incompress.Close_files();
    }


    // DIAGNOSTICS
    else if (global.program.iter_or_diag == "DIAGNOSTICS") {
        string filename;
        
        FluidVF  U(global.field.diss_coefficients[0], global.field.hyper_diss_coefficients[0], global.field.hyper_diss_exponents[0], global.force.U_switch, "U");
        FluidVF  B(global.field.diss_coefficients[1], global.field.hyper_diss_coefficients[1], global.field.hyper_diss_exponents[1], global.force.W_switch, "B");

        FluidVF helicalU("helicalU");
        FluidVF helicalW("helicalW");

        
        Pressure P;

        fluidIO_incompress.Init_energy_transfer();
        fluidIO_incompress.Read_init_cond(U, B);
        
        int i=0;
        while (i< global.io.diagnostic_procedures.size()) {
            switch (global.io.diagnostic_procedures[i])  {
            
                case (0) : {
                    filename = "/out/glob.d";
                    filename = global.io.data_dir+ filename;
                    fluidIO_incompress.global_file.open(filename.c_str());
                    if (!fluidIO_incompress.global_file.is_open())
                        cout << "UNABLE TO OPEN FILE global_file (glob.d) " << endl;
                    fluidIO_incompress.Output_global(U, B);
                    fluidIO_incompress.Close_files();
                    break;
                }
                    
				case (1) : {
                    filename = "/out/spectrum.d";
                    filename = global.io.data_dir+ filename;
                    fluidIO_incompress.spectrum_file.open(filename.c_str());
                    fluidIO_incompress.Output_shell_spectrum(U, B, helicalU, helicalW);
                    fluidIO_incompress.Close_files();
                    break;
                }
                    
				case (2) : {
                    filename = "/out/ring_spectrum.d";
                    filename = global.io.data_dir+ filename;
                    fluidIO_incompress.ring_spectrum_file.open(filename.c_str());
                    fluidIO_incompress.Output_ring_spectrum(U, B);
                    fluidIO_incompress.Close_files();
                    break;
                }
                    
                    
				case (3) : {
                    filename = "/out/cyl_ring_spectrum.d";
                    filename = global.io.data_dir+ filename;
                    fluidIO_incompress.cylindrical_ring_spectrum_file.open(filename.c_str());
                    fluidIO_incompress.Output_cylindrical_ring_spectrum(U, B);
                    fluidIO_incompress.Close_files();
                    break;
                }
                    
				case (4) : {
                    filename = "/out/flux.d";
                    filename = global.io.data_dir+ filename;
                    fluidIO_incompress.flux_file.open(filename.c_str());
                    fluidIO_incompress.Output_flux(U, B, P, helicalU, helicalW);
                    fluidIO_incompress.Close_files();
                    break;
                }
                    
                    // force reqd for force-feed calculations
				case (5) : {
                    filename = "/out/shell_to_shell.d";
                    filename = global.io.data_dir+ filename;
                    fluidIO_incompress.shell_to_shell_file.open(filename.c_str());
                    fluidIO_incompress.Output_shell_to_shell(U, B, P, helicalU, helicalW);
                    fluidIO_incompress.Close_files();
                    break;
                }
                    
                case (6) : {
                    filename = "/out/ring_to_ring.d";
                    filename = global.io.data_dir+ filename;
                    fluidIO_incompress.ring_to_ring_file.open(filename.c_str());
                    fluidIO_incompress.Output_ring_to_ring(U, B, P, helicalU, helicalW);
                    fluidIO_incompress.Close_files();
                    break;
                }
                    
                    
				case (7) : {
                    filename = "/out/cylindrical_ring_to_ring.d";
                    filename = global.io.data_dir+ filename;
                    fluidIO_incompress.cylindrical_ring_to_ring_file.open(filename.c_str());
                    fluidIO_incompress.Output_cylindrical_ring_to_ring(U, B, P);
                    fluidIO_incompress.Close_files();
                    break;
                }
                    //		case (7) : fluidIO_incompress.Output_structure_fn(U);  break;
                    //		case (8) : fluidIO_incompress.Output_planar_structure_fn(U);  break;
				case (10) : {
                    filename = "/out/field_k_out_"+To_string(my_id)+".d";
                    filename = global.io.data_dir+ filename;
                    fluidIO_incompress.field_k_out_file.open(filename.c_str());
                    fluidIO_incompress.Output_field_k(U, B);
                    fluidIO_incompress.Close_files();
                    break;
                }
                    
				case (11) : {
                    filename = "/out/field_r_out_"+To_string(my_id)+".d";
                    filename = global.io.data_dir+ filename;
                    fluidIO_incompress.field_r_out_file.open(filename.c_str());
                    fluidIO_incompress.Output_field_r(U, B);
                    fluidIO_incompress.Close_files();
                    break;
                }
                    
				case (13) : {
                    U.Inverse_transform();
                    B.Inverse_transform();
                    fluidIO_incompress.Output_real_field(U, B);
                    break;
                }
				case (14) : {
                    filename = "/out/field_out_reduced.d";
                    filename = global.io.data_dir+ filename;
                    fluidIO_incompress.field_out_reduced_file.open(filename.c_str());
                    fluidIO_incompress.Output_reduced_complex_field(U, B);
                    fluidIO_incompress.Close_files();
                    break;
                }
            }

                i++;
        }	

    }

    return(1);

} 


//********************************** End of Ifluid_main.cc ************************************	



