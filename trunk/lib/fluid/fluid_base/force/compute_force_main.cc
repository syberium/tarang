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

/*! \file  compute_force_main.cc
 * 
 * @brief  Compute force and put it in F.
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Dec. 2008
 *
 * @bug  No known bugs
 */

#include "FORCE.h"

//*********************************************************************************************


void FORCE::Compute_force(FluidVF& U)
{
	switch (global.force.field_procedure) {		
		case (0) : Compute_force_decay(U);	break;		
		case (1) : Compute_force_const_energy_helicity_supply(U); break;	
		case (2) : Compute_force_const_energy_helicity(U); break;
		case (3) : Compute_force_given_modes(U); break;
		case (4) : Compute_force_Taylor_Green(U); break;
		case (5) : Compute_force_ABC(U); break;
		case (6) : Compute_force_using_random_energy_helicity_spectrum(U); break;
		case (21) : Compute_force_Coriolis(U); break;
		case (22) : Compute_force_Keplerian(U); break;
		case (23) : Compute_force_Keplerian_SB(U); break;

		case (24) : Compute_force_hydro_helical_decompose(U); break;

		case (101): Compute_force_Liquid_metal(U); break;
		case (102): Compute_force_Kolmogorov_flow(U); break;
		case (103): Compute_force_Liquid_metal_const_energy_supply(U); break;
		case (104): Compute_force_Ekman_friction(U); break;
		case (105): Compute_force_Ekman_friction_const_energy_supply(U); break;	
		case (420): Compute_force_pressure_grad(U); break;
	}
	
	int force_proc = global.force.field_procedure;
	if (U.force_switch) {
		if ((force_proc==1) || (force_proc==2)) {
			U.Satisfy_strong_reality_condition_force_field();
			
			if ((basis_type == "SFF") || (basis_type == "SSF") || (basis_type == "SSS"))
				universal->Zero_modes(U.Force1, U.Force2, U.Force3);
		}
	}
}

//
//

void FORCE::Compute_force(FluidVF& U, FluidSF& T)
{
	switch (global.force.field_procedure) {		
		case (0) : Compute_force_decay(U, T);	break;		
		case (1) : Compute_force_const_energy_helicity_supply(U, T); break;
		case (2) : Compute_force_const_energy_helicity(U, T); break;		
		case (3) : Compute_force_given_modes(U, T); break;
		case (4) : Compute_force_Taylor_Green(U, T); break;
		case (5) : Compute_force_ABC(U, T); break;
		case (6) : Compute_force_using_random_energy_helicity_spectrum(U, T); break;
		case (11): Compute_force_const_energy_helicity_supply(U); break;
		case (21) : Compute_force_Coriolis(U, T); break;	
		case (51) : Compute_force_RBC(U, T); break;	
		case (52) : Compute_force_RBC_rotation(U, T); break;
		case (53) : Compute_force_stratified_random(U, T); break;
		case (101): Compute_force_Liquid_metal(U, T); break;
	}
	
	int force_proc = global.force.field_procedure;
	if (U.force_switch) {
		if ((force_proc==1) || (force_proc==2)) {
			U.Satisfy_strong_reality_condition_force_field();
			
			if ((basis_type == "SFF") || (basis_type == "SSF") || (basis_type == "SSS"))
				universal->Zero_modes(U.Force1, U.Force2, U.Force3);
		}
	}
	
	if (T.force_switch) {
		if ((force_proc==1) || (force_proc==2)) {
			T.Satisfy_strong_reality_condition_force_field();
			
			if ((basis_type == "SFF") || (basis_type == "SSF") || (basis_type == "SSS"))
				universal->Zero_modes(T.Force);
		}
	}
}



void FORCE::Compute_force(FluidVF& U, FluidSF& T1, FluidSF& T2)
{
	switch (global.force.field_procedure) {		
		case (0) : Compute_force_decay(U, T1, T2);	break;		
		case (51) : Compute_force_MRBC(U, T1, T2); break;	
	}
}

//
//

void FORCE::Compute_force(FluidVF& U, FluidVF& W)
{
	switch (global.force.field_procedure) {		
		case (0)   : Compute_force_decay(U, W);	break;				
		case (1)   : Compute_force_const_energy_helicity_supply(U, W); break;	
		case (2)   : Compute_force_const_energy_helicity(U, W); break;	
		case (3)   : Compute_force_given_modes(U, W); break;
		case (4)   : Compute_force_Taylor_Green(U, W); break;
		case (5)   : Compute_force_ABC(U, W); break;
		case (6)   : Compute_force_using_random_energy_helicity_spectrum(U, W); break;
		case (21)  : Compute_force_Coriolis(U, W); break;
		case (22)  : Compute_force_Keplerian(U, W); break;
		case (23)  : Compute_force_Keplerian_SB(U, W); break;
		case (24)  : Compute_force_crosshelicity_supply(U, W); break;
		
		case (25)  : Main_random_direction_forcing(U, W); break;
		case (26)  : Compute_random_crosshelical_forcing(U, W); break;
		case (261) : Compute_random_crosshelical_forcing_addition(U, W); break;
		
		case (27)  : Compute_force_MHD_helical_decompose(U,W); break;
		
		case (271)  : Compute_force_MHD_helical_decompose_shell_distribution(U,W); break;

		case (2715)  : Compute_random_crosshelical_forcing_addition_plus_helical_decompose_kinetic_helicity(U,W); break;
		case (2716)  : Compute_random_crosshelical_forcing_addition_plus_helical_decompose_magnetic_helicity(U,W); break;
		case (2717)  : Compute_random_crosshelical_forcing_addition_plus_helical_decompose_kinetic_magnetic_helicity(U,W); break;

		case (272)  : Compute_force_MHD_helical_decompose_shell_distribution_differ_scales(U,W); break;

		case (28)  : Compute_force_MHD_helical_decompose_with_crosshel(U,W); break;
		case (101) : Compute_force_DYNAMO_SIX_MODE(U, W); break;
	}
	
	int force_proc = global.force.field_procedure;
	if (U.force_switch) {
		if ((force_proc==1) || (force_proc==2)) {
			U.Satisfy_strong_reality_condition_force_field();
			
			if ((basis_type == "SFF") || (basis_type == "SSF") || (basis_type == "SSS"))
				universal->Zero_modes(U.Force1, U.Force2, U.Force3);
		}
	}
	
	if (W.force_switch) {
		if ((force_proc==1) || (force_proc==2)) {
			W.Satisfy_strong_reality_condition_force_field();
			
			if ((basis_type == "SFF") || (basis_type == "SSF") || (basis_type == "SSS"))
				universal->Zero_modes(W.Force1, W.Force2, W.Force3);
		}
	}
}

//
//

void FORCE::Compute_force(FluidVF& U, FluidVF& W, FluidSF& T)
{
	switch (global.force.field_procedure) {		
		case (0) : Compute_force_decay(U, W, T);	break;					
		case (1) : Compute_force_const_energy_helicity_supply(U, W, T); break;	
		case (2) : Compute_force_const_energy_helicity(U, W, T); break;
		case (3) : Compute_force_given_modes(U, W, T); break;
		case (4) : Compute_force_Taylor_Green(U, W, T); break;
		case (5) : Compute_force_ABC(U, W, T); break;
		case (6) : Compute_force_using_random_energy_helicity_spectrum(U, W, T); break;
		case (21) : Compute_force_Coriolis(U, W, T); break;
		case (51) : Compute_force_RBC(U, W,T); break;
		case (52) : Compute_force_RBC_rotation(U, W, T); break;
		case (55) : Compute_force_astro(U, W, T); break;
	}
	
	int force_proc = global.force.field_procedure;
	if (U.force_switch) {
		if ((force_proc==1) || (force_proc==2)) {
			U.Satisfy_strong_reality_condition_force_field();
			
			if ((basis_type == "SFF") || (basis_type == "SSF") || (basis_type == "SSS"))
				universal->Zero_modes(U.Force1, U.Force2, U.Force3);
		}
	}
	
	if (W.force_switch) {
		if ((force_proc==1) || (force_proc==2)) {
			W.Satisfy_strong_reality_condition_force_field();
			
			if ((basis_type == "SFF") || (basis_type == "SSF") || (basis_type == "SSS"))
				universal->Zero_modes(W.Force1, W.Force2, W.Force3);
		}
	}
	
	if (T.force_switch) {
		if ((force_proc==1) || (force_proc==2)) {
			T.Satisfy_strong_reality_condition_force_field();
			
			if ((basis_type == "SFF") || (basis_type == "SSF") || (basis_type == "SSS"))
				universal->Zero_modes(T.Force);
		}
	}
}



void FORCE::Compute_force(FluidVF& U, FluidVF& W, FluidSF& T, FluidSF& C)
{
	switch (global.force.field_procedure) {
		case (55) : Compute_force_astro(U, W, T, C); break;
	}
}



//  GP

void FORCE::Compute_force(FluidSF& T)
{
	
}


//*******************************  End of compute_force_main.cc *******************************



