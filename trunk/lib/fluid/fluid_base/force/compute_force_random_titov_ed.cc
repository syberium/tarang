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

/*! \file  compute_force_random_titov_ed.cc
 * 
 * @brief Compute force when ek and hk supply rate is given
 *
  * @note 3D;   F(k) = alpha * V(k) + beta(k) * Omega(k)
 *				alpha and beta determined from the supply rates
 *
 * @note:   Satisfy reality condition is critical here for kz=0 and N[3]/2 planes.   
 *			Do not remove this function.
 *
 * @author  V. V. Titov
 * @version 4.0 MPI
 * @date Jan. 2018
 *
 * @bug sk=1, -1 needs to be handled separately.
 */

#include "FORCE.h"

extern Uniform<Real> SPECrand;

//*********************************************************************************************

void FORCE::Compute_force_random_direction(FluidVF& U, FluidVF& W)
{
		
	Real inner_radius = global.force.double_para(0);
	Real outer_radius = global.force.double_para(1);
	Real Au = global.force.double_para(2);
	Real Bu = global.force.force_beta;				
	Real Cu = global.force.force_gamma;				

	Real inner_radiusM = global.force.double_para(5);
	Real outer_radiusM = global.force.double_para(6);
	Real Ab = global.force.double_para(7);
	Real Bb = global.force.mag_force_beta;
	Real Cb = global.force.mag_force_gamma;
	
	if (global.force.update_time + global.force.delta_t < global.time.now) {
		global.force.update_time = global.time.now;
		// if (global.mpi.my_id == global.mpi.master_id) cout<<"DEBUG!!! Force updated at "<<global.force.update_time<<endl;
		if (U.force_switch) {
			Kinetic_force_random_direction_assign(U, W, inner_radius, outer_radius, Au, Bu, Cu);
		}

		if (W.force_switch) {
			Magnetic_force_random_direction_assign(U, W, inner_radiusM, outer_radiusM, Ab, Bb, Cb);
		}
	}


}

void FORCE::Kinetic_force_random_direction_assign(FluidVF& U, FluidVF& W, Real inner_radius, Real outer_radius, Real Au, Real Bu, Real Cu){
// from adapting force (case 24)	
	Real helicity_supply = Bu;
	Real crosshelicity_supply = Cu;
	Real helicity_supply_per_mode, crosshelicity_supply_per_mode;
//end of from adapting force (case 24)
	int kx_max, ky_max, kz_max, kx_min, ky_min, kz_min;
	int nf;
// service part
	kx_max = (int) ceil(outer_radius/kfactor[1]);
	
	if (Ny > 1)
		ky_max = (int) ceil(outer_radius/kfactor[2]);
	else
		ky_max = 0;	
	
	kz_max = (int) ceil(outer_radius/kfactor[3]);
	
	kx_min = ky_min = kz_min = 0;
	
	if (basis_type == "FFF" || basis_type == "FFFW")
		kx_min = -kx_max;
	
	if ((basis_type == "FFF" || basis_type == "FFFW") || (basis_type == "SFF"))
		ky_min = -ky_max;

	nf = universal->Get_number_modes_in_shell(inner_radius, outer_radius);
// Exclude modes with 0 components
	int lx, ly, lz;
	Real Kmag;
	for (int kx = kx_min; kx <= kx_max; kx++)
	for (int ky = ky_min; ky <= ky_max; ky++)  
	for (int kz = 0; kz <= kz_max; kz++) {
		if ((kx==0) || (ky==0) || (kz==0)) 
		if (universal->Probe_in_me(kx,ky,kz))  {
			Kmag = sqrt(pow2(kx)+pow2(ky)+pow2(kz));
			if ((Kmag > inner_radius) && (Kmag <= outer_radius)) {
				nf-=(kz==0?1:2);
			}
		}
	}
// end of service part	
	TinyVector<Complex,3> e, k, ExK, f;
	Real norm, alpha, beta, gamma;

	alpha =  sqrt( 2 * Au / (global.force.delta_t * nf) );
// from adapting force (case 24)
	helicity_supply_per_mode = helicity_supply / nf;
	crosshelicity_supply_per_mode = crosshelicity_supply / nf;
	Real modal_energy, modal_vorticity, modal_mag_energy;
	TinyVector<Complex,3> localW, localC, temp;
	Real epsilon = 0.01;
// end of "from adapting force (case 24)
	for (int kx = kx_min; kx <= kx_max; kx++)
	for (int ky = ky_min; ky <= ky_max; ky++)  
	for (int kz = 0; kz <= kz_max; kz++) {
		if ((kx!=0) && (ky!=0) && (kz!=0))
		if (universal->Probe_in_me(kx,ky,kz))  {
			lx = universal->Get_lx(kx);
			ly = universal->Get_ly(ky);
			lz = universal->Get_lz(kz);
				
			Kmag = universal->Kmagnitude(lx, ly, lz);
			if ((Kmag > inner_radius) && (Kmag <= outer_radius)) {
				e = (2*SPECrand.random()-1) + I*(2*SPECrand.random()-1),
					(2*SPECrand.random()-1) + I*(2*SPECrand.random()-1),
					(2*SPECrand.random()-1) + I*(2*SPECrand.random()-1);
				k = kx, ky, kz;

				ExK = mycross(e,k);

				norm = mynorm(ExK);
// from adapting force (case 24)
				U.cvf.Compute_Modal_vorticity(lx, ly, lz, localW);
				localC = universal->Get_local_spectral_field(lx, ly, lz, W.cvf.V1, W.cvf.V2, W.cvf.V3);
				
				modal_energy = 2 * U.cvf.Modal_energy(lx, ly, lz) + epsilon;
				modal_vorticity = Vsqr(localW) + epsilon;
				modal_mag_energy = 2 * W.cvf.Modal_energy(lx, ly, lz) + epsilon;

				//kinetic force
				beta = helicity_supply_per_mode / sqrt(modal_energy * modal_vorticity);
					
				gamma = crosshelicity_supply_per_mode / sqrt(modal_energy * modal_mag_energy);	
// end of from adapting force (case 24)

				// First realisation
				// f = alpha * ExK / norm; 

				f = alpha * ExK / norm + beta * localW + gamma * localC;

				universal->Assign_local_spectral_field(lx, ly, lz, U.Force1, U.Force2, U.Force3, f);

				temp = alpha * ExK / norm;
				localW *= beta;
				localC *= gamma;
				// universal->Assign_local_spectral_field(lx, ly, lz, global.temp_array.ForceA1, global.temp_array.ForceA2, global.temp_array.ForceA3, temp);
				// universal->Assign_local_spectral_field(lx, ly, lz, global.temp_array.ForceB1, global.temp_array.ForceB2, global.temp_array.ForceB3, localW);
				// universal->Assign_local_spectral_field(lx, ly, lz, global.temp_array.ForceC1, global.temp_array.ForceC2, global.temp_array.ForceC3, localC);
			}
		}

	}
	
}

void FORCE::Magnetic_force_random_direction_assign(FluidVF& U, FluidVF& W, Real inner_radius, Real outer_radius, Real Ab, Real Bb, Real Cb){
// from adapting force (case 24)
	Real helicity_supply = Bb;
	Real crosshelicity_supply = Cb;
	Real helicity_supply_per_mode, crosshelicity_supply_per_mode;
//end of from adapting force (case 24)
	int kx_max, ky_max, kz_max, kx_min, ky_min, kz_min;
	int nf;
// service part
	kx_max = (int) ceil(outer_radius/kfactor[1]);
	
	if (Ny > 1)
		ky_max = (int) ceil(outer_radius/kfactor[2]);
	else
		ky_max = 0;	
	
	kz_max = (int) ceil(outer_radius/kfactor[3]);
	
	kx_min = ky_min = kz_min = 0;
	
	if (basis_type == "FFF" || basis_type == "FFFW")
		kx_min = -kx_max;
	
	if ((basis_type == "FFF" || basis_type == "FFFW") || (basis_type == "SFF"))
		ky_min = -ky_max;

	nf = universal->Get_number_modes_in_shell(inner_radius, outer_radius);
// Exclude modes with 0 components
	int lx, ly, lz;
	Real Kmag;
	for (int kx = kx_min; kx <= kx_max; kx++)
	for (int ky = ky_min; ky <= ky_max; ky++)  
	for (int kz = 0; kz <= kz_max; kz++) {
		if ((kx==0) || (ky==0) || (kz==0)) 
		if (universal->Probe_in_me(kx,ky,kz))  {
			Kmag = sqrt(pow2(kx)+pow2(ky)+pow2(kz));
			if ((Kmag > inner_radius) && (Kmag <= outer_radius)) {
				nf-=(kz==0?1:2);
			}
		}
	}
// end of service part	
	TinyVector<Complex,3> e, k, ExK, f;
	Real norm, alpha, beta, gamma, k2;

	alpha =  sqrt( 2 * Ab / (global.force.delta_t * nf) );
// from adapting force (case 24)
	helicity_supply_per_mode = helicity_supply / nf;
	crosshelicity_supply_per_mode = crosshelicity_supply / nf;
	Real modal_energy, modal_potential, modal_kin_energy;
	TinyVector<Complex,3> localW, localC, temp;
	Real epsilon = 0.01;
// end of "from adapting force (case 24)
	for (int kx = kx_min; kx <= kx_max; kx++)
	for (int ky = ky_min; ky <= ky_max; ky++)  
	for (int kz = 0; kz <= kz_max; kz++) {
		if ((kx!=0) && (ky!=0) && (kz!=0))
		if (universal->Probe_in_me(kx,ky,kz))  {
			lx = universal->Get_lx(kx);
			ly = universal->Get_ly(ky);
			lz = universal->Get_lz(kz);
				
			Kmag = universal->Kmagnitude(lx, ly, lz);
			if ((Kmag > inner_radius) && (Kmag <= outer_radius)) {
				e = (2*SPECrand.random()-1) + I*(2*SPECrand.random()-1),
					(2*SPECrand.random()-1) + I*(2*SPECrand.random()-1),
					(2*SPECrand.random()-1) + I*(2*SPECrand.random()-1);
				k = kx, ky, kz;

				ExK = mycross(e,k);

				norm = mynorm(ExK);
// from adapting force (case 24)
				W.cvf.Compute_Modal_vorticity(lx, ly, lz, localW);
				localW/=(Kmag==0?1:Kmag*Kmag);
				localC = universal->Get_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3);

				modal_energy = 2 * W.cvf.Modal_energy(lx, ly, lz) + epsilon;
				modal_potential = Vsqr(localW) + epsilon;
				modal_kin_energy = 2 * U.cvf.Modal_energy(lx, ly, lz) + epsilon;

				beta = helicity_supply_per_mode / sqrt(modal_energy * modal_potential);

				gamma = crosshelicity_supply_per_mode / sqrt(modal_energy * modal_kin_energy);

				temp = alpha * ExK / norm + beta * localW + gamma * localC;

				k2=Kmag*Kmag;

			    f(0) = (1-kx*kx/k2)*temp(0) -     kx*ky/k2*temp(1) -     kx*kz/k2*temp(2);
			    f(1) =    -ky*kx/k2*temp(0) + (1-ky*ky/k2)*temp(1) -     ky*kz/k2*temp(2);
			    f(2) =    -kz*kx/k2*temp(0) -     kz*ky/k2*temp(1) + (1-kz*kz/k2)*temp(2);

// end of from adapting force (case 24)

				// First realisation
				// f = alpha * ExK / norm; 

				universal->Assign_local_spectral_field(lx, ly, lz, W.Force1, W.Force2, W.Force3, f);
			}
		}

	}
}
//*****************************************************************
//Crosshelical random force
void FORCE::Compute_random_crosshelical_forcing(FluidVF& U, FluidVF& W) {
	Real inner_radius = global.force.double_para(0);
	Real outer_radius = global.force.double_para(1);
	Real Au = global.force.double_para(2);
	Real Bb = global.force.double_para(3);				
	Real Cu = global.force.double_para(4);
	Real Cb = global.force.double_para(5);

	if (global.force.update_time + global.force.delta_t < global.time.now) {
		global.force.update_time = global.time.now;

		if (U.force_switch && W.force_switch) {
			Kinetic_and_magnetic_random_crosshelical_forcing_assign(U, W, inner_radius, outer_radius, Au, Bb, Cu, Cb);

		}
	}

}

void FORCE::Kinetic_and_magnetic_random_crosshelical_forcing_assign(FluidVF& U, FluidVF& W, Real inner_radius, Real outer_radius, Real Au, Real Bb, Real Cu, Real Cb){ 
	int kx_max, ky_max, kz_max, kx_min, ky_min, kz_min;
	int nf;
// service part
	kx_max = (int) ceil(outer_radius/kfactor[1]);
	
	if (Ny > 1)
		ky_max = (int) ceil(outer_radius/kfactor[2]);
	else
		ky_max = 0;	
	
	kz_max = (int) ceil(outer_radius/kfactor[3]);
	
	kx_min = ky_min = kz_min = 0;
	
	if (basis_type == "FFF" || basis_type == "FFFW")
		kx_min = -kx_max;
	
	if ((basis_type == "FFF" || basis_type == "FFFW") || (basis_type == "SFF"))
		ky_min = -ky_max;

	nf = universal->Get_number_modes_in_shell(inner_radius, outer_radius);
// Exclude modes with 0 components
	int lx, ly, lz;
	Real Kmag;
	for (int kx = kx_min; kx <= kx_max; kx++)
	for (int ky = ky_min; ky <= ky_max; ky++)  
	for (int kz = 0; kz <= kz_max; kz++) {
		if ((kx==0) || (ky==0) || (kz==0)) 
		if (universal->Probe_in_me(kx,ky,kz))  {
			Kmag = sqrt(pow2(kx)+pow2(ky)+pow2(kz));
			if ((Kmag > inner_radius) && (Kmag <= outer_radius)) {
				nf-=(kz==0?1:2);
			}
		}
	}
// end of service part		
	TinyVector<Complex,3> eu, eb, ec, k, e1, e2, e3, fu, fb;
	Real alpha, beta, gammaU, gammaB, k2;

	alpha  = sqrt( 2 * sqrt(2.0) * Au / (global.force.delta_t * nf) );
	beta  = sqrt( 2 * sqrt(2.0) * Bb / (global.force.delta_t * nf) );
	gammaU = sqrt( 2 * sqrt(2.0) * Cu / (global.force.delta_t * nf) );
	gammaB = sqrt( 2 * sqrt(2.0) * Cb / (global.force.delta_t * nf) );

	TinyVector<Complex,3> temp;

	for (int kx = kx_min; kx <= kx_max; kx++)
	for (int ky = ky_min; ky <= ky_max; ky++)  
	for (int kz = 0; kz <= kz_max; kz++) {
		if ((kx!=0) && (ky!=0) && (kz!=0))
		if (universal->Probe_in_me(kx,ky,kz))  {
			lx = universal->Get_lx(kx);
			ly = universal->Get_ly(ky);
			lz = universal->Get_lz(kz);
				
			Kmag = universal->Kmagnitude(lx, ly, lz);
			if ((Kmag > inner_radius) && (Kmag <= outer_radius)) {
// preparing vectors
				e1 = (2*SPECrand.random()-1) + I*(2*SPECrand.random()-1),
					 (2*SPECrand.random()-1) + I*(2*SPECrand.random()-1),
					 (2*SPECrand.random()-1) + I*(2*SPECrand.random()-1);

				e2 = (2*SPECrand.random()-1) + I*(2*SPECrand.random()-1),
					 (2*SPECrand.random()-1) + I*(2*SPECrand.random()-1),
					 (2*SPECrand.random()-1) + I*(2*SPECrand.random()-1);

				e3 = (2*SPECrand.random()-1) + I*(2*SPECrand.random()-1),
					 (2*SPECrand.random()-1) + I*(2*SPECrand.random()-1),
					 (2*SPECrand.random()-1) + I*(2*SPECrand.random()-1);

				k = kx, ky, kz;

				eu = mycross(k,e1);
				eu /= mynorm(eu);

				ec = mycross(k,e2); 
				ec /= mynorm(ec);

				eb = mycross(k,e3);
				eb /= mynorm(eb);
//kinetic force
				fu = alpha * eu + gammaU * ec;			
//magnetic force
				temp = beta * eb + gammaB * ec;
//div filter
				k2=Kmag*Kmag;

			    fb(0) = (1-kx*kx/k2)*temp(0) -     kx*ky/k2*temp(1) -     kx*kz/k2*temp(2);
			    fb(1) =    -ky*kx/k2*temp(0) + (1-ky*ky/k2)*temp(1) -     ky*kz/k2*temp(2);
			    fb(2) =    -kz*kx/k2*temp(0) -     kz*ky/k2*temp(1) + (1-kz*kz/k2)*temp(2);
//assign forces
				universal->Assign_local_spectral_field(lx, ly, lz, U.Force1, U.Force2, U.Force3, fu);
				universal->Assign_local_spectral_field(lx, ly, lz, W.Force1, W.Force2, W.Force3, fb);

				// if (my_id == master_id) 
				// 	cout<<"Crosshelical force DEBUG:\n"
				// 		<<"[Au, Bb, Cu, Cb] = ["<<Au<<", "<<Bb<<", "<<Cu<<", "<<Cb<<"]\n"
				// 		<<"[alpha, beta, gammaU, gammaB] = ["<<alpha<<", "<<beta<<", "<<gammaU<<", "<<gammaB<<"]\n"
				// 		<<"----------------------------------------------;"<<endl<<endl;
			}
		}
	}
}

//*****************************************************************
//Crosshelical random force which adding directly to the fields without assign force fields
void FORCE::Compute_random_crosshelical_forcing_addition(FluidVF& U, FluidVF& W) {
	Real inner_radius = global.force.double_para(0);
	Real outer_radius = global.force.double_para(1);
	Real Au = global.force.double_para(2);
	Real Bb = global.force.double_para(3);				
	Real Cu = global.force.double_para(4);
	Real Cb = global.force.double_para(5);


	if ((!global.force.force_U_lock) && (!global.force.force_B_lock)) {
		global.force.force_U_lock = global.force.force_B_lock = true;

			if (U.force_switch && W.force_switch) {
				Kinetic_and_magnetic_random_crosshelical_forcing_addition(U, W, inner_radius, outer_radius, Au, Bb, Cu, Cb);
			}
	}

}

void FORCE::Kinetic_and_magnetic_random_crosshelical_forcing_addition(FluidVF& U, FluidVF& W, Real inner_radius, Real outer_radius, Real Au, Real Bb, Real Cu, Real Cb){ 
	int kx_max, ky_max, kz_max, kx_min, ky_min, kz_min;
	int nf;
// service part
	kx_max = (int) ceil(outer_radius/kfactor[1]);
	
	if (Ny > 1)
		ky_max = (int) ceil(outer_radius/kfactor[2]);
	else
		ky_max = 0;	
	
	kz_max = (int) ceil(outer_radius/kfactor[3]);
	
	kx_min = ky_min = kz_min = 0;
	
	if (basis_type == "FFF" || basis_type == "FFFW")
		kx_min = -kx_max;
	
	if ((basis_type == "FFF" || basis_type == "FFFW") || (basis_type == "SFF"))
		ky_min = -ky_max;

//	nf = universal->Get_number_modes_in_shell(inner_radius, outer_radius);
// Calculate modes without 0 components
	nf=0;
	int lx, ly, lz;
	Real Kmag;
	for (int kx = kx_min; kx <= kx_max; kx++)
	for (int ky = ky_min; ky <= ky_max; ky++)  
	for (int kz = 0; kz <= kz_max; kz++) {
		if ((kx!=0) && (ky!=0) && (kz!=0)) {
			Kmag = sqrt(pow2(kx)+pow2(ky)+pow2(kz));
			if ((Kmag > inner_radius) && (Kmag <= outer_radius)) {
				nf+=(kz==0?1:2);
			}
		}
	}
// end of service part	
// for output purposes and compute eps_j_i rates
	Real in[6], out[6];
// end

	TinyVector<Complex,3> eu, eb, ec, k, e1, e2, e3, fu, fb, vectorU, vectorB;
	Real alpha, beta, gammaU, gammaB, k2;
	Real eps_E_a, eps_E_b, eps_E_g, eps_E_d;
	Real eps_C_a, eps_C_b, eps_C_g, eps_C_d;

	eps_E_a = eps_E_b = eps_E_g = eps_E_d = eps_C_a = eps_C_b = eps_C_g = eps_C_d = 0.0;

	alpha  = sqrt( 2.0 * Au / (global.time.dt * nf) );
	beta  = sqrt( 2.0 * Bb / (global.time.dt * nf) );
	gammaU = sqrt( 2.0 * Cu / (global.time.dt * nf) );
	gammaB = sqrt( 2.0 * Cb / (global.time.dt * nf) );

	TinyVector<Complex,3> temp;

	for (int kx = kx_min; kx <= kx_max; kx++)
	for (int ky = ky_min; ky <= ky_max; ky++)  
	for (int kz = 0; kz <= kz_max; kz++) {
		if ((kx!=0) && (ky!=0) && (kz!=0))
		if (universal->Probe_in_me(kx,ky,kz))  {
			lx = universal->Get_lx(kx);
			ly = universal->Get_ly(ky);
			lz = universal->Get_lz(kz);
				
			Kmag = universal->Kmagnitude(lx, ly, lz);
			if ((Kmag > inner_radius) && (Kmag <= outer_radius)) {
// preparing vectors
				e1 = (2*SPECrand.random()-1) + I*(2*SPECrand.random()-1),
					 (2*SPECrand.random()-1) + I*(2*SPECrand.random()-1),
					 (2*SPECrand.random()-1) + I*(2*SPECrand.random()-1);

				e2 = (2*SPECrand.random()-1) + I*(2*SPECrand.random()-1),
					 (2*SPECrand.random()-1) + I*(2*SPECrand.random()-1),
					 (2*SPECrand.random()-1) + I*(2*SPECrand.random()-1);

				e3 = (2*SPECrand.random()-1) + I*(2*SPECrand.random()-1),
					 (2*SPECrand.random()-1) + I*(2*SPECrand.random()-1),
					 (2*SPECrand.random()-1) + I*(2*SPECrand.random()-1);

				k = kx, ky, kz;

				eu = mycross(k,e1);
				eu /= mynorm(eu);

				ec = mycross(k,e2); 
				ec /= mynorm(ec);

				eb = mycross(k,e3);
				eb /= mynorm(eb);
//kinetic force
				fu = alpha * eu + gammaU * ec;			
//magnetic force
				fb = beta * eb + gammaB * ec;

				fu*=global.time.dt;
				fb*=global.time.dt;
//div filter
				// k2=Kmag*Kmag;

			 //    fb(0) = (1-kx*kx/k2)*temp(0) -     kx*ky/k2*temp(1) -     kx*kz/k2*temp(2);
			 //    fb(1) =    -ky*kx/k2*temp(0) + (1-ky*ky/k2)*temp(1) -     ky*kz/k2*temp(2);
			 //    fb(2) =    -kz*kx/k2*temp(0) -     kz*ky/k2*temp(1) + (1-kz*kz/k2)*temp(2);

//assign forces
				universal->Add_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3, fu);
				universal->Add_local_spectral_field(lx, ly, lz, W.cvf.V1, W.cvf.V2, W.cvf.V3, fb);

				//for output and controlling purposes
				vectorU = universal->Get_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3);
				vectorB = universal->Get_local_spectral_field(lx, ly, lz, W.cvf.V1, W.cvf.V2, W.cvf.V3);

				temp = alpha * eu * global.time.dt; 	eps_E_a+=mydot(vectorU,temp); eps_C_a+=mydot(vectorB,temp);
				temp = beta * eb * global.time.dt; 		eps_E_b+=mydot(vectorB,temp); eps_C_b+=mydot(vectorU,temp);

				temp = gammaU * ec * global.time.dt; 	eps_E_g+=mydot(vectorU,temp); eps_C_g+=mydot(vectorB,temp); 
				temp = gammaB * ec * global.time.dt; 	eps_E_d+=mydot(vectorB,temp); eps_C_d+=mydot(vectorU,temp); 
			}
		}
	}
	in[0] = eps_E_a; 
	in[1] = eps_E_b;
	in[2] = eps_E_g;
	in[3] = eps_E_d;
	in[4] = eps_C_a; 
	in[5] = eps_C_b;
	in[6] = eps_C_g;
	in[7] = eps_C_d;


	MPI_Reduce(&in[0], &out[0], 8, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);

	global.force.A1 += out[0];
	global.force.A2 += out[1];
	global.force.A3 += out[2];

	global.force.B1 += out[3];
	global.force.B2 += out[4];
	global.force.B3 += out[5];

	global.force.C1 += out[6];
	global.force.C2 += out[7];
}





//*****************************************************************
//Kinetic and magnetic random helical forces (from mahyd 2016 with modifications)
void FORCE::Compute_random_helical_forcing_Stepanov_Titov_mod(FluidVF& U, FluidVF& W) {
	Real inner_radius = global.force.double_para(0);
	Real outer_radius = global.force.double_para(1);
	Real Au = global.force.double_para(2);
	Real Bu = global.force.double_para(3);				
	Real Cu = global.force.double_para(4);

	Real inner_radiusM = global.force.double_para(5);
	Real outer_radiusM = global.force.double_para(6);
	Real Ab = global.force.double_para(7);
	Real Bb = global.force.double_para(8);				
	Real Cb = global.force.double_para(9);


	if (global.force.update_time + global.force.delta_t < global.time.now) {
		global.force.update_time = global.time.now;

		if (U.force_switch) {
			Kinetic_random_helical_forcing_Stepanov_Titov_assign_mod(U, W, inner_radius, outer_radius, Au, Bu);
		}
	}
	if (W.force_switch) {
		Magnetic_random_helical_forcing_Stepanov_Titov_assign_mod(U, W, inner_radiusM, outer_radiusM, Ab, Bb);
	}
}

void FORCE::Kinetic_random_helical_forcing_Stepanov_Titov_assign_mod(FluidVF& U, FluidVF& W, Real inner_radius, Real outer_radius, Real Au, Real Bu){
	int kx_max, ky_max, kz_max, kx_min, ky_min, kz_min;
	int nf;
// service part
	kx_max = (int) ceil(outer_radius/kfactor[1]);
	
	if (Ny > 1)
		ky_max = (int) ceil(outer_radius/kfactor[2]);
	else
		ky_max = 0;	
	
	kz_max = (int) ceil(outer_radius/kfactor[3]);
	
	kx_min = ky_min = kz_min = 0;
	
	if (basis_type == "FFF" || basis_type == "FFFW")
		kx_min = -kx_max;
	
	if ((basis_type == "FFF" || basis_type == "FFFW") || (basis_type == "SFF"))
		ky_min = -ky_max;

	nf = universal->Get_number_modes_in_shell(inner_radius, outer_radius);
// Exclude modes with 0 components
	int lx, ly, lz;
	Real Kmag;
	for (int kx = kx_min; kx <= kx_max; kx++)
	for (int ky = ky_min; ky <= ky_max; ky++)  
	for (int kz = 0; kz <= kz_max; kz++) {
		if ((kx==0) || (ky==0) || (kz==0)) 
		if (universal->Probe_in_me(kx,ky,kz))  {
			Kmag = sqrt(pow2(kx)+pow2(ky)+pow2(kz));
			if ((Kmag > inner_radius) && (Kmag <= outer_radius)) {
				nf-=(kz==0?1:2);
			}
		}
	}
// end of service part	
	TinyVector<Real, 3> k, e1,KxE;
	TinyVector<Complex, 3> f, localV, localW;
	Real norm, k2;
	Real alpha;

	alpha  = sqrt( 2 * Au / (global.force.delta_t * nf) );

	for (int kx = kx_min; kx <= kx_max; kx++)
	for (int ky = ky_min; ky <= ky_max; ky++)  
	for (int kz = 0; kz <= kz_max; kz++) {
		if ((kx!=0) && (ky!=0) && (kz!=0))
		if (universal->Probe_in_me(kx,ky,kz))  {
			lx = universal->Get_lx(kx);
			ly = universal->Get_ly(ky);
			lz = universal->Get_lz(kz);
				
			Kmag = universal->Kmagnitude(lx, ly, lz);
			if ((Kmag > inner_radius) && (Kmag <= outer_radius)) {
					k = kx, ky, kz;
					e1 = (2*SPECrand.random()-1), (2*SPECrand.random()-1), (2*SPECrand.random()-1);
					e1/=mynorm(e1);

					KxE = mycross(k,e1);

					f = alpha * (Bu * mycross(k,KxE) - I*Kmag*KxE);

					k2 = pow(Kmag,2.0);

					if (Bu == 1.0) {
						norm = sqrt(2.0)*k2*sqrt(1-pow(mydot(k,e1),2)/k2);
					}
					else if (Bu == 0.0) {
						norm = 1.0*k2*sqrt(1-pow(mydot(k,e1),2)/k2);
					}

					f/=norm;


					// localV = universal->Get_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3);
					// U.cvf.Compute_Modal_vorticity(lx, ly, lz, localW);
					// cout<<"Force debug:"<<endl
					// 	<<"k = " << k << endl
					// 	<<"e1 = "<< e1 << endl
					// 	<<"KxE = "<< KxE << endl
					// 	<<"f = "<< f<<endl
					// 	<<"Au = "<< Au <<", Bu = " << Bu << ", nf = "<< nf <<endl
					// 	<<"norm(f) = " << mynorm(f) <<", norm(v) = " << mynorm(localV) <<endl
					// 	<<"f.u* = "<<mydot(f, localV) <<endl
					//  	<<"f.w* = "<<mydot(f, localW) <<endl
					// 	<<"*********End! Thank you for cooperation...."<<endl;

					universal->Assign_local_spectral_field(lx, ly, lz, U.Force1, U.Force2, U.Force3, f);
				}
			}
		}

}

void FORCE::Magnetic_random_helical_forcing_Stepanov_Titov_assign_mod(FluidVF& U, FluidVF& W, Real inner_radius, Real outer_radius, Real Ab, Real Bb){
	int kx_max, ky_max, kz_max, kx_min, ky_min, kz_min;
	int nf;
// service part
	kx_max = (int) ceil(outer_radius/kfactor[1]);
	
	if (Ny > 1)
		ky_max = (int) ceil(outer_radius/kfactor[2]);
	else
		ky_max = 0;	
	
	kz_max = (int) ceil(outer_radius/kfactor[3]);
	
	kx_min = ky_min = kz_min = 0;
	
	if (basis_type == "FFF" || basis_type == "FFFW")
		kx_min = -kx_max;
	
	if ((basis_type == "FFF" || basis_type == "FFFW") || (basis_type == "SFF"))
		ky_min = -ky_max;

	nf = universal->Get_number_modes_in_shell(inner_radius, outer_radius);
// Exclude modes with 0 components
	int lx, ly, lz;
	Real Kmag;
	for (int kx = kx_min; kx <= kx_max; kx++)
	for (int ky = ky_min; ky <= ky_max; ky++)  
	for (int kz = 0; kz <= kz_max; kz++) {
		if ((kx==0) || (ky==0) || (kz==0)) 
		if (universal->Probe_in_me(kx,ky,kz))  {
			Kmag = sqrt(pow2(kx)+pow2(ky)+pow2(kz));
			if ((Kmag > inner_radius) && (Kmag <= outer_radius)) {
				nf-=(kz==0?1:2);
			}
		}
	}
// end of service part	
	Real Dh, alpha, aScal, bScal, tmp, epsilon;
	TinyVector<Real, 3> k, e1, e2, KxE1, KxE2;
	TinyVector<Complex, 3> f, localB, localW, pVect, kcomp;

	Real relative_helicity_corr;
	Real Nb;

	// alpha = Ab / nf;
	// alpha /= (inner_radius + outer_radius) / 2.0;
	alpha = Ab;
	epsilon = Bb;

	for (int kx = kx_min; kx <= kx_max; kx++)
	for (int ky = ky_min; ky <= ky_max; ky++)  
	for (int kz = 0; kz <= kz_max; kz++) {
		if ((kx!=0) && (ky!=0) && (kz!=0))

		if (universal->Probe_in_me(kx,ky,kz))  {
			lx = universal->Get_lx(kx);
			ly = universal->Get_ly(ky);
			lz = universal->Get_lz(kz);
				
			Kmag = universal->Kmagnitude(lx, ly, lz);
			if ((Kmag > inner_radius) && (Kmag <= outer_radius)) {
					k = kx, ky, kz;
					kcomp = kx, ky, kz;

					e1 = (2*SPECrand.random()-1), (2*SPECrand.random()-1), (2*SPECrand.random()-1);
					e1/=mynorm(e1);
					
					e2 = (2*SPECrand.random()-1), (2*SPECrand.random()-1), (2*SPECrand.random()-1);
					e2/=mynorm(e2);

					localB = universal->Get_local_spectral_field(lx, ly, lz, W.cvf.V1, W.cvf.V2, W.cvf.V3);
					W.cvf.Compute_Modal_vorticity(lx, ly, lz, localW);
					localW/=(Kmag==0?1:Kmag*Kmag);
					
					KxE1 = mycross(k,e1);
					KxE2 = mycross(k,e2);

					aScal = mydot(imag(localB),KxE2);
					bScal = mydot(real(localB),KxE1);

					pVect = -aScal*e1 + I*bScal*e2;

					f = mycross(kcomp, pVect);

					// Dh = aScal*mydot(imag(localB),e1) + bScal*mydot(real(localB),e2);

					// Nb = mydot(f,mycross(-I*kcomp,localB)) / (Vsqr(kcomp)*Vsqr(localB));

					//relative_helicity_corr =  Kmag * mydot(localB, localW) / (2 * W.cvf.Modal_energy(lx, ly, lz));
					relative_helicity_corr =  W.cvf.Modal_energy(lx, ly, lz) / Kmag;

					// cout<<"B = "<<localB<<endl
					// 	<<"A = "<<localW<<endl
					// 	<<"e1 = "<<e1<<endl
					// 	<<"e2 = "<<e2<<endl
					// 	<<"k = "<<k<<endl;

					// cout<<"k Hb/Eb = "<< relative_helicity_corr<<endl
					// 	<<"Max Hb = "<< W.cvf.Modal_energy(lx, ly, lz) / Kmag<<endl
					// 	<<"f.b* = "<<mydot(f, localB) <<endl
					//  	<<"f.a* = "<<mydot(f, localW) <<endl
					//  	<<endl
					// 	<<"F1.b* = "<<mydot(alpha*f/Dh, localB) <<endl
					//  	<<"F1.a* = "<<mydot(alpha*f/Dh, localW) <<endl
					//  	<<"norm D = "<<Dh<<endl
					// 	<<"F2.b* = "<<mydot(alpha*f*Nb, localB) <<endl
					//  	<<"F2.a* = "<<mydot(alpha*f*Nb, localW) <<endl
					//  	<<"norm Nb = "<<Nb<<endl
					//   	<<"----------------------------------------------\n\n";

					Dh = mydot(f,localW);
					//epsilon = 0.1;
					tmp = pow(mynorm(k),3.0) * pow(mynorm(localW),2.0)* epsilon;

					Nb = (fabs(Dh) >= tmp?1.0/Dh:0);

					f*=alpha*Nb*relative_helicity_corr;

					// cout<<"alpha = "<<alpha<<endl
					// 	<<"nf = "<<nf<<endl
					// 	<<"Dh = "<<Dh<<endl
					// 	<<"tmp = "<<tmp<<endl
					// 	<<"Nb = "<<Nb<<endl
					// 	<<"f.b* = "<<mydot(f, localB) <<endl
					// 	<<"f.a* = "<<mydot(f, localW) <<endl
					// 	<<"A.B* = "<<mydot(localB, localW) <<endl
					// 	<<"HbRel = "<<total_helicity<<endl
					// 	<<"---------------------------------"<<endl;


					// if (mydot(f,localB) >= pow(10,-8)) {
					// 	cout<<"*****************************************************"<<endl
					// 		<<"***              ALARM!!!                         ***"<<endl
					// 		<<"*****************************************************"<<endl
					// 		<<"k = "<<k<<endl
					// 		<<"F.B* = "<<mydot(f,localB)<<endl
					// 		<<"F.A* = "<<mydot(f,localW)<<endl
					// 		<<"*****************************************************"<<endl
					// 		<<"***              END OF ALARM                     ***"<<endl
					// 		<<"*****************************************************"<<endl;
					// }

					universal->Assign_local_spectral_field(lx, ly, lz, W.Force1, W.Force2, W.Force3, f);
				}
			}
		}

}

//*****************************************************************
//Kinetic and magnetic random helical forces (from mahyd 2016 original)
void FORCE::Compute_random_helical_forcing_Stepanov_Titov(FluidVF& U, FluidVF& W) {
	Real inner_radius = global.force.double_para(0);
	Real outer_radius = global.force.double_para(1);
	Real Au = global.force.double_para(2);
	Real Bu = global.force.double_para(3);				
	Real Cu = global.force.double_para(4);

	Real inner_radiusM = global.force.double_para(5);
	Real outer_radiusM = global.force.double_para(6);
	Real Ab = global.force.double_para(7);
	Real Bb = global.force.double_para(8);				
	Real Cb = global.force.double_para(9);

	if (U.force_switch) {
		Kinetic_random_helical_forcing_Stepanov_Titov_assign(U, W, inner_radius, outer_radius, Au, Bu);
	}
	if (W.force_switch) {
		Magnetic_random_helical_forcing_Stepanov_Titov_assign(U, W, inner_radiusM, outer_radiusM, Ab, Bb);
	}
}

void FORCE::Kinetic_random_helical_forcing_Stepanov_Titov_assign(FluidVF& U, FluidVF& W, Real inner_radius, Real outer_radius, Real Au, Real Bu){
	
	// global.force.force_U_lock = true;
	
	int kx_max, ky_max, kz_max, kx_min, ky_min, kz_min;
// service part
	kx_max = (int) ceil(outer_radius/kfactor[1]);
	
	if (Ny > 1)
		ky_max = (int) ceil(outer_radius/kfactor[2]);
	else
		ky_max = 0;	
	
	kz_max = (int) ceil(outer_radius/kfactor[3]);
	
	kx_min = ky_min = kz_min = 0;
	
	if (basis_type == "FFF" || basis_type == "FFFW")
		kx_min = -kx_max;
	
	if ((basis_type == "FFF" || basis_type == "FFFW") || (basis_type == "SFF"))
		ky_min = -ky_max;

	int lx, ly, lz;
	Real Kmag;

	int kx, ky, kz;
	int kk[3];

	TinyVector<Real, 3> k, e1, KxE;
	TinyVector<Complex, 3> f, localV;
	Real D, norm, epsilon;	

	if (!global.force.force_U_lock) {
		if (my_id == master_id) do {
			kk[0] = (int)((2*SPECrand.random()-1) * kx_max);
			kk[1] = (int)((2*SPECrand.random()-1) * ky_max);
			kk[2] = (int)(SPECrand.random() * kz_max);

			if ((kk[0]==0) || (kk[1]==0) || (kk[2]==0)) continue;
			Kmag = sqrt(kk[0]*kk[0] + kk[1]*kk[1] +kk[2]*kk[2]);

		} while (fabs(Kmag - outer_radius) > 0.5);

		MPI_Bcast(kk, 3, MPI_INT, master_id, MPI_COMM_WORLD);
		kx = kk[0]; ky = kk[1]; kz = kk[2];

		global.force.Ku = kx, ky, kz;
		global.force.force_U_lock = true;

//		cout<<"generated new K"<<endl;

	}		
		
	epsilon = pow(10.0, -2.0);
	k = global.force.Ku * 1.0;
	// if (my_id == master_id) cout<< "force started!"<<endl;

	for (kx = kx_min; kx <= kx_max; kx++)
	for (ky = ky_min; ky <= ky_max; ky++)  
	for (kz = 0; kz <= kz_max; kz++) {
		if (universal->Probe_in_me(kx,ky,kz))  {
			lx = universal->Get_lx(kx);
			ly = universal->Get_ly(ky);
			lz = universal->Get_lz(kz);
			if ((kx==k(0)) && (ky==k(1)) && (kz==k(2))) {
				Kmag = universal->Kmagnitude(lx, ly, lz);
				localV = universal->Get_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3);	

				e1 = (2*SPECrand.random()-1), (2*SPECrand.random()-1), (2*SPECrand.random()-1);
				e1/=mynorm(e1);

				KxE = mycross(k,e1);

				f = Bu * mycross(k,KxE) - I*Kmag*KxE;

				D = mydot(f,localV);

				norm = Au / (D<0?D-epsilon:D+epsilon);

				// cout<<"DEBUG: "<<endl
				// 	<<"k = "<< k << endl
				// 	<<"eps = " <<epsilon << endl
				// 	<<"f = "<< f << endl
				// 	<<"localV = "<< localV << endl
				// 	<<"D = "<< D << endl
				// 	<<"norm = "<< norm << endl
				// 	<<"--------------------------------------------\n";
				// 	;

				f*=norm;

				// cout<<"my_id is "<<my_id<<endl;
				// cout<<"k = "<< k << endl;
				// cout<<"norm = "<<norm<<endl;
				// cout<<"f.U = "<< mydot(f,localV) << endl <<endl
				//   	;
				// cout<<"Kinetic force finished"<<endl<<"------------------------------------------\n";
			}
			else {
				f = 0.0, 0.0, 0.0;
			}
	
			universal->Assign_local_spectral_field(lx, ly, lz, U.Force1, U.Force2, U.Force3, f);
		}
	}


}

void FORCE::Magnetic_random_helical_forcing_Stepanov_Titov_assign(FluidVF& U, FluidVF& W, Real inner_radius, Real outer_radius, Real Ab, Real Bb){
	int kx_max, ky_max, kz_max, kx_min, ky_min, kz_min;
	int nf;
// service part
	kx_max = (int) ceil(outer_radius/kfactor[1]);
	
	if (Ny > 1)
		ky_max = (int) ceil(outer_radius/kfactor[2]);
	else
		ky_max = 0;	
	
	kz_max = (int) ceil(outer_radius/kfactor[3]);
	
	kx_min = ky_min = kz_min = 0;
	
	if (basis_type == "FFF" || basis_type == "FFFW")
		kx_min = -kx_max;
	
	if ((basis_type == "FFF" || basis_type == "FFFW") || (basis_type == "SFF"))
		ky_min = -ky_max;

	nf = universal->Get_number_modes_in_shell(inner_radius, outer_radius);

	int lx, ly, lz;
	int kx, ky, kz;
	int kk[3];
	Real Kmag;

	if (!global.force.force_B_lock) {
		if (my_id == master_id) do {
			kk[0] = (int)((2*SPECrand.random()-1) * kx_max);
			kk[1] = (int)((2*SPECrand.random()-1) * ky_max);
			kk[2] = (int)(SPECrand.random() * kz_max);

			if ((kk[0]==0) || (kk[1]==0) || (kk[2]==0)) continue;
			Kmag = sqrt(kk[0]*kk[0] + kk[1]*kk[1] +kk[2]*kk[2]);

		} while (fabs(Kmag - outer_radius) > 0.5);

		MPI_Bcast(kk, 3, MPI_INT, master_id, MPI_COMM_WORLD);
		kx = kk[0]; ky = kk[1]; kz = kk[2];

		global.force.Kb = kx, ky, kz;
		global.force.force_B_lock = true;
	}

	Real Dh, alpha, aScal, bScal, tmp, epsilon;
	TinyVector<Real, 3> k, e1, e2, KxE1, KxE2;
	TinyVector<Complex, 3> f, localB, localW, pVect, kcomp;

	k = global.force.Kb * 1.0;

	Real relative_helicity_corr;
	Real Nb;

	alpha = Ab;
	epsilon = Bb;

	for (kx = kx_min; kx <= kx_max; kx++)
	for (ky = ky_min; ky <= ky_max; ky++)  
	for (kz = 0; kz <= kz_max; kz++) {
		if (universal->Probe_in_me(kx,ky,kz))  {
			lx = universal->Get_lx(kx);
			ly = universal->Get_ly(ky);
			lz = universal->Get_lz(kz);
				
			Kmag = universal->Kmagnitude(lx, ly, lz);
				if ((kx==k(0)) && (ky==k(1)) && (kz==k(2))) {
					// k = kx, ky, kz;
					kcomp = kx, ky, kz;

					e1 = (2*SPECrand.random()-1), (2*SPECrand.random()-1), (2*SPECrand.random()-1);
					e1/=mynorm(e1);
					
					e2 = (2*SPECrand.random()-1), (2*SPECrand.random()-1), (2*SPECrand.random()-1);
					e2/=mynorm(e2);

					localB = universal->Get_local_spectral_field(lx, ly, lz, W.cvf.V1, W.cvf.V2, W.cvf.V3);
					W.cvf.Compute_Modal_vorticity(lx, ly, lz, localW);
					//localW/=(Kmag==0?1:Kmag*Kmag);
					
					KxE1 = mycross(k,e1);
					KxE2 = mycross(k,e2);

					aScal = mydot(imag(localB),KxE2);
					bScal = mydot(real(localB),KxE1);

					pVect = -aScal*e1 + I*bScal*e2;

					f = mycross(kcomp, pVect);

					relative_helicity_corr =  W.cvf.Modal_energy(lx, ly, lz) / Kmag;

					Dh = mydot(f,localW);
					//epsilon = 0.1;
					tmp = pow(mynorm(k),3.0) * pow(mynorm(localW),2.0)* epsilon;

					Nb = (fabs(Dh) >= tmp?1.0/Dh:0);

					// Nb = (fabs(Nb)>pow(10.0,6)?1e6:Nb);

					// cout<<"f = "<<f<<endl
					// 	<<"Hb = "<<relative_helicity_corr<<endl
					// 	<<"Dh = "<<Dh<<endl
					// 	<<"tmp = "<<tmp<<endl
					// 	<<"Nb = "<<Nb<<endl
					// 	<<"alpha = "<<alpha<<endl
					// 	;

					f*=alpha*Nb*relative_helicity_corr;
					//f*=alpha*Nb;

					// cout<<"f = "<<f<<endl
					// 	<<"f.b* = "<<mydot(f, localB) <<endl
					// 	<<"f.a* = "<<mydot(f, localW) <<endl
					// 	<<"A.B* = "<<mydot(localB, localW) <<endl
					// 	<<"---------------------------------"<<endl;
					if (mydot(f,localB) >= pow(10,-8)) {
						cout<<"*****************************************************"<<endl
							<<"***              ALARM!!!                         ***"<<endl
							<<"*****************************************************"<<endl
							<<"k = "<<k<<endl
							<<"F.B* = "<<mydot(f,localB)<<endl
							<<"F.A* = "<<mydot(f,localW)<<endl
							<<"*****************************************************"<<endl
							<<"***              END OF ALARM                     ***"<<endl
							<<"*****************************************************"<<endl;
					}
				} else {
					f = 0.0, 0.0, 0.0;
				}
				universal->Assign_local_spectral_field(lx, ly, lz, W.Force1, W.Force2, W.Force3, f);
			}
		}

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////NEW RANDOM FORCING WRAPPER AND ANOTHER CODE///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void FORCE::Main_random_direction_forcing(FluidVF& U, FluidVF& W){
	Real inner_radius = global.force.double_para(0);
	Real outer_radius = global.force.double_para(1);
	Real Au = global.force.double_para(2);
	Real Bu = global.force.double_para(3);				
	Real Cu = global.force.double_para(4);

	Real inner_radiusM = global.force.double_para(5);
	Real outer_radiusM = global.force.double_para(6);
	Real Ab = global.force.double_para(7);
	Real Bb = global.force.double_para(8);				
	Real Cb = global.force.double_para(9);
// Later I will add switch or something like this


	if ((!global.force.force_U_lock) && (!global.force.force_B_lock)) {
		if (U.force_switch) {
			global.force.force_U_lock = global.force.force_B_lock = true;
			Compute_simple_kinetic_random_forcing_addition_shell_distr(U, W, inner_radius, outer_radius, Au);
		}
	}
}


void FORCE::Compute_simple_kinetic_random_forcing_addition_shell_distr(FluidVF& U, FluidVF& W, Real inner_radius, Real outer_radius, Real Au){
	int kx_max, ky_max, kz_max, kx_min, ky_min, kz_min;
	int nf;
// service part
	kx_max = (int) ceil(outer_radius/kfactor[1]);
	
	if (Ny > 1)
		ky_max = (int) ceil(outer_radius/kfactor[2]);
	else
		ky_max = 0;	
	
	kz_max = (int) ceil(outer_radius/kfactor[3]);
	
	kx_min = ky_min = kz_min = 0;
	
	if (basis_type == "FFF" || basis_type == "FFFW")
		kx_min = -kx_max;
	
	if ((basis_type == "FFF" || basis_type == "FFFW") || (basis_type == "SFF"))
		ky_min = -ky_max;

//	nf = universal->Get_number_modes_in_shell(inner_radius, outer_radius);
// Calculate modes without 0 components
	nf=0;
	int lx, ly, lz;
	Real Kmag;
	for (int kx = kx_min; kx <= kx_max; kx++)
	for (int ky = ky_min; ky <= ky_max; ky++)  
	for (int kz = 0; kz <= kz_max; kz++) {
		Kmag = sqrt(pow2(kx)+pow2(ky)+pow2(kz));
		if ((Kmag > inner_radius) && (Kmag <= outer_radius)) {
			nf+=(kz==0?1:2);
		}
	}
// end of service part	
// for output purposes and compute eps_j_i rates
	Real in[9], out[9];
// end

	TinyVector<Complex,3> eu, k, e1, fu, vectorU, vectorB, vectorW;
	Real alpha, k2;
	Real eps_Eu, eps_Hu, eps_Hc;

	eps_Eu = eps_Hu = eps_Hc =0;

	alpha  = sqrt( 2.0 * Au / (global.time.dt * nf) );

	TinyVector<Complex,3> temp;

	for (int kx = kx_min; kx <= kx_max; kx++)
	for (int ky = ky_min; ky <= ky_max; ky++)  
	for (int kz = 0; kz <= kz_max; kz++) {
		// if ((kx!=0) && (ky!=0) && (kz!=0))
		if (universal->Probe_in_me(kx,ky,kz))  {
			lx = universal->Get_lx(kx);
			ly = universal->Get_ly(ky);
			lz = universal->Get_lz(kz);
				
			Kmag = universal->Kmagnitude(lx, ly, lz);
			if ((Kmag > inner_radius) && (Kmag <= outer_radius)) {
// preparing vectors
				e1 = (2*SPECrand.random()-1) + I*(2*SPECrand.random()-1),
					 (2*SPECrand.random()-1) + I*(2*SPECrand.random()-1),
					 (2*SPECrand.random()-1) + I*(2*SPECrand.random()-1);

				k = kx, ky, kz;

				eu = mycross(k,e1);
				eu /= mynorm(eu);
//kinetic force
				fu = alpha * eu;			

				fu*=global.time.dt;
//assign forces
				universal->Add_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3, fu);

				//for output and controlling purposes
				vectorU = universal->Get_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3);
				U.cvf.Compute_Modal_vorticity(lx, ly, lz, vectorW);
				vectorB = universal->Get_local_spectral_field(lx, ly, lz, W.cvf.V1, W.cvf.V2, W.cvf.V3);

				temp = alpha * eu * global.time.dt;

				eps_Eu+=mydot(vectorU,temp); 
				eps_Hu+=mydot(vectorW,temp);
				eps_Hc+=mydot(vectorB,temp);

				temp = alpha * eu;

				universal->Add_local_spectral_field(lx, ly, lz, U.Force1, U.Force2, U.Force3, temp);
			}
		}
	}
	global.force.empty_force = true;
	in[0] = eps_Eu; 
	in[1] = eps_Hu;
	in[2] = eps_Hc;

	MPI_Reduce(&in[0], &out[0], 3, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);

	global.force.A1 += out[0];
	global.force.B1 += out[1];
	global.force.C1 += out[2];
}


//***********************  End of compute_force_random_titov_ed.cc *******************************
