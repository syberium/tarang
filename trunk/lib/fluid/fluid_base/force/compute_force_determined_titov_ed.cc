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

/*! \file  compute_force_determined_titov_ed.cc
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
 * @date Oct. 2018
 *
  */

#include "FORCE.h"

extern Uniform<Real> SPECrand;

//*********************************************************************************************

void FORCE::Compute_force_hydro_helical_decompose(FluidVF& U) {
	Real inner_radius = global.force.double_para(0);
	Real outer_radius = global.force.double_para(1);
	Real Au = global.force.double_para(2);
	Real Bu = global.force.double_para(3);			
	Real Cu = global.force.double_para(4);
	
	if ((!global.force.force_U_lock)) {
		global.force.force_U_lock = true;

			if (U.force_switch) {
				Compute_force_hydro_helical_decompose_addition(U, inner_radius, outer_radius, Au, Bu, Cu, 0);
			}
	}
}

void FORCE::Compute_force_hydro_helical_decompose_addition(FluidVF& U, Real inner_radius, Real outer_radius, Real Au, Real Bu, Real Cu, int MHD_switch) {
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

//	nf = universal->Get_number_modes_in_shell(inner_radius, outer_radius);
// Calculate modes without 0 components

//OLD STYLE

	int nf = 0;
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
// for output purposes and compute eps_i rates
	Real in[2], out[2];
	TinyVector<Complex,3> dU, w1;
	Real eps_E, eps_Hk;
// end

	TinyVector<Complex,3> k, z, hP, hM, U0, U1;
	Complex u0P, u0M, u1P, u1M;
	Real phiP, phiM;

	Real dE, dH, dEl, dHl;

	dE = Au * global.time.dt / nf;
	dH = Bu * global.time.dt / nf;

	eps_E = eps_Hk = 0.0;

	//some auxiliary variables
	TinyVector<Complex,3> ZxK;
	Real limit, protector;

	for (int kx = kx_min; kx <= kx_max; kx++)
	for (int ky = ky_min; ky <= ky_max; ky++)  
	for (int kz = 0; kz <= kz_max; kz++) {
		if (universal->Probe_in_me(kx,ky,kz))  {
			lx = universal->Get_lx(kx);
			ly = universal->Get_ly(ky);
			lz = universal->Get_lz(kz);
				
			Kmag = universal->Kmagnitude(lx, ly, lz);
			if ((Kmag > inner_radius) && (Kmag <= outer_radius)) {
// preparing vectors
				k = kx, ky, kz;
				do {
					z = 2*SPECrand.random()-1,
						2*SPECrand.random()-1,
						2*SPECrand.random()-1;
					ZxK = mycross(z,k);
				} while (mynorm(ZxK) < pow(10.0,-3));

				hP = mycross(ZxK,k)/( mynorm(k) * mynorm(ZxK) ) + I* ZxK/mynorm(ZxK);
				hM = mycross(ZxK,k)/( mynorm(k) * mynorm(ZxK) ) - I* ZxK/mynorm(ZxK);

				U0 = universal->Get_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3);

				u0P = (+I*mydot_complex(mycross(k,z),U0) - mynorm(k) * mydot_complex(z,U0) ) / (2*mynorm(ZxK));
				u0M = (-I*mydot_complex(mycross(k,z),U0) - mynorm(k) * mydot_complex(z,U0) ) / (2*mynorm(ZxK));

				phiP = atan2(imag(u0P),real(u0P));
				phiM = atan2(imag(u0M),real(u0M));

				if (MHD_switch == 0) {
					phiP += MY_PI/4.0*(2*SPECrand.random()-1);
					phiM += MY_PI/4.0*(2*SPECrand.random()-1);
				}

				// PROTECTOR part. Parameter limitation 
				// using only for helical part

				//dE, if negative, has to be \le of u0P and u0M simultaneously

				if (dE < 0) {
					// cout<< "was dE = "<< dE <<endl;
					limit = - min(2.0 * Vsqr(u0P), 2.0 * Vsqr(u0M));
					dEl = dE < limit? limit : dE;
					// cout<< "now dE = "<< dEl <<endl<<endl;
				} else {
					dEl = dE;
				}
				//dH
				// if (dH < 0) {
				// 	// cout<< "was dH = "<< dH <<endl;
				// 	limit = - (dEl + 2.0 * Vsqr(u0P));
				// 	if (MHD_switch == 0) limit*=mynorm(k); else limit/=mynorm(k);
				// 	dHl = dH < limit ? limit : dH;
				// 	// cout<< "now dH = "<< dHl <<endl<<endl;
				// } else {
				// 	// cout<< "was dH = "<< dH <<endl;
				// 	limit = (dEl + 2.0 * Vsqr(u0M));
				// 	if (MHD_switch == 0) limit*=mynorm(k); else limit/=mynorm(k);
				// 	dHl = dH > limit ? limit : dH;
				// 	// cout<< "now dH = "<< dHl <<endl<<endl;
				// }

				if (dH < 0) {
					limit = - (dEl + 2.0 * Vsqr(u0P));
					if (MHD_switch == 0) limit*=mynorm(k); else limit/=mynorm(k);
					dHl = limit*Cu * global.time.dt;
				} else {
					limit = (dEl + 2.0 * Vsqr(u0M));
					if (MHD_switch == 0) limit*=mynorm(k); else limit/=mynorm(k);
					dHl = limit*Cu * global.time.dt;
				}

				if (MHD_switch == 0) {
					protector = (dEl + dHl/mynorm(k))/2.0 + Vsqr(u0P);
					u1P = exp(I*phiP)*sqrt( protector<0?0:protector );

					protector = (dEl - dHl/mynorm(k))/2.0 + Vsqr(u0M);
					u1M = exp(I*phiM)*sqrt( protector<0?0:protector );
				} else {
					protector = (dEl + dHl*mynorm(k))/2.0 + Vsqr(u0P);
					u1P = exp(I*phiP)*sqrt( protector<0?0:protector );

					protector = (dEl - dHl*mynorm(k))/2.0 + Vsqr(u0M);
					u1M = exp(I*phiM)*sqrt( protector<0?0:protector );
				}

				U1 = u1P*hP + u1M*hM;

				universal->Assign_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3, U1);
				// calculations for output 
				w1 = I * mycross(k,U1);
				if (MHD_switch != 0) w1/=pow(mynorm(k),2.0);

				dU = U1 - U0;

				eps_E += mydot(U1, dU);
				eps_Hk += mydot(w1, dU);

				dU/=global.time.dt;
				universal->Assign_local_spectral_field(lx, ly, lz, U.Force1, U.Force2, U.Force3, dU);
				global.force.empty_force = true;
			}
		}
	}

	in[0] = 2.0 * eps_E; 
	in[1] = 2.0 * eps_Hk;

	MPI_Reduce(&in[0], &out[0], 2, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);

	if(MHD_switch == 0) {
		global.force.A1 += out[0];
		global.force.B1 += out[1];
	} else {
		global.force.A2 += out[0];
		global.force.B2 += out[1];
	}
}

void FORCE::Compute_force_hydro_helical_decompose_addition_shell_distribution(FluidVF& U, Real inner_radius, Real outer_radius, Real Au, Real Bu, Real Cu, int MHD_switch) {
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


// Calculate modes without 0 components

//NEW STYLE

	int shellsNum = outer_radius - inner_radius;
	Array<int,1> nfS(shellsNum);
	Array<Real,1> koef(shellsNum), dES(shellsNum), dHS(shellsNum);
	int lx, ly, lz;
	Real Kmag;

	nfS = 0;
	for (int kx = kx_min; kx <= kx_max; kx++)
	for (int ky = ky_min; ky <= ky_max; ky++)  
	for (int kz = 0; kz <= kz_max; kz++) {
		// if ((kx!=0) && (ky!=0) && (kz!=0)) {
			Kmag = sqrt(pow2(kx)+pow2(ky)+pow2(kz));
			if ((Kmag > inner_radius) && (Kmag <= outer_radius)) {
				nfS(floor(Kmag - inner_radius - 0.00001 ))+= (kz==0?1:2);
			}
		// }
	}

	// cout<< "nfS = " << nfS <<endl;

// end of service part	
// for output purposes and compute eps_i rates
	Real in[2], out[2];
	TinyVector<Complex,3> dU, w1;
	Real eps_E, eps_Hk;
// end

	TinyVector<Complex,3> k, z, hP, hM, U0, U1;
	Complex u0P, u0M, u1P, u1M;
	Real phiP, phiM;

	Real dE, dH, dEl, dHl, index, norm;

	dE = Au * global.time.dt;
	dH = Bu * global.time.dt;

	eps_E = eps_Hk = 0.0;

	for (int i=0; i<shellsNum; i++) koef(i) =2.0 * (shellsNum - 0.5 - i);
	norm = sum(koef);

	// cout<< "norm = " << norm << endl
	// 	<< "dE = " << dE << endl
	// 	<< "dH = " << dH << endl;

	dES = dE * koef / norm / nfS;
	dHS = dH * koef / norm / nfS;

	//some auxiliary variables
	TinyVector<Complex,3> ZxK;
	Real limit, limit_coefficient, protector;
	Real small_eps2 = 1e-06;

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
				index = floor(Kmag - inner_radius - 0.00001 );
				k = kx, ky, kz;
				do {
					z = 2*SPECrand.random()-1,
						2*SPECrand.random()-1,
						2*SPECrand.random()-1;
					ZxK = mycross(z,k);
				} while (mynorm(ZxK) < pow(10.0,-3));

				hP = mycross(ZxK,k)/( mynorm(k) * mynorm(ZxK) ) + I* ZxK/mynorm(ZxK);
				hM = mycross(ZxK,k)/( mynorm(k) * mynorm(ZxK) ) - I* ZxK/mynorm(ZxK);

				U0 = universal->Get_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3);

				u0P = (+I*mydot_complex(mycross(k,z),U0) - mynorm(k) * mydot_complex(z,U0) ) / (2*mynorm(ZxK));
				u0M = (-I*mydot_complex(mycross(k,z),U0) - mynorm(k) * mydot_complex(z,U0) ) / (2*mynorm(ZxK));

				phiP = atan2(imag(u0P),real(u0P));
				phiM = atan2(imag(u0M),real(u0M));

				// PROTECTOR part. Parameter limitation 
				// using only for helical part
				// if (MHD_switch != 0) limit_coefficient = 0.10; else limit_coefficient = 0.20;
				limit_coefficient = 1.0;
				//dE, if negative, has to be \le of u0P and u0M simultaneously

				if (dES(index) < 0) {
					limit = - min(2.0 * Vsqr(u0P), 2.0 * Vsqr(u0M));
					dEl = dES(index) < limit? limit : dES(index);
				} else {
					dEl = dES(index);
				}
				//dH
				if (dHS(index) < 0) {
					limit = - (dEl + 2.0 * Vsqr(u0P));
					if (MHD_switch == 0) limit*=mynorm(k); else limit/=mynorm(k);
					dHl = dHS(index) < limit ? limit : dHS(index);
				} else {
					limit = (dEl + 2.0 * Vsqr(u0M));
					if (MHD_switch == 0) limit*=mynorm(k); else limit/=mynorm(k);
					dHl = dHS(index) > limit ? limit : dHS(index);
				}

				//Clipping by TAS! Bless God this good guy :-)
				if (Vsqr(u0M)<small_eps2) { dHl = 0.0; }

				if (MHD_switch == 0) {
					protector = (dEl + dHl/mynorm(k))/2.0 + Vsqr(u0P);
					u1P = exp(I*phiP)*sqrt( protector<0?0:protector );

					protector = (dEl - dHl/mynorm(k))/2.0 + Vsqr(u0M);
					u1M = exp(I*phiM)*sqrt( protector<0?0:protector );
				} else {
					protector = (dEl + dHl*mynorm(k))/2.0 + Vsqr(u0P);
					u1P = exp(I*phiP)*sqrt( protector<0?0:protector );

					protector = (dEl - dHl*mynorm(k))/2.0 + Vsqr(u0M);
					u1M = exp(I*phiM)*sqrt( protector<0?0:protector );
				}

				U1 = u1P*hP + u1M*hM;

				universal->Assign_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3, U1);
				// calculations for output 
				w1 = I * mycross(k,U1);
				if (MHD_switch != 0) w1/=pow(mynorm(k),2.0);

				dU = U1 - U0;

				eps_E += mydot(U1, dU);
				eps_Hk += mydot(w1, dU);

				dU/=global.time.dt;
				universal->Add_local_spectral_field(lx, ly, lz, U.Force1, U.Force2, U.Force3, dU);
				global.force.empty_force = true;
			}
		}
	}

	in[0] = 2.0 * eps_E; 
	in[1] = 2.0 * eps_Hk;

	MPI_Reduce(&in[0], &out[0], 2, MPI_Real, MPI_SUM, master_id, MPI_COMM_WORLD);

	if(MHD_switch == 0) {
		global.force.A1 += out[0];
		global.force.B1 += out[1];
	} else {
		global.force.A2 += out[0];
		global.force.B2 += out[1];
	}
}


//*****************************************************************************************************

void FORCE::Compute_force_MHD_helical_decompose(FluidVF& U,FluidVF& W) {

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

	if (!global.force.force_U_lock) {
		global.force.force_U_lock = true;

			if (U.force_switch) {
				Compute_force_hydro_helical_decompose_addition(U, inner_radius, outer_radius, Au, Bu, Cu, 0);
			}
	}

	if (!global.force.force_B_lock) {
		global.force.force_B_lock = true;

			if (W.force_switch) {
				Compute_force_hydro_helical_decompose_addition(W, inner_radiusM, outer_radiusM, Ab, Bb, Cb, 1);
				if(Ab == 0) {
					// Compute_force_hydro_helical_decompose_addition(W, 1.0, 2.0, 0, -Bb, Cb, 1);
				}
			}
	}

}


void FORCE::Compute_force_MHD_helical_decompose_shell_distribution(FluidVF& U,FluidVF& W) {

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

	if (!global.force.force_U_lock) {
		global.force.force_U_lock = true;

			if (U.force_switch) {
				Compute_force_hydro_helical_decompose_addition_shell_distribution(U, inner_radius, outer_radius, Au, Bu, Cu, 0);
			}
	}

	if (!global.force.force_B_lock) {
		global.force.force_B_lock = true;

			if (W.force_switch) {
				Compute_force_hydro_helical_decompose_addition_shell_distribution(W, inner_radiusM, outer_radiusM, Ab, Bb, Cb, 1);
				if(Ab == 0) {
					Compute_force_hydro_helical_decompose_addition_shell_distribution(W, 1.0, 2.0, 0, -Bb, Cb, 1);
				}
			}
	}

}


void FORCE::Compute_force_MHD_helical_decompose_shell_distribution_differ_scales(FluidVF& U,FluidVF& W) {

	Real inner_radiusE = global.force.double_para(0);
	Real outer_radiusE = global.force.double_para(1);
	Real Au = global.force.double_para(2);

	Real inner_radiusH = global.force.double_para(3);
	Real outer_radiusH = global.force.double_para(4);
	Real Bu = global.force.double_para(5);

	Real inner_radiusEM = global.force.double_para(6);
	Real outer_radiusEM = global.force.double_para(7);
	Real Ab = global.force.double_para(8);

	Real inner_radiusHM = global.force.double_para(9);
	Real outer_radiusHM = global.force.double_para(10);
	Real Bb = global.force.double_para(11);


	if (!global.force.force_U_lock) {
		global.force.force_U_lock = true;
			if (U.force_switch) {
				Compute_force_hydro_helical_decompose_addition_shell_distribution(U, inner_radiusE, outer_radiusE, Au, 0, 1, 0);
				Compute_force_hydro_helical_decompose_addition_shell_distribution(U, inner_radiusH, outer_radiusH, 0, Bu, 1, 0);
			}
	}

	if (!global.force.force_B_lock) {
		global.force.force_B_lock = true;
			if (W.force_switch) {
				Compute_force_hydro_helical_decompose_addition_shell_distribution(W, inner_radiusEM, outer_radiusEM, Ab, 0, 1, 1);
				Compute_force_hydro_helical_decompose_addition_shell_distribution(W, inner_radiusHM, outer_radiusHM, 0, Bb, 1, 1);
			}
	}

}


//*****************************************************************************************************

void FORCE::Compute_force_MHD_helical_decompose_with_crosshel(FluidVF& U,FluidVF& W) {

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

	if (!global.force.force_U_lock) {
		global.force.force_U_lock = true;

			if (U.force_switch) {
				Compute_force_helical_decompose_crosshel_addition(U, inner_radius, outer_radius, Au, Bu, Cu, 0, W);
			}
	}

	if (!global.force.force_B_lock) {
		global.force.force_B_lock = true;

			if (W.force_switch) {
				Compute_force_helical_decompose_crosshel_addition(W, inner_radiusM, outer_radiusM, Ab, Bb, Cb, 1, U);
			}
	}

}

void FORCE::Compute_force_helical_decompose_crosshel_addition(FluidVF& U, Real inner_radius, Real outer_radius, Real Au, Real Bu, Real Cu, int MHD_switch, FluidVF& B) {
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
// for output purposes and compute eps_i rates
	Real in[2], out[2];
	TinyVector<Complex,3> dU;
	Real eps_E, eps_Hk;
// end

	TinyVector<Complex,3> k, z, hP, hM, U0, U1, B0;
	Complex u0P, u0M, u1P, u1M, b0P, b0M, b1P, b1M;
	Real phiP, phiM;

	Real dE, dH;

	dE = Au * global.time.dt / nf;
	dH = Bu * global.time.dt / nf;

	eps_E = eps_Hk = 0.0;

	//some auxiliary variables
	TinyVector<Complex,3> ZxK;
	Real protector;

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
				k = kx, ky, kz;
				do {
					z = 2*SPECrand.random()-1,
						2*SPECrand.random()-1,
						2*SPECrand.random()-1;
					ZxK = mycross(z,k);
				} while (mynorm(ZxK) < pow(10.0,-3));

				hP = mycross(ZxK,k)/( mynorm(k) * mynorm(ZxK) ) + I* ZxK/mynorm(ZxK);
				hM = mycross(ZxK,k)/( mynorm(k) * mynorm(ZxK) ) - I* ZxK/mynorm(ZxK);

				U0 = universal->Get_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3);

				u0P = (+I*mydot_complex(mycross(k,z),U0) - mynorm(k) * mydot_complex(z,U0) ) / (2*mynorm(ZxK));
				u0M = (-I*mydot_complex(mycross(k,z),U0) - mynorm(k) * mydot_complex(z,U0) ) / (2*mynorm(ZxK));

				phiP = atan2(imag(u0P),real(u0P));
				phiM = atan2(imag(u0M),real(u0M));

				// NOW PROTECTOR WORKS ONLY FOR CASE dE>0, dH>0

				protector = MHD_switch==0?
						((dE - dH/mynorm(k))/2.0 + Vsqr(u0M)):
						((dE - dH*mynorm(k))/2.0 + Vsqr(u0M));
				
				if (protector < 0) {
					dH = MHD_switch==0?
						(dE + 2.0*Vsqr(u0M)) * mynorm(k) * 0.995:
						(dE + 2.0*Vsqr(u0M)) / mynorm(k) * 0.995;
					// protector = MHD_switch==0?
					// 		((dE - dH/mynorm(k))/2.0 + Vsqr(u0M)):
					// 		((dE - dH*mynorm(k))/2.0 + Vsqr(u0M));
					// if (MHD_switch !=0 ) cout<<"now = "<< protector<<endl<<endl;
				}

				if (MHD_switch == 0) {
					u1P = exp(I*phiP)*sqrt( (dE + dH/mynorm(k))/2.0 + Vsqr(u0P) );
					u1M = exp(I*phiM)*sqrt( (dE - dH/mynorm(k))/2.0 + Vsqr(u0M) );
				} else {
					u1P = exp(I*phiP)*sqrt( (dE + dH*mynorm(k))/2.0 + Vsqr(u0P) );
					u1M = exp(I*phiM)*sqrt( (dE - dH*mynorm(k))/2.0 + Vsqr(u0M) );
				}

				U1 = u1P*hP + u1M*hM;

				// cout
				//	<<"k = "<< k <<endl
				// 	<<"z = "<< z <<endl
				// 	<<"ZxK = "<< ZxK <<endl
				// 	<<"norm(ZxK) = "<< mynorm(ZxK) <<endl<<endl
					
					// <<"hP = "<< hP <<endl
					// <<"hM = "<< hM <<endl<<endl
					
					// <<"U0 = "<< U0 <<endl<<endl
					
					// <<"u0P = "<< u0P <<endl
					// <<"u0M = "<< u0M <<endl<<endl
					
					// <<"phiP = "<< phiP <<endl
					// <<"phiM = "<< phiM <<endl<<endl
					
					// <<"dE = " << dE <<endl
					// <<"dH = " << dH <<endl<<endl

					// <<"u1P = "<< u1P <<endl
					// <<"u1M = "<< u1M <<endl<<endl

					// <<"U1 = "<< U1 <<endl
					// <<"---------------------------------"<<endl<<endl;

				universal->Assign_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3, U1);
				// calculations for output 
				dU = U1 - U0;
			}
		}
	}
}




void FORCE::Compute_random_crosshelical_forcing_addition_plus_helical_decompose_kinetic_helicity(FluidVF& U,FluidVF& W){
	
	Real inner_radius = global.force.double_para(0);
	Real outer_radius = global.force.double_para(1);
	// Real Au = global.force.double_para(2);
	// Real Bu = global.force.double_para(3);			
	// Real Cu = global.force.double_para(4);			

	Real inner_radiusM = global.force.double_para(5);
	Real outer_radiusM = global.force.double_para(6);
	// Real Ab = global.force.double_para(7);
	// Real Bb = global.force.double_para(8);			
	// Real Cb = global.force.double_para(9);			

	if (U.force_switch && W.force_switch) {
		Kinetic_and_magnetic_random_crosshelical_forcing_addition(U, W, inner_radius, outer_radius, .5, .5, .25, .25);
	}

	if (!global.force.force_U_lock) {
		global.force.force_U_lock = true;
			if (U.force_switch) {
				Compute_force_hydro_helical_decompose_addition_shell_distribution(U, inner_radius, outer_radius, 0, 1, 1, 0);
			}
	}

	if (!global.force.force_B_lock) {
		global.force.force_B_lock = true;
	}
}


void FORCE::Compute_random_crosshelical_forcing_addition_plus_helical_decompose_magnetic_helicity(FluidVF& U,FluidVF& W){
	
	Real inner_radius = global.force.double_para(0);
	Real outer_radius = global.force.double_para(1);
	// Real Au = global.force.double_para(2);
	// Real Bu = global.force.double_para(3);			
	// Real Cu = global.force.double_para(4);			

	Real inner_radiusM = global.force.double_para(5);
	Real outer_radiusM = global.force.double_para(6);
	// Real Ab = global.force.double_para(7);
	// Real Bb = global.force.double_para(8);			
	// Real Cb = global.force.double_para(9);			

	if (U.force_switch && W.force_switch) {
		Kinetic_and_magnetic_random_crosshelical_forcing_addition(U, W, inner_radius, outer_radius, .5, .5, .25, .25);
	}

	if (!global.force.force_B_lock) {
		global.force.force_B_lock = true;
			if (W.force_switch) {
				Compute_force_hydro_helical_decompose_addition_shell_distribution(W, inner_radiusM, outer_radiusM, 0., .1, 1, 1);
			}
	}
}

void FORCE::Compute_random_crosshelical_forcing_addition_plus_helical_decompose_kinetic_magnetic_helicity(FluidVF& U,FluidVF& W){
	
	Real inner_radius = global.force.double_para(0);
	Real outer_radius = global.force.double_para(1);
	// Real Au = global.force.double_para(2);
	// Real Bu = global.force.double_para(3);			
	// Real Cu = global.force.double_para(4);			

	Real inner_radiusM = global.force.double_para(5);
	Real outer_radiusM = global.force.double_para(6);
	// Real Ab = global.force.double_para(7);
	// Real Bb = global.force.double_para(8);			
	// Real Cb = global.force.double_para(9);			

	if (U.force_switch && W.force_switch) {
		Kinetic_and_magnetic_random_crosshelical_forcing_addition(U, W, 1, 3, .5, .5, .25, .25);
	}

	if (!global.force.force_U_lock) {
		global.force.force_U_lock = true;
		if (U.force_switch) {
			Compute_force_hydro_helical_decompose_addition_shell_distribution(U, inner_radius, outer_radius, 0, 1, 1, 0);
		}
	}

	if (!global.force.force_B_lock) {
		global.force.force_B_lock = true;
		if (W.force_switch) {
			Compute_force_hydro_helical_decompose_addition_shell_distribution(W, inner_radiusM, outer_radiusM, 0, .1, 1, 1);
		}
	}
}