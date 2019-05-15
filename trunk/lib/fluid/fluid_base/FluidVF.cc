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

/*! \file  Cvf.cc
 * 
 * @brief  Class declaration of Cvf, a Complex Vector Field 
 *
 * @sa Cvf.h
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Sept 2008
 *
 * @bugs No known bugs
 */

#include "FluidVF.h"

					
/**********************************************************************************************

   Class constructor: Allocates for Vi, Ek etc.; Initializes the arrays.
   
**********************************************************************************************/

FluidVF::FluidVF(string field_name): PlainFluidVF(field_name) {
	nlin1.resize(shape_complex_array);
    nlin2.resize(shape_complex_array);
	nlin3.resize(shape_complex_array);
	
	
	// Memory allocation if force_switch is on
	if (force_switch) {
        Force1.resize(shape_complex_array);
        Force2.resize(shape_complex_array);
        Force3.resize(shape_complex_array);
	}
}

FluidVF::FluidVF
(
	Real dissipation_coefficient, 
	Real hyper_dissipation_coefficient, 
	int hyper_dissipation_exponent,
	bool force_switch,
	string field_name
 ) : PlainFluidVF(field_name)
{			

	this->dissipation_coefficient = dissipation_coefficient;
	this->hyper_dissipation_coefficient = hyper_dissipation_coefficient;
	this->hyper_dissipation_exponent = hyper_dissipation_exponent;
	
	if (hyper_dissipation_exponent == 0)
		hyper_dissipation_switch = false;
	
	this->force_switch = force_switch;
    
    nlin1.resize(shape_complex_array);
    nlin2.resize(shape_complex_array);
	nlin3.resize(shape_complex_array);
	
	
	// Memory allocation if force_switch is on
	if (force_switch) {
        Force1.resize(shape_complex_array);
        Force2.resize(shape_complex_array);
        Force3.resize(shape_complex_array);
	}
}
	

	//********************************************************************************************* 


void FluidVF::Compute_divergence_nlin(Array<Complex,3> div)
{
    Real total_abs_div;  // not reqd for this function.
    universal->Compute_divergence(nlin1, nlin2, nlin3, div, "nlin", total_abs_div, false);
    // PS: global.temp_array.X is used in Compute_divergence_field.. so
    // be careful.
    // Used in pressure computation.
}

//********************************************************************************************* 


/** @brief Compute divergence of VF \f$ \vec{V} \f$.
 *
 *	@note	Divergence is stored in CVF of class NLIN.  We use the same class to store
 *			pressure also.  Hence divergence function should be called in the beginning 
 *			or end of loop when *F is free.
 *
 *  @return  \f$ *F = \mathcal{F}(D_i V_i) \f$. 
 */
void FluidVF::Compute_divergence_field(Array<Complex,3> div, Real &total_abs_div, bool print_switch)
{
    universal->Compute_divergence(cvf.V1, cvf.V2, cvf.V3, div, "field", total_abs_div, print_switch);
    // PS: global.temp_array.X is used in Compute_divergence_field.. so
    // be careful.  We call div using global.temp_array.X2
}
	
//*********************************************************************************************	

Real	FluidVF::Get_mag_V0()
{
	Real modV0;
	
	if (my_id == master_id)
		modV0 = sqrt( pow2(real((cvf.V1)(0,0,0)))  +  pow2(real((cvf.V2)(0,0,0)))  +  pow2(real((cvf.V3)(0,0,0))) ) ;	

	MPI_Bcast( &modV0, 1, MPI_Real, master_id, MPI_COMM_WORLD);

	return modV0;														
}	
	
	
//*********************************************************************************************	

void FluidVF::Satisfy_strong_reality_condition_field()
{
	universal->Satisfy_strong_reality_condition_in_Array(cvf.V1);
	universal->Satisfy_strong_reality_condition_in_Array(cvf.V2);
	universal->Satisfy_strong_reality_condition_in_Array(cvf.V3);
}


//*********************************************************************************************	

void FluidVF::Satisfy_strong_reality_condition_force_field()
{

	universal->Satisfy_strong_reality_condition_in_Array(Force1);
	universal->Satisfy_strong_reality_condition_in_Array(Force2);
	universal->Satisfy_strong_reality_condition_in_Array(Force3);

}

//*********************************************************************************************

void FluidVF::Satisfy_weak_reality_condition_field()
{
	universal->Satisfy_weak_reality_condition_in_Array(cvf.V1);
	universal->Satisfy_weak_reality_condition_in_Array(cvf.V2);
	universal->Satisfy_weak_reality_condition_in_Array(cvf.V3);
}


//*********************************************************************************************

void FluidVF::Satisfy_weak_reality_condition_force_field()
{
    
	universal->Satisfy_weak_reality_condition_in_Array(Force1);
	universal->Satisfy_weak_reality_condition_in_Array(Force2);
	universal->Satisfy_weak_reality_condition_in_Array(Force3);
    
}
	

//******************************************************************************
void FluidVF::Test_reality_condition_field()
{
	universal->Test_reality_condition_in_Array(cvf.V1);
	universal->Test_reality_condition_in_Array(cvf.V2);
	universal->Test_reality_condition_in_Array(cvf.V3);
}
	

void FluidVF::Test_reality_condition_force_field()
{
    
	universal->Test_reality_condition_in_Array(Force1);
	universal->Test_reality_condition_in_Array(Force2);
	universal->Test_reality_condition_in_Array(Force3);
    
}

//******************************************************************************
void FluidVF::Dealias_force_field()
{
	if (global.program.alias_option == "DEALIAS") {
		universal->Dealias(Force1);
		universal->Dealias(Force2);
		universal->Dealias(Force3);
	}
}

//*********************************************************************************************	

void FluidVF::Zero_Prandtl_number_compute_temperature(FluidSF& T)
{
	T.csf.F = cvf.V1;
	T.csf.Divide_ksqr();
}


void FluidVF::Infinite_Prandtl_number_compute_velocity(FluidSF& T)
{
	//global.temp_array.X = T.csf.F;
	global.program.sincostr_switch = sincostr_switch_Vx; // cvf.V1 and T.csf.F has same basis
	universal->Xderiv(T.csf.F, global.temp_array.X);
	universal->Array_divide_ksqr(global.temp_array.X);  //X contains -Pressure
	
	// Compute Vy = -Dy(P)/(K^2); Vz = -Dz(P)/(K^2); Vx = (-Dx(P)+theta)/(K^2)
	// Adjusting the sign..
	if (N[2] > 1) {
		global.program.sincostr_switch = global.program.sincostr_switch_divergence;
		universal->Yderiv(global.temp_array.X, cvf.V2); 
		universal->Array_divide_ksqr(cvf.V2);
	}
	
	global.program.sincostr_switch = global.program.sincostr_switch_divergence;
	universal->Zderiv(global.temp_array.X, cvf.V3);
	universal->Array_divide_ksqr(cvf.V3);
	
	global.program.sincostr_switch = global.program.sincostr_switch_divergence;
	universal->Xderiv(global.temp_array.X, cvf.V1);   // V1 = Dx(P)
	cvf.V1 = cvf.V1 + T.csf.F;
	universal->Array_divide_ksqr(cvf.V1);
	
	if (global.PHYSICS.Uscaling == "ULARGE") {
		cvf.V1 = cvf.V1*sqrt(global.PHYSICS.Rayleigh);
		cvf.V2 = cvf.V2*sqrt(global.PHYSICS.Rayleigh);
		cvf.V3 = cvf.V3*sqrt(global.PHYSICS.Rayleigh);
	}
    
    else if (global.PHYSICS.Uscaling == "USMALL") {
		cvf.V1 = cvf.V1*global.PHYSICS.Rayleigh;
		cvf.V2 = cvf.V2*global.PHYSICS.Rayleigh;
		cvf.V3 = cvf.V3*global.PHYSICS.Rayleigh;
	}
    
}


/**********************************************************************************************
 Copy_field_to(W):	  W.V <- V
 Copy_field_from(W):   V <-W. V
***********************************************************************************************/

void FluidVF::Copy_field_to(PlainCVF& W)
{	
	W.V1 = cvf.V1;  W.V2 = cvf.V2;	W.V3 = cvf.V3;
} 

void FluidVF::Copy_field_to(CVF& W)
{	
	W.V1 = cvf.V1;  W.V2 = cvf.V2;	W.V3 = cvf.V3;
}

void FluidVF::Copy_field_from(PlainCVF& W)
{	
	cvf.V1 = W.V1;  cvf.V2 = W.V2;	cvf.V3 = W.V3;
}
  

void FluidVF::Copy_field_from(CVF& W)
{	
	cvf.V1 = W.V1;  cvf.V2 = W.V2;	cvf.V3 = W.V3;
}


/**********************************************************************************************
 
Function 1:  V = V + nlin*dt;  
Function 2: for CVF: Y = Y+ factor*dt*U.nlin; 
 
***********************************************************************************************/

void FluidVF::Add_nlin_factor_dt(Real factor)
{	
	if (abs(factor) > MYEPS2) {
		cvf.V1 += (factor*global.time.dt)*(nlin1);
		cvf.V2 += (factor*global.time.dt)*(nlin2);
		cvf.V3 += (factor*global.time.dt)*(nlin3);
	}
	
}


void FluidVF::Add_field_nlin_factor_dt(PlainCVF& Y, Real factor)
{	
	if (abs(factor) > MYEPS2) {
		Y.V1 += (factor*global.time.dt)* (nlin1);         
		Y.V2 += (factor*global.time.dt)* (nlin2);
		Y.V3 += (factor*global.time.dt)* (nlin3);
	}
}

//*********************************************************************************************

/** @brief Multiply vector field by \f$ \exp(-\nu K^2 dt) \f$ or 
 *			\f$ \exp(-\nu K^2 dt) + \exp(-\nu_h)  K^4 dt) \f$
 * 
 *  @param dt
 *
 * @return \f$ V(k) = V(k) * \exp(-\nu K^2 dt) \f$.
 * @return when hyper_dissipation_switch == 1, 
 *			\f$ V(k) = V(k) * \exp(-\nu K^2 dt) + \exp(-\nu_h)  K^4 dt) \f$
 */
void FluidVF::Mult_field_exp_ksqr_dt(Real a)
{
	if (abs(a) > MYEPS) {
		if (!hyper_dissipation_switch) {
			universal->Array_exp_ksqr(cvf.V1, -dissipation_coefficient*a*global.time.dt);
			universal->Array_exp_ksqr(cvf.V2, -dissipation_coefficient*a*global.time.dt);
			universal->Array_exp_ksqr(cvf.V3, -dissipation_coefficient*a*global.time.dt);
		}
	
		else {
			universal->Array_exp_ksqr(cvf.V1, -dissipation_coefficient*a*global.time.dt, -hyper_dissipation_coefficient*a*global.time.dt, hyper_dissipation_exponent);
			
			universal->Array_exp_ksqr(cvf.V2, -dissipation_coefficient*a*global.time.dt, -hyper_dissipation_coefficient*a*global.time.dt, hyper_dissipation_exponent);
			
			universal->Array_exp_ksqr(cvf.V3, -dissipation_coefficient*a*global.time.dt, -hyper_dissipation_coefficient*a*global.time.dt, hyper_dissipation_exponent);
		}
	}
	
}



//*********************************************************************************************
/** @brief Multiply nonlinear field by \f$ \exp(-\nu K^2 dt) \f$ or 
 *			\f$ \exp(-\nu K^2 dt) + \exp(-\nu_h)  K^4 dt) \f$
 * 
 *  @param dt
 *
 * @return \f$ N(k) = N(k) * \exp(-\nu K^2 dt) \f$.
 * @return when hyper_dissipation_switch == 1, 
 *			\f$ N(k) = N(k) * \exp(-\nu K^2 dt) + \exp(-\nu_h)  K^4 dt) \f$
 */

void FluidVF::Mult_nlin_exp_ksqr_dt(Real a)
{
	if (abs(a) > MYEPS) {
		if (!hyper_dissipation_switch) {
			universal->Array_exp_ksqr(nlin1, -dissipation_coefficient*a*global.time.dt);
			universal->Array_exp_ksqr(nlin2, -dissipation_coefficient*a*global.time.dt);
			universal->Array_exp_ksqr(nlin3, -dissipation_coefficient*a*global.time.dt);
		}
	
		else {
			universal->Array_exp_ksqr(nlin1, -dissipation_coefficient*a*global.time.dt, -hyper_dissipation_coefficient*a*global.time.dt, hyper_dissipation_exponent);
			
			universal->Array_exp_ksqr(nlin2, -dissipation_coefficient*a*global.time.dt, -hyper_dissipation_coefficient*a*global.time.dt, hyper_dissipation_exponent);
			
			universal->Array_exp_ksqr(nlin3, -dissipation_coefficient*a*global.time.dt, -hyper_dissipation_coefficient*a*global.time.dt, hyper_dissipation_exponent);
		}
	}
}



//******************************************************************************

/** @brief D=3: Add complex conjugate at \f$ \vec{k} \f$ 
 *				once \f$ \vec{V}(kx, ky, kz) \f$ is given.
 *  
 *	@bug  Add_complex_conj does not work across the processors.
 * @note  For ky = 0 and ky = N(2)/2 planes only.
 * 
 * @return  FFF: \f$ \vec{V}(-kx, -ky, kz) = \vec{V}^*(kx, ky, kz)\f$
 * @return  SCFT: For (ky != 0) -> \f$ \vec{V}(kx, -ky, kz) = \vec{V}^*(kx, ky, kz)\f$.
 * @return  SCFT: For (ky = 0) -> \f$ imag(V_y(kx, 0, kz)) = imag(V_z(kx, 0, kz)) = 0\f$,
 *					and \f$ V_x(kx, 0, kz) = 0 \f$.
 *				
 */
void FluidVF::Add_complex_conj(int kx, int ky, int kz, Complex Vx, Complex Vy, Complex Vz)			
{	
	TinyVector<Complex,3> localV;
	
		// On kz=0 or kz=N[3]/2 plane
	if ((basis_type == "FFF" || basis_type == "FFFW") && (kz == 0)) {
		localV = conj(Vx), conj(Vy), conj(Vz);
		universal->Assign_spectral_field(-kx, -ky, kz, cvf.V1, cvf.V2, cvf.V3, localV);
		if (master)
			cout << "Complex-conj(V) added for k = (" << -kx << "," << -ky << "," << kz << ")" << endl;	
	}
	
	else if ((basis_type == "SFF") && (kz == 0)) {          
		localV = conj(Vx), conj(Vy), conj(Vz);
		universal->Assign_spectral_field(kx, -ky, kz, cvf.V1, cvf.V2, cvf.V3, localV);
		if (master)
			cout << "Complex-conj(V) added for k = (" << kx << "," << -ky << "," << kz << ")" << endl;	
	}
}

void FluidVF::Assign_field_and_comp_conj(int kx, int ky, int kz, Complex Vx, Complex Vy, Complex Vz)
{	
    
    TinyVector<Complex,3> localV(Vx, Vy, Vz);

	universal->Assign_spectral_field(kx, ky, kz, cvf.V1, cvf.V2, cvf.V3, localV);	// in appropriate proc.
	
	Add_complex_conj(kx, ky, kz, Vx, Vy, Vz);	
}

/// 3D: Generate random vector at (kx,ky,kz) with range = rand_range

void FluidVF::Assign_random_complex_vector(int kx, int ky, int kz, Real rand_range)
{
	Complex  Vx, Vy, Vz;
	
	Vx = Complex( 2*rand_range*(SPECrand.random()-0.5), 2*rand_range*(SPECrand.random()-0.5) );
	
	Vy = Complex( 2*rand_range*(SPECrand.random()-0.5), 2*rand_range*(SPECrand.random()-0.5) );
    
   // universal->Last_component(kx, ky, kz, Vx, Vy, Vz);
	
	Assign_field_and_comp_conj(kx, ky, kz, Vx, Vy, Vz);
}

void FluidVF::Assign_random_real_vector(int kx, int ky, int kz, Real rand_range)
{
	Real Vx, Vy, Vz;
	
	Vx = 2*rand_range*(SPECrand.random()-0.5);
	
	Vy = 2*rand_range*(SPECrand.random()-0.5);
	
	universal->Last_component(kx, ky, kz, Vx, Vy, Vz);
	
	//Assign_field(kx, ky, kz, Vx, Vy, Vz);
}


void FluidVF::Assign_field(int kx, int ky, int kz, Real Vx, Real Vy, Real Vz)
{	
    
    TinyVector<Real,3> localV(Vx, Vy, Vz);
	
	universal->Assign_spectral_field(kx, ky, kz, cvf.V1, cvf.V2, cvf.V3, localV);	// in appropriate proc.
}

//*********************************************************************************************
/// 3D: Return Tk = Real(-nlin(k). conj(V(k))  for vector V
Real FluidVF::Get_Tk(int kx, int ky, int kz)
{	
	TinyVector<Complex,3> localV_complex, localnlin_complex;
	TinyVector<Real,3> localV_real, localnlin_real;
	
	if (universal->Probe_in_me(kx,ky,kz)) {
		if (basis_type == "SSS") {
			localV_real = real(universal->Get_spectral_field(kx, ky, kz, cvf.V1, cvf.V2, cvf.V3));
			localnlin_real = real(universal->Get_spectral_field(kx, ky, kz, nlin1, nlin2, nlin3));
			return dot(localV_real, localnlin_real);
		}
		
		else {
			localV_complex = universal->Get_spectral_field(kx, ky, kz, cvf.V1, cvf.V2, cvf.V3);
			localnlin_complex = universal->Get_spectral_field(kx, ky, kz, nlin1, nlin2, nlin3);
			return mydot(localV_complex, localnlin_complex);
		}
	}

	return 0;
}


//*********************************************************************************************

void FluidVF::Get_local_max_real_space(Real local_ur[])
{	
	local_ur[0] = max(abs(rvf.V1r));
	local_ur[2] = max(abs(rvf.V3r));
		
	if (Ny > 1) 
		local_ur[1] = max(abs(rvf.V2r));
	else
		local_ur[1] = 0.0;
}

//*********************************************************************************************
/*
Real FluidVF::Get_dt()
{
	Real local_ux,local_uy,local_uz;
	Real ux, uy, uz;
	
	if (global.program.dt_option == 0) {
		return global.time.dt_fixed;
	}
	
	else if (global.program.dt_option == 1) {
		Get_local_max_real_space(local_ux,local_uy,local_uz);
		
		MPI_Reduce(&local_ux, &ux, 1, MPI_Real, MPI_MAX, master_id, MPI_COMM_WORLD);
		MPI_Reduce(&local_uy, &uy, 1, MPI_Real, MPI_MAX, master_id, MPI_COMM_WORLD);
		MPI_Reduce(&local_uz, &uz, 1, MPI_Real, MPI_MAX, master_id, MPI_COMM_WORLD);
		
		Real sum_ubydx = ux/global.field.Delta_x[1] +  uz/global.field.Delta_x[3];
		
		if (Ny > 1)
			sum_ubydx += uy/global.field.Delta_x[2];
		
		Real dt = global.time.Courant_no/sum_ubydx;
		
		MPI_Bcast(&dt, 1, MPI_DOUBLE, master_id, MPI_COMM_WORLD); 
		
		return min(dt, global.time.dt_fixed);
	}
}	
*/



Real FluidVF::Get_dt()
{
	Real local_ur[3];
	Real ur[3];
	
	if (global.program.dt_option == 0) {
		return global.time.dt_fixed;
	}
	
	else if (global.program.dt_option == 1) {
		
		Get_local_max_real_space(local_ur);
		
		MPI_Allreduce(local_ur, ur, 3, MPI_Real, MPI_MAX, MPI_COMM_WORLD);
		
		Real sum_ubydx = ur[0]/global.field.Delta_x[1] +  ur[2]/global.field.Delta_x[3];
		
		if (Ny > 1)
			sum_ubydx += ur[1]/global.field.Delta_x[2];
		
		Real dt = global.time.Courant_no/sum_ubydx;
		
		return min(dt, global.time.dt_fixed);
	}

	return -1;
}	


Real FluidVF::Get_dt(FluidSF& T)
{
	return Get_dt();
}


// Alfven waves move with the speed of local mean magnetic field.. So we look for Bmax
// in the real space.
Real FluidVF::Get_dt(FluidVF& W)
{
	
	
	if (global.program.dt_option == 0) {
		return global.time.dt_fixed;
	}
	
	else if (global.program.dt_option == 1) {
		
		Real local_ur[3], local_wr[3];
		
		Real local_uwr[6];
		Real uwr[6];
		
		Get_local_max_real_space(local_ur);
		W.Get_local_max_real_space(local_wr);
		
		local_uwr[0] = local_ur[0];
		local_uwr[1] = local_ur[1];
		local_uwr[2] = local_ur[2];
		local_uwr[3] = local_wr[0]; // w fields
		local_uwr[4] = local_wr[1];
		local_uwr[5] = local_wr[2];
		
		MPI_Allreduce(local_uwr, uwr, 6, MPI_Real, MPI_MAX, MPI_COMM_WORLD);
		
		// for fluid
		Real sum_ubydx = (uwr[0]+uwr[3])/global.field.Delta_x[1] +  (uwr[2]+uwr[5])/global.field.Delta_x[3];
		
		if (Ny > 1)
			sum_ubydx += (uwr[1]+uwr[4])/global.field.Delta_x[2];
		
		Real dt = global.time.Courant_no/sum_ubydx;
		
		return min(dt, global.time.dt_fixed);
	}

	return -1;
}


Real FluidVF::Get_dt(FluidVF& W, FluidSF& T)
{
	return Get_dt(W);
}


void FluidVF::Make_incomplete_shells_zero() {
    int kx_max, ky_max, kz_max, kx_min, ky_min, kz_min;
    int nf;
    int lx, ly, lz;
    Real Kmag, kmagMax;

    kmagMax = universal->Max_radius_inside();
// service part
    kx_max = (int) ceil(kmagMax/kfactor[1]);
    
    if (Ny > 1)
        ky_max = (int) ceil(kmagMax/kfactor[2]);
    else
        ky_max = 0; 
    
    kz_max = (int) ceil(kmagMax/kfactor[3]);
    
    kx_min = ky_min = kz_min = 0;
    
    if (basis_type == "FFF" || basis_type == "FFFW")
        kx_min = -kx_max;
    
    if ((basis_type == "FFF" || basis_type == "FFFW") || (basis_type == "SFF"))
        ky_min = -ky_max;

// end of service part  
    TinyVector<Complex,3> f;
    f = 0.0;

    for (int kx = kx_min; kx <= kx_max; kx++)
    for (int ky = ky_min; ky <= ky_max; ky++)  
    for (int kz = 0; kz <= kz_max; kz++) {
        if (universal->Probe_in_me(kx,ky,kz))  {
            lx = universal->Get_lx(kx);
            ly = universal->Get_ly(ky);
            lz = universal->Get_lz(kz);
                
            Kmag = universal->Kmagnitude(lx, ly, lz);
            if ((Kmag > kmagMax)) {
                universal->Assign_local_spectral_field(lx, ly, lz, cvf.V1, cvf.V2, cvf.V3, f);
            }
        }
    }
}

void FluidVF::Filter_field_and_write_real(Real k0, int dir, int type, Real kf, string field_name) {
//dir = 1 - procedure leaves only larger than k0 modes
//dir = -1 - smaller
    CVF complexF("buffer");

    Filter_field_and_write_real(k0, dir, type, kf, field_name, complexF);

}

void FluidVF::Filter_field_and_write_real(Real k0, int dir, int type, Real kf, string field_name,CVF complexF) {
//dir = 1 - procedure leaves only larger than k0 modes
//dir = -1 - smaller
	int kx, ky, kz;
    TinyVector<Real, 3> k, filter;
    TinyVector<Complex, 3> result, fff;

    for (int lx=0; lx<global.field.maxlx; lx++)
    for (int ly=0; ly<global.field.maxly; ly++)
    for (int lz=0; lz<global.field.maxlz; lz++) {
                    
        kx = universal->Get_kx(lx);
        ky = universal->Get_ky(ly);
        kz = universal->Get_kz(lz);

        k = kx, ky, kz;

        fff = universal->Get_local_spectral_field(lx, ly, lz, cvf.V1, cvf.V2, cvf.V3);
        
        if(type == 3) {
			filter = 1.0, 1.0, 1.0;
			filter *= exp(-fabs(mynorm(k) - k0)/kf);
        } else {
	        if (dir*mynorm(k) >= dir*k0) {
	            filter = 1.0, 1.0, 1.0;
	        } else {
	        	if (type == 1) { 
	        	    filter = 0.0, 0.0, 0.0;
				}
				if (type == 2) {
					filter = 1.0, 1.0, 1.0;
	                filter *= exp( (mynorm(k) - k0)/kf );
				}
			}
		}
        result = filter * fff;
        universal->Assign_local_spectral_field(lx, ly, lz, complexF.V1, complexF.V2, complexF.V3, result);
    }             
    RVF realF(field_name);
    realF.Inverse_transform(complexF);
    realF.Write_real_field();
}

void FluidVF::Helical_decomposition(FluidVF& plus, FluidVF& minus){
	int kx, ky, kz;
   	TinyVector<Complex,3> k, z, U, Uplus, Uminus;
   	TinyVector<Complex,3> hP, hM;
	Complex uP, uM;
   	TinyVector<Complex,3> ZxK;
   	TinyVector<Complex,3> hr, hi;
   	TinyVector<Complex,3> temp;

    for (int lx=0; lx<global.field.maxlx; lx++)
    for (int ly=0; ly<global.field.maxly; ly++)
    for (int lz=0; lz<global.field.maxlz; lz++) {
                    
        kx = universal->Get_kx(lx);
        ky = universal->Get_ky(ly);
        kz = universal->Get_kz(lz);

        k = kx, ky, kz;

        if (mynorm(k)==0) {
        	Uplus = 0.0, 0.0, 0.0;
        	Uminus = 0.0, 0.0, 0.0;

	        universal->Assign_local_spectral_field(lx, ly, lz, plus.cvf.V1, plus.cvf.V2, plus.cvf.V3, Uplus);
	        universal->Assign_local_spectral_field(lx, ly, lz, minus.cvf.V1, minus.cvf.V2, minus.cvf.V3, Uminus);

        	continue;
        }

		do {
			z = 2*SPECrand.random()-1,
				2*SPECrand.random()-1,
				2*SPECrand.random()-1;
			ZxK = mycross(z,k);
		} while (mynorm(ZxK) < pow(10.0,-3));

    	// z = -0.6233810140363079, 0.6889268375950501, -0.04868331309267637;

		U = universal->Get_local_spectral_field(lx, ly, lz, cvf.V1, cvf.V2, cvf.V3);

		hr = mycross(ZxK,k)/( mynorm(k) * mynorm(ZxK) );
		hi = I*ZxK/mynorm(ZxK);
      
        Uplus = mydot_complex(U, hr-hi) * (hr+hi)  / 2.0;
        Uminus = mydot_complex(U, hr+hi) * (hr-hi) / 2.0;
        
        universal->Assign_local_spectral_field(lx, ly, lz, plus.cvf.V1, plus.cvf.V2, plus.cvf.V3, Uplus);
        universal->Assign_local_spectral_field(lx, ly, lz, minus.cvf.V1, minus.cvf.V2, minus.cvf.V3, Uminus);
    }  
}


void FluidVF::Write_real_field(string field_name){
	RVF realF(field_name);
	realF.Inverse_transform(cvf);
    realF.Write_real_field();
}


Real FluidVF::Get_dt_new(FluidVF& W)
{
	if (global.program.dt_option == 0) {
		return global.time.dt_fixed;
	}
	
	else if (global.program.dt_option == 1) {
		
		Real local_ur[3], local_wr[3];
		
		Real local_uwr[6];
		Real uwr[6];
		
		Get_local_max_real_space(local_ur);
		W.Get_local_max_real_space(local_wr);
		
		local_uwr[0] = local_ur[0];
		local_uwr[1] = local_ur[1];
		local_uwr[2] = local_ur[2];
		local_uwr[3] = local_wr[0]; // w fields
		local_uwr[4] = local_wr[1];
		local_uwr[5] = local_wr[2];
		
		MPI_Allreduce(local_uwr, uwr, 6, MPI_Real, MPI_MAX, MPI_COMM_WORLD);
		
		// for fluid
		Real sum_ubydx = (uwr[0]+uwr[3])/global.field.Delta_x[1] +  (uwr[2]+uwr[5])/global.field.Delta_x[3];
		
		if (Ny > 1)
			sum_ubydx += (uwr[1]+uwr[4])/global.field.Delta_x[2];
		
		Real dt = global.time.Courant_no/sum_ubydx;

		return min(global.time.dt * 1.1, dt);
		
		//return min(dt, global.time.dt_fixed);
	}

	return -1;
}


//*********************************************************************************************



//************************ END of CVF class Definitions ***************************************
