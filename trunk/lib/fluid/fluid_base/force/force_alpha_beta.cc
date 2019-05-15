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

/*! \file  force_alpha_beta.cc
 * 
 * @brief Compute Force given alpha and beta
 *
 * @note 2D:   F(k) = alpha * V(k)
 * @note 3D;   F(k) = alpha * V(k) + beta(k) * Omega(k)
 *
 * @sa void IncVF::Put_random_vector_add_conj(int i1, int ly, Real amp)
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Feb 2009
 *
 * @bug   No known bugs
 */


#include "FORCE.h"


//*********************************************************************************************
/** @brief Assign a force vector at (kx, ky, kz) given alpha and beta,
 *					and its complex conj at \f$ -\vec{K} \f$. 
 * 
 *  @param kx, ky
 *	@param alpha, beta
 *
 *  @return  F(k) = alpha * V(k)  + beta(k) * Omega(k)
 */
void FORCE::Const_energy_supply_alpha_beta(FluidVF& U, int lx, int ly, int lz, Real alpha, Real beta, bool add_flag)
{
	if (basis_type == "SSS") { // Hk =0
		TinyVector<Real,3> localV_real, localForce_real;
		
		localV_real = real(universal->Get_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3));
		localForce_real = alpha*localV_real;
		
		if (!add_flag)
			universal->Assign_local_spectral_field(lx, ly, lz, U.Force1, U.Force2, U.Force3, localForce_real);
		else
			universal->Add_local_spectral_field(lx, ly, lz, U.Force1, U.Force2, U.Force3, localForce_real);
	}
		
	else {
		static int index=0;
		TinyVector<Complex,3> localV, localForce, vorticity;
		
		localV = universal->Get_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3);
		U.cvf.Compute_Modal_vorticity(lx, ly, lz, vorticity);

		localForce = alpha*localV + beta*vorticity;
		
		if (!add_flag)
			universal->Assign_local_spectral_field(lx, ly, lz, U.Force1, U.Force2, U.Force3, localForce);
		else 
			universal->Add_local_spectral_field(lx, ly, lz, U.Force1, U.Force2, U.Force3, localForce);
	}
}

//*********************************************************************************************
/** @brief Assign a force vector at (kx, ky, kz) given alpha and beta,
 *					and its complex conj at \f$ -\vec{K} \f$. 
 * 
 *  @param kx, ky
 *	@param alpha, beta
 *
 *  @return  F(k) = alpha * V(k)  + beta(k) * Omega(k)
 */
void FORCE::Const_energy_supply_alpha_beta(FluidVF& U, FluidVF& W, int lx, int ly, int lz, Real alpha, Real beta, bool add_flag)
{
	if (basis_type == "SSS") {  // Hk =0
	/*
  		TinyVector<Real,3> localV_real, localForce_real;
		
		localV_real = real(universal->Get_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3));
		localForce_real = alpha*localV_real;
		
		if (!add_flag)
			universal->Assign_local_spectral_field(lx, ly, lz, U.Force1, U.Force2, U.Force3, localForce_real);
		else
			universal->Add_local_spectral_field(lx, ly, lz, U.Force1, U.Force2, U.Force3, localForce_real);
	*/
		cerr << "FORCE::Const_energy_supply_alpha_beta(FluidVF& U, FluidVF& W, int lx, int ly, int lz, Real alpha, Real beta, bool add_flag): Not implemented for SSS basis" << endl;

	} else {
		static int index=0;
		TinyVector<Complex,3> localV, localForce, localW;
		Real phiX, phiY, phiZ, div;

		localV = universal->Get_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3);
		//localW = universal->Get_local_spectral_field(lx, ly, lz, W.cvf.V1, W.cvf.V2, W.cvf.V3);

		U.cvf.Compute_Modal_vorticity(lx, ly, lz, localW);

		if (global.force.field_procedure == 24) {
			beta = 0.1/(sqrt(W.cvf.Modal_energy(lx, ly, lz)) + pow(10,-8));
			localW = universal->Get_local_spectral_field(lx, ly, lz, W.cvf.V1, W.cvf.V2, W.cvf.V3);
		}

		localForce = alpha*localV + beta*localW;

/*		localW = universal->Get_local_spectral_field(lx, ly, lz, W.cvf.V1, W.cvf.V2, W.cvf.V3);		
		div = 0.1;
		localForce += div*localW;

		if (global.force.field_procedure == 24) {
			localW = universal->Get_local_spectral_field(lx, ly, lz, W.cvf.V1, W.cvf.V2, W.cvf.V3);

			phiX = 2 * M_PI * SPECrand.random();
			phiY = 2 * M_PI * SPECrand.random();
			phiZ = 2 * M_PI * SPECrand.random();

			localV(0) = cos(phiX) + I*sin(phiX);
			localV(1) = cos(phiY) + I*sin(phiY);
			localV(2) = cos(phiZ) + I*sin(phiZ);

			localForce = alpha * localV + beta * localW;

		}
		
*/		if (!add_flag)
			universal->Assign_local_spectral_field(lx, ly, lz, U.Force1, U.Force2, U.Force3, localForce);
		else 
			universal->Add_local_spectral_field(lx, ly, lz, U.Force1, U.Force2, U.Force3, localForce);
	}
}

//*********************************************************************************************
/** @brief Assign a force vector at (kx, ky, kz) given alpha and beta,
 *					and its complex conj at \f$ -\vec{K} \f$. 
 * 
 *  @param kx, ky
 *	@param alpha, beta
 *
 *  @return  F(k) = alpha * V(k)  + beta(k) * Omega(k)
 */
void FORCE::Const_magnetic_energy_supply_alpha_beta(FluidVF& U, FluidVF& W, int lx, int ly, int lz, Real alpha, Real beta, bool add_flag)
{
	if (basis_type == "SSS") {  // Hk =0
	/*
  		TinyVector<Real,3> localV_real, localForce_real;
		
		localV_real = real(universal->Get_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3));
		localForce_real = alpha*localV_real;
		
		if (!add_flag)
			universal->Assign_local_spectral_field(lx, ly, lz, U.Force1, U.Force2, U.Force3, localForce_real);
		else
			universal->Add_local_spectral_field(lx, ly, lz, U.Force1, U.Force2, U.Force3, localForce_real);
	*/
		cerr << "FORCE::Const_energy_supply_alpha_beta(FluidVF& U, FluidVF& W, int lx, int ly, int lz, Real alpha, Real beta, bool add_flag): Not implemented for SSS basis" << endl;

	} else {
		static int index=0;
		TinyVector<Complex,3> localV, localForce, localW;
		Real phiX, phiY, phiZ, div;

		localV = universal->Get_local_spectral_field(lx, ly, lz, W.cvf.V1, W.cvf.V2, W.cvf.V3);
		//localW = universal->Get_local_spectral_field(lx, ly, lz, W.cvf.V1, W.cvf.V2, W.cvf.V3);

		W.cvf.Compute_Modal_vorticity(lx, ly, lz, localW);

		// if (global.force.field_procedure == 24) {
		// 	beta = 0.1/(sqrt(W.cvf.Modal_energy(lx, ly, lz)) + pow(10,-8));
		// 	localW = universal->Get_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3);
		// }

		localForce = alpha*localV + beta*localW;
		if (!add_flag)
			universal->Assign_local_spectral_field(lx, ly, lz, W.Force1, W.Force2, W.Force3, localForce);
		else 
			universal->Add_local_spectral_field(lx, ly, lz, W.Force1, W.Force2, W.Force3, localForce);
	}
}

//*********************************************************************************************
/** @brief Assign a force vector at (kx, ky, kz) given alpha and beta,
 *					and its complex conj at \f$ -\vec{K} \f$. 
 * 
 *  @param kx, ky
 *	@param alpha, beta
 *
 *  @return  F(k) = alpha * V(k)  + beta(k) * Omega(k)
 */
void FORCE::Const_energy_supply_alpha_beta_gamma(FluidVF& U, FluidVF& W, int lx, int ly, int lz, Real alpha, Real beta, Real gamma, bool add_flag)
{
	if (basis_type == "SSS") {  // Hk =0
	/*
  		TinyVector<Real,3> localV_real, localForce_real;
		
		localV_real = real(universal->Get_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3));
		localForce_real = alpha*localV_real;
		
		if (!add_flag)
			universal->Assign_local_spectral_field(lx, ly, lz, U.Force1, U.Force2, U.Force3, localForce_real);
		else
			universal->Add_local_spectral_field(lx, ly, lz, U.Force1, U.Force2, U.Force3, localForce_real);
	*/
		cerr << "FORCE::Const_energy_supply_alpha_beta_gamma(FluidVF& U, FluidVF& W, int lx, int ly, int lz, Real alpha, Real beta, Real gamma, bool add_flag): Not implemented for SSS basis" << endl;

	} else {
		static int index=0;
		TinyVector<Complex,3> localV, localForce, localW, localC;
		Real phiX, phiY, phiZ, div;

		localV = universal->Get_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3);

		U.cvf.Compute_Modal_vorticity(lx, ly, lz, localW);

		localC = universal->Get_local_spectral_field(lx, ly, lz, W.cvf.V1, W.cvf.V2, W.cvf.V3);

		localForce = alpha*localV + beta*localW + gamma*localC;

		if (!add_flag) {
			universal->Assign_local_spectral_field(lx, ly, lz, U.Force1, U.Force2, U.Force3, localForce);
			////////////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////// ONLY FOR ADAPTIVE FORCE//////////////////////////////////////
			// localV*=alpha;
			// localW*=beta;
			// localC*=gamma;
			// universal->Assign_local_spectral_field(lx, ly, lz, global.temp_array.ForceA1, global.temp_array.ForceA2, global.temp_array.ForceA3, localV);
			// universal->Assign_local_spectral_field(lx, ly, lz, global.temp_array.ForceB1, global.temp_array.ForceB2, global.temp_array.ForceB3, localW);
			// universal->Assign_local_spectral_field(lx, ly, lz, global.temp_array.ForceC1, global.temp_array.ForceC2, global.temp_array.ForceC3, localC);
		}
		else {
			universal->Add_local_spectral_field(lx, ly, lz, U.Force1, U.Force2, U.Force3, localForce);
		}
	}
}

//*********************************************************************************************
/** @brief Assign a force vector at (kx, ky, kz) given alpha and beta,
 *					and its complex conj at \f$ -\vec{K} \f$. 
 * 
 *  @param kx, ky
 *	@param alpha, beta
 *
 *  @return  F(k) = alpha * V(k)  + beta(k) * Omega(k)
 */
void FORCE::Const_magnetic_energy_supply_alpha_beta_gamma(FluidVF& U, FluidVF& W, int lx, int ly, int lz, Real alpha, Real beta, Real gamma, bool add_flag)
{
	if (basis_type == "SSS") {  // Hk =0
	/*
  		TinyVector<Real,3> localV_real, localForce_real;
		
		localV_real = real(universal->Get_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3));
		localForce_real = alpha*localV_real;
		
		if (!add_flag)
			universal->Assign_local_spectral_field(lx, ly, lz, U.Force1, U.Force2, U.Force3, localForce_real);
		else
			universal->Add_local_spectral_field(lx, ly, lz, U.Force1, U.Force2, U.Force3, localForce_real);
	*/
		cerr << "FORCE::Const_magnetic_energy_supply_alpha_beta_gamma(FluidVF& U, FluidVF& W, int lx, int ly, int lz, Real alpha, Real beta, Real gamma, bool add_flag): Not implemented for SSS basis" << endl;

	} else {
		static int index=0;
		TinyVector<Complex,3> localV, localForce, localW, localC, temp;
		Real Kmag, div, kx, ky, kz, k2;
		Complex tempK;

		kx = universal->Get_kx(lx);
		ky = universal->Get_ky(ly);
		kz = universal->Get_kz(lz);

		Kmag = universal->Kmagnitude(lx, ly, lz);
		k2=Kmag*Kmag;
		div = (Kmag==0?1:Kmag*Kmag);

		localV = universal->Get_local_spectral_field(lx, ly, lz, W.cvf.V1, W.cvf.V2, W.cvf.V3);
				
		W.cvf.Compute_Modal_vorticity(lx, ly, lz, localW); //Here we get J

		localW /= div; //now we have A
		
		localC = universal->Get_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3);

		temp = alpha*localV + beta*localW + gamma*localC;

	    localForce(0) = (1-kx*kx/k2)*temp(0) -     kx*ky/k2*temp(1) -     kx*kz/k2*temp(2);
	    localForce(1) =    -ky*kx/k2*temp(0) + (1-ky*ky/k2)*temp(1) -     ky*kz/k2*temp(2);
	    localForce(2) =    -kz*kx/k2*temp(0) -     kz*ky/k2*temp(1) + (1-kz*kz/k2)*temp(2);


	    // tempK = kx*localForce(0) + ky*localForce(1) + kz*localForce(2) ;
	    // cout<< " k.f(k) = " << tempK<<endl;

		if (!add_flag)
			universal->Assign_local_spectral_field(lx, ly, lz, W.Force1, W.Force2, W.Force3, localForce);
		else 
			universal->Add_local_spectral_field(lx, ly, lz, W.Force1, W.Force2, W.Force3, localForce);
	}
}

//*********************************************************************************************

/** @brief Assign a force vector at (kx, ky, kz) given alpha,
 *					and its complex conj at \f$ -\vec{K} \f$. 
 * 
 *  @param kx, ky
 *	@param alpha
 *
 *  @return  F(k) = alpha * V(k).
 *  @note in SCFT basis, conversion to FFF and back to SCFT implies multiplication
 *			of (-I) and (I).  Hence the conversion remains unchanged.
 */
void FORCE::Const_energy_supply_alpha(FluidSF& T, int lx, int ly, int lz, Real alpha, bool add_flag)
{
	
	if (basis_type == "SSS")  { // Hk =0
		Real localG_real, localForce_real;
		
		localG_real = real(universal->Get_local_spectral_field(lx, ly, lz, T.csf.F));
		localForce_real = alpha*localG_real;
		
		if (!add_flag)
			universal->Assign_local_spectral_field(lx, ly, lz, T.Force, localForce_real);
		else
			universal->Add_local_spectral_field(lx, ly, lz, T.Force, localForce_real);
	}
	
	else {
		Complex localG, localForce;
		
		localG = universal->Get_local_spectral_field(lx, ly, lz, T.csf.F);
		
		localForce = alpha*localG;
		
		if (!add_flag)
			universal->Assign_local_spectral_field(lx, ly, lz, T.Force, localForce);
		else
			universal->Add_local_spectral_field(lx, ly, lz, T.Force, localForce);
	}
}


//*********************************************************************************************
/** @brief Assign a force vector at (kx, ky, kz) given alpha and beta,
 *					and its complex conj at \f$ -\vec{K} \f$. 
 * 
 *  @param kx, ky
 *	@param alpha, beta
 *
 *  @return  V(k) = alpha * V(k)  + beta(k) * Omega(k)
 */
void FORCE::Const_energy_alpha_beta(FluidVF& U, int lx, int ly, int lz, Real alpha, Real beta, bool add_flag)
{
	if (basis_type == "SSS")  { // Hk =0
		TinyVector<Real,3> localV_real;
		
		localV_real = real(universal->Get_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3));
		localV_real = alpha*localV_real;
		
		if (!add_flag)
			universal->Assign_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3, localV_real);
		else
			universal->Add_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3, localV_real);
	}
	
	else {
		TinyVector<Complex,3> localV, vorticity;
		
		localV = universal->Get_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3);
		U.cvf.Compute_Modal_vorticity(lx, ly, lz, vorticity);
		
		localV = alpha*localV + beta*vorticity;
		
		if (!add_flag)
			universal->Assign_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3, localV);
		else
			universal->Add_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3, localV);
	}
}

	//*********************************************************************************************

/** @brief Assign a force vector at (kx, ky, kz) given alpha,
 *					and its complex conj at \f$ -\vec{K} \f$. 
 * 
 *  @param kx, ky
 *	@param alpha
 *
 *  @return  F(k) = alpha * V(k).
 *  @note in SCFT basis, conversion to FFF and back to SCFT implies multiplication
 *			of (-I) and (I).  Hence the conversion remains unchanged.
 */
void FORCE::Const_energy_alpha(FluidSF& T, int lx, int ly, int lz, Real alpha, bool add_flag)
{
	
	if (basis_type == "SSS") { // Hk =0
		Real localG_real;
		
		localG_real  = real(universal->Get_local_spectral_field(lx, ly, lz, T.csf.F));
		localG_real = alpha*localG_real;
		
		if (!add_flag)
			universal->Assign_local_spectral_field(lx, ly, lz, T.csf.F, localG_real);
		else
			universal->Add_local_spectral_field(lx, ly, lz, T.csf.F, localG_real);
	}
	
	else {
		Complex localG;
		
		localG = universal->Get_local_spectral_field(lx, ly, lz, T.csf.F);
		
		localG = alpha*localG;
		if (!add_flag)
			universal->Assign_local_spectral_field(lx, ly, lz, T.csf.F, localG);
		else
			universal->Add_local_spectral_field(lx, ly, lz, T.csf.F, localG);
	}
}

//********************************** force_alpha_beta.cc **************************************

	
