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

/*! \file  EnergyTr.h
 * 
 * @brief  Class declaration of energy transfer vars and functions. 
 * 
 *	We omputes energy transfers from region A1 of Giver field G, to the receiver field R
 *	with the help of helper field H.
 *
 *	The Giver field is filled in k-space appropriately (region A1).  <BR>
 *	We do the same for the receiver field (region A2). <BR>
 * 
 *	Definitions:
 *		\f$ T^{G,A1}_{R,A2} = Sign \sum_{A2} R_i^*(\vec{K})  N_i(\vec{K}) \f$ <BR>
 *		where the nonlinear term 
 *		\f$ \mathcal{N} = \mathcal{F}( (\vec{H}.\nabla) \vec{G}^{A1} ) \f$. <BR>
 *
 *	For u-to-u:	Sign = -; H, G, R = \f$ \vec{V} \f$. <BR>
 *	For b-to-b:	Sign = -; H = \f$ \vec{V} \f$; G, R = \f$ \vec{B} \f$. <BR> 
 *	For u-to-b: Sign = +; G = \f$ \vec{V} \f$; H, R = \f$ \vec{B} \f$. <BR>
 *  For b-to-u: Sign = +; R = \f$ \vec{V} \f$; G, H = \f$ \vec{B} \f$. <BR>
 *	
 *	For \f$\psi\f$ to \f$ \psi \f$: Sign = -; H = \f$ \vec{V} \f$, and G,R= \f$ \psi \f$.
 *
 *	
 *	Shell radii 0,r(1),r(2), .., r(n-1), r(noshells)=Infty. <BR>
 *	Shell(n) = (r(n-1),r(n)]. <BR> 
 *	There are n+1 shells with Shell(0) being the origin. <BR>
 *	The last shell contains modes beyond the largest sphere that fits inside the box.
 *
 *	For rings, we divide angles into sectors, with angles being 
 *		\f$ 0, \theta_{max}/n, 2*\theta_{max}/n, .., \theta_{max} \f$. <BR>
 *	thetamax can be either \f$ \pi \f$ or \f$ pi/2 \f$.
 *	Sector(n) = \f$ (\theta(n-1), \theta(n)] \f$.  <BR>
 *	Sector(1) contains \f$ \theta=0 \f$ points as well.
 *	Sector(0) does not contain any point.
 *
 *	For cylindrical rings, we divide the anisotropic directions into noslabs slabs.
 *	Slab(n) = (H(n-1), H(n)] with H(n)=0, .., kmax. 
 *
 *	ET_shell_radii_sector_array() -- entries <BR>
 *		Shell: 0, r(1), r(2), .. r(n-1); Inf skipped.  <BR>
 *		Ring-shell: 0, r(1),  r(2), .. r(n-1); Inf skipped.  <BR>
 *		Ring-sector: 0, s(1), s(2), ... s(n).  <BR>
 *		cylinder-ring-shell: 0, r(1), r(2), ..., r(n-1); Inf skipped.<BR>
 *
 *	Shell_to_shell_transfer is computed for shells 1:noshells. Note that the last shell
 *	contains modes outside the largest sphere.
 *
 *	For ring-to-ring transfer we ignore the last shell since the ring-to-ring transfer
 *	takes too much time.  <BR>
 *	For ring-to-ring transfer it is better to use-defined r(i).
 *   
 * Skpq commented out...
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Dec 2008
 *
 * @bug   No known bugs
 */


//*********************************************************************************************

#ifndef _EnergyTr
#define _EnergyTr

#include "fluid_base.h"
#include "Pressure.h"
#include "Nlin.h"

//*********************************************************************************************

class EnergyTr
{ 
 public:
	
	PlainFluidVF Giver;  // Gr_i, G_i

public:


	//! for energy transfer to be passed Shell_mult_all
	Array<Real,1>  temp_shell_tr;				
	
	//! for energy transfer to be passed Shell_mult_all
	Array<Real,2>  temp_ring_tr;					
	
	//! for energy transfer to be passed Shell_mult_all
	Array<Real,2>  temp_cylindrical_ring_tr;

	// S(k,p,q) calculations
	
/*	int				no_of_Skpq_triads;
	
	Array<int,2>	*triad_array;					// grid index e.g., 1,2,...
	Array<Real,1>		*Sself_array;
	Array<Real,1>		*S_to_VF_array;
	Array<Real,1>		*S_SF_array; */

	// For FLUIDS
	// Fluid fluxes
	Array<Real,1> flux_self;	
	Array<Real,1> flux_hk;

	
	// Isotropic shell-to-shell
	Array<Real,2>  shelltoshell_self;
	Array<Real,2>  shelltoshell_hk;
	Array<Real,2> shelltoshell_hk_U_to_helicalU;
	Array<Real,2> shelltoshell_hk_helicalU_to_helicalU;
	
	// Ring to ring
	Array<Real,4>  ring_to_ring_self;
	Array<Real,4> ring_to_ring_U_to_helicalU;
	Array<Real,4> ring_to_ring_helicalU_to_helicalU;
	
	// Cylinderical ring to ring
	Array<Real,4>  cylindrical_ring_to_ring_self;
	
	
	// FOR A SCALAR
	// flux
	Array<Real,1>  flux_SF;						//  T< to T>
	
	// shell-to-shell
	Array<Real,2>  shelltoshell_SF;				// theta->theta
	
	// Ring to ring
	Array<Real,4>  ring_to_ring_SF;				// theta->theta
	
	// Cylinderical ring to ring
	Array<Real,4>  cylindrical_ring_to_ring_SF;				// theta->theta
	
	
	// FOR MHD and Enstrophy calculation
	// fluxes
    Array<Real,1>  flux_VF;
	Array<Real,1>  flux_VF_Win_Wout;				// from Win to Wout
	Array<Real,1>  flux_VF_Uin_Wout;				// from Uin to Wout
	Array<Real,1>  flux_VF_Uin_Uout;                // U giver U receiver
	Array<Real,1>  flux_VF_Bin_Wout;                // B giver W receiver
	Array<Real,1>  flux_VF_Bin_Uout;                // B giver U receiver
	Array<Real,1>  flux_VF_Jin_Uout;                // J giver U receiver

	//Helical decomposition fluxes for HD
	Array<Real,1>  flux_U_grad_Uplus_Uplus;                // U+ giver U+ receiver U  mediator
	Array<Real,1>  flux_U_grad_Uplus_Uminus;               // U+ giver U- receiver U  mediator
	Array<Real,1>  flux_U_grad_Uminus_Uplus;               // U- giver U+ receiver U  mediator
	Array<Real,1>  flux_U_grad_Uminus_Uminus;              // U- giver U- receiver U  mediator
	Array<Real,1>  flux_Uplus_grad_Uplus_Uplus;          // U+ giver U+ receiver U+ mediator
	Array<Real,1>  flux_Uplus_grad_Uplus_Uminus;         // U+ giver U- receiver U+ mediator
	Array<Real,1>  flux_Uminus_grad_Uminus_Uminus;       // U- giver U- receiver U- mediator
	Array<Real,1>  flux_Uminus_grad_Uminus_Uplus;        // U- giver U+ receiver U- mediator
	//Helical decomposition fluxes for MHD	
	Array<Real,1>  flux_B_grad_Bplus_Uplus;                // U+ giver U+ receiver U  mediator
	Array<Real,1>  flux_B_grad_Bplus_Uminus;               // U+ giver U- receiver U  mediator
	Array<Real,1>  flux_B_grad_Bminus_Uplus;               // U- giver U+ receiver U  mediator
	Array<Real,1>  flux_B_grad_Bminus_Uminus;              // U- giver U- receiver U  mediator
	Array<Real,1>  flux_Bplus_grad_Bplus_Uplus;          // U+ giver U+ receiver U+ mediator
	Array<Real,1>  flux_Bplus_grad_Bplus_Uminus;         // U+ giver U- receiver U+ mediator
	Array<Real,1>  flux_Bminus_grad_Bminus_Uminus;       // U- giver U- receiver U- mediator
	Array<Real,1>  flux_Bminus_grad_Bminus_Uplus;        // U- giver U+ receiver U- mediator

	Array<Real,1>  flux_U_grad_Bplus_Bplus;                // U+ giver U+ receiver U  mediator
	Array<Real,1>  flux_U_grad_Bplus_Bminus;               // U+ giver U- receiver U  mediator
	Array<Real,1>  flux_U_grad_Bminus_Bplus;               // U- giver U+ receiver U  mediator
	Array<Real,1>  flux_U_grad_Bminus_Bminus;              // U- giver U- receiver U  mediator
	Array<Real,1>  flux_Uplus_grad_Bplus_Bplus;          // U+ giver U+ receiver U+ mediator
	Array<Real,1>  flux_Uplus_grad_Bplus_Bminus;         // U+ giver U- receiver U+ mediator
	Array<Real,1>  flux_Uminus_grad_Bminus_Bminus;       // U- giver U- receiver U- mediator
	Array<Real,1>  flux_Uminus_grad_Bminus_Bplus;        // U- giver U+ receiver U- mediator

	Array<Real,1>  flux_B_grad_Uplus_Bplus;                // U+ giver U+ receiver U  mediator
	Array<Real,1>  flux_B_grad_Uplus_Bminus;               // U+ giver U- receiver U  mediator
	Array<Real,1>  flux_B_grad_Uminus_Bplus;               // U- giver U+ receiver U  mediator
	Array<Real,1>  flux_B_grad_Uminus_Bminus;              // U- giver U- receiver U  mediator
	Array<Real,1>  flux_Bplus_grad_Uplus_Bplus;          // U+ giver U+ receiver U+ mediator
	Array<Real,1>  flux_Bplus_grad_Uplus_Bminus;         // U+ giver U- receiver U+ mediator
	Array<Real,1>  flux_Bminus_grad_Uminus_Bminus;       // U- giver U- receiver U- mediator
	Array<Real,1>  flux_Bminus_grad_Uminus_Bplus;        // U- giver U+ receiver U- mediator

	// For magnetic helicity fluxes
	Array<Real,1>  flux_VF_Bin_Aout;				// from Win to Aout
	Array<Real,1>  flux_VF_Uin_Aout_1;				// from Uin to Aout
	Array<Real,1>  flux_VF_Ain_Bout;                // A giver W receiver
	Array<Real,1>  flux_VF_Uin_Aout_2;                // U giver A receiver
    Array<Real,1>  flux_HM;
	

	// FOR MHD
    Array<Real,1>  flux_VF_U;
    Array<Real,1>  flux_VF_B;
	Array<Real,1>  flux_VF_Uin_Win;				// from Uin to Win
	Array<Real,1>  flux_VF_Win_Uout;				// from Win to Uout
	Array<Real,1>  flux_VF_Uout_Wout;			// from Wout to Uout
	Array<Real,1>  flux_Elsasser_plus;			// From Zplus_in to Zplus_out
	Array<Real,1>  flux_Elsasser_minus;			// From minus_in to Zplus_minus
	
	Array<Real,1>  flux_Whk;
	
	// shell-to-shell
	Array<Real,2>  shelltoshell_VF_WtoW;			// from W to W
	Array<Real,2>  shelltoshell_VF_UtoW;
  	Array<Real,2>  shelltoshell_VF_WtoU;
  	Array<Real,2>  shelltoshell_VF_UtoU;
  
    Array<Real,2>  shelltoshell_VF_BtoW;
    Array<Real,2>  shelltoshell_VF_JtoU;
    Array<Real,2>  shelltoshell_VF_BtoU;
  
	Array<Real,2>  shelltoshell_Elsasser_plus;	// from Zplus to Zplus
	Array<Real,2>  shelltoshell_Elsasser_minus;	// from Zminus to Zminus
	
	Array<Real,2>  shelltoshell_Whk;
	
  
	// Ring to ring
	Array<Real,4>  ring_to_ring_VF_WtoW;			// from W to W
	Array<Real,4>  ring_to_ring_VF_UtoW;			// from U to W
  	Array<Real,4>  ring_to_ring_VF_WtoU;
  	Array<Real,4>  ring_to_ring_VF_UtoU;
    Array<Real,4>  ring_to_ring_VF_BtoW;			// from U to W
    Array<Real,4>  ring_to_ring_VF_JtoU;
    Array<Real,4>  ring_to_ring_VF_BtoU;
	Array<Real,4>  ring_to_ring_Elsasser_plus;	// from Zplus to Zplus
	Array<Real,4>  ring_to_ring_Elsasser_minus;	// from Zminus to Zminus
	
	// Cylinderical ring to ring
	Array<Real,4>  cylindrical_ring_to_ring_VF_WtoW;		// from W to W
	Array<Real,4>  cylindrical_ring_to_ring_VF_UtoW;		// from U to W
	Array<Real,4>  cylindrical_ring_to_ring_Elsasser_plus;	// from Zplus to Zplus
	Array<Real,4>  cylindrical_ring_to_ring_Elsasser_minus;	// from Zminus to Zminus
	
	// B0 effect
	//! Energy transfer from b to u due to B0.
	//! \f$ \sum (\vec{B0} \cdot \vec{K}) \Im(\vec{u}(\vec{K'}) \vec{b}^*(\vec{K})) \f$
	Array<Real,1>	energy_tr_shell_B0;
	Array<Real,2>	energy_tr_ring_B0;
	Array<Real,2>	energy_tr_cylindrical_ring_B0;	
	
	// Force*V  or Force*T
	Array<Real,1> sphere_force_x_field;
	Array<Real,2> ring_force_x_field;
	Array<Real,2> cylindrical_ring_force_x_field;
	
public:
	EnergyTr();
	
	void Compute_flux(FluidVF &U);
	void Compute_flux(FluidVF& U, FluidSF& T);
	void Compute_flux_scalar(FluidVF& U, FluidSF& T);
	void Compute_flux_RBC(FluidVF& U, FluidSF& T);
	void Compute_flux(FluidVF& U, FluidVF& W);
	void Compute_flux(FluidVF& U, FluidVF& W, FluidSF& T);
    
    void Compute_flux_Vpll(FluidVF& U, Pressure& P);
    void Compute_flux_Vpll(FluidVF& U, FluidVF& W, Pressure& P);
    void Compute_shell_tr_Vpll(FluidVF& U, Pressure& P);
    void Compute_shell_tr_Vpll(FluidVF& U, FluidVF& W, Pressure& P);
    void Compute_ring_tr_Vpll(FluidVF& U, Pressure& P);
    void Compute_ring_tr_Vpll(FluidVF& U, FluidVF& W, Pressure& P);
    void Compute_cylindrical_ring_tr_Vpll(FluidVF& U, Pressure& P);
    void Compute_cylindrical_ring_tr_Vpll(FluidVF& U, FluidVF& W, Pressure& P);
    
    void Compute_flux_helical_decomposition(FluidVF &U);
    void Compute_flux_helical_decomposition(FluidVF &U, FluidVF &W);
	
	void Compute_kinetic_helicity_flux(FluidVF& U);
	void Compute_enstrophy_flux(FluidVF& U);

    void Compute_kinetic_helicity_flux(FluidVF& U, FluidVF& helicalU);
    void Compute_kinetic_helicity_flux_old(FluidVF& U, FluidVF& helicalU);
    void Compute_kinetic_helicity_flux(FluidVF& U, FluidVF& W, FluidVF& helicalU, FluidVF& helicalW);
    void Compute_magnetic_helicity_flux(FluidVF& U, FluidVF& W);
	void Compute_magnetic_helicity_flux(FluidVF& U, FluidVF& W, FluidVF& helicalW);
	void Compute_magnetic_enstrophy_flux(FluidVF& U, FluidVF& W);
	
	// Shell-to-shell
	void Compute_shell_tr(FluidVF &U);
	void Compute_shell_tr(FluidVF &U, FluidSF& T);
	void Compute_shell_tr_scalar(FluidVF &U, FluidSF& T);
	void Compute_shell_tr_RBC(FluidVF &U, FluidSF& T);
	void Compute_shell_tr(FluidVF &U, FluidVF& W);
	void Compute_shell_tr(FluidVF &U, FluidVF& W, FluidSF& T);
	
	
	void Compute_kinetic_helicity_shell_tr(FluidVF &U);
  	void Compute_kinetic_helicity_shell_tr(FluidVF &U, FluidVF &helicalU);
    void Compute_kinetic_helicity_shell_tr(FluidVF& U, FluidVF& W, FluidVF& helicalU, FluidVF& helicalW);
	void Compute_magnetic_helicity_shell_tr(FluidVF &U, FluidVF& W);
	void Compute_enstrophy_shell_tr(FluidVF &U, FluidVF &helicalU);
	void Compute_magnetic_enstrophy_shell_tr(FluidVF &U, FluidVF& W);
	
		// Ring_to_ring energy transfer
	void Compute_ring_tr(FluidVF &U);
	void Compute_enstrophy_ring_tr(FluidVF& U, FluidVF& helicalU);
    void Compute_kinetic_helicity_ring_tr(FluidVF& U,  FluidVF& helicalU);
	void Compute_ring_tr(FluidVF &U, FluidSF& T);
	void Compute_ring_tr_scalar(FluidVF &U, FluidSF& T);
	void Compute_ring_tr_RBC(FluidVF &U, FluidSF& T);
	void Compute_ring_tr(FluidVF &U, FluidVF& W);
    void Compute_kinetic_helicity_ring_tr(FluidVF& U, FluidVF& W, FluidVF& helicalU, FluidVF& helicalW);
	void Compute_ring_tr(FluidVF &U, FluidVF& W, FluidSF& T);
	
	
		// Cylinderical ring_to_ring energy transfer
	void Compute_cylindrical_ring_tr(FluidVF &U);
	void Compute_cylindrical_ring_tr(FluidVF &U, FluidSF& T);
	void Compute_cylindrical_ring_tr_scalar(FluidVF &U, FluidSF& T);
	void Compute_cylindrical_ring_tr_RBC(FluidVF &U, FluidSF& T);
	void Compute_cylindrical_ring_tr(FluidVF &U, FluidVF& W);
	void Compute_cylindrical_ring_tr(FluidVF &U, FluidVF& W, FluidSF& T);
	
	void Compute_shell_ET_B0(FluidVF &U, FluidVF& W);
	void Compute_ring_ET_B0(FluidVF &U, FluidVF& W);
	void Compute_cylindrical_ring_ET_B0(FluidVF &U, FluidVF& W);
	
	
		// force*V or force*T
	void Power_supply_within_sphere(FluidVF& U);
	void Power_supply_within_sphere(FluidSF& T);
	
	void Power_supply_ring(FluidVF& U);
	void Power_supply_ring(FluidSF& T);
	
	void Power_supply_cylindrical_ring(FluidVF& U);
	void Power_supply_cylindrical_ring(FluidSF& T);
	
	void Fill_in_sphere(int sphere_index, FluidVF& U);
	void Fill_in_sphere(int sphere_index, FluidSF& T);
    void Fill_in_sphere_Vpll(int sphere_index, FluidVF& U);
	
	void Fill_out_sphere(int sphere_index, FluidVF& U);						
	void Fill_out_sphere(int sphere_index, FluidSF& T);
    void Fill_out_sphere_Vpll(int sphere_index, FluidVF& U);
	
	void Fill_shell(int shell_from_index, FluidVF& U);						
	void Fill_shell(int shell_from_index, FluidSF& T);
    void Fill_shell_Vpll(int shell_from_index, FluidVF& U);
	
	void Fill_ring(int shell_from_index, int sector_from_index, FluidVF& U);
	void Fill_ring(int shell_from_index, int sector_from_index, FluidSF& T);
	void Fill_ring_Vpll(int ring_shell_from_index, int sector_from_index, FluidVF& U);
    ;
    
	void Fill_cylindrical_ring(int cylindrical_shell_from_index, int slab_from_index, FluidVF& U);
	void Fill_cylindrical_ring(int cylindrical_shell_from_index, int slab_from_index, FluidSF& T);
    void Fill_cylindrical_ring_Vpll(int cylindrical_shell_from_index, int slab_from_index, FluidVF& U);
	
	Real Prod_out_sphere_nlinV(int sphere_index,FluidVF& U, FluidVF& W);
	Real Prod_out_sphere_nlinT(int sphere_index, FluidVF& U, FluidSF& T);
    Real Prod_out_sphere_nlin_Vpll(int sphere_index, FluidVF& U, FluidVF& W);
	
	Real Prod_in_sphere_nlinV(int sphere_index, FluidVF& U, FluidVF& W);
	Real Prod_in_sphere_nlinT(int sphere_index, FluidVF& U, FluidSF& T);
	Real Prod_in_sphere_nlin_Vpll(int sphere_index, FluidVF& U, FluidVF& W);
	
	Real Prod_out_sphere_nlin_vorticity(int sphere_index, FluidVF& U, FluidVF& W);
	Real Prod_out_sphere_nlin_vector_potential(int sphere_index, FluidVF& U, FluidVF& W);

	Array<Real,1> Prod_by_spheres_nlinV(FluidVF& U, FluidVF& W);
};

#endif

//***************************   End Class decl of EnergyTr   **********************************


 

