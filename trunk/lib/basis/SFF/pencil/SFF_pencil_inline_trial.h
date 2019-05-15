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


/*! \file four_inline.h 
 * 
 * @brief Inline functions to compute array indices \f$ \vec{i} \f$ given wavenumber  
 * \f$ \vec{k} \f$ and viceversa.
 * Also contains other useful functions like Kmagnitude, Max radius etc.
 *
 * Grid wavenumber \f$ \vec{k} \f$ is integer, that can be computed using the grid index 
 * \f$ \vec{i} \f$.
 * Actual wavenumber \f$ \vec{K} \f$ is computed from the grid wavenumber \f$ \vec{k} \f$ 
 *		using \f$ K_i = k_i * f_i \f$
 * where \f$ f_i \f$ is the kfactor[i]. 
 *
 * For shell/ring energy spectrum \f$ E(k) \f$, the wavenumber is either actual wavenumber 
 *	or gridwavenumber depending on the switch WAVENOACTUAL or WAVENOGRID.   
 *	Actual wavenumber is preferred.
 * 
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date	August 2008
 * @bug		No known bugs
 *
 */ 
 
#include "SFF_pencil.h"

//*********************************************************************************************

/*! @brief	Get grid waveno kx given first local array index lx.
 * 
 *  if i1 <=Ni/2, ki=i1; else kx=i1-Ni.  <BR>
 *	i1= local_N1_start + lx.
 * 
 * \param lx  first local index of an array

 */
inline int SFF_PENCIL::Get_kx(int lx) {return  (lx <= Nx/2) ? lx : (lx-Nx); } 
 

/*! @brief	Get local array index lx given grid waveno kx.
 * 
 *  If kx>=0, i1=kx; else i1=kx+N1.  <BR>
 *	i1= local_N1_start + lx.
 * 
 * \param	kx  grid wavenumber along x
 * \return	lx  local array index along x
 */
inline int SFF_PENCIL::Get_lx(int kx) {return (kx >= 0) ? kx : (kx + Nx); } 


/*! @brief	Get array index ix given grid waveno kx.
 * 
 *  If kx>=0, i1=kx; else i1=kx+N1.  <BR>
 * 
 * \param	kx  grid wavenumber along x
 * \return	i1  array index along x
 */
inline int SFF_PENCIL::Get_ix(int kx) {return  (kx >= 0) ? kx : (kx + Nx);  }		


/*! @brief	Get grid waveno ky given first local array index ly.
 * 
 *  If ly<=N2/2, k2=ly; else ky=ly-N2.
 * 
 * \param ly  second local index of an array
 * \return kx corresponding to lx
 */
inline int SFF_PENCIL::Get_ky(int ly) {return  ((local_Ny_start + ly) <= Ny/2) ? (local_Ny_start + ly) : (local_Ny_start + ly-Ny); }

/*! @brief	Get local array index ly given grid waveno ky.
 * 
 *  If ky>=0, ly=ky; else ly=ky+N2.
 * 
 * \param	ky  grid wavenumber along y
 * \return	ly  local array index along y
 */
inline int SFF_PENCIL::Get_ly(int ky) { return  (ky >= 0) ? (ky-local_Ny_start) : (ky + Ny-local_Ny_start);  } 

inline int SFF_PENCIL::Get_iy(int ky) { return  (ky >= 0) ? ky : (ky + Ny); }

inline int SFF_PENCIL::Get_kz(int lz)  {return local_Nz_start + lz;} 

inline int SFF_PENCIL::Get_lz(int kz)  {return kz - local_Nz_start;}

inline int SFF_PENCIL::Get_iz(int kz)  {return kz;}

inline bool SFF_PENCIL::Probe_in_me(int kx, int ky, int kz) 
{  
    int ly = Get_ly(ky);
    int lz = Get_lz(kz);
    
	return ( ((ly >= 0) && (ly < local_Ny_vert)) && ((lz >= 0) && (lz < local_Nz_hor)) );
}


inline Complex SFF_PENCIL::Get_spectral_field(int kx, int ky, int kz, Array<Complex,3> A)
{ 
    int lx = Get_lx(kx);
    int ly = Get_ly(ky);
    int lz = Get_lz(kz);
    
    if ( ((ly >= 0) && (ly < local_Ny_vert)) && ((lz >= 0) && (lz < local_Nz_hor)) )
        return A(ly, lz, lx);

}

inline TinyVector<Complex,3> SFF_PENCIL::Get_spectral_field(int kx, int ky, int kz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az)
{
	
	int lx = Get_lx(kx);
    int ly = Get_ly(ky);
    int lz = Get_lz(kz);
    
    if ( ((ly >= 0) && (ly < local_Ny_vert)) && ((lz >= 0) && (lz < local_Nz_hor)) )
        return TinyVector<Complex,3>(Ax(ly, lz, lx), Ay(ly, lz, lx), Az(ly, lz, lx));
}


//  Assign
inline void SFF_PENCIL::Assign_spectral_field(int kx, int ky, int kz, Array<Complex,3> A, Complex field)
{

	int lx = Get_lx(kx);
    int ly = Get_ly(ky);
    int lz = Get_lz(kz);
    
    if ( ((ly >= 0) && (ly < local_Ny_vert)) && ((lz >= 0) && (lz < local_Nz_hor)) )
        A(ly, lz, lx) = field;

}

inline void SFF_PENCIL::Assign_spectral_field(int kx, int ky, int kz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, TinyVector<Complex,3> V)
{
	
	int lx = Get_lx(kx);
    int ly = Get_ly(ky);
    int lz = Get_lz(kz);
    
    if ( ((ly >= 0) && (ly < local_Ny_vert)) && ((lz >= 0) && (lz < local_Nz_hor)) ) {
        Ax(ly, lz, lx) = V(0);
        Ay(ly, lz, lx) = V(1);
        Az(ly, lz, lx) = V(2);
    }

}

inline void SFF_PENCIL::Assign_spectral_field(int kx, int ky, int kz, Array<Complex,3> A,Real field)
{ 
	cout << "MYERROR: SFF_PENCIL::Assign_spectral_field(); Use complex data type " << endl; 
}

inline void SFF_PENCIL::Assign_spectral_field(int kx, int ky, int kz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, TinyVector<Real,3> V)
{
	
	cout << "MYERROR: SFF_PENCIL::Assign_spectral_field(); Use complex data type " << endl;
}

inline void SFF_PENCIL::Add_spectral_field(int kx, int ky, int kz, Array<Complex,3> A, Complex field)
{

	int lx = Get_lx(kx);
    int ly = Get_ly(ky);
    int lz = Get_lz(kz);
    
    if ( ((ly >= 0) && (ly < local_Ny_vert)) && ((lz >= 0) && (lz < local_Nz_hor)) )
        A(ly, lz, lx) += field;

}

inline void SFF_PENCIL::Add_spectral_field(int kx, int ky, int kz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, TinyVector<Complex,3> V)
{
	
	int lx = Get_lx(kx);
    int ly = Get_ly(ky);
    int lz = Get_lz(kz);
    
    if ( ((ly >= 0) && (ly < local_Ny_vert)) && ((lz >= 0) && (lz < local_Nz_hor)) ) {
        Ax(ly, lz, lx) += V(0);
        Ay(ly, lz, lx) += V(1);
        Az(ly, lz, lx) += V(2);
    }

}

inline void SFF_PENCIL::Add_spectral_field(int kx, int ky, int kz, Array<Complex,3> A,Real field)
{ 
	cout << "MYERROR: SFF_PENCIL::Assign_spectral_field(); Use complex data type " << endl; 
}

inline void SFF_PENCIL::Add_spectral_field(int kx, int ky, int kz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, TinyVector<Real,3> V)
{
	
	cout << "MYERROR: SFF_PENCIL::Assign_spectral_field(); Use complex data type " << endl;
}




// Local field given local lx,ly,lz

inline Complex SFF_PENCIL::Get_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> A)
{ 
    return A(ly, lz, lx);
}

inline TinyVector<Complex,3> SFF_PENCIL::Get_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az)
{
	
	return TinyVector<Complex,3>(Ax(ly, lz, lx), Ay(ly, lz, lx), Az(ly, lz, lx));
}


	//  Assign
inline void SFF_PENCIL::Assign_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> A, Complex field)
{ 
    if ( ((ly >= 0) && (ly < local_Ny_vert)) && ((lz >= 0) && (lz < local_Nz_hor)) )
        A(ly, lz, lx) = field;
}

inline void SFF_PENCIL::Assign_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, TinyVector<Complex,3> V)
{
    if ( ((ly >= 0) && (ly < local_Ny_vert)) && ((lz >= 0) && (lz < local_Nz_hor)) ) {
        Ax(ly, lz, lx) = V(0);
        Ay(ly, lz, lx) = V(1);
        Az(ly, lz, lx) = V(2);
    }
}

inline void SFF_PENCIL::Assign_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> A,Real field)
{ 
	cout << "MYERROR: SFF_PENCIL::Assign_local_spectral_field(); Use complex data type " << endl; 
}

inline void SFF_PENCIL::Assign_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, TinyVector<Real,3> V)
{
	
	cout << "MYERROR: SFF_PENCIL::Assign_local_spectral_field(); Use complex data type " << endl;
}



inline void SFF_PENCIL::Add_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> A, Complex field)
{ 
    if ( ((ly >= 0) && (ly < local_Ny_vert)) && ((lz >= 0) && (lz < local_Nz_hor)) )
        A(ly, lz, lx) += field;
}

inline void SFF_PENCIL::Add_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, TinyVector<Complex,3> V)
{
    if ( ((ly >= 0) && (ly < local_Ny_vert)) && ((lz >= 0) && (lz < local_Nz_hor)) ) {
        Ax(ly, lz, lx) += V(0);
        Ay(ly, lz, lx) += V(1);
        Az(ly, lz, lx) += V(2);
    }
}

inline void SFF_PENCIL::Add_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> A,Real field)
{ 
	cout << "MYERROR: SFF_PENCIL::Assign_local_spectral_field(); Use complex data type " << endl; 
}

inline void SFF_PENCIL::Add_local_spectral_field(int lx, int ly, int lz, Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, TinyVector<Real,3> V)
{
	
	cout << "MYERROR: SFF_PENCIL::Assign_local_spectral_field(); Use complex data type " << endl;
}


// REAL - SPACE FNS
inline int SFF_PENCIL::Get_lx_real_space(int rx)  {return  rx - my_vert_pcoord*local_Nx_vert;}

inline int SFF_PENCIL::Get_ly_real_space(int ry)  {return  ry - my_hor_pcoord*local_Ny_hor; }	

inline int SFF_PENCIL::Get_lz_real_space(int rz) {return rz;}


inline int SFF_PENCIL::Get_rx_real_space(int lx)  {return  lx + my_vert_pcoord*local_Nx_vert;}

inline int SFF_PENCIL::Get_ry_real_space(int ly)  {return  ly + my_hor_pcoord*local_Ny_hor; }		

inline int SFF_PENCIL::Get_rz_real_space(int lz) {return lz;}



inline bool SFF_PENCIL::Probe_in_me_real_space(int rx, int ry, int rz) 
{
    int lx = Get_lx_real_space(rx);
    int ly = Get_ly_real_space(ry);
    
	return ( ((lx >= 0) && (lx < local_Nx_vert)) && ((ly >= 0) && (ly < local_Ny_hor)) );
}


// lz is the coord of the complex array
inline Real SFF_PENCIL::Get_real_field(int rx, int ry, int rz, Array<Real,3> A)
{	
    int lx = Get_lx_real_space(rx);
    int ly = Get_ly_real_space(ry);
    int lz = Get_lz_real_space(rz);
    
    if ( ((lx >= 0) && (lx < local_Nx_vert)) && ((ly >= 0) && (ly < local_Ny_hor)) )
        return (A(ly, lz, lx));
}

inline TinyVector<Real,3> SFF_PENCIL::Get_real_field(int rx, int ry, int rz, Array<Real,3> Ax, Array<Real,3> Ay, Array<Real,3> Az)
{
    
	int lx = Get_lx_real_space(rx);
    int ly = Get_ly_real_space(ry);
	int lz = Get_lz_real_space(rz);
    
    if ( ((lx >= 0) && (lx < local_Nx_vert)) && ((ly >= 0) && (ly < local_Ny_hor)) )
		return TinyVector<Real,3>(Ax(ly, lz, lx), Ay(ly, lz, lx), Az(ly, lz, lx));
	
}


// lz is the coord of the complex array
inline void SFF_PENCIL::Assign_real_field(int rx, int ry, int rz, Array<Real,3> A, Real field)
{	
    int lx = Get_lx_real_space(rx);
    int ly = Get_ly_real_space(ry);
	int lz = Get_lz_real_space(rz);
    
    if ( ((lx >= 0) && (lx < local_Nx_vert)) && ((ly >= 0) && (ly < local_Ny_hor)) )
		A(ly, lz, lx) = field;
}

inline void SFF_PENCIL::Assign_real_field(int rx, int ry, int rz, Array<Real,3> Ax, Array<Real,3> Ay, Array<Real,3> Az, TinyVector<Real,3> V)
{
    
	int lx = Get_lx_real_space(rx);
    int ly = Get_ly_real_space(ry);
	int lz = Get_lz_real_space(rz);
    
    if ( ((lx >= 0) && (lx < local_Nx_vert)) && ((ly >= 0) && (ly < local_Ny_hor)) ) {
        Ax(ly, lz, lx) = V(0);
		Ay(ly, lz, lx) = V(1);
		Az(ly, lz, lx) = V(2);
    }    
	
}





/**********************************************************************************************

	Compute Wavenumber components for grid index (i1,i2).

***********************************************************************************************/

// Real K
inline void SFF_PENCIL::Wavenumber(int lx, int ly, int lz, TinyVector<Real,3> & K)
{
	 K = Get_kx(lx)*kfactor[1],  Get_ky(ly)*kfactor[2], Get_kz(lz)*kfactor[3];
}


// Complex K; The imaginary part is zero.  Written to use cross function of blitz.
// Omega = cross(V,K).
inline void SFF_PENCIL::Wavenumber(int lx, int ly, int lz, TinyVector<Complex,3> & K)
{
	 K = Complex(Get_kx(lx)*kfactor[1], 0.0), Complex(Get_ky(ly)*kfactor[2], 0.0), Complex(Get_kx(lz)*kfactor[3], 0.0);
}




/**********************************************************************************************

				 If wavenos computed using actual wavenumber:  
						globalvar_waveno_switch = 0;  
						Ki = Kfactor[i]*grid[i]
				 
				 If wavenos computed using grid wavenumber:  Ki = grid[i]
						globalvar_waveno_switch = 1;  
						Ki = grid[i]

***********************************************************************************************/


///  WAVENOACTUAL: \f$ K = \sqrt{K_x^2 + K_y^2 + K_z^2} \f$
inline Real SFF_PENCIL::Kmagnitude(int lx, int ly, int lz)
{
	if	(global.field.waveno_switch)
		return sqrt(  pow2(Get_kx(lx)*kfactor[1])+ pow2(Get_ky(ly)*kfactor[2]) + pow2(Get_kz(lz)*kfactor[3]) );
	
	else 
		return sqrt(  pow2(lx) + pow2(Get_ky(ly)) + pow2(Get_kz(lz)) ); 
}

/// WAVENOACTUAL -- Radius of the smallest sphere that contains the wavenumber K box Ni's. <BR>
/// The range of ki = [-Ni/2+1 : Ni/2].
inline int SFF_PENCIL::Min_radius_outside() 
{
	if	(global.field.waveno_switch)
		return (int) ceil(sqrt(  pow2((Nx/2)*kfactor[1])+ pow2(Ny/2*kfactor[2]) + pow2(Nz/2*kfactor[3]) ));
	
	else
		return (int) ceil(sqrt(  pow2(Nx/2) + pow2(Ny/2) + pow2(Nz/2) ));
}


/// WAVENOACTUAL -- Radius of the largest sphere that fits inside the wavenumber K box Ni's. <BR>
/// The range of ki = [-Ni/2+1 : Ni/2].
inline int SFF_PENCIL::Max_radius_inside() 
{
	
	int ans = 1;
	Real Kmag;
	
	if	(global.field.waveno_switch)	{
        Kmag = min( (Nx/2)*kfactor[1], (Ny/2)*kfactor[2]);
        Kmag = min(Kmag, (Nz/2)*kfactor[3]); 

		ans = ((int) Kmag);
	}
	
	else  {
        ans = min(Nx/2, Ny/2);
        ans = min(ans, Nz/2);
    }
	
	return ans;
}


/*! \brief WAVENOACTUAL -- Returns the approximate number of modes in a wavenumber K semishell 
 *							(upper hemisphere) of radius "radius".
 * 
 *  We divide the area in Fourier space by volume of unit lattice \f$ \Pi_i f_i \f$, 
 *	where \f$ f_i \f$ is the factor[i].
 *  In 1D, The area in NOT divided by kfactor.
 *
 * For shell/ring energy spectrum \f$ E(k) \f$, the wavenumber is either actual wavenumber 
 * or gridwavenumber depending on the switch
 * WAVENOACTUAL or WAVENOGRID.   Actual wavenumber is preferred.
 *
 * \param  radius
 * \return The number of modes in a hemispheric shell of radius. 
 *			In 3D, it is hemisphere sphere with (kz>=0).
 */
inline Real SFF_PENCIL::Approx_number_modes_in_shell(int radius)
{	
    if	(global.field.waveno_switch)
        return (4*M_PI*radius*radius)/(kfactor[1]*kfactor[2]*kfactor[3]);	

    else 
       return (4*M_PI*radius*radius);
}


//*********************************************************************************************



/*! \brief Returns multiplication factor for computing enregy spectrum etc. 
 * 
 *  \f$ E(k) = |\vect{u}(\vect{k})|^2 /2 \f$, but in simulation we
 *  double the energy of most of the modes because complex conjugates modes with 
 * -ky are not stored in the simulation.  
 * The modes on the x-axis are not doubled because their c.c. are already counted.
 * 
 * \param  lx, ly
 * \return Multiplication factor to \f$ |\vect{u}(\vect{k})|^2  \f$ for computing 
 *		energy spectrum etc. factor = 1 implies that the modal energy is already doubled.
 */
 
inline Real SFF_PENCIL::Multiplicity_factor(int lx, int ly, int lz)
{
	Real factor;

	int kx = lx;
	int ky = Get_ky(ly);
    int kz = Get_kz(lz);
	
	if (kz > 0)
		factor = 1.0;
		
	else
		factor = 0.5;
		
	if ( (kx == Nx/2) || ((ky == Ny/2)&& (Ny > 1)) )
		return 2*factor;		// for kx = -Nx/2 or ky = -Ny/2 ; 
								// Ignoring corner which would have factor 4
	else
		return factor;
}



/// Modal energy -- \f$ E(k) = |\vect{u}(\vect{k})|^2 /2 \f$
inline Real SFF_PENCIL::Modal_energy(int lx, int ly, int lz, Array<Complex,3> A)
{
	return pow2(abs(A(lz, ly, lx)))/2;	
}

/**********************************************************************************************

	Compute Modal helicity for grid index (i1, i2, i3) = K. [Vr x Vi].

***********************************************************************************************/

inline Real SFF_PENCIL::Get_Modal_helicity
(
	int lx, int ly, int lz,
	Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az
)
{

	TinyVector<Real,3> Vreal, Vimag, VrcrossVi,  K;
	
	Vreal = real(Ax(ly, lz, lx)), real(Ay(ly, lz, lx)), real(Az(ly, lz, lx));
	Vimag = imag(Ax(ly, lz, lx)), imag(Ay(ly, lz, lx)), imag(Az(ly, lz, lx));
		
	VrcrossVi = cross(Vreal, Vimag);
	Wavenumber(lx, ly, lz,  K);
	
	return (dot( K, VrcrossVi));	

}



/**********************************************************************************************=

	Compute Modal Vorticity for grid index (i1, i2, i3) = I K x V(k).

***********************************************************************************************/


inline void SFF_PENCIL::Compute_Modal_vorticity
(
	int lx, int ly, int lz, 
	Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, 
	TinyVector<Complex,3> &vorticity
)
{

	TinyVector<Real,3>  K;
	TinyVector<Complex,3> Vi;
	
	Vi = Ax(ly, lz, lx), Ay(ly, lz, lx), Az(ly, lz, lx);

	Wavenumber(lx, ly, lz,  K);
	
	vorticity(0) = I *  ( K(1)*Vi(2) -  K(2)*Vi(1));
	vorticity(1) = I *  ( K(2)*Vi(0) -  K(0)*Vi(2));
	vorticity(2) = I *  ( K(0)*Vi(1) -  K(1)*Vi(0));
}




inline void SFF_PENCIL::Compute_Modal_vorticity_y_component
(
 int lx, int ly, int lz, 
 Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az, 
 Complex &vort_y
 )
{
	
	TinyVector<Real,3>  K;
	TinyVector<Complex,3> Vi;
	
	Vi = Ax(ly, lz, lx), 0, Az(ly, lz, lx);
	//Vi(1) set to zero to save time..
	
	Wavenumber(lx, ly, lz,  K);
	
	vort_y = I *  ( K(2)*Vi(0) -  K(0)*Vi(2));
}




//*********************************************************************************************	

			
/// Anisotropic axis along x1: for anisotropic energy spectrum and energy 
///		transfer calculations,  \f$ K_{||} = K_1 \f$			
inline Real SFF_PENCIL::AnisKpll(int lx, int ly, int lz)
{	
	if (global.field.anisotropy_dirn == 1)
		return (Get_kx(lx) * kfactor[1]); 
	
	else if (global.field.anisotropy_dirn == 2)
		return (Get_ky(ly) * kfactor[2]);
		
	else if (global.field.anisotropy_dirn == 3)
		return (Get_kz(lz) * kfactor[3]);
	
	else
		return 0;		// for -Wall
		
}

/// Anisotropic axis along x1: for anisotropic energy spectrum and energy 
///		transfer calculations,  \f$ K_\perp =\sqrt{K_2^2 + K_3^2} \f$			
inline Real SFF_PENCIL::AnisKperp(int lx, int ly, int lz)
{
	if (global.field.anisotropy_dirn == 1)
		return sqrt( pow2(Get_ky(ly) * kfactor[2]) + pow2(Get_kz(lz)* kfactor[3]) ); 
	
	else if (global.field.anisotropy_dirn == 2)
		return sqrt( pow2(Get_kx(lx) * kfactor[1]) + pow2(Get_kz(lz)* kfactor[3]) );
		
	else if (global.field.anisotropy_dirn == 3)
		return sqrt( pow2(Get_kx(lx) * kfactor[1]) 
						+ pow2(Get_ky(ly) * kfactor[2]) );
	else
		return 0;		// for -Wall
}

/// Anisotropic axis along x1: for anisotropic energy spectrum and energy 
///		transfer calculations,  horizontal direction 1, \f$ K_{h1} = K_2 \f$										
inline Real SFF_PENCIL::AnisKh1(int lx, int ly, int lz)
{	
	if (global.field.anisotropy_dirn == 1)
		return (Get_ky(ly) * kfactor[2]);  
	
	else if (global.field.anisotropy_dirn == 2)
		return (Get_kz(lz) * kfactor[3]);
		
	else if (global.field.anisotropy_dirn == 3)
		return (Get_kx(lx) * kfactor[1]);
	
	else
		return 0;		// for -Wall
}

/// Anisotropic axis along x1: for anisotropic energy spectrum and energy 
///		transfer calculations,  horizontal direction 2, \f$ K_{h2} = K_3 \f$				
inline Real SFF_PENCIL::AnisKh2(int lx, int ly, int lz)
{	
	if (global.field.anisotropy_dirn == 1)
		return (Get_kz(lz) * kfactor[3]);  
	
	else if (global.field.anisotropy_dirn == 2)
		return (Get_kx(lx)  * kfactor[1]);
		
	else if (global.field.anisotropy_dirn == 3)
		return (Get_ky(ly) * kfactor[2]); 
	
	else
		return 0;		// for -Wall
}
			
			
/// Cylindrical: Anis_min_Kpll
inline Real SFF_PENCIL::Anis_min_Kpll() 
{ 
	return 0.0;
}
				
/// Cylindrical: Anis_max_Kpll
inline Real SFF_PENCIL::Anis_max_Kpll() 
{ 
	Real maxkpll = 0;
		
    if (global.field.anisotropy_dirn == 1)
        maxkpll = ((Nx/2) * kfactor[1]); 
    
    else if (global.field.anisotropy_dirn == 2)
        maxkpll = ((Ny/2) * kfactor[2]); 
        
    else if (global.field.anisotropy_dirn == 3)
        maxkpll = ((Nz/2) * kfactor[3]);
		
	return maxkpll;
}

/// 3D Cylindrical: Anis_max_Krho_radius_inside the wavenumber box
inline int SFF_PENCIL::Anis_max_Krho_radius_inside() 	
{
	Real Kmag = 0.0;
	
    if (global.field.anisotropy_dirn == 1)
        Kmag = min( (Ny/2)*kfactor[2], (Nz/2)*kfactor[3] );
    
    else if (global.field.anisotropy_dirn == 2)
        Kmag = min( (Nx/2)*kfactor[1], (Nz/2)*kfactor[3] );  
    
    else if (global.field.anisotropy_dirn == 3)
        Kmag = min( (Nx/2)*kfactor[1], (Ny/2)*kfactor[2] ); 

    return ((int) Kmag);
}

// Max polar angle
inline Real SFF_PENCIL::Get_max_polar_angle() 
{ 
	
	return M_PI/2;
	
}			

//*********************************************************************************************	


/*! \brief Returns the angle K vector makes with the anisotropic axis 
 * 
 * The range of angle is \f$ [0:\pi] \f$.
 *
 * \param  lx, ly, lz (3D)
 * \return \f$ \tan^{-1}(K_{\perp}/K_{||}) \f$.
 * \return \f$ \pi/2 \f$ if \f$ K_{||} = 0 \f$.
 */	
inline Real SFF_PENCIL::AnisKvect_polar_angle(int lx, int ly, int lz)
{
	Real kkpll, kkperp;
	
	kkpll = AnisKpll(lx, ly, lz);
	kkperp = AnisKperp(lx, ly, lz);
	
	return Get_polar_angle(kkperp, kkpll);
}


/*! \brief 3D: Returns the azimutal angle.
 * 
 * The range of angle is \f$ [0:\pi] \f$.
 *
 * \param  lx, ly, lz (3D)
 * \return \f$ \tan^{-1}(Ky}/Kx \f$.
 * \return \f$ \pi/2 \f$ if \f$ K_{||} = 0 \f$.
 */	
inline Real SFF_PENCIL::AnisKvect_azimuthal_angle(int lx, int ly, int lz)
{
	
	Real kkh1 = AnisKh1(lx, ly, lz);
	Real kkh2 = AnisKh2(lx, ly, lz);
	
	return Get_azimuthal_angle(kkh1, kkh2);
}			


//********************************** END of four_inline.h *************************************


