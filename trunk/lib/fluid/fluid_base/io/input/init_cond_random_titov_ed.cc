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

/*! \file  init_cond_energy.ccinit_cond_double_para
 * 
 * @brief Construct energy spectrum given parameter values.  Given spectrum construct
 *			random initial condition (phases).
 *
 * @note 	Given  Model energy spectrum for initial condition
 *		Sk(k) = a k^4 exp(-b k^1.1) / (k^4 + q^4)^(2.8/12)
 *		with q = 1.5, b = 0.02
 *
 * @note:   Satisfy reality condition is critical here for kz=0 and N[3]/2 planes.   
 *			Do not remove this function.
 *
 *	Notation:  (Ki =) kki = ki * kfactor[i]
 *
 * @author  M. K. Verma
 * @version 4.0
 * @date Feb 2009
 *
 * @bug  No known bugs
 */


#include "FluidIO.h"


extern Uniform<Real> SPECrand;

//*********************************************************************************************

void  FluidIO::Init_cond_random_noise(FluidVF& U) {
    // Kinetic and Magnetic fields

    Real EuA = global.io.double_para(0), EuP = global.io.double_para(1); 

    Real kin_energy = EuA * pow(10.0, EuP);
    Real norm, Kmag;

    int kx, ky, kz;
    int px, py, pz;
    TinyVector<Complex,3> e, e2, k, u;

    for (int lx=0; lx<global.field.maxlx; lx++)
    for (int ly=0; ly<global.field.maxly; ly++)
    for (int lz=0; lz<global.field.maxlz; lz++) {
                    
        kx = universal->Get_kx(lx);
        ky = universal->Get_ky(ly);
        kz = universal->Get_kz(lz);

        k = kx, ky, kz;

        Kmag = universal->Kmagnitude(lx, ly, lz);

        Kmag=(Kmag==0?1:Kmag);
                    
        e = (2*SPECrand.random()-1) + I*(2*SPECrand.random()-1),
            (2*SPECrand.random()-1) + I*(2*SPECrand.random()-1),
            (2*SPECrand.random()-1) + I*(2*SPECrand.random()-1);

        e2 = mycross(e,k);

        norm = mynorm(e2);
        norm = (norm==0?1:norm);
        e2/=norm;

        u = 2 * pow(kin_energy,1.0/2.0) * e2 / Kmag;

        universal->Assign_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3, u);

    }  // of for loop
}

void  FluidIO::Init_cond_random_noise(FluidVF& U, FluidVF& W) {
    // Kinetic and Magnetic fields

	Real EuA = global.io.double_para(0), EuP = global.io.double_para(1); 
	Real EbA = global.io.double_para(2), EbP = global.io.double_para(3);

    Real kin_energy = EuA * pow(10.0, EuP);
    Real mag_energy = EbA * pow(10.0, EbP);
    Real norm, Kmag;

    int kx, ky, kz;
    int px, py, pz;
    TinyVector<Complex,3> e, e2, k, b, u;

    for (int lx=0; lx<global.field.maxlx; lx++)
    for (int ly=0; ly<global.field.maxly; ly++)
    for (int lz=0; lz<global.field.maxlz; lz++) {
                    
        kx = universal->Get_kx(lx);
        ky = universal->Get_ky(ly);
        kz = universal->Get_kz(lz);

        k = kx, ky, kz;

        Kmag = universal->Kmagnitude(lx, ly, lz);

        Kmag=(Kmag==0?1:Kmag);
                    
        e = (2*SPECrand.random()-1) + I*(2*SPECrand.random()-1),
            (2*SPECrand.random()-1) + I*(2*SPECrand.random()-1),
            (2*SPECrand.random()-1) + I*(2*SPECrand.random()-1);

        e2 = mycross(e,k);

        norm = mynorm(e2);
        norm = (norm==0?1:norm);
        e2/=norm;

        u = 2 * pow(kin_energy,1.0/2.0) * e2 / Kmag;

        e = (2*SPECrand.random()-1) + I*(2*SPECrand.random()-1),
            (2*SPECrand.random()-1) + I*(2*SPECrand.random()-1),
            (2*SPECrand.random()-1) + I*(2*SPECrand.random()-1);

        e2 = mycross(e,k);

        norm = mynorm(e2);
        norm = (norm==0?1:norm);
        e2/=norm;

        b = 2 * pow(mag_energy,1.0/2.0) * e2 / Kmag;

        universal->Assign_local_spectral_field(lx, ly, lz, U.cvf.V1, U.cvf.V2, U.cvf.V3, u);
        universal->Assign_local_spectral_field(lx, ly, lz, W.cvf.V1, W.cvf.V2, W.cvf.V3, b);

    }  // of for loop
}



//********************************** init_cond_random_titov_ed.cc **************************************



  
