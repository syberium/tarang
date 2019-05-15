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



/*! \file  Rvf.cc
 * 
 * @brief  Class declaration of Rvf, a Real Vector Field 
 *
 * @sa	Rvf.h
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 *
 * @bug  No known bugs
 */
 
 

#include "Rvf.h"
#include "Cvf.h"

//*********************************************************************************************

RVF::RVF(string field_name)
{
	this->field_name=field_name;
	
	V1r.resize(shape_real_array);
	V2r.resize(shape_real_array);
	V3r.resize(shape_real_array);

}

void RVF::Set_field_name(string f) {
	field_name = f;
}

/**********************************************************************************************

				Forward transform: Inplace 
   
**********************************************************************************************/
/**********************************************************************************************

			Forward_transform(*Vir) = Vi 
				Vir unchanged.
   
**********************************************************************************************/


void RVF::Forward_transform(CVF cvf)
{
	global.program.sincostr_switch = sincostr_switch_Vx;
	universal->Forward_transform(V1r, cvf.V1);
	
	if (!global.program.two_dimension) {  // for 3d and 2.5d
		global.program.sincostr_switch = sincostr_switch_Vy;
		universal->Forward_transform(V2r, cvf.V2);
	}
	
	global.program.sincostr_switch = sincostr_switch_Vz;
	universal->Forward_transform(V3r, cvf.V3);
}


/**********************************************************************************************

		Inplace Inverse Fourier transform.
   
**********************************************************************************************/

/**********************************************************************************************

			Inverse_transform(Vi) = *Vir 
				Keeping Vi unchanged....
				temp = N1 x N2 x N3
   
**********************************************************************************************/

void RVF::Inverse_transform(CVF cvf)
{

	global.program.sincostr_switch = sincostr_switch_Vx;
	universal->Inverse_transform(cvf.V1, V1r);

	if (!global.program.two_dimension) { // for 3d and 2.5d
		global.program.sincostr_switch = sincostr_switch_Vy;
		universal->Inverse_transform(cvf.V2, V2r);

	}
	
	global.program.sincostr_switch = sincostr_switch_Vz;
	universal->Inverse_transform(cvf.V3, V3r);
}

void RVF::Read_real_field()
{
	V1r = 0;
	V2r = 0;
	V3r = 0;

	universal->Read(V1r, universal->H5_real, field_name+".V1r");

	if (!global.program.two_dimension)
		universal->Read(V2r, universal->H5_real, field_name+".V2r");
	
	universal->Read(V3r, universal->H5_real, field_name+".V3r");
}

// field_kind = Ur or Wr
void RVF::Write_real_field()
{
	string folder_name="real_" + To_string(global.time.now);

	universal->Write(V1r, universal->H5_real, "w", folder_name, field_name+".V1r");

	if (!global.program.two_dimension)
		universal->Write(V2r, universal->H5_real, "w", folder_name, field_name+".V2r");
	
	universal->Write(V3r, universal->H5_real, "w", folder_name, field_name+".V3r");
}


void RVF::Write_real_field_slice(unsigned int slice_file_counter)
{
	string file_name;

	for (vector<h5::Plan>::size_type i=0; i<universal->H5_slices.size(); i++) {
		file_name = "slice_" + To_string(i) + "_" + To_string(slice_file_counter);

		universal->Write(V1r, universal->H5_slices[i], "a", "slices", file_name, "V1r");

		if (!global.program.two_dimension)
			universal->Write(V2r, universal->H5_slices[i], "a", "slices", file_name, "V2r");

		universal->Write(V3r, universal->H5_slices[i], "a", "slices", file_name, "V3r");
	}
}

//**************************  End of RVF class definitions ************************************




