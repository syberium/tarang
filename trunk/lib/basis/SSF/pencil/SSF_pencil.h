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


/*! \file  universal_basic.h
 * 
 * @brief  Universal functions based on basis fns (MPI)
 *
 * @author  M. K. Verma
 * @version 4.0 MPI
 * @date Sept 2008
 * @bug	No known bugs
 */


#ifndef _SSF_PENCIL_H
#define _SSF_PENCIL_H

#include "def_vars.h"
#include "basicfn_inline.h"
#include "fftk.h"
#include "universal.h"
#include "ArrayOps.h"

//*********************************************************************************************	


class SSF_PENCIL:public Universal
{
public:
	
	SSF_PENCIL();
	
	#include "universal_fn_names.h"
	void Zero_modes(Array<Complex,3> Ax, Array<Complex,3> Ay, Array<Complex,3> Az);
	void Zero_modes(Array<Complex,3> F);
	
	void Assign_sub_array(Range y_range, Range z_range, Range x_range, Array<Complex,3> A, Complex value);
};

#endif


//******************************** End of scft-slab.h  **************************************


