
//
// ginac_aux_functions.h
//
//
//
//  Copyright (C) 2008 Warren Weckesser
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License, Version 2, as
//  published by the Free Software Foundation.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with this program.  The file LICENSE-GPL2 in the VFGEN source directory
//  contains a copy. You may also obtain a copy by writing to the Free Software
//  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
//


#ifndef GINAC_AUX_FUNCTIONS_H_INCLUDED
#define GINAC_AUX_FUNCTIONS_H_INCLUDED

#include <ginac/ginac.h>

#include "ginac_declare_funcs.h"



GiNaC::ex iterated_subs(GiNaC::ex f, GiNaC::lst e);
GiNaC::ex delay_vars(const GiNaC::ex& e, const GiNaC::ex& lag, const GiNaC::lst& vars);
GiNaC::ex delay_transform(const GiNaC::ex& f, const GiNaC::lst& vars);

#endif
