
//
// ginac_aux_functions.cpp
//
//
// Utilities for working with GiNaC data.
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

#include <ginac/ginac.h>

#include "ginac_declare_funcs.h"


using namespace GiNaC;

ex iterated_subs(ex f, lst e)
    {
    ex g = f;
    int n = e.nops();
    for (int i = n-1; i >= 0; --i)
        {
        g = g.subs(e[i]);
        }
    return g;
    }

//
// Replace each var in e with delay(var,lag).
//
ex delay_vars(const ex& e, const ex& lag, const lst& vars)
    {
    ex g = e;
    for (lst::const_iterator i = vars.begin(); i != vars.end(); ++i)
        {
        g = g.subs(*i == delay(*i,lag));
        }
    return g;
    }

//
// Transform f so that the first operand in any delay function is
// a variable (not an expression).  For example, if
//    f = 1 + delay(a+b*x*sin(y),tau) - delay(x^2,2*delta)
// and vars = {x,y}, this function will return
//    1+ a+b*delay(x,tau)*sin(delay(y,tau)) - delay(x,2*delta)^2
//
ex delay_transform(const ex& f, const lst& vars)
    {
    if (f.has( delay(wild(1),wild(2)) ) )
        {
        exset dlist;
        f.find( delay(wild(1),wild(2)) , dlist);
        ex g = f;
        for (exset::const_iterator i = dlist.begin(); i != dlist.end(); ++i)
            {
            g = g.subs( *i == delay_vars(i->op(0),i->op(1),vars) );
            }
        return g;
        }
    else
        return f;
    }
