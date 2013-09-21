
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

#include <cassert>
#include <ginac/ginac.h>

#include "ginac_declare_funcs.h"


using namespace GiNaC;

ex iterated_subs(ex f, lst e)
{
    ex g = f;
    int n = e.nops();
    for (int i = n-1; i >= 0; --i) {
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
    for (lst::const_iterator i = vars.begin(); i != vars.end(); ++i) {
        g = g.subs(*i == delay(*i, lag));
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
ex delay_transform(const ex& f, const lst& vars) {
    if (f.has( delay(wild(1),wild(2)) ) ) {
        exset dlist;
        f.find( delay(wild(1),wild(2)) , dlist);
        ex g = f;
        for (exset::const_iterator i = dlist.begin(); i != dlist.end(); ++i) {
            g = g.subs( *i == delay_vars(i->op(0),i->op(1),vars) );
        }
        return g;
    }
    else {
        return f;
    }
}


//
// Merge the expressions in exprs and formulas into a single expression.
// They are "merged" using the equality relation.  That is, if
// exprs = {e1, e2, e3} and formulas = {f1, f2}, the result is
// e1 == (e2 == (e3 == (f1 == f2)))
// Using == as the operator is just for convenience.  Any binary operation
// that does not occur in the expressions could have been used.
//
ex to_nested_tuple(const lst& exprs, const lst& formulas)
{
    lst all;
    int ne = exprs.nops();
    int nf = formulas.nops();

    for (lst::const_iterator iter = exprs.begin(); iter != exprs.end(); ++iter) {
        all.append(*iter);
    }
    for (lst::const_iterator iter = formulas.begin(); iter != formulas.end(); ++iter) {
        all.append(*iter);
    }

    int n = all.nops();
    assert(n != 0);

    if (n == 1) {
        return all.op(0);
    }
    ex t = all.op(n-2) == all.op(n-1);
    for (int k = n - 3; k >= 0; --k) {
        t = all.op(k) == t;
    }
    return t;
}


//
// Converts an expression created by `to_nested_tuple` into a list (GiNaC:lst)
// of expressions.
//

lst to_list(const ex& expr)
{
    lst exprs;
    if (!is_a<relational>(expr)) {
        exprs = lst(expr);
        return exprs;
    }
    ex e = expr;
    while (is_a<relational>(e)) {
        exprs.append(e.op(0));
        e = e.op(1);
    }
    exprs.append(e);
    return exprs;
}
