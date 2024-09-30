//
//  vf_r.cpp
//
//  This file defines the VectorField::PrintRode and
//  VectorField::PrintRdede methods.
//
//
//  Copyright (C) 2013 Warren Weckesser
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
#include <fstream>
#include <sstream>
#include <string>
#include <ginac/ginac.h>

#include "vf.h"
#include "vf_utils.h"

using namespace std;
using namespace GiNaC;


static void
vf_header_comment(ofstream& fout, string name, string indvar)
{
    fout << "#" << endl;
    fout << "# " << name << "(" << indvar << ", state, parameters)" << endl;
    fout << "#" << endl;
    fout << "# The vector field function" << endl;
    fout << "#" << endl;
    return;
}

static void
jac_header_comment(ofstream& fout, string name, string indvar)
{
    fout << "#" << endl;
    fout << "# " << name << "_jac(" << indvar << ", state, parameters)" << endl;
    fout << "#" << endl;
    fout << "# The jacobian function" << endl;
    fout << "#" << endl;
    return;
}

static void
generate_demo_common(ofstream& out, string &vf_filename, bool HasPi,
                     lst& conname_list, lst& convalue_list,
                     lst& parname_list, lst& pardefval_list,
                     lst& varname_list, lst& vardefic_list)
{
    int nc = conname_list.nops();
    int nv = varname_list.nops();
    int np = parname_list.nops();

    out << "#" << endl;
    PrintVFGENComment(out, "# ");
    out << "#" << endl;
    out << endl;
    out << "library(deSolve)" << endl;
    out << endl;
    out << "# Load the vector field definition and the jacobian." << endl;
    out << "source(\"" << vf_filename << "\")" << endl;
    out << endl;
    if (nc > 0) {
        out << "# Constants.  The names of these constants can be used in the x_mdialog." << endl;
    }
    if (HasPi) {
        out << "Pi = pi\n";
    }
    AssignNameValueLists(out, "", conname_list, "<-", convalue_list, "");
    out << endl;
    out << "# Default parameter values" << endl;
    AssignNameValueLists(out, "", parname_list, "<-", pardefval_list, "");
    out << endl;
    if (np > 0) {
        out << "# Parameters list\n";
        out << "parameters = c(\n";
        for (int i = 0; i < np; ++i) {
            out << "    " << parname_list[i] << " = " << parname_list[i];
            if (i < np - 1) {
                out << ",";
            }
            out << "\n";
        }
        out << ")\n\n";
    }
    out << "# Initial conditions\n";
    out << "state = c(\n";
    for (int i = 0; i < nv; ++i) {
        out << "    " << varname_list[i] << " = " << vardefic_list[i];
        if (i < nv - 1) {
            out << ",";
        }
        out << "\n";
    }
    out << ")\n";
    out << endl;
    out << "# Time values\n";
    out << "times = seq(0, 10, by = 0.02)\n";
    out << endl;
}

//
// PrintRode -- The R Code Generator
//    Generate code for the 'ode' function of the 'deSolve' package.
//

void VectorField::PrintRode(map<string, string> options)
{
    symbol t(IndependentVariable);
    int nv, na, nf;

    nv = varname_list.nops();
    na = exprname_list.nops();
    nf = funcname_list.nops();

    //
    // The file name is the vector field name with the extension ".R"
    //
    string vf_filename = Name() + ".R";
    ofstream fout;
    fout.open(vf_filename.c_str());
    fout << left;
    //
    //  Create the vector field function.
    //
    fout << "#" << endl;
    fout << "# " << vf_filename << endl;
    fout << "#" << endl;
    fout << "# R vector field functions for: " << Name() << endl;
    fout << "#" << endl;
    PrintVFGENComment(fout, "# ");
    fout << "#" << endl;
    fout << endl;
    fout << endl;
    vf_header_comment(fout, Name(), IndependentVariable);

    fout << Name() << " <- function(" << t << ", state, parameters) {\n";
    string pname = Name() + "_parameters";

    if (HasPi) {
        fout << "    Pi <- pi\n";
    }
    AssignNameValueLists(fout, "    ", conname_list, "<-", convalue_list, "");
    GetFromVector(fout, "    ", varname_list, "<-", "state", "[]", 1, "");
    GetFromVector(fout, "    ", parname_list, "<-", "parameters", "[]", 1, "");
    fout << endl;
    for (int i = 0; i < na; ++i) {
        fout << "    " << exprname_list[i] << " <- " << exprformula_list[i] << endl;
    }
    fout << "    vf_ <- vector(len = " << nv << ")" << endl;
    for (int i = 0; i < nv; ++i) {
        fout << "    vf_[" << (i+1) << "]" << " = " << varvecfield_list[i] << endl;
    }
    fout << "    return(list(vf_))\n";
    fout << "}" << endl;
    fout << endl;
    //
    //  Create the Jacobian function.
    //
    jac_header_comment(fout, Name(), IndependentVariable);

    fout << Name() << "_jac <- function(" << t << ", state, parameters) {\n";
    if (HasPi) {
        fout << "    Pi <- pi\n";
    }
    AssignNameValueLists(fout, "    ", conname_list, "<-", convalue_list, "");
    GetFromVector(fout,"    ",varname_list, "<-", "state", "[]", 1, "");
    GetFromVector(fout,"    ",parname_list, "<-", "parameters", "[]", 1, "");
    fout << "    jac_ = matrix(nrow = " << nv << ", ncol = " << nv << ")" << endl;
    for (int i = 0; i < nv; ++i) {
        ex f = iterated_subs(varvecfield_list[i], expreqn_list);
        for (int j = 0; j < nv; ++j) {
            symbol v = ex_to<symbol>(varname_list[j]);
            ex df = f.diff(v);
            fout << "    jac_[" << (i+1) << "," << (j+1) << "]" << " = " << f.diff(v) << endl;
        }
    }
    fout << "    return(jac_)\n";
    fout << "}" << endl;
    fout << endl;

    if (options["func"] == "yes") {
        //
        //  Create the user-defined functions.
        //  Two R functions are created for each user-defined function.
        //  One has the signature
        //      <vectorfieldname>_<functionname>(t, state, parameters)
        //  and the other has the signature
        //      <vectorfieldname>_sol_<functionname>(sol, parameters)
        //
        for (int n = 0; n < nf; ++n) {
            fout << endl;
            symbol fn = ex_to<symbol>(funcname_list[n]);
            string funcname = Name() + "_sol_" + fn.get_name();
            fout << "#" << endl;
            fout << "# " << funcname << endl;
            fout << "#" << endl;
            fout << "# This function implements the user-defined function \"" << fn.get_name() << "\"" << endl;
            fout << "#" << endl;
            fout << funcname << " <- function(sol, parameters) {\n";
            if (HasPi) {
                fout << "    Pi <- pi\n";
            }
            AssignNameValueLists(fout, "    ", conname_list, "<-", convalue_list, "");
            fout << "    " << t << " <- sol[, 1]\n";
            GetFromVector2(fout, "    ", varname_list, "<-", "sol", "[, ", "]", 2, "");
            GetFromVector(fout, "    ", parname_list, "<-", "parameters", "[]", 1, "");
            for (int i = 0; i < na; ++i) {
                fout << "    " << exprname_list[i] << " <- " << exprformula_list[i] << endl;
            }
            fout << "    r_ <- " << funcformula_list[n] << endl;
            fout << "    return(r_)" << endl;
            fout << "}" << endl;
            fout << endl;

            string funcname2 = Name() + "_" + fn.get_name();
            fout << "#" << endl;
            fout << "# " << funcname2 << endl;
            fout << "#" << endl;
            fout << "# This function implements the user-defined function \"" << fn.get_name() << "\"" << endl;
            fout << "#" << endl;
            fout << funcname2 << " <- function(" << t << ", state, parameters) {\n";
            if (HasPi) {
                fout << "    Pi <- pi\n";
            }
            AssignNameValueLists(fout, "    ", conname_list, "<-", convalue_list, "");
            GetFromVector(fout, "    ", varname_list, "<-", "state", "[]", 1, "");
            GetFromVector(fout, "    ", parname_list, "<-", "parameters", "[]", 1, "");
            for (int i = 0; i < na; ++i) {
                fout << "    " << exprname_list[i] << " <- " << exprformula_list[i] << endl;
            }
            fout << "    r_ <- " << funcformula_list[n] << endl;
            fout << "    return(r_)" << endl;
            fout << "}" << endl;
        }
    }
    fout.close();

    if (options["demo"] == "yes") {
        //
        //  Create a demonstration script.
        //
        string script_filename = Name() + "_demo.R";
        fout.open(script_filename.c_str());
        fout << left;

        fout << "#" << endl;
        fout << "# " << script_filename << endl;
        fout << "#" << endl;
        fout << "# R demonstration script that uses the vector field" << endl;
        fout << "# defined in " << vf_filename << endl;
        fout << "#" << endl;
        generate_demo_common(fout, vf_filename, HasPi,
                             conname_list, convalue_list,
                             parname_list, pardefval_list,
                             varname_list, vardefic_list);
        fout << "# Call the ODE solver" << endl;
        fout << "sol = ode(y = state, times = times, func = " << Name() << ", parms = parameters,\n";
        fout << "          jactype = \"fullusr\", jacfunc = " << Name() << "_jac,\n";
        fout << "          atol = 1e-12, rtol = 1e-8)\n";
        fout << endl;
        fout << "# Plot the solution" << endl;
        fout << "par(mfcol = c(" << nv  << ", 1))" << endl;
        fout << "t <- sol[, \"time\"]" << endl;
        for (int i = 0; i < nv; ++i) {
            fout << "plot(t, sol[, \"" << varname_list[i] << "\"], type = \"l\", col = \"blue\",\n";
            fout << "     xlab = \"" << t << "\", ylab = \"" << varname_list[i] << "\")\n";
        }
        fout.close();
    }
}


//
//  On return, lags is a GiNac::lst, where each element is a GiNac::lst
//  of length 4 containing {lagsym, variable_index + 1, var, lag_time}
//
void VectorField::convert_delay_to_lagvalue(ex& f, lst &lags)
{
    symbol t(IndependentVariable);
    exset dlist;
    f.find(delay(wild(1), wild(2)), dlist);
    for (exset::const_iterator iter = dlist.begin(); iter != dlist.end(); ++iter) {
        ex delayfunc = *iter;
        ex delayexpr = delayfunc.op(0);
        lst vars = FindVarsInEx(delayexpr);
        ex del = delayfunc.op(1);
        for (lst::const_iterator iter = vars.begin(); iter != vars.end(); ++iter) {
            ostringstream os;
            lst tmp;
            os << lags.nops() + 1;
            symbol lagsym("lag" + os.str());
            int vindex = FindVar(ex_to<symbol>(*iter));
            delayexpr = delayexpr.subs(*iter == lagsym);
            tmp = {lagsym, vindex + 1, *iter, del};
            lags.append(tmp);
        }
        f = f.subs(delayfunc == delayexpr);
    }
}


static void generate_lag_assignment(ofstream& fout, const lst& lag,
        symbol& indvar, ex& history_formula)
{
    fout << "    if (" << indvar << " < " << lag.op(3) << ") {\n";
    fout << "        " << lag.op(0) << " = " << history_formula.subs(indvar == indvar - lag.op(3)) << endl;
    fout << "    }\n";
    fout << "    else {\n";
    fout << "        " << lag.op(0) << " = lagvalue(" << indvar - lag.op(3) <<  ", " <<  lag.op(1) << ")\n";
    fout << "    }\n";
}

//
// PrintRdede -- The R Code Generator
//     Generate code for use with the 'dede' function of the deSolve package.
//

void VectorField::PrintRdede(map<string, string> options)
{
    symbol t(IndependentVariable);
    int nc, nv, na, nf;

    nc = conname_list.nops();
    nv = varname_list.nops();
    na = exprname_list.nops();
    nf = funcname_list.nops();

    vector<lst> references;

    assert((na + nv) > 0);

    ex tuple_expr = to_nested_tuple(exprformula_list, varvecfield_list);
    lst tmp;
    convert_delay_to_lagvalue(tuple_expr, tmp);
    lst new_list = to_list(tuple_expr);
    exprformula_list = lst();
    for (int k = 0; k < na; ++k) {
        exprformula_list.append(new_list.op(k));
    }
    varvecfield_list = lst();
    for (int k = na; k < na + nv; ++k) {
        varvecfield_list.append(new_list.op(k));
    }
    int k = 0;
    for (lst::const_iterator iter = new_list.begin(); iter != new_list.end(); ++iter) {
        lst ref;
        ex e = *iter;
        for (lst::const_iterator p = tmp.begin(); p != tmp.end(); ++p) {
            lst lag_data = ex_to<lst>(*p);
            if (e.has(lag_data.op(0))) {
                ref.append(lag_data);
            }
        }
        ++k;
        references.push_back(ref);
    }

    //
    // The file name is the vector field name with the extension ".R"
    //
    string vf_filename = Name() + ".R";
    ofstream fout;
    fout.open(vf_filename.c_str());
    fout << left;

    fout << "#" << endl;
    fout << "# " << vf_filename << endl;
    fout << "#" << endl;
    fout << "# R vector field functions for: " << Name() << endl;
    fout << "#" << endl;
    PrintVFGENComment(fout, "# ");
    fout << "#" << endl;
    fout << endl;
    fout << endl;

    //
    //  Create the vector field function.
    //
    vf_header_comment(fout, Name(), IndependentVariable);

    fout << Name() << " <- function(" << t << ", state, parameters) {\n";
    string pname = Name() + "_parameters";

    if (HasPi) {
        fout << "    Pi <- pi\n";
    }
    for (int i = 0; i < nc; ++i) {
        fout << "        " << conname_list[i] << " <- " << convalue_list[i] << endl;
    }
    GetFromVector(fout, "    ", varname_list, "<-", "state", "[]", 1, "");
    GetFromVector(fout, "    ", parname_list, "<-", "parameters", "[]", 1, "");
    fout << endl;

    ex generated = 0;
    for (int i = 0; i < na; ++i) {
        for (lst::const_iterator iter = references[i].begin(); iter != references[i].end(); ++iter) {
            lst e = ex_to<lst>(*iter);
            if (!generated.has(e.op(0))) {
                generate_lag_assignment(fout, e, IndVar,
                    vardefhist_list[e.op(1)-1]);
                generated = generated + e.op(0);
            }
        }

        fout << "    " << exprname_list[i] << " <- " << exprformula_list[i] << endl;
    }

    for (int i = na; i < na + nv; ++i) {
        for (lst::const_iterator iter = references[i].begin(); iter != references[i].end(); ++iter) {
            lst e = ex_to<lst>(*iter);
            if (!generated.has(e.op(0))) {
                generate_lag_assignment(fout, e, IndVar,
                    vardefhist_list[e.op(1)-1]);
                generated = generated + e.op(0);
            }
        }
    }

    fout << "    vf_ <- vector(len = " << nv << ")" << endl;
    for (int i = 0; i < nv; ++i) {
        fout << "    vf_[" << (i+1) << "]" << " = " << varvecfield_list[i] << endl;
    }
    fout << "    return(list(vf_))\n";
    fout << "}" << endl;
    fout << endl;

    if (options["func"] == "yes") {
        //
        //  Create the user-defined functions.
        //  Two R functions are created for each user-defined function.
        //  One has the signature
        //      <vectorfieldname>_<functionname>(t, state, parameters)
        //  and the other has the signature
        //      <vectorfieldname>_sol_<functionname>(sol, parameters)
        //
        for (int n = 0; n < nf; ++n) {
            fout << endl;
            symbol fn = ex_to<symbol>(funcname_list[n]);
            string funcname = Name() + "_sol_" + fn.get_name();
            fout << "#" << endl;
            fout << "# " << funcname << endl;
            fout << "#" << endl;
            fout << "# This function implements the user-defined function \"" << fn.get_name() << "\"" << endl;
            fout << "#" << endl;
            fout << funcname << " <- function(sol, parameters) {\n";
            if (HasPi) {
                fout << "    Pi <- pi\n";
            }
            AssignNameValueLists(fout, "    ", conname_list, "<-", convalue_list, "");
            fout << "    " << t << " <- sol[, 1]\n";
            GetFromVector2(fout, "    ", varname_list, "<-", "sol", "[, ", "]", 2, "");
            GetFromVector(fout, "    ", parname_list, "<-", "parameters", "[]", 1, "");
            for (int i = 0; i < na; ++i) {
                fout << "    " << exprname_list[i] << " <- " << exprformula_list[i] << endl;
            }
            fout << "    r_ <- " << funcformula_list[n] << endl;
            fout << "    return(r_)" << endl;
            fout << "}" << endl;
            fout << endl;

            string funcname2 = Name() + "_" + fn.get_name();
            fout << "#" << endl;
            fout << "# " << funcname2 << endl;
            fout << "#" << endl;
            fout << "# This function implements the user-defined function \"" << fn.get_name() << "\"" << endl;
            fout << "#" << endl;
            fout << funcname2 << " <- function(" << t << ", state, parameters) {\n";
            if (HasPi) {
                fout << "    Pi <- pi;\n";
            }
            for (int i = 0; i < nc; ++i) {
                fout << "    " << conname_list[i] << " <- " << convalue_list[i] << ";" << endl;
            }
            GetFromVector(fout, "    ", varname_list, "<-", "state", "[]", 1, "");
            GetFromVector(fout, "    ", parname_list, "<-", "parameters", "[]", 1, "");
            for (int i = 0; i < na; ++i) {
                fout << "    " << exprname_list[i] << " <- " << exprformula_list[i] << endl;
            }
            fout << "    r_ <- " << funcformula_list[n] << endl;
            fout << "    return(r_)" << endl;
            fout << "}" << endl;
        }
    }
    fout.close();

    if (options["demo"] == "yes") {
        //
        //  Create a demonstration script.
        //
        string script_filename = Name() + "_demo.R";
        fout.open(script_filename.c_str());
        fout << left;

        fout << "#" << endl;
        fout << "# " << script_filename << endl;
        fout << "#" << endl;
        fout << "# R demonstration script that uses the vector field" << endl;
        fout << "# defined in " << vf_filename << endl;
        fout << "#" << endl;
        generate_demo_common(fout, vf_filename, HasPi,
                             conname_list, convalue_list,
                             parname_list, pardefval_list,
                             varname_list, vardefic_list);
        fout << "# Call the DDE solver" << endl;
        fout << "sol = dede(y = state, times = times, func = " << Name() << ", parms = parameters,\n";
        fout << "           atol = 1e-12, rtol = 1e-8)\n";
        fout << endl;
        fout << "# Plot the solution" << endl;
        fout << "par(mfcol = c(" << nv  << ", 1))" << endl;
        fout << "t <- sol[, \"time\"]" << endl;
        for (int i = 0; i < nv; ++i) {
            fout << "plot(t, sol[, \"" << varname_list[i] << "\"], type = \"l\", col = \"blue\",\n";
            fout << "     xlab = \"" << t << "\", ylab = \"" << varname_list[i] << "\")\n";
        }
        fout.close();
    }
}
