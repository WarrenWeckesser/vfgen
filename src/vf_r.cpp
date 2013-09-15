//
//  vf_r.cpp
//
//  This file defines the VectorField::PrintR method.
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


#include <fstream>
#include <string>
#include <ginac/ginac.h>

#include "vf.h"
#include "vf_utils.h"

using namespace std;
using namespace GiNaC;

//
// PrintR -- The R Code Generator.
//

void VectorField::PrintR(map<string,string> options)
{
    symbol t(IndependentVariable);
    int nc, np, nv, na, nf;

    nc = conname_list.nops();
    nv = varname_list.nops();
    np = parname_list.nops();
    na = exprname_list.nops();
    nf = funcname_list.nops();

    //
    // The file name is the vector field name with the extension ".R"
    //
    string vf_filename = Name()+".R";
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
    fout << "#" << endl;
    fout << "# " << Name() << "_vf" << endl;
    fout << "#" << endl;
    fout << "# The vector field " << endl;
    fout << "#" << endl;
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
    for (int i = 0; i < na; ++i) {
        fout << "    " << exprname_list[i] << " <- " << exprformula_list[i] << ";" << endl;
    }
    fout << "    vf_ <- vector(len = " << nv << ")" << endl;
    for (int i = 0; i < nv; ++i) {
        fout << "    vf_[" << (i+1) << "]" << " = " << varvecfield_list[i] << ";" << endl;
    }
    fout << "    return(list(vf_))\n";
    fout << "}" << endl;
    fout << endl;
    //
    //  Create the Jacobian function.
    //
    fout << "#" << endl;
    fout << "# " << Name() << "_jac" << endl;
    fout << "#" << endl;
    fout << "# The Jacobian of the vector field" << endl;
    fout << "#" << endl;
    fout << Name() << "_jac <- function(" << t << ", state, parameters) {\n";
    if (HasPi) {
        fout << "    Pi <- pi;\n";
    }
    for (int i = 0; i < nc; ++i) {
        fout << "    " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
    }
    GetFromVector(fout,"    ",varname_list, "<-", "state", "[]", 1, "");
    GetFromVector(fout,"    ",parname_list, "<-", "parameters", "[]", 1, "");
    fout << "    jac_ = matrix(nrow = " << nv << ", ncol = " << nv << ")" << endl;
    for (int i = 0; i < nv; ++i) {
        ex f = iterated_subs(varvecfield_list[i],expreqn_list);
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
                fout << "    Pi <- pi;\n";
            }
            for (int i = 0; i < nc; ++i) {
                fout << "    " << conname_list[i] << " <- " << convalue_list[i] << ";" << endl;
            }
            fout << "    " << t << " <- sol[, 1]\n";
            GetFromVector2(fout, "    ", varname_list, "<-", "sol", "[, ", "]", 2, "");
            GetFromVector(fout, "    ", parname_list, "<-", "parameters", "[]", 1, "");
            for (int i = 0; i < na; ++i) {
                fout << "    " << exprname_list[i] << " <- " << exprformula_list[i] << endl;
            }
            fout << "    r_ <- " << funcformula_list[n] << ";" << endl;
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
            fout << "    r_ <- " << funcformula_list[n] << ";" << endl;
            fout << "    return(r_)" << endl;
            fout << "}" << endl;
        }
    }
    fout.close();

    if (options["demo"] == "yes") {
        //
        //  Create a demonstration script.
        //
        string script_filename = Name()+"_demo.R";
        fout.open(script_filename.c_str());
        fout << left;

        fout << "#" << endl;
        fout << "# " << script_filename << endl;
        fout << "#" << endl;
        fout << "# R demonstration script that uses the vector field" << endl;
        fout << "# defined in " << vf_filename << endl;
        fout << "#" << endl;
        PrintVFGENComment(fout,"# ");
        fout << "#" << endl;
        fout << endl;
        fout << "library(deSolve)" << endl;
        fout << endl;
        fout << "# Load the vector field definition and the jacobian." << endl;
        fout << "source(\"" << vf_filename << "\")" << endl;
        fout << endl;
        if (nc > 0) {
            fout << "# Constants.  The names of these constants can be used in the x_mdialog." << endl;
        }
        if (HasPi) {
            fout << "Pi = pi;\n";
        }
        for (int i = 0; i < nc; ++i) {
            fout << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
        fout << endl;
        if (np > 0) {
            fout << "# --- Parameters ---\n";
            fout << "parameters = c(\n";
            for (int i = 0; i < np; ++i) {
                fout << "    " << parname_list[i] << " = " << pardefval_list[i];
                if (i < np - 1) {
                    fout << ",";
                }
                fout << "\n";
            }
            fout << ")\n\n";
        }
        fout << "# --- Initial conditions ---\n";
        fout << "state = c(\n";
        for (int i = 0; i < nv; ++i) {
            fout << "    " << varname_list[i] << " = " << vardefic_list[i];
            if (i < nv - 1) {
                fout << ",";
            }
            fout << "\n";
        }
        fout << ")\n\n";
        fout << endl;
        fout << "# --- Time values ---\n";
        fout << "times = seq(0, 10, by = 0.02)\n";
        fout << endl;
        fout << "# --- Call the ODE solver ---" << endl;
        fout << "sol = ode(y = state, times = times, func = " << Name() << ", parms = parameters,\n";
        fout << "          jactype = \"fullusr\", jacfunc = " << Name() << "_jac,\n";
        fout << "          atol = 1e-8, rtol = 1e-6)\n";
        fout << endl;
        fout << "# --- Plot the solution ---" << endl;
        fout << "par(mfcol = c(" << nv  << ", 1))" << endl;
        fout << "t <- sol[, \"time\"]" << endl;
        for (int i = 0; i < nv; ++i) {
            fout << "plot(t, sol[, \"" << varname_list[i] << "\"], type = \"l\", col = \"blue\",\n";
            fout << "     xlab = \"" << t << "\", ylab = \"" << varname_list[i] << "\")\n";
        }
        fout.close();
    }
}
