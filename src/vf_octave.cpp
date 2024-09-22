//
//  vf_octave.cpp
//
//  This file defines the VectorField::PrintOctave method.
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


#include <fstream>
#include <string>
#include <ginac/ginac.h>

#include "vf.h"
#include "vf_utils.h"

using namespace std;
using namespace GiNaC;

//
// PrintOctave -- The Octave Code Generator.
//

void VectorField::PrintOctave(map<string, string> options)
{
    symbol t(IndependentVariable);
    int nc, np, nv, nf;
    string pname;

    if (options["parstyle"] == "") {
        options["parstyle"] = "global";
    }
    if ((options["parstyle"] != "global") &&
        (options["parstyle"] != "list") &&
        (options["parstyle"] != "vector")) {
        cout << "Invalid value '" << options["parstyle"] << "' for the octave parstyle option." << endl;
        exit(-1);
    }
    if (options["demo"] == "") {
        options["demo"] = "no";
    }
    if ((options["demo"] != "no") &&
        (options["demo"] != "yes")) {
        cout << "Invalid value '" << options["demo"] << "' for the octave demo option." << endl;
        exit(-1);
    }

    nc = conname_list.nops();
    nv = varname_list.nops();
    np = parname_list.nops();
    nf = funcname_list.nops();

    //
    // The file name is the vector field name with the extension ".m"
    //
    string vf_filename = Name() + ".m";
    ofstream fout;
    fout.open(vf_filename.c_str());
    fout << left;
    //
    //  Create the vector field function.
    //
    fout << "#" << endl;
    fout << "# " << vf_filename << endl;
    fout << "#" << endl;
    fout << "# Octave vector field functions for: " << Name() << endl;
    fout << "#" << endl;
    PrintVFGENComment(fout, "# ");
    fout << "#" << endl;
    fout << endl;
    fout << "# The next line, containing just the number 1, is not a mistake!\n";
    fout << "# It allows the file to define more than one function.\n";
    fout << "1;\n";
    fout << endl;
    fout << "#" << endl;
    fout << "# " << Name() << "_vf" << endl;
    fout << "#" << endl;
    fout << "# The vector field " << endl;
    fout << "#" << endl;
    fout << "function vf_ = " << Name() << "_vf(x_, " << t;
    if (options["parstyle"] == "global") {
        pname = Name() + "_parameters";
    }
    else {
        pname = "p_";
    }
    if (np > 0) {
        if (options["parstyle"] == "vector") {
            fout << ", p_";
        }
        else if (options["parstyle"] == "list") {
            fout << ", ";
            PrintNameList(fout, parname_list);
        }
    }
    fout << ")" << endl;
    if ((np > 0) && (options["parstyle"] == "global")) {
        fout << "    global " << pname << ";\n";
    }
    if (HasPi) {
        fout << "    Pi = pi;\n";
    }
    AssignNameValueLists(fout, "    ", conname_list, "=", convalue_list, ";");
    GetFromVector(fout, "    ", varname_list, "=", "x_", "()", 1, ";");
    if (options["parstyle"] == "global") {
        GetFromVector(fout, "    ", parname_list, "=", pname.c_str(), "()", 1, ";");
    }
    else if (options["parstyle"] == "vector") {
        GetFromVector(fout, "    ", parname_list, "=", "p_", "()", 1, ";");
    }
    AssignNameValueLists(fout, "    ", exprname_list, "=", exprformula_list, ";");
    fout << "    vf_ = zeros(" << nv << ", 1);" << endl;
    for (int i = 0; i < nv; ++i) {
        fout << "    vf_(" << (i+1) << ")" << " = " << varvecfield_list[i] << ";" << endl;
    }
    fout << "endfunction" << endl;
    fout << endl;
    //
    //  Create the Jacobian function.
    //
    fout << "#" << endl;
    fout << "# " << Name() << "_jac" << endl;
    fout << "#" << endl;
    fout << "# The Jacobian of the vector field" << endl;
    fout << "#" << endl;
    fout << "function jac_ = " << Name() << "_jac(x_, " << t;
    if (np > 0) {
        if (options["parstyle"] == "vector") {
            fout << ", p_";
        }
        else if (options["parstyle"] == "list") {
            fout << ", ";
            PrintNameList(fout, parname_list);
        }
    }
    fout << ")" << endl;
    if ((np > 0) && (options["parstyle"] == "global")) {
        fout << "    global " << pname << ";\n";
    }
    if (HasPi) {
        fout << "    Pi = pi;\n";
    }
    AssignNameValueLists(fout, "    ", conname_list, "=", convalue_list, ";");
    GetFromVector(fout, "    ", varname_list, "=", "x_", "()", 1, ";");
    if (options["parstyle"] == "global") {
        GetFromVector(fout, "    ", parname_list, "=", pname.c_str(), "()", 1, ";");
    }
    else if (options["parstyle"] == "vector") {
        GetFromVector(fout, "    ", parname_list, "=", "p_", "()", 1, ";");
    }
    fout << "    jac_ = zeros(" << nv << ", " << nv << ");" << endl;
    for (int i = 0; i < nv; ++i) {
        ex f = iterated_subs(varvecfield_list[i], expreqn_list);
        for (int j = 0; j < nv; ++j) {
            symbol v = ex_to<symbol>(varname_list[j]);
            ex df = f.diff(v);
            if (df != 0) {
                fout << "    jac_(" << (i+1) << ", " << (j+1) << ")" << " = " << f.diff(v) << ";" << endl;
            }
        }
    }
    fout << "endfunction" << endl;
    fout << endl;

    if (options["func"] == "yes") {
        //
        //  Create the user-defined functions.
        //
        for (int n = 0; n < nf; ++n) {
            symbol fn = ex_to<symbol>(funcname_list[n]);
            string funcname = Name() + "_" + fn.get_name();
            fout << "#" << endl;
            fout << "# " << funcname << endl;
            fout << "#" << endl;
            fout << "# This function implements the user-defined function \"" << fn.get_name() << "\"" << endl;
            fout << "#" << endl;
            fout << "function r_ = " << funcname << "(x_, " << t;
            if (np > 0) {
                if ((options["parstyle"] == "vector") || (options["parstyle"] == "global")) {
                    fout << ", p_";
                }
                else if (options["parstyle"] == "list") {
                    fout << ", ";
                    PrintNameList(fout, parname_list);
                }
            }
            fout << ")" << endl;

            if (HasPi) {
                fout << "    Pi = pi;\n";
            }
            AssignNameValueLists(fout, "    ", conname_list, "=", convalue_list, ";");
            GetFromVector(fout, "    ", varname_list, "=", "x_", "()", 1, ";");
            if ((options["parstyle"] == "vector") || (options["parstyle"] == "global")) {
                GetFromVector(fout, "    ", parname_list, "=", "p_", "()", 1, ";");
            }
            AssignNameValueLists(fout, "    ", exprname_list, "=", exprformula_list, ";");
            fout << "    r_ = " << funcformula_list[n] << ";" << endl;
            fout << "endfunction" << endl;
        }
    }
    fout.close();

    if (options["demo"] == "yes") {
        //
        //  Create a demonstration script.
        //
        char colors[] = "bgrcm";
        string script_filename = Name() + "_demo.m";
        fout.open(script_filename.c_str());
        fout << left;

        fout << "#" << endl;
        fout << "# " << script_filename << endl;
        fout << "#" << endl;
        fout << "# Octave demonstration script that uses the vector field" << endl;
        fout << "# defined in " << vf_filename << endl;
        fout << "#" << endl;
        PrintVFGENComment(fout, "# ");
        fout << "#" << endl;
        fout << endl;
        if (np > 0) {
            if (options["parstyle"] == "global") {
                fout << "# Global vector for the parameters.\n";
                fout << "global " << pname << ";\n";
                fout << endl;
            }
        }
        fout << "# Load the vector field definition and the jacobian." << endl;
        fout << "source \"" << vf_filename << "\";" << endl;
        fout << endl;
        // Output the constants; the initial conditions in
        // the x_mdialog can be expressed in terms of
        // these constants.
        if (nc > 0) {
            fout << "# Constants.  The names of these constants can be used in the x_mdialog." << endl;
        }
        if (HasPi) {
            fout << "Pi = pi;\n";
        }
        AssignNameValueLists(fout, "", conname_list, "=", convalue_list, ";");
        fout << endl;
        if (np > 0) {
            fout << "# --- Parameters ---\n";
            if ((options["parstyle"] == "global") || (options["parstyle"] == "vector")) {
                for (int i = 0; i < np; ++i) {
                    fout << pname << "(" << i+1 << ") = " << pardefval_list[i]
                         << ";   # " << parname_list[i] << endl;
                }
            }
            else {
                // options["parstyle"] == "list"
                for (int i = 0; i < np; ++i) {
                    fout << parname_list[i] << " = " << pardefval_list[i] << ";" << endl;
                }
            }
            fout << endl;
        }

        fout << "# --- Initial conditions ---\n";
        SetVectorFromNames(fout, "", "x_", vardefic_list, "()", 1, ";");
        fout << endl;
        fout << "# --- Time values ---\n";
        fout << "t0_ = 0.0;" << endl;
        fout << "t1_ = 10.0;" << endl;
        fout << "numpoints = 201;" << endl;
        fout << "t_ = linspace(t0_, t1_, numpoints);\n";
        fout << endl;
        fout << "# --- Solver error tolerances ---\n";
        fout << "lsode_options(\"relative tolerance\", 1e-6);\n";
        fout << "lsode_options(\"absolute tolerance\", 1e-8);\n";
        fout << endl;
        fout << "# --- Call the ODE solver ---" << endl;
        if ((np == 0) || (options["parstyle"] == "global")) {
            fout << "xsol_ = lsode({\"" << Name() << "_vf\", \"" << Name() << "_jac\"}, x_, t_);\n";
        }
        else {
            string parargs;
            string name_vf = Name() + "_vf";
            string name_jac = Name() + "_jac";
            if (options["parstyle"] == "vector") {
                parargs = pname;
            } else {
                // options["parstyle"] == "list"
                parargs = Parameters[0]->Name();
                for (int i = 1; i < np; ++i) {
                    parargs += ", " + Parameters[i]->Name();
                }
            }
            fout << "xsol_ = lsode({@(x,t) " << name_vf << "(x, t, " << parargs << ")," <<
                                  " @(x,t) " << name_jac << "(x, t, " << parargs << ")}, x_, t_);\n";
        }
        fout << endl;
        fout << "# --- Plot the solution ---" << endl;
        fout << "args = argv();\n";
        fout << "fig = figure;\n";
        fout << "if length(args) == 1\n";
        fout << "    set(fig, \"visible\", \"off\")\n";
        fout << "end\n";
        for (int i = 0; i < nv; ++i) {
            int k = i % strlen(colors);
            fout << "plot(t_, xsol_(1:end, " << i+1 << "), \""
                 << colors[k] << ";" << varname_list[i] << ";\", "
                 << "\"linewidth\", 1)\n";
            if (i == 0) {
                fout << "hold on\n";
            }
        }
        fout << "grid on\n";
        fout << "xlabel(\"t\")\n";
        fout << "if length(args) == 1\n";
        fout << "    print(args{1})\n";
        fout << "end\n";
        fout.close();
    }
}
