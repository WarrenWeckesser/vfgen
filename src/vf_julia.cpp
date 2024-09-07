
//
//  vf_julia.cpp
//
//  This file defines the VectorField::Julia method.
//
//
//  Copyright (C) 2024 Warren Weckesser
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
// PrintJulia -- The Julia Code Generator.
//

void VectorField::PrintJulia(map<string,string> options)
{
    int nv, np;   // currently na and nf are not used.

    nv = varname_list.nops();
    np = parname_list.nops();
    // na = exprname_list.nops();
    // nf = funcname_list.nops();

    string filename = Name() + ".jl";
    ofstream fout;
    fout.open(filename.c_str());
    fout << csrc << left;

    //
    //  Print Python file header information.
    //
    fout << "#" << endl;
    fout << "# " << filename << endl;
    fout << "#" << endl;
    fout << "# Function for the vector field named: " << Name() << endl;
    fout << "#\n";
    PrintVFGENComment(fout,"# ");
    fout << "#\n";

    //
    //  Print the vector field function.
    //
    fout << endl;
    fout << "function " << Name() << "!(du_, u_, p_, " << IndependentVariable << ")" << endl;
    if (HasPi) {
        fout << "    Pi = pi\n";
    }
    AssignNameValueLists(fout, "    ", conname_list, "=", convalue_list, "");
    GetFromVector(fout, "    ", varname_list, "=", "u_", "[]", 1, "");
    fout << endl;
    if (np > 0) {
        GetFromVector(fout, "    ", parname_list, "=", "p_", "[]", 1, "");
        fout << endl;
    }
    AssignNameValueLists(fout, "    ", exprname_list, "=", exprformula_list, "");
    for (int i = 0; i < nv; ++i) {
        fout << "    du_[" << i+1 << "] = " << varvecfield_list[i] << endl;
    }
    fout << "end" << endl;
    fout.close();

    if (options["demo"] == "yes") {
        //
        //  Create a self-contained ODE solver for this vector field.
        //

        string tfilename = Name()+"_demo.jl";
        ofstream tout;
        tout.open(tfilename.c_str());
        tout << csrc << left;
        tout << "#" << endl;
        tout << "# " << tfilename << endl;
        tout << "#" << endl;
        tout << "# Demo program for the vector field named: " << Name() << endl;
        tout << "#" << endl;
        PrintVFGENComment(tout,"# ");
        tout << "#\n" ;
        tout << endl;
        tout << "using DifferentialEquations" << endl;
        tout << "using Plots" << endl;
        tout << endl;
        tout << "include(\"" << Name() << ".jl\")" << endl;
        tout << endl << endl;

        if (HasPi) {
            tout << "Pi = pi\n";
        }
        AssignNameValueLists(tout, "", conname_list, "=", convalue_list, "");
        tout << "u0 = [";
        PrintList(tout, vardefic_list);
        tout << "]\n" ;
        if (np > 0) {
            tout << "p = [" ;
            PrintList(tout, pardefval_list);
            tout << "]\n" ;
        }
        tout << "tspan = (0.0, 20.0)" << endl;
        tout << endl;
        tout << "prob = ODEProblem(" << Name() << "!, u0, tspan";
        if (np > 0) {
            tout << ", p";
        }
        tout << ")" << endl; 
        tout << "sol = solve(prob)" << endl;
        tout << endl;
        tout << "ylabels = ";
        ListOfStrings(tout, varname_list, "[]", " ");
        tout << endl;
        tout << "plot(sol, ";
        tout <<      "xaxis=\"" << IndependentVariable << "\", ";
        tout <<      "ylabel=ylabels, ";
        tout <<      "layout=(" << nv << ", 1), ";
        tout <<      "legend=false)" << endl;
        tout.close();
    }
}
