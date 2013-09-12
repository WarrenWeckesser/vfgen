//
//  vf_scilab.cpp
//
//  This file defines the VectorField::PrintScilab method.
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
// PrintScilab -- The Scilab Code Generator.
//

void VectorField::PrintScilab(map<string,string> options)
    {
    symbol t(IndependentVariable);
    int nc, np, nv, na, nf;

    nc = conname_list.nops();
    nv = varname_list.nops();
    np = parname_list.nops();
    na = exprname_list.nops();
    nf = funcname_list.nops();

    //
    // The file name is the vector field name with the extension ".sci"
    //
    string vf_filename = Name()+".sci";
    ofstream fout;
    fout.open(vf_filename.c_str());
    fout << left;
    //
    //  Create the vector field function.
    //
    fout << "//" << endl;
    fout << "// " << vf_filename << endl;
    fout << "//" << endl;
    fout << "// Scilab vector field functions for: " << Name() << endl;
    fout << "//" << endl;
    PrintVFGENComment(fout,"// ");
    fout << "//" << endl;
    fout << endl;
    fout << "//" << endl;
    fout << "// " << Name() << "_vf" << endl;
    fout << "//" << endl;
    fout << "// The vector field " << endl;
    fout << "//" << endl;
    fout << "function vf_ = " << Name() << "_vf(" << t << ",x_";
    if (np > 0)
        {
        fout << ",";
        if (options["parstyle"] == "list")
            PrintNameList(fout,parname_list);
        else
            fout << "p_";
        }
    fout << ")" << endl;
    if (HasPi)
        {
        fout << "    Pi = %pi;\n";
        }
    for (int i = 0; i < nc; ++i)
        {
        fout << "    " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    GetFromVector(fout,"    ",varname_list,"x_","()",1,";");
    if (options["parstyle"] != "list")
        GetFromVector(fout,"    ",parname_list,"p_","()",1,";");
    for (int i = 0; i < na; ++i)
        {
        fout << "    " << exprname_list[i] << " = " << exprformula_list[i] << ";" << endl;
        }
    fout << "    vf_ = zeros(" << nv << ",1);" << endl;
    for (int i = 0; i < nv; ++i)
        {
        fout << "    vf_(" << (i+1) << ")" << " = " << varvecfield_list[i] << ";" << endl;
        }
    fout << "endfunction" << endl;
    fout << endl;
    //
    //  Create the Jacobian function.
    //
    fout << "//" << endl;
    fout << "// " << Name() << "_jac" << endl;
    fout << "//" << endl;
    fout << "// The Jacobian of the vector field" << endl;
    fout << "//" << endl;
    fout << "function jac_ = " << Name() << "_jac(" << t << ",x_";
    if (np > 0)
        {
        fout << ",";
        if (options["parstyle"] == "list")
            PrintNameList(fout,parname_list);
        else
            fout << "p_";
        }
    fout << ")" << endl;
    if (HasPi)
        {
        fout << "    Pi = %pi;\n";
        }
    for (int i = 0; i < nc; ++i)
        {
        fout << "    " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    GetFromVector(fout,"    ",varname_list,"x_","()",1,";");
    if (options["parstyle"] != "list")
        GetFromVector(fout,"    ",parname_list,"p_","()",1,";");
    fout << "    jac_ = zeros(" << nv << "," << nv << ");" << endl;
    for (int i = 0; i < nv; ++i)
        {
        ex f = iterated_subs(varvecfield_list[i],expreqn_list);
        for (int j = 0; j < nv; ++j)
            {
            symbol v = ex_to<symbol>(varname_list[j]);
            ex df = f.diff(v);
            if (df != 0)
                fout << "    jac_(" << (i+1) << "," << (j+1) << ")" << " = " << f.diff(v) << ";" << endl;
            }
        }
    fout << "endfunction" << endl;
    fout << endl;

    if (options["func"] == "yes")
        {
        //
        //  Create the user-defined functions.
        //
        for (int n = 0; n < nf; ++n)
            {
            symbol fn = ex_to<symbol>(funcname_list[n]);
            string funcname = Name() + "_" + fn.get_name();
            fout << "//" << endl;
            fout << "// " << funcname << endl;
            fout << "//" << endl;
            fout << "// This function implements the user-defined function \"" << fn.get_name() << "\"" << endl;
            fout << "//" << endl;
            fout << "function r_ = " << funcname << "(" << t << ",x_";
            if (np > 0)
                {
                fout << ",";
                if (options["parstyle"] == "list")
                    PrintNameList(fout,parname_list);
                else
                    fout << "p_";
                }
            fout << ")" << endl;
            if (HasPi)
                {
                fout << "    Pi = %pi;\n";
                }
            for (int i = 0; i < nc; ++i)
                {
                fout << "    " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
                }
            GetFromVector(fout,"    ",varname_list,"x_","()",1,";");
            if (options["parstyle"] != "list")
                GetFromVector(fout,"    ",parname_list,"p_","()",1,";");
            for (int i = 0; i < na; ++i)
                {
                fout << "    " << exprname_list[i] << " = " << exprformula_list[i] << ";" << endl;
                }
            fout << "    r_ = " << funcformula_list[n] << ";" << endl;
            fout << "endfunction" << endl;
            }
        }
    fout.close();

    if (options["demo"] == "yes")
        {
        //
        //  Create a demonstration script.
        //
        string script_filename = Name()+"_demo.sce";
        fout.open(script_filename.c_str());
        fout << left;

        fout << "//" << endl;
        fout << "// " << script_filename << endl;
        fout << "//" << endl;
        fout << "// Scilab demonstration script that uses the vector field" << endl;
        fout << "// defined in " << vf_filename << endl;
        fout << "//" << endl;
        PrintVFGENComment(fout,"// ");
        fout << "//" << endl;
        fout << endl;
        fout << "// Load the vector field definition and the jacobian." << endl;
        fout << "exec('" << vf_filename << "');" << endl;
        fout << endl;
        // Output the constants; the initial conditions in
        // the x_mdialog can be expressed in terms of
        // these constants.
        if (nc > 0)
            fout << "// Constants.  The names of these constants can be used in the x_mdialog." << endl;
        if (HasPi)
            {
            fout << "Pi = %pi;\n";
            }
        for (int i = 0; i < nc; ++i)
            {
            fout << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
            }
        fout << endl;
        fout << "// Create data for an x_mdialog." << endl;
        fout << "tstr = 'Enter initial conditions, parameters, stop time, and number of samples:';" << endl;
        fout << "field_names = [";
        for (int i = 0; i < nv; ++i)
            {
            fout << "'" << varname_list[i] << "';";
            }
        for (int i = 0; i < np; ++i)
            {
            fout << "'" << parname_list[i] << "';";
            }
        fout << "'Stop Time';'Num Samples'];" << endl;
        fout << "default_field_values = [";
        for (int i = 0; i < nv; ++i)
            {
            fout << "'" << vardefic_list[i] << "';";
            }
        for (int i = 0; i < np; ++i)
            {
            fout << "'" << pardefval_list[i] << "';";
            }
        fout << "'10.0';'201'];" << endl;
        fout << "t0 = 0.0;" << endl;
        fout << "field_values = x_mdialog(tstr,field_names, default_field_values);" << endl;
        fout << "while (field_values ~= [])" << endl;
        fout << "    // Pull the data from the x_mdialog values. " << endl;
        fout << "    real_values = evstr(field_values);" << endl;
        fout << "    x0 = real_values(1:" << nv << ");" << endl;
        if (np > 0)
            {
            if (options["parstyle"] != "list")
                fout << "    params = real_values(" << nv+1 << ":" << nv+np << ");" << endl;
            else
                GetFromVector(fout,"    ",parname_list,"real_values","()",nv+1,";");
            }
        fout << "    tfinal = real_values(" << nv+np+1 << ");" << endl;
        fout << "    nsamples = real_values(" << nv+np+2 << ");" << endl;
        fout << "    tsamples = linspace(t0,tfinal,nsamples);" << endl;
        fout << endl;
        // if (options["parstyle"] != "vector")
        //    GetFromVector(fout,"    ",parname_list,"params","()",1,";");
        fout << "    // Call the ODE solver." << endl;
        fout << "    sol = ode(x0,t0,tsamples,list(" << Name() << "_vf";
        if (np > 0)
            {
            fout << ",";
            if (options["parstyle"] == "list")
                PrintNameList(fout,parname_list);
            else
                fout << "params";
            }
        fout << "),list(" << Name() << "_jac";
        if (np > 0)
            {
            fout << ",";
            if (options["parstyle"] == "list")
                PrintNameList(fout,parname_list);
            else
                fout << "params";
            }
        fout << "));" << endl;

        fout << endl;
        fout << "    // Plot the solution." << endl;
        fout << "    n = size(sol,2);" << endl;
        fout << "    clf;" << endl;

        for (int i = 0; i < nv; ++i)
            {
            fout << "    subplot(" << nv << ",1," << i+1 << ");\n";
            fout << "    plot(tsamples(1:n),sol(" << i+1 << ",:));\n";
            fout << "    ax = gca();\n";
            fout << "    ax.x_label.text = 't';\n";
            fout << "    ax.y_label.text = '" << varname_list[i] << "';\n";
            }
        fout << "    // Get another set of data from the user." << endl;
        fout << "    field_values = x_mdialog(tstr,field_names, field_values);" << endl;
        fout << "end;" << endl;
        fout.close();
        }
    }

