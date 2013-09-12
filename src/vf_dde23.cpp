
//
//  vf_dde23.cpp
//
//  This file defines the VectorField::PrintDDE23 method.
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

#include <fstream>
#include <sstream>
#include <string>
#include <cassert>
#include <ginac/ginac.h> 

#include "vf.h"
#include "vf_utils.h"

using namespace std;
using namespace GiNaC;

//
// DDE23_ConvertDelaysToZlags(ex& f)
//
// This function converts each subexpression of the form delay(delayexpr,del)
// in f to an expression in terms of Zlags_(i,j).
//

void VectorField::DDE23_ConvertDelaysToZlags(ex& f)
    {
    exset dlist;
    f.find(delay(wild(1),wild(2)),dlist);
    for (exset::const_iterator iter = dlist.begin(); iter != dlist.end(); ++iter)
        {
        ex delayfunc = *iter;
        ex delayexpr = delayfunc.op(0);
        lst vars = FindVarsInEx(delayexpr);
        ex del = delayfunc.op(1);
        int dindex = FindDelay(del);
        assert(dindex != -1);
        for (lst::const_iterator iter = vars.begin(); iter != vars.end(); ++iter)
            {
            int vindex = FindVar(ex_to<symbol>(*iter));
            delayexpr = delayexpr.subs(*iter == Zlags_(vindex+1,dindex+1));
            }
        f = f.subs(delayfunc == delayexpr);
        }
    }

//
// PrintDDE23 -- The DDE23 Matlab function code generator.
//

void VectorField::PrintDDE23(map<string,string> options)
    {
    unsigned nc, np, nv, na, nf;

    nc = conname_list.nops();
    nv = varname_list.nops();
    np = parname_list.nops();
    na = exprname_list.nops();
    nf = funcname_list.nops();

    //
    //  Create the vector field function.
    //
    string vf_filename = Name()+"_dde23.m";
    ofstream fout;
    fout.open(vf_filename.c_str());
    fout << left;

    fout << "%" << endl;
    fout << "% " << vf_filename << endl;
    fout << "%" << endl;
    fout << "% MATLAB vector field function for: " << Name() << endl;
    fout << "%" << endl;
    PrintVFGENComment(fout,"% ");
    fout << "%" << endl;
    fout << "%" << endl;
    // Include a list of the lags in the comments.
    fout << "% The lags are: {";
    for (unsigned k = 0; k < Delays.size(); ++k)
        {
        fout << Delays[k];
        if (k < Delays.size()-1)
            fout << ", "; 
        }
    fout << "}" << endl;
    // Function definition starts here.
    fout << "%" << endl;
    fout << "function vf_ = " << Name() << "_dde23(" << IndependentVariable << ",x_,Zlags_";
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
        fout << "    Pi = pi;\n";
        }
    //
    // Constants...
    //
    for (unsigned i = 0; i < nc; ++i)
        {
        fout << "    " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    fout << "    % State variables\n";
    GetFromVector(fout,"    ",varname_list,"x_","()",1,";");
    //
    // Parameters...
    //    If parstyle is vector, get the actual parameter variables from
    //    the vector argument.
    //
    if (options["parstyle"] != "list")
        {
        if (np > 0)
            {
            fout << "    % Parameters\n";
            GetFromVector(fout,"    ",parname_list,"p_","()",1,";");
            fout << endl;
            }
        }
    //
    // The following code assumes that the delays are single parameters,
    // and not mathematical expressions.
    //
    // Expressions...
    //
    for (unsigned i = 0; i < na; ++i)
        {
        ex f = exprformula_list[i];
        if (f.has(delay(wild(1),wild(2))))
            DDE23_ConvertDelaysToZlags(f);
        fout << "    " << exprname_list[i] << " = " << f << ";" << endl;
        }
    //
    // StateVariables...
    //
    fout << "    vf_ = zeros(" << nv << ",1);" << endl;
    for (unsigned i = 0; i < nv; ++i)
        {
        ex f = varvecfield_list[i];
        if (f.has(delay(wild(1),wild(2))))
            DDE23_ConvertDelaysToZlags(f);
        fout << "    vf_(" << (i+1) << ")" << " = " << f << ";" << endl;
        }
    fout << endl;
    fout.close();

    if (options["demo"] == "yes")
        {
        //
        // Create the demo function.
        //
        string filename = Name() + "_dde23_demo.m";
        ofstream fout;
        fout.open(filename.c_str());
        fout << left;

        fout << "%" << endl;
        fout << "% " << filename << endl;
        fout << "%" << endl;
        fout << "% MATLAB demo script for the vector field: " << Name() << endl;
        fout << "%" << endl;
        PrintVFGENComment(fout,"% ");
        fout << "%" << endl;
        fout << "%" << endl;
        fout << endl;

        fout << "function " + Name() + "_dde23_demo(stoptime)\n";
 
        for (unsigned k = 0; k < conname_list.nops(); ++k)
            {
            fout << "    " << conname_list[k] << " = " << convalue_list[k] << ";\n";
            }

        string parg = "";
        if (np > 0)
            {
            for (unsigned k = 0; k < parname_list.nops(); ++k)
                {
                fout << "    " << parname_list[k] << " = " << pardefval_list[k] << ";\n"; 
                }
            if (options["parstyle"] != "list")
                {
                // Create the code that creates the vector of parameters
                fout << "    p_ = zeros(" << np << ",1);\n";
                for (unsigned k = 0; k < np; ++k)
                    { 
                    fout << "    p_(" << k+1 << ") = " << parname_list[k] << ";\n";
                    }
                fout << endl;
                parg = ",p_";
                }
            else
                {
                ostringstream os;
                for (unsigned k = 0; k < np; ++k)
                    os << "," << parname_list[k];
                parg = os.str();
                }
            }

        fout << "    x = zeros(" << varname_list.nops() << ",1);\n";
        for (unsigned k = 0; k < varname_list.nops(); ++k)
            {
            fout << "    x(" << k+1 << ") = " << vardefic_list[k] << ";\n";
            }
        fout << endl;
        fout << "    lags = [";
        for (unsigned k = 0; k < Delays.size(); ++k)
            {
            fout << Delays[k];
            if (k < Delays.size()-1)
                fout << ", ";
            }
        fout << "];\n";
        fout << endl;
        fout << "    x0_ = [";
        for (unsigned k = 0; k < vardefic_list.nops(); ++k)
            {
            fout << vardefic_list[k];
            if (k < vardefic_list.nops()-1)
                fout << ", ";
            }
        fout << "];\n";
        fout << endl;
        fout << "    opts = ddeset('reltol',1e-8,'abstol',1e-11,'InitialY',x0_);\n";
        string pre1 = "(t_,y_,Z_)";
        string pre2 = "(t_)";
        if (np == 0)
            {
            pre1 = "";
            pre2 = "";
            }
        fout << "    sol = dde23(@" << pre1 << Name() << "_dde23";
        if (np > 0)
            fout << "(t_,y_,Z_" << parg << ")";
        fout << ",lags,@" << pre2 << Name() << "_history";
        if (np > 0)
            fout << "(t_" << parg << ")";
        fout << ",[0 stoptime],opts);\n";
        fout << endl;
        fout << "    num_plot_samples = 500;\n";
        fout << "    tint = linspace(0,stoptime,num_plot_samples);\n";
        fout << "    xint = deval(sol,tint);\n";
        fout << "    clf\n";
        fout << "    plot(tint,xint,'linewidth',2);\n";
        fout << "    grid on\n";
        fout << "    xlabel('" << IndependentVariable << "');\n";
        fout << "    legend(";
        for (unsigned i = 0; i < nv; ++i)
            {
            fout << "'" << varname_list[i] << "',";
            }
        fout << "'Location','Best')\n";
        fout << "end\n";
        fout << endl;
        fout << "function x_ = " + Name() + "_history(t" << parg << ")\n";
        for (unsigned i = 0; i < nc; ++i)
            {
            fout << "    " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
            }
        if (options["parstyle"] != "list")
            GetFromVector(fout,"    ",parname_list,"p_","()",1,";");
        for (unsigned k = 0; k < vardefhist_list.nops(); ++k)
            {
            fout << "    x_(" << k+1 << ") = " << vardefhist_list[k] << ";\n";
            }
        fout << "end\n";
        fout.close();
        }
    }
