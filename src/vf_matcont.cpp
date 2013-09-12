
//
//  vf_matcon.cpp
//
//  This file defines the VectorField::PrintMATCONT method.
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
// PrintMATCONT -- The MATCONT Code Generator.
//

void VectorField::PrintMATCONT(map<string,string> options)
    {
    int nc, np, nv, na, nf;

    nc = conname_list.nops();
    nv = varname_list.nops();
    np = parname_list.nops();
    na = exprname_list.nops();
    nf = funcname_list.nops();

    string sep = "%--------------------------------------------------";

    //
    //  Print the main function.
    //
    string filename = Name()+".m";
    ofstream fout;
    fout.open(filename.c_str());
    fout << left;

    fout << "%" << endl;
    fout << "% " << filename << endl;
    fout << "%" << endl;
    fout << "% MATLAB file to be used with MATCONT." << endl;
    fout << "%" << endl;
    PrintVFGENComment(fout,"% ");
    fout << "%" << endl;
    fout << "%" << endl;
    fout << "function out = " << Name() << endl;
    fout << "    out{1} = @" << Name() << "_init;" << endl;
    fout << "    out{2} = @" << Name() << "_vf;"    << endl;
    fout << "    out{3} = @" << Name() << "_jac;"   << endl;
    fout << "    out{4} = @" << Name() << "_jacp;"  << endl;
    fout << "    out{5} = @" << Name() << "_hess;"  << endl;
    fout << "    out{6} = @" << Name() << "_hessp;" << endl;
    fout << "    out{7} = @" << Name() << "_der3;"  << endl;
    fout << "    out{8} = [];" << endl;
    fout << "    out{9} = [];" << endl;
    for (int i = 0; i < nf; ++i)
        {
        fout << "    out{" << 10+i << "} = @" << Name() << "_" << funcname_list[i] << ";" << endl;
        }
    fout << endl;

    //
    //  Print the init function.
    //
    fout << sep << endl;
    fout << endl;
    fout << "function [tspan, y0, options] = " << Name() << "_init\n";
    fout << endl;
    fout << "    tspan = [0; 10];\n";
    fout << endl;
    fout << "    % Get the constants and default parameters, because the default\n";
    fout << "    % initial conditions can depend on them.\n";
    if (HasPi)
        {
        fout << "    Pi = pi;\n";
        }
    for (int i = 0; i < nc; ++i)
        {
        fout << "    " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    for (int i = 0; i < np; ++i)
        fout << "    " << parname_list[i] << " = " << pardefval_list[i] << ";\n" ;
    fout << "    % Set the default initial point.\n";
    fout << "    y0 = zeros(" << nv << ",1);\n";
    for (int i = 0; i < nv; ++i)
        {
        fout << "    y0(" << i+1 << ") = " << vardefic_list[i] << ";\n";
        }
    fout << endl;
    fout << "    handles = feval(@" << Name() << ");\n";
    fout << "    options = odeset('Jacobian',handles(3),'JacobianP',handles(4), ...\n";
    fout << "                     'Hessians',handles(5),'HessiansP',handles(6), ...\n";
    fout << "                     'Der3',handles(7) );\n";
    fout << endl;
    //
    //  Print the vector field function.
    //
    fout << sep << endl;
    fout << "%" << endl;
    fout << "% The vector field" << endl;
    fout << "%" << endl;
    fout << "function vf_ = " << Name() << "_vf(" << IndependentVariable << ",x_,";
    PrintNameList(fout,parname_list);
    fout << ")" << endl;
    for (int i = 0; i < nc; ++i)
        {
        fout << "    " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    GetFromVector(fout,"    ",varname_list,"x_","()",1,";");
    fout << endl;
    for (int i = 0; i < na; ++i)
        {
        fout << "    " << exprname_list[i] << " = " << exprformula_list[i] << ";" << endl;
        }
    if (na > 0)
        fout << endl;
    fout << "    vf_ = zeros(" << nv << ",1);" << endl;
    for (int i = 0; i < nv; ++i)
        {
        fout << "    vf_(" << (i+1) << ")" << " = " << varvecfield_list[i] << ";" << endl;
        }
    fout << endl;
    //
    //  Print the Jacobian with respect to the variables.
    //
    fout << sep << endl;
    fout << "%" << endl;
    fout << "% The Jacobian of the vector field with respect to the variables" << endl;
    fout << "%" << endl;
    fout << "function jac_ = " << Name() << "_jac(" << IndependentVariable << ",x_,";
    PrintNameList(fout,parname_list);
    fout << ")" << endl;
    for (int i = 0; i < nc; ++i)
        {
        fout << "    " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    GetFromVector(fout,"    ",varname_list,"x_","()",1,";");
    fout << endl;
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
    fout << endl;
    //
    //  Print the Jacobian with respect to the parameters.
    //
    fout << sep << endl;
    fout << "%" << endl;
    fout << "% The Jacobian of the vector field with respect to the parameters" << endl;
    fout << "%" << endl;
    fout << "function jacp_ = " << Name() << "_jacp(" << IndependentVariable << ",x_,";
    PrintNameList(fout,parname_list);
    fout << ")" << endl;
    for (int i = 0; i < nc; ++i)
        {
        fout << "    " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    GetFromVector(fout,"    ",varname_list,"x_","()",1,";");
    fout << endl;
    fout << "    jacp_ = zeros(" << nv << "," << np << ");" << endl;
    for (int i = 0; i < nv; ++i)
        {
        ex f = iterated_subs(varvecfield_list[i],expreqn_list);
        for (int j = 0; j < np; ++j)
            {
            symbol p = ex_to<symbol>(parname_list[j]);
            ex df = f.diff(p);
            if (df != 0)
                fout << "    jacp_(" << (i+1) << "," << (j+1) << ")" << " = " << f.diff(p) << ";" << endl;
            }
        }
    fout << endl;
    //
    //  Print the Hessians function.
    //
    fout << sep << endl;
    fout << "%" << endl;
    fout << "% The Hessians function." << endl;
    fout << "% This function returns a 3D matrix." << endl;
    fout << "% hess_(n,:,:) is the Hessian of the n-th component of the vector field. That is," << endl;
    fout << "% hess_(n,i,j) is the second partial derivative of the n-th component" << endl;
    fout << "% of the vector field, taken with respect to the i-th and j-th variables." << endl;
    fout << "%" << endl;
    fout << "function hess_ = " << Name() << "_hess(" << IndependentVariable << ",x_,";
    PrintNameList(fout,parname_list);
    fout << ")" << endl;
    for (int i = 0; i < nc; ++i)
        {
        fout << "    " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    GetFromVector(fout,"    ",varname_list,"x_","()",1,";");
    fout << endl;
    fout << "    hess_ = zeros(" << nv << "," << nv << "," << nv << ");" << endl;
    for (int n = 0; n < nv; ++n)
        {
        fout << endl;
        // Get the n-th component of the vector field.
        ex f = iterated_subs(varvecfield_list[n],expreqn_list);
        for (int i = 0; i < nv; ++i)
            {
            // Get the i-th variable
            symbol x_i = ex_to<symbol>(varname_list[i]);
            // Differentiate f once
            ex df_i = f.diff(x_i);
            for (int j = i; j < nv; ++j)
                {
                // Get the j-th variable
                symbol x_j = ex_to<symbol>(varname_list[j]);
                // Differentiate again
                ex df_ij = df_i.diff(x_j);
                if (df_ij != 0)
                    { 
                    fout << "    hess_(" << (n+1) << "," << (i+1) << "," << (j+1) << ")" << " = " << df_ij << ";" << endl;
                    if (j > i)
                        fout << "    hess_(" << (n+1) << "," << (j+1) << "," << (i+1) << ") = hess_(" << (n+1) << "," << (i+1) << "," << (j+1) << ");" << endl;
                    }
                }
            }
        }
    fout << endl;
    //
    //  Print the Hessians with respect to parameters.
    //
    fout << sep << endl;
    fout << "%" << endl;
    fout << "% The Hessians with respect to the parameters." << endl;
    fout << "% This function returns a 3D matrix." << endl;
    fout << "% hessp_(n,i,j) is the second partial derivative of the n-th component" << endl;
    fout << "% of the vector field, taken with respect to the i-th variable and" << endl;
    fout << "% the j-th parameter." << endl;
    fout << "%" << endl;
    fout << "function hessp_ = " << Name() << "_hessp(" << IndependentVariable << ",x_,";
    PrintNameList(fout,parname_list);
    fout << ")" << endl;
    for (int i = 0; i < nc; ++i)
        {
        fout << "    " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    GetFromVector(fout,"    ",varname_list,"x_","()",1,";");
    fout << endl;
    fout << "    hessp_ = zeros(" << nv << "," << nv << "," << np << ");" << endl;
    for (int n = 0; n < nv; ++n)
        {
        fout << endl;
        // Get the n-th component of the vector field.
        ex f = iterated_subs(varvecfield_list[n],expreqn_list);
        for (int i = 0; i < nv; ++i)
            {
            // Get the i-th variable
            symbol x_i = ex_to<symbol>(varname_list[i]);
            // Differentiate f once
            ex df_i = f.diff(x_i);
            for (int j = 0; j < np; ++j)
                {
                // Get the j-th parameter
                symbol p_j = ex_to<symbol>(parname_list[j]);
                // Differentiate again
                ex df_ij = df_i.diff(p_j);
                if (df_ij != 0)
                    fout << "    hessp_(" << (n+1) << "," << (i+1) << "," << (j+1) << ")" << " = " << df_ij << ";" << endl;
                }
            }
        }
    fout << endl;
    //
    //  Print the function that computes the third derivatives
    //  with respect to the variables.
    //
    fout << sep << endl;
    fout << "%" << endl;
    fout << "% Third derivatives of the vector field." << endl;
    fout << "% This function returns a 4D matrix." << endl;
    fout << "% der3_(n,i,j,k) is the third partial derivative of the n-th component" << endl;
    fout << "% of the vector field, taken with respect to the i-th, j-th and k-th variables." << endl;
    fout << "%" << endl;
    fout << "function der3_ = " << Name() << "_der3(" << IndependentVariable << ",x_,";
    PrintNameList(fout,parname_list);
    fout << ")" << endl;
    for (int i = 0; i < nc; ++i)
        {
        fout << "    " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    GetFromVector(fout,"    ",varname_list,"x_","()",1,";");
    fout << endl;
    fout << "    der3_ = zeros(" << nv << "," << nv << "," << nv << "," << nv << ");" << endl;
    for (int n = 0; n < nv; ++n)
        {
        fout << endl;
        // Get the n-th component of the vector field.
        ex f = iterated_subs(varvecfield_list[n],expreqn_list);
        for (int i = 0; i < nv; ++i)
            {
            // Get the i-th variable
            symbol x_i = ex_to<symbol>(varname_list[i]);
            // Differentiate f once
            ex df_i = f.diff(x_i);
            for (int j = i; j < nv; ++j)
                {
                // Get the j-th variable
                symbol x_j = ex_to<symbol>(varname_list[j]);
                // Differentiate again
                ex df_ij = df_i.diff(x_j); 
                for (int k = j; k < nv; ++k)
                    {
                    // Get the k-th variable
                    symbol x_k = ex_to<symbol>(varname_list[k]);
                    // Take the third derivative
                    ex df_ijk = df_ij.diff(x_k);
                    if (df_ijk != 0)
                        {
                        fout << "    der3_(" << (n+1) << "," << (i+1) << "," << (j+1) << "," << (k+1) << ") = " << df_ijk << ";" << endl;
                        if (j < k)
                            fout << "    der3_(" << (n+1) << "," << (i+1) << "," << (k+1) << "," << (j+1) << ") = der3_(" << (n+1) << "," << (i+1) << "," << (j+1) << "," << (k+1) << ");" << endl;
                        if (i < k)
                            fout << "    der3_(" << (n+1) << "," << (k+1) << "," << (j+1) << "," << (i+1) << ") = der3_(" << (n+1) << "," << (i+1) << "," << (j+1) << "," << (k+1) << ");" << endl;
                        if (i < j)
                            fout << "    der3_(" << (n+1) << "," << (j+1) << "," << (i+1) << "," << (k+1) << ") = der3_(" << (n+1) << "," << (i+1) << "," << (j+1) << "," << (k+1) << ");" << endl;
                        }
                    }
                }
            }
        }
    fout << endl;

    //
    //  Create the user-defined functions.
    //
    for (int n = 0; n < nf; ++n)
        {
        symbol fn = ex_to<symbol>(funcname_list[n]);
        string funcname = fn.get_name();
        fout << left;
        fout << sep << endl;
        fout << endl;
        fout << "function r_ = " << Name() << "_" << funcname << "(" << IndependentVariable << ",x_,";
        PrintNameList(fout,parname_list);
        fout << ")" << endl;
        if (HasPi)
            {
            fout << "    Pi = pi;\n";
            }
        for (int i = 0; i < nc; ++i)
            {
            fout << "    " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
            }
        GetFromVector(fout,"    ",varname_list,"x_","()",1,";");
        for (int i = 0; i < na; ++i)
            {
            fout << "    " << exprname_list[i] << " = " << exprformula_list[i] << ";" << endl;
            }
        fout << "    r_ = " << funcformula_list[n] << ";" << endl;
        }

    fout.close();

    }
