
//
//  vf_dstool.cpp
//
//  This file defines the VectorField::PrintDSTOOL method.
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
#include <string>
#include <ginac/ginac.h>
 
#include "vf.h"
#include "vf_utils.h"

using namespace std;
using namespace GiNaC;

//
// PrintDSTOOL -- The DSTOOL Code Generator.
//

void VectorField::PrintDSTool(void)
    {
    int nc, np, nv, na, nf;

    nc = conname_list.nops();
    nv = varname_list.nops();
    np = parname_list.nops();
    na = exprname_list.nops();
    nf = funcname_list.nops();

    string filename = Name()+"_def.c";
    ofstream fout;
    fout.open(filename.c_str());
    fout << csrc << left;

    //
    //  Print the DSTOOL C file header information.
    //
    fout << "/*" << endl;
    fout << " *  " << filename << endl;
    fout << " *" << endl;
    fout << " *  DSTOOL C file for the vector field named: " << Name() << endl;
    fout << " *" << endl;
    PrintVFGENComment(fout," *  ");
    fout << " */" << endl;
    fout << endl;
    fout << "#include <model_headers.h>" << endl;
    //
    // I'm not sure if math.h is needed; for now, I'll keep it in there.
    //
    fout << "#include <math.h>" << endl;
    fout << endl;
    //
    //  Print the DSTOOL C vector field function.
    //
    fout << "/*" << endl;
    fout << " *  The vector field." << endl;
    fout << " */" << endl;
    fout << endl;
    //
    //  Note: Following the instructions in the user manual, the function
    //  will be declared with type int.  However, it does not return anything,
    //  so it should really be void.
    //
    fout << "int " << Name() << "(double *f_, double *Y_, double *p_)" << endl;
    fout << "    {" << endl;
    if (HasPi)
        {
        fout << "    const double Pi = M_PI;\n";
        }
    for (int i = 0; i < nc; ++i)
        {
        fout << "    const double " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    CDeclare_double(fout,varname_list);
    CDeclare_double(fout,parname_list);
    CDeclare_double(fout,exprname_list);
    fout << endl;
    GetFromVector(fout,"    ",varname_list,"Y_","[]",0,";");
    fout << endl;
    GetFromVector(fout,"    ",parname_list,"p_","[]",0,";");
    fout << endl;
    fout << "    double " << IndependentVariable << " = Y_[" << nv << "];" << endl;
    fout << endl;
    for (int i = 0; i < na; ++i)
        {
        fout << "    " << exprname_list[i] << " = " << exprformula_list[i] << ";" << endl;
        }
    fout << endl;
    for (int i = 0; i < nv; ++i)
        {
        fout << "    f_[" << i << "]" << " = " << varvecfield_list[i] << ";" << endl;
        }
    fout << "    }" << endl;
    fout << endl;
    //
    // Print the DSTOOL C function jacobian
    //
    fout << "/*" << endl;
    fout << " *  The Jacobian." << endl;
    fout << " */" << endl;
    fout << endl;
    fout << "int " << Name() << "_jac(double **jac_, double *Y_, double *p_)" << endl;
    fout << "    {" << endl;
    if (HasPi)
        {
        fout << "    const double Pi = M_PI;\n";
        }
    for (int i = 0; i < nc; ++i)
        {
        fout << "    const double " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    CDeclare_double(fout,varname_list);
    CDeclare_double(fout,parname_list);
    fout << endl;
    GetFromVector(fout,"    ",varname_list,"Y_","[]",0,";");
    fout << endl;
    GetFromVector(fout,"    ",parname_list,"p_","[]",0,";");
    fout << endl;
    fout << "    double " << IndependentVariable << " = Y_[" << nv << "];" << endl;
    fout << endl;
    for (int i = 0; i < nv; ++i)
        {
        ex f = varvecfield_list[i].subs(expreqn_list);
        // fout << "    /*  Row " << (i+1) << " of the Jacobian:      */"<< endl;
        for (int j = 0; j < nv; ++j)
            {
            symbol v = ex_to<symbol>(varname_list[j]);
            fout << "    jac_[" << i << "][" << j << "] = " << f.diff(v) << ";" << endl;
            }
        }
    fout << "    }" << endl;
    fout << endl;
    //
    // Print the DSTOOL C function dfdt
    //
    fout << "/*" << endl;
    fout << " *  The derivative with respect to the independent variable." << endl;
    fout << " */" << endl;
    fout << endl;
    fout << "int " << Name() << "_dfdt(double *dfdt_, double *Y_, double *p_)" << endl;
    fout << "    {" << endl;
    if (HasPi)
        {
        fout << "    const double Pi = M_PI;\n";
        }
    for (int i = 0; i < nc; ++i)
        {
        fout << "    const double " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    CDeclare_double(fout,varname_list);
    CDeclare_double(fout,parname_list);
    fout << endl;
    GetFromVector(fout,"    ",varname_list,"Y_","[]",0,";");
    fout << endl;
    GetFromVector(fout,"    ",parname_list,"p_","[]",0,";");
    fout << endl;
    fout << "    double " << IndependentVariable << " = Y_[" << nv << "];" << endl;
    fout << endl;
    symbol t(IndependentVariable);
    for (int i = 0; i < nv; ++i)
        {
        ex f = varvecfield_list[i].subs(expreqn_list);
        fout << "    dfdt_[" << i << "] = " << f.diff(t) << ";" << endl;
        }
    fout << "    }" << endl;
    fout << endl;

    //
    // Print the function that computes the Jacobian with respect
    // to the parameters.
    //
    fout << "/*" << endl;
    fout << " *  The Jacobian with respect to the parameters." << endl;
    fout << " */" << endl;
    fout << endl;
    fout << "int " << Name() << "_dfdp(double **dfdp_, double *Y_, double *p_)" << endl;
    fout << "    {" << endl;
    if (HasPi)
        {
        fout << "    const double Pi = M_PI;\n";
        }
    for (int i = 0; i < nc; ++i)
        {
        fout << "    const double " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    CDeclare_double(fout,varname_list);
    CDeclare_double(fout,parname_list);
    fout << endl;
    GetFromVector(fout,"    ",varname_list,"Y_","[]",0,";");
    fout << endl;
    GetFromVector(fout,"    ",parname_list,"p_","[]",0,";");
    fout << endl;
    for (int i = 0; i < nv; ++i)
        {
        ex f = iterated_subs(varvecfield_list[i],expreqn_list);
        for (int j = 0; j < np; ++j)
            {
            symbol p = ex_to<symbol>(parname_list[j]);
            fout << "    dfdp_[" << i << "][" << j << "] = " << f.diff(p) << ";" << endl;
            }
        }
    fout << "    }" << endl;
    fout << endl;
    //
    // Print the user-defined aux. functions.
    //
    fout << "/*" << endl;
    fout << " *  User-defined aux functions" << endl;
    fout << " */" << endl;
    fout << endl;
    fout << "int " << Name() << "_aux(double *f_, double *Y_, double *p_)" << endl;
    fout << "    {" << endl;
    if (HasPi)
        {
        fout << "    const double Pi = M_PI;\n";
        }
    for (int i = 0; i < nc; ++i)
        {
        fout << "    const double " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    CDeclare_double(fout,varname_list);
    CDeclare_double(fout,parname_list);
    fout << endl;
    GetFromVector(fout,"    ",varname_list,"Y_","[]",0,";");
    fout << endl;
    GetFromVector(fout,"    ",parname_list,"p_","[]",0,";");
    fout << endl;
    fout << "    double " << IndependentVariable << " = Y_[" << nv << "];" << endl;
    fout << endl;
    for (int i = 0; i < na; ++i)
        {
        fout << "    " << exprname_list[i] << " = " << exprformula_list[i] << ";" << endl;
        }
    for (int n = 0; n < nf; ++n)
        {
        fout << endl;
        fout << "    /*" << endl;
        fout << "     *  Aux function: " << funcname_list[n] << endl;
        fout << "     */" << endl;
        fout << "    f_[" << n << "] = " << funcformula_list[n] << ";" << endl;
        }
    fout << "    }" << endl;
    fout << endl;
    //
    //  Print the user_init function.
    //
    fout << "/*" << endl;
    fout << " *  The initialization function" << endl;
    fout << " */" << endl;
    fout << endl;
    fout << "int " << Name() << "_init()" << endl;
    fout << "    {" << endl;
    fout << "    /* ------------ Define the dynamical system in this segment ---------- */" << endl;
    fout << endl;
    if (HasPi)
        {
        fout << "    const double Pi = M_PI;\n";
        }
    for (int i = 0; i < nc; ++i)
        {
        fout << "    const double " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    fout << "    int             n_varb = " << nv << ";" << endl;
    fout << "    static char     *variable_names[] = {";
    for (int n = 0; n < nv; ++n)
        {
        if (n > 0)
            {
            fout << ",";
            }
        fout << "\"" << varname_list[n] << "\"";
        }
    fout << "};" << endl;
    fout << "    static double   variables[] = {";
    for (int n = 0; n < nv; ++n)
        {
        if (n > 0)
            {
            fout << ",";
            }
        fout << "0.0" ;
        }
    fout << "};" << endl;
    //
    // The VectorField object does not have max or min data,
    // so we'll fill in this data with -10 and 10 for variables
    // and parameters, and -100 and 100 for functions.
    //
    fout << "    static double   variable_min[] = {";
    for (int n = 0; n < nv; ++n)
        {
        if (n > 0)
            {
            fout << ",";
            }
        fout << "-10.0";
        }
    fout << "};" << endl;
    fout << "    static double   variable_max[] = {";
    for (int n = 0; n < nv; ++n)
        {
        if (n > 0)
            {
            fout << ",";
            }
        fout << "10.0";
        }
    fout << "};" << endl;
    fout << "    static char     *indep_varb_name = \"" << IndependentVariable << "\";" << endl;
    fout << "    static double   indep_varb_min = 0.0;" << endl;
    fout << "    static double   indep_varb_max = 10000.0;" << endl;
    fout << endl;
    fout << "    int             n_param = " << np << ";" << endl;
    fout << "    static char     *parameter_names[] = {";
    for (int n = 0; n < np; ++n)
        {
        if (n > 0)
            {
            fout << ",";
            }
        fout << "\"" << parname_list[n] << "\"";
        }
    fout << "};" << endl;
    fout << "    static double   parameters[] = {";
    for (int n = 0; n < np; ++n)
        {
        if (n > 0)
            {
            fout << ",";
            }
        fout << pardefval_list[n];
        }
    fout << "};" << endl;
    fout << "    static double   parameter_min[] = {";
    for (int n = 0; n < np; ++n)
        {
        if (n > 0)
            {
            fout << ",";
            }
        fout << "-10.0";
        }
    fout << "};" << endl;
    fout << "    static double   parameter_max[] = {";
    for (int n = 0; n < np; ++n)
        {
        if (n > 0)
            {
            fout << ",";
            }
        fout << "10.0";
        }
    fout << "};" << endl;
    fout << endl;
    fout << "    int             n_funct = " << nf << ";" << endl;
    fout << "    static char     *funct_names[] = {";
    for (int n = 0; n < nf; ++n)
        {
        if (n > 0)
            {
            fout << ",";
            }
        fout << "\"" << funcname_list[n] << "\"";
        }
    fout << "};" << endl;
    fout << "    static double   *funct_min[] = {";
    for (int n = 0; n < nf; ++n)
        {
        if (n > 0)
            {
            fout << ",";
            }
        fout << "-100.0";
        }
    fout << "};" << endl;
    fout << "    static double   *funct_max[] = {";
    for (int n = 0; n < nf; ++n)
        {
        if (n > 0)
            {
            fout << ",";
            }
        fout << "100.0";
        }
    fout << "};" << endl;
    fout << endl;
    //
    // If any StateVariable is periodic, the manifold_type is PERIODIC.
    //
    bool periodic = false;
    for (int n = 0; n < nv & !periodic; ++n)
        periodic = periodic | StateVariables[n]->IsPeriodic();
    fout << "    int             manifold_type = ";
    if (periodic)
        fout << "PERIODIC;" << endl;
    else
        fout << "EUCLIDEAN;" << endl;
    fout << "    static int      periodic_varb[] = {";
    for (int n = 0; n < nv; ++n)
        {
        if (n > 0)
            {
            fout << ",";
            }
        if (StateVariables[n]->IsPeriodic())
            fout << "TRUE";
        else
            fout << "FALSE";
        }
    fout << "};" << endl;
    fout << "    static double   period_start[] = {";
    for (int n = 0; n < nv; ++n)
        {
        if (n > 0)
            {
            fout << ",";
            }
        if (StateVariables[n]->IsPeriodic())
            fout << StateVariables[n]->PeriodicFrom();
        else
            fout << "0.0";
        }
    fout << "};" << endl;
    fout << "    static double   period_end[] = {";
    for (int n = 0; n < nv; ++n)
        {
        if (n > 0)
            {
            fout << ",";
            }
        if (StateVariables[n]->IsPeriodic())
            fout << StateVariables[n]->PeriodicTo();
        else
            fout << "0.0";
        }
    fout << "};" << endl;
    fout << endl;
    fout << "    int             mapping_toggle = FALSE;" << endl;
    fout << "    int             inverse_flag   = FALSE;" << endl;
    fout << endl;
    fout << "    int             (*def_name)() = " << Name() << ";" << endl;
    fout << "    int             (*jac_name)() = " << Name() << "_jac;" << endl;
    fout << "    int             (*aux_func_name)() = " << Name() << "_aux;" << endl;
    fout << "    int             (*inv_name)() = NULL;" << endl;
    fout << "    int             (*dfdt_name)() = " << Name() << "_dfdt;" << endl;
    fout << "    int             (*dfdparam_name)() = " << Name() << "_dfdp;" << endl;
    fout << endl;
    fout << "    c_filename = __FILE__;" << endl;
    fout << endl;
    fout << "    /* ------------ End of dynamical system definition          ---------- */" << endl;
    fout << "    #include <ds_define.c>" << endl;
    fout << "    return 0;" << endl;
    fout << "    }" << endl;
    fout.close();
    }

