
//
//  vf_adolc.cpp
//
//  This file defines the VectorField::ADOLC method.
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
// PrintADOLC -- The ADOLC Code Generator
//

void VectorField::PrintADOLC(map<string,string> options)
    {

    int nc, np, nv, na, nf;

    symbol t(IndependentVariable);
    nc = conname_list.nops();
    nv = varname_list.nops();
    np = parname_list.nops();
    na = exprname_list.nops();
    nf = funcname_list.nops();

    string filename = Name()+"_adolc.cpp";
    ofstream fout;
    fout.open(filename.c_str());
    fout << csrc << left;

    string pfilename = Name()+"_adolc.h";
    ofstream pout;
    pout.open(pfilename.c_str());
    pout << csrc << left;

    //
    //  Print C file header information.
    //
    fout << "//" << endl;
    fout << "//  " << filename << endl;
    fout << "//" << endl;
    fout << "//  ADOLC C++ file for the vector field named: " << Name() << endl;
    fout << "//" << endl;
    PrintVFGENComment(fout,"//  ");
    fout << "//" << endl;
    fout << endl;

    pout << "//" << endl;
    pout << "//  " << pfilename << endl;
    pout << "//" << endl;
    pout << "//  ADOLC C++ header file for the functions defined in " << filename << endl;
    pout << "//" << endl;
    PrintVFGENComment(pout,"//  ");
    pout << "//" << endl;
    pout << endl;

    fout << "#include <math.h>" << endl;
    fout << "#include \"adolc/adolc.h\"\n";
    fout << endl;
    //
    //  Print the vector field function.
    //
    fout << "//" << endl;
    fout << "//  The vector field." << endl;
    fout << "//" << endl;
    fout << endl;
    fout << "void " << Name() << "_vf(short int tag, double *y_, double *f_, double *params_)" << endl;
    pout << "void " << Name() << "_vf(short int, double *, double *, double *);" << endl;
    fout << "    {" << endl;
    fout << "    adouble ay_[" << nv << "];" << endl;
    fout << "    adouble af_[" << nv << "];" << endl;
    fout << endl;
    if (HasPi)
        {
        fout << "    const adouble Pi = M_PI;\n";
        }
    for (int i = 0; i < nc; ++i)
        {
        fout << "    const double " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    CDeclare(fout, "adouble", varname_list);
    CDeclare(fout, "adouble", parname_list);
    CDeclare(fout, "adouble", exprname_list);
    fout << endl;
    fout << "    trace_on(tag);" << endl;
    fout << "    for (int i = 0; i < " << nv << "; i++)" << endl;
        {
        fout << "        ay_[i] <<= y_[i];" << endl;
        }
    GetFromVector(fout,"    ", varname_list, "ay_", "[]", 0, ";");
    fout << endl;
    GetFromVector(fout,"    ", parname_list, "params_", "[]", 0, ";");
    fout << endl;
    for (int i = 0; i < na; ++i)
        {
        fout << "    " << exprname_list[i] << " = " << exprformula_list[i] << ";" << endl;
        }
    if (na > 0)
        fout << endl;
    for (int i = 0; i < nv; ++i)
        {
        fout << "    af_[" << i << "] = " << varvecfield_list[i] << ";" << endl;
        }
    fout << "    for (int i = 0; i < " << nv << "; i++)" << endl;
        {
        fout << "        af_[i] >>= f_[i];" << endl;
        }        
    fout << "    trace_off(tag);" << endl;
    fout << "    }" << endl;
    fout << endl;
    }

