
//
//  vf_latex.cpp
//
//  This file defines the VectorField::PrintLatex method.
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
// PrintLatex -- The Latex Code Generator.
//

void VectorField::PrintLatex(map<string,string> options)
    {
    int nc, np, nv, na, nf;
    int i;

    nc = conname_list.nops();
    nv = varname_list.nops();
    np = parname_list.nops();
    na = exprname_list.nops();
    nf = funcname_list.nops();

    string filename = Name()+".tex";
    ofstream fout;
    fout.open(filename.c_str());
    fout << "%" << endl;
    fout << "% " << filename << endl;
    fout << "%" << endl;
    fout << "% Latex file for the vector field: " << Name() << endl;
    fout << "%" << endl;
    PrintVFGENComment(fout,"% ");
    fout << "%" << endl;
    fout << "% It is assumed that this Latex fragment will be\n";
    fout << "% included in a larger document.  The \"split\" environment\n";
    fout << "% is used, which is defined in the AMS package \"amsmath\".\n";
    fout << "% The main Latex document must include the line\n";
    fout << "%    \\usepackage{amsmath}\n";
    fout << "%\n";
/*
    //
    //  Print the constants.
    //
    for (i = 0; i < nc; ++i)
        {
        fout << "number " << conname_list[i] << "=" << convalue_list[i] << endl;
        }
    fout << "%" << endl;
    //
    //  Print the parameters.
    //
    for (i = 0; i < np; ++i)
        {
        fout << "par " << parname_list[i] << "=" << pardefval_list[i] << endl;
        }
    fout << "%" << endl;
    //
    //  Print the intermediate/auxiliary formulas.
    //
    for (i = 0; i < na; ++i)
        {
        fout << exprname_list[i] << "=" << exprformula_list[i]  << endl;
        }
*/
    //
    //  Print the vector field expressions.
    //
    fout << "%" << endl;
    fout << "% The vector field\n";
    fout << "%" << endl;
    fout << GiNaC::latex ;
    fout << "\\begin{equation}\n";
    fout << "\\begin{split}\n";
    for (i = 0; i < nv; ++i)
        {
        fout << "  \\frac{d " << varname_list[i] << "}{d " << IndependentVariable << "} & =" << varvecfield_list[i] << "\\\\" << endl;
        }
    fout << "\\end{split}\n";
    fout << "\\label{eqn:" << Name() << "}\n";
    fout << "\\end{equation}\n";
    fout << "%" << endl;
    fout << dflt;
/*
    for (int i = 0; i < nf; ++i)
        {
        fout << "aux " << funcname_list[i] << "=" << funcformula_list[i] << endl;
        }
    fout << "#" << endl;
    fout << "done" << endl;
*/
    fout.close();
    }
