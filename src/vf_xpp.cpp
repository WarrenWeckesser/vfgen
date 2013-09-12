
//
//  vf_xpp.cpp
//
//  This file defines the VectorField::PrintXPP method.
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
#include <cmath>

#include "vf.h"
#include "vf_utils.h"

using namespace std;
using namespace GiNaC;

//
// PrintXPP -- The XPP Code Generator.
//

void VectorField::PrintXPP(map<string,string> options)
    {
    int nc, np, nv, na, nf;
    int i;

    nc = conname_list.nops();
    nv = varname_list.nops();
    np = parname_list.nops();
    na = exprname_list.nops();
    nf = funcname_list.nops();

    string filename = Name()+".ode";
    ofstream fout;
    fout.open(filename.c_str());
    fout << "#" << endl;
    fout << "# " << filename << endl;
    fout << "#" << endl;
    fout << "# XPP ODE file for the vector field: " << Name() << endl;
    fout << "#" << endl;
    PrintVFGENComment(fout,"# ");
    fout << "#" << endl;
    //
    //  Print the constants.
    //
    for (i = 0; i < nc; ++i)
        {
        ex val = convalue_list[i];
        for (int j = i-1; j >= 0; j--)
            val = val.subs(conname_list[j]==convalue_list[j]);
        val = val.subs(Pi==M_PI);
        fout << "number " << conname_list[i] << "=" << val << endl;
        }
    fout << "#" << endl;
    //
    //  Print the parameters.
    //
    for (i = 0; i < np; ++i)
        {
        ex val = pardefval_list[i];
        for (int j = i-1; j >= 0; j--)
            val = val.subs(parname_list[j]==pardefval_list[j]);
        for (int j = nc-1; j >= 0; j--)
            val = val.subs(conname_list[j]==convalue_list[j]);
        val = val.subs(Pi==M_PI);            
        fout << "par " << parname_list[i] << "=" << val << endl;
        }
    fout << "#" << endl;
    //
    //  Print the intermediate/auxiliary formulas.
    //
    for (i = 0; i < na; ++i)
        {
        fout << exprname_list[i] << "=" << delay_transform(exprformula_list[i],varname_list)  << endl;
        }
    //
    //  Print the vector field expressions.
    //
    fout << "#" << endl;
    fout << "# The vector field and initial conditions" << endl;
    fout << "#" << endl;
    for (i = 0; i < nv; ++i)
        {
        fout << varname_list[i] << "'=" << delay_transform(varvecfield_list[i],varname_list) << endl;
        //
        // Replace Pi, Constants and Parameters in the initial condition expression.
        // Work backwards through the Parameters, then the Constants, then Pi.
        //
        ex ic = vardefic_list[i];
        for (int j = np-1; j >= 0; j--)
            ic = ic.subs(parname_list[j]==pardefval_list[j]);
        for (int j = nc-1; j >= 0; j--)
            ic = ic.subs(conname_list[j]==convalue_list[j]);
        ic = ic.subs(Pi==M_PI);        
        fout << "init " << varname_list[i] << "=" << ic << endl;
        if (IsDelay)
            fout << varname_list[i] << "(0)=" << vardefhist_list[i] << endl;
        }
    if (nf > 0)
        fout << "#" << endl;
    for (int i = 0; i < nf; ++i)
        {
        fout << "aux " << funcname_list[i] << "=" << funcformula_list[i] << endl;
        }
    fout << "#" << endl;
    if (options["extra"] != "")
        {
        string s = options["extra"];
        char delim = ';';
        string::size_type k = 0;
        while ((k = s.find(delim,k)) != string::npos)
            {
            s[k]='\n';
            }
        if (s[s.length()-1] != '\n')
            s = s + '\n';
        fout << "#--- Lines provided by the user with the \"extra\" option ---\n";
        fout << s;
        fout << "#----------------------------------------------------------\n";
        }
    fout << "done" << endl;
    fout.close();
    }
