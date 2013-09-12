
//
//  vf_evf.cpp
//
//  This file defines the VectorField::PrintEVF method.
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
#include <ginac/ginac.h> 

#include "vf.h"
#include "vf_utils.h"

using namespace std;
using namespace GiNaC;

int findpar(string parname, vector<Parameter *> pars)
    {
    for (unsigned i = 0; i < pars.size(); ++i)
        {
        if (pars[i]->Name() == parname)
            return (int) i;
        }
    return -1;
    }

//
// PrintEVF -- The EVF Code Generator.
//
// This function outputs a new vector field file to stdout.
// The new vector field is the original vector field extended
// with its variational equation.
//

void VectorField::PrintEVF(map<string,string> options)
    {
    //                                                            
    // This function adds variables to the vector field.
    // The vector field for the new variables is given by the
    // variational equations.
    // The code should probably check for a partial derivative being
    // zero and not bother generating the corresponding output.  Or it
    // could create a ginac ex to hold the expression until all derivatives
    // have been added to it, then create the string in the new state variable.
    // Ginac would automatically drop the zero terms.  I'll get around to
    // this some time...
    //
    // Note that this function MODIFIES THE OBJECT!!!
    // It then calls PrintXML to output the extended vector field to the
    // standard output in XML format.
    //
    int kpar = -1;
    symbol p;
    
    if (options.find("par") != options.end())
        {
        // cerr << "The option par=" << options["par"] << " has been given.\n";
        kpar = findpar(options["par"],Parameters);
        if (kpar == -1)
            {
            cerr << "Error: Unknown parameter \"" << options["par"] << "\"\n";
            cerr << "The parameters are: ";
            for (unsigned j = 0; j < parname_list.nops(); ++j)
                {
                if (j != 0)
                    cerr << ", ";
                cerr << parname_list[j];
                }
            cerr << endl;
            exit(-1);
            }
        }
    if (kpar != -1)
        {
        // cerr << "parname_list[" << kpar << "] = " << parname_list[kpar] << endl;
        p = ex_to<symbol>(parname_list[kpar]);
        }
    ostringstream oss;
    int nv = varname_list.nops();

    for (int i = 0; i < nv; ++i)
        {
        oss.str("");
        ex f = iterated_subs(varvecfield_list[i],expreqn_list);
        int num_terms_output = 0;
        for (int j = 0; j < nv; ++j)
            {
            symbol v = ex_to<symbol>(varname_list[j]);
            ex df = f.diff(v);
            if (df != 0)
                {
                if (num_terms_output > 0)
                    oss << " + ";
                if (df == 1)
                    oss << "d" << varname_list[j];
                else
                    oss << "(" << df << ")*d" << varname_list[j];
                num_terms_output = num_terms_output + 1;
                }
            }
        if (kpar != -1)
            {
            ex dfdp = f.diff(p);
            if (dfdp != 0)
                {
                if (num_terms_output > 0)
                    oss << " + ";
                oss << "(" << dfdp << ")";
                num_terms_output = num_terms_output + 1;
                }
            }
        if (num_terms_output == 0)
            oss << "0";
        StateVariable *var = new StateVariable("d"+StateVariables[i]->Name());
        var->Formula(oss.str());
        AddStateVariable(var);
        }
    Name(Name()+"_evf");
    PrintXML("evf");
    }
