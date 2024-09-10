
//
//  vf_scipy.cpp
//
//  This file defines the VectorField::SciPy method.
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

static void PrintArgDescription(ofstream &fout,lst varname_list, lst parname_list,
                                bool tfirst)
{
    fout << "    Arguments:\n";
    if (tfirst) {
        fout << "        t_ :  time\n";
    }
    fout << "        y_ :  vector of the state variables:\n";
    for (unsigned i = 0; i < varname_list.nops(); ++i) {
        fout << "                  y_[" << i << "] is " << varname_list[i] << endl;
    }
    if (!tfirst) {
        fout << "        t_ :  time\n";
    }
    if (parname_list.nops() > 0) {
        fout << "        p_ :  vector of the parameters\n";
        for (unsigned i = 0; i < parname_list.nops(); ++i) {
            fout << "                  p_[" << i << "] is " << parname_list[i] << endl;
        }
    }
}

//
// PrintSciPy -- The SciPy Code Generator.
//

void VectorField::PrintSciPy(map<string,string> options)
{
    int nv, np, na, nf;

    nv = varname_list.nops();
    np = parname_list.nops();
    na = exprname_list.nops();
    nf = funcname_list.nops();

    string filename = Name() + ".py";
    ofstream fout;
    fout.open(filename.c_str());
    fout << csrc << left;

    //
    //  Print Python file header information.
    //
    fout << "#" << endl;
    fout << "# " << filename << endl;
    fout << "#" << endl;
    fout << "# Python file for the vector field named: " << Name() << endl;
    fout << "#\n";
    fout << endl;
    fout << "\"\"\"\n";
    fout << "This module implements the vector field name \"" << Name() << "\" as" << endl;
    fout << "the function vectorfield().  The Jacobian of the vector field\n";
    fout << "is computed by jacobian().  These functions can be used with" << endl;
    fout << "the SciPy odeint function." << endl;
    if (options["func"] == "yes" && nf > 0) {
        fout << "This module also defines the function";
        if (nf > 1) {
            fout << "s";
        }
        for (int i = 0; i < nf; ++i) {
            if (i == 0) {
                fout << " ";
            }
            else if (nf > 1) {
                if (i == nf-1) {
                    fout << " and ";
                }
                else {
                    fout << ", ";
                }
            }
            fout << funcname_list[i] << "()";
        }
        fout << ".\n";
    }
    fout << endl;
    fout << "For example:\n";
    fout << endl;
    fout << "    from scipy.integrate import odeint\n";
    fout << "    import " << Name() << endl;
    fout << endl;
    if (np > 0) {
        fout << "    params = [";
        PrintNameList(fout,parname_list);
        fout << "]   # Assume the parameters have been set elsewhere\n";
    }
    fout << "    t = [i/10.0 for i in range(0, 101)]\n";
    fout << "    ic = ";
    for (int i = 0; i < nv; ++i) {
        fout << ((i == 0) ? "[1.0" : ",0.0");
    }
    fout << "]\n";
    fout << "    sol = odeint(" << Name() << ".vectorfield, ic, t,";
    if (options["tfirst"] == "yes") {
        fout << " tfirst=True,";
    }
    if (np > 0) {
        fout << " args=(params,),";
    }
    fout << " Dfun=" << Name() << ".jacobian)\n";
    fout << endl;
    PrintVFGENComment(fout,"");
    fout << endl;
    fout << "\"\"\"\n";
    fout << endl;
    fout << "from math import *" << endl;
    fout << "import numpy as np" << endl;
    fout << endl;
    //
    //  Print the vector field function.
    //
    fout << endl;
    fout << "def vectorfield(";
    if (options["tfirst"] == "yes") {
        fout << "t_, y_";
    }
    else {
        fout << "y_, t_";
    }

    if (np > 0) {
        fout << ", p_";
    }
    fout << "):" << endl;
    fout << "    \"\"\"\n";
    fout << "    The vector field function for the vector field \"" << Name() << "\"\n";
    PrintArgDescription(fout, varname_list, parname_list, options["tfirst"] == "yes");
    fout << "    \"\"\"\n";
    // fout << endl;
    if (HasPi) {
        fout << "    Pi = pi\n";
    }
    AssignNameValueLists(fout, "    ", conname_list, "=", convalue_list, "");
    GetFromVector(fout, "    ", varname_list, "=", "y_", "[]", 0, "");
    fout << endl;
    if (np > 0) {
        GetFromVector(fout, "    ", parname_list, "=", "p_", "[]", 0, "");
        fout << endl;
    }
    AssignNameValueLists(fout, "    ", exprname_list, "=", exprformula_list, "");
    if (na > 0) {
        fout << endl;
    }
    fout << "    f_ = np.zeros((" << nv << ",))" << endl;
    for (int i = 0; i < nv; ++i) {
        fout << "    f_[" << i << "] = " << varvecfield_list[i] << endl;
    }
    fout << endl;

    fout << "    return f_" << endl;
    fout << endl;

    //
    // Print the Jacobian function.
    //
    fout << endl;
    fout << "def jacobian(";
    if (options["tfirst"] == "yes") {
        fout << "t_, y_";
    }
    else {
        fout << "y_, t_";
    }
    if (np > 0) {
        fout << ", p_";
    }
    fout << "):" << endl;
    
    fout << "    \"\"\"\n";
    fout << "    The Jacobian of the vector field \"" << Name() << "\"\n";
    PrintArgDescription(fout, varname_list, parname_list, options["tfirst"] == "yes");
    fout << "    \"\"\"\n";
    if (HasPi) {
        fout << "    Pi = pi\n";
    }
    AssignNameValueLists(fout, "    ", conname_list, "=", convalue_list, "");
    GetFromVector(fout, "    ", varname_list, "=", "y_", "[]", 0, "");
    if (np > 0) {
        GetFromVector(fout, "    ", parname_list, "=", "p_", "[]", 0, "");
    }
    fout << endl;
    fout << "    # Create the Jacobian matrix:" << endl; 
    fout << "    jac_ = np.zeros((" << nv << ", " << nv << "))" << endl;
    for (int i = 0; i < nv; ++i) {
        ex f = iterated_subs(varvecfield_list[i],expreqn_list);
        for (int j = 0; j < nv; ++j) {
            symbol v = ex_to<symbol>(varname_list[j]);
            ex df = f.diff(v);
            if (df != 0) {
                fout << "    jac_[" << i << ", " << j << "] = " << df << endl;
            }
        }
    }
    fout << "    return jac_\n";

    if (options["func"] == "yes") {
        //
        // Print the user-defined functions.
        //
        for (int n = 0; n < nf; ++n) {
            fout << endl;
            fout << "#" << endl;
            fout << "# User function: " << funcname_list[n] << endl;
            fout << "#" << endl;
            fout << endl;
            fout << "def " << funcname_list[n];
            if (options["tfirst"] == "yes") {
                fout << "(t_, y_";
            }
            else {
                fout << "(y_, t_";
            }
            if (np > 0) {
                fout << ", p_";
            }
            fout << "):" << endl;
            fout << "    \"\"\"\n";
            fout << "    The user-defined function \"" << funcname_list[n] << "\" for the vector field \"" << Name() << "\"\n";
            PrintArgDescription(fout, varname_list, parname_list, options["tfirst"] == "yes");
            fout << "    \"\"\"\n";
            if (HasPi) {
                fout << "    Pi = pi\n";
            }
            AssignNameValueLists(fout, "    ", conname_list, "=", convalue_list, "");
            GetFromVector(fout, "    ", varname_list, "=", "y_", "[]", 0, "");
            fout << endl;
            if (np > 0) {
                GetFromVector(fout, "    ", parname_list, "=", "p_", "[]", 0, "");
                fout << endl;
            }
            AssignNameValueLists(fout, "    ", exprname_list, "=", exprformula_list, "");
            if (na > 0) {
                fout << endl;
            }
            fout << "    return " << funcformula_list[n] << "" << endl;
        }
    }

    fout.close();

    if (options["demo"] == "yes") {
        //
        //  Create a self-contained ODE solver for this vector field
        //  that allows the user to give the initial conditions,
        //  parameter values, and some solver control parameters
        //  on the command line.
        //

        string tfilename = Name()+"_demo.py";
        ofstream tout;
        tout.open(tfilename.c_str());
        tout << csrc << left;
        tout << "#" << endl;
        tout << "# " << tfilename << endl;
        tout << "#\n" ;
        tout << "#" << endl;
        tout << "# Python ODE solver for the vector field named: " << Name() << endl;
        tout << "#" << endl;
        tout << "# This script uses the scipy odeint function to solve the differential equations." << endl;
        tout << "#" << endl;
        PrintVFGENComment(tout,"# ");
        tout << "#\n" ;
        tout << endl;
        tout << "import sys" << endl;
        if (HasPi) {
            tout << "from math import pi" << endl;
        }
        tout << "import numpy as np" << endl;
        tout << "from scipy.integrate import odeint" << endl;
        tout << "import " << Name() << endl;
        tout << endl << endl;

        tout << "def use():" << endl;
        tout << "    print('use: ',sys.argv[0],' [options]')" << endl;
        tout << "    print('options:')" << endl;
        tout << "    print('    -h    Print this message.')" << endl;
        tout << "    for i in range(N_):" << endl;
        tout << "        print('    '+varnames_[i]+'=<initial_condition>  Default value is ', def_y_[i])" << endl;
        if (np > 0) {
            tout << "    for i in range(P_):" << endl;
            tout << "        print('    '+parnames_[i]+'=<parameter_value>    Default value is ', def_p_[i])" << endl;
        }
        tout << "    print('    abserr=<absolute_error_tolerance>   Default value is ', abserr)" << endl;
        tout << "    print('    relerr=<relative_error_tolerance>   Default value is ', relerr)" << endl;
        tout << "    print('    stoptime=<stop_time>                Default value is ', stoptime)" << endl;
        tout << "    print('    numpoints=<number_of_output_points> Default value is ', numpoints)" << endl;
        tout << endl;

        tout << "#" << endl;
        tout << "# Main script begins here..." << endl;
        tout << "#" << endl;

        if (HasPi) {
            tout << "Pi = pi\n";
        }
        AssignNameValueLists(tout, "", conname_list, "=", convalue_list, "");
        tout << "N_ = " << nv << "\n" ;
        if (np > 0) {
            tout << "P_ = " << np << "\n" ;
        }
        tout << "# Default values for the initial conditions";
        if (np > 0) {
            tout << ", parameters";
        }
        tout << " and solver parameters\n";
        tout << "def_y_ = [";
        for (int i = 0; i < nv; ++i) {
            tout << vardefic_list[i];
            if (i != nv-1) {
                tout << ", " ;
            }
        }
        tout << "]\n" ;
        if (np > 0) {
            tout << "def_p_ = [" ;
            for (int i = 0; i < np; ++i) {
                tout << pardefval_list[i] ;
                if (i != np-1) {
                    tout << ", " ;
                }
            }
            tout << "]\n" ;
        }
        tout << "abserr = 1.0e-8\n";
        tout << "relerr = 1.0e-6\n";
        tout << "stoptime = 10.0\n";
        tout << "numpoints = 250\n";
        MakePythonListOfStrings(tout,"varnames_",varname_list,"");
        if (np > 0) {
            MakePythonListOfStrings(tout,"parnames_",parname_list,"");
        }
        tout << endl;
        tout << "# Create a dict of all the options that can be given on the command line.\n";
        tout << "# Set the values to the default value of option.\n";
        tout << "options = {}" << endl;
        tout << "for i in range(N_):" << endl;
        tout << "    options[varnames_[i]] = def_y_[i]" << endl;
        if (np > 0) {
            tout << "for i in range(P_):" << endl;
            tout << "    options[parnames_[i]] = def_p_[i]" << endl;
        }
        tout << "options['abserr'] = abserr" << endl;
        tout << "options['relerr'] = relerr" << endl;
        tout << "options['stoptime'] = stoptime" << endl;
        tout << "options['numpoints'] = numpoints" << endl;
        tout << endl;
        tout << "# Process the command line arguments.\n";
        tout << "for a in sys.argv[1:]:" << endl;
        tout << "    if a == '-h' or a == '-help' or a == '--help' or a == 'help':" << endl;
        tout << "        use()" << endl;
        tout << "        sys.exit()" << endl;
        tout << "    eqloc = a.find('=')" << endl;
        tout << "    if (eqloc == -1):" << endl;
        tout << "        print('Invalid argument (missing =): ', a)" << endl;
        tout << "        use()" << endl;
        tout << "        sys.exit()" << endl;
        tout << "    else:" << endl;
        tout << "        var = a[0:eqloc]" << endl;
        tout << "        val = a[eqloc+1:]" << endl;
        tout << "        if var in options:" << endl;
        tout << "            options[var] = float(val)" << endl;
        tout << "        else:" << endl;
        tout << "            print('Unknown argument: ', a)" << endl;
        tout << "            use()" << endl;
        tout << "            sys.exit()" << endl;
        tout << endl;
        tout << "# Get the values for the initial conditions";
        if (np > 0) {
            tout << " and parameters";
        }
        tout << " from the options dict.\n";
        tout << "y_ = []" << endl;
        tout << "for i in range(N_):" << endl;
        tout << "    y_ = y_ + [options[varnames_[i]]]" << endl;
        if (np > 0) {
            tout << "p_ = []" << endl;
            tout << "for i in range(P_):" << endl;
            tout << "    p_ = p_ + [options[parnames_[i]]]" << endl;
        }
        tout << endl;
        tout << "# Create the time samples for the output of the ODE solver." << endl;
        tout << "tfinal = options['stoptime']" << endl;
        tout << "N = int(options['numpoints'])" << endl;
        tout << "if N < 2:\n";
        tout << "    print('The number of points must be at least 2.')\n";
        tout << "    sys.exit()\n";
        tout << "t = np.linspace(0, tfinal, N)\n";
        tout << endl;
        tout << "# Call the ODE solver.\n";
        tout << "ysol = odeint(" << Name() << ".vectorfield, y_, t";
        if (options["tfirst"] == "yes") {
            tout << ", tfirst=True";
        }
        if (np > 0) {
            tout << ", args=(p_,)";
        }
        tout << ",\n              Dfun=" << Name() << ".jacobian, atol=options['abserr'], rtol=options['relerr'])\n";
        tout << endl;
        tout << "# Print the solution.\n";
        tout << "for t1, y1 in zip(t, ysol):\n";
        tout << "    print(t1, ";
        for (int i = 0; i < nv; ++i) {
            tout << "y1[" << i << "], ";
        }
        tout << "end=' ')\n";
        tout << "    print(";
        for (int i = 0; i < nf; ++i) {
            tout << Name() << "." << funcname_list[i];
            if (options["tfirst"] == "yes") {
                tout << "(t1, y1";
            }
            else {
                tout << "(y1, t1";
            }
            if (np > 0) {
                tout << ", p_";
            }
            tout << ")";
            if (i < nf-1) {
                tout << ",";
            }
        }
        tout << ")\n";
        tout.close();
    }
}
