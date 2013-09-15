
//
//  vf_pygsl.cpp
//
//  This file defines the VectorField::PrintPyGSL method.
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
// PrintPyGSL -- The PyGSL Code Generator.
//

void VectorField::PrintPyGSL(map<string,string> options)
{
    int nc, nv, np, na, nf;

    nc = conname_list.nops();
    nv = varname_list.nops();
    np = parname_list.nops();
    na = exprname_list.nops();
    nf = funcname_list.nops();

    string filename = Name()+".py";
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
    fout << "# The functions defined here can be used with the ODEIV routines of PyGSL." << endl;
    fout << "#" << endl;
    PrintVFGENComment(fout,"# ");
    fout << "#" << endl;
    fout << endl;
    fout << "from math import *" << endl;
    fout << "import numpy" << endl;
    fout << endl;
    //
    //  Print the vector field function.
    //
    fout << "#" << endl;
    fout << "# The vector field." << endl;
    fout << "#" << endl;
    fout << endl;
    fout << "def " << "vectorfield(" << IndependentVariable << ",x_,args):" << endl;
    fout << "    \"\"\"\n";
    fout << "    The vector field function for the vector field \"" << Name() << "\"\n";
    fout << "    \"\"\"\n";
    if (HasPi) {
        fout << "    Pi = numpy.pi\n";
    }
    for (int i = 0; i < nc; ++i) {
        fout << "    " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
    }
    GetFromVector(fout, "    ", varname_list, "=", "x_", "[]", 0, "");
    GetFromVector(fout, "    ", parname_list, "=", "args", "[]", 0, "");
    for (int i = 0; i < na; ++i) {
        fout << "    " << exprname_list[i] << " = " << exprformula_list[i] << endl;
    }
    fout << endl;
    fout << "    f_ = numpy.zeros((" << nv << ",))" << endl;
    for (int i = 0; i < nv; ++i) {
        fout << "    f_[" << i << "] = " << varvecfield_list[i] << endl;
    }
    fout << endl;
    fout << "    return f_" << endl;
    fout << endl;

    //
    // Print the Jacobian function.
    //
    fout << "#" << endl;
    fout << "#  The Jacobian." << endl;
    fout << "#" << endl;
    fout << endl;
    fout << "def " << "jacobian(t_, y_, args):" << endl;
    fout << "    \"\"\"\n";
    fout << "    The Jacobian of the vector field \"" << Name() << "\"\n";
    fout << "    \"\"\"\n";
    if (HasPi) {
        fout << "    Pi = numpy.pi\n";
    }
    for (int i = 0; i < nc; ++i) {
        fout << "    " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
    }
    GetFromVector(fout, "    ", varname_list, "=", "y_", "[]", 0, "");
    GetFromVector(fout, "    ", parname_list, "=", "args", "[]", 0, "");
    for (int i = 0; i < na; ++i) {
        fout << "    " << exprname_list[i] << " = " << exprformula_list[i] << endl;
    }
    fout << endl;
    fout << "    # Create the Jacobian matrix, initialized with zeros." << endl; 
    fout << "    jac_ = numpy.zeros((" << nv << "," << nv << "))" << endl; 
    for (int i = 0; i < nv; ++i) {
        ex f = iterated_subs(varvecfield_list[i],expreqn_list);
        for (int j = 0; j < nv; ++j) {
            symbol v = ex_to<symbol>(varname_list[j]);
            ex df = f.diff(v);
            if (df != 0) {
                fout << "    jac_[" << i << "," << j << "] = " << df << endl;
            }
        }
    }
    fout << endl;
    fout << "    dfdt_ = numpy.zeros((2,),dtype=numpy.float)" << endl;
    symbol t(IndependentVariable);
    for (int i = 0; i < nv; ++i) {
        ex f = iterated_subs(varvecfield_list[i],expreqn_list);
        ex df = f.diff(t);
        if (df != 0) {
            fout << "    dfdt_[" << i << "] = " << df << endl;
        }
    }
    fout << endl;
    fout << "    return jac_,dfdt_" << endl;

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
            fout << "def " << funcname_list[n] << "(t_, y_, args):" << endl;
            fout << "    \"\"\"\n";
            fout << "    The user-defined function \"" << funcname_list[n] << "\" for the vector field \"" << Name() << "\"\n";
            fout << "    \"\"\"\n";
            if (HasPi) {
                fout << "    Pi = numpy.pi\n";
            }
            for (int i = 0; i < nc; ++i) {
                fout << "    " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
            }
            GetFromVector(fout, "    ", varname_list, "=", "y_", "[]", 0, "");
            GetFromVector(fout, "    ", parname_list, "=", "args", "[]", 0, "");
            for (int i = 0; i < na; ++i) {
                fout << "    " << exprname_list[i] << " = " << exprformula_list[i] << "" << endl;
            }
            fout << endl;
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
        tout << "# This script uses the PyGSL odeiv functions to solve the differential equations." << endl;
        tout << "#" << endl;
        PrintVFGENComment(tout,"# ");
        tout << "#\n" ;
        tout << endl;
        tout << "import sys" << endl;
        if (HasPi) {
            tout << "from math import pi as Pi" << endl;
        }
        tout << "from pygsl import odeiv" << endl;
        tout << "import " << Name() << endl;
        tout << endl;

        tout << "def use():" << endl;
        tout << "    print 'use: ',sys.argv[0],' [options]'" << endl;
        tout << "    print 'options:'" << endl;
        tout << "    print '    -h    Print this message.'" << endl;
        tout << "    for vname in varnames_:" << endl;
        tout << "        print '    '+vname+'=<initial_condition>  Default value is ',def_options[vname]" << endl;
        tout << "    for pname in parnames_:" << endl;
        tout << "        print '    '+pname+'=<parameter_value>    Default value is ',def_options[pname]" << endl;
        tout << "    print '    abserr=<absolute_error_tolerance>  Default value is ',def_options['abserr']" << endl;
        tout << "    print '    relerr=<relative_error_tolerance>  Default value is ',def_options['relerr']" << endl;
        tout << "    print '    stoptime=<stop_time>               Default value is ',def_options['stoptime']" << endl;
        tout << endl;

        tout << "#" << endl;
        tout << "# Main script begins here..." << endl;
        tout << "#" << endl;

        // if (HasPi) {
        //     tout << "Pi = " << numstr << ".pi\n";
        // }
        for (int i = 0; i < nc; ++i) {
            tout << conname_list[i] << " = " << convalue_list[i] << endl;
        }
        tout << "N_ = " << nv << "\n" ;

        MakePythonListOfStrings(tout,"varnames_",varname_list,"");
        MakePythonListOfStrings(tout,"parnames_",parname_list,"");
        tout << "solver_param_names_ = [\"abserr\",\"relerr\",\"stoptime\"]\n" ;
        tout << endl;
        tout << "#" << endl;
        tout << "# Set the default values of everything in options dict." << endl;
        tout << "#" << endl;
        tout << "options = {}" << endl;
        tout << "# Default initial conditions" << endl;
        for (int i = 0; i < nv; ++i) {
            tout << "options['" << varname_list[i] << "'] = " << vardefic_list[i] << endl;
        }
        tout << "# Default vector field parameter values" << endl;
        for (int i = 0; i < np; ++i) {
            tout << "options['" << parname_list[i] << "'] = " << pardefval_list[i] << endl;
        }
        tout << "# Default ODE solver parameters:" << endl;
        tout << "options['abserr'] = 1.0e-8" << endl;
        tout << "options['relerr'] = 1.0e-6" << endl;
        tout << "options['stoptime'] = 10.0" << endl;
        tout << endl;
        tout << "# Make a copy of the default options so use() can show the default values" << endl;
        tout << "def_options = options.copy()" << endl;
        tout << endl;
        tout << "#" << endl;
        tout << "# Process the command-line arguments" << endl;
        tout << "#" << endl;
        tout << "for a in sys.argv[1:]:" << endl;
        tout << "    if a == '-h' or a == '--help' or a == '-help' or a == 'help':" << endl;
        tout << "        use()" << endl;
        tout << "        sys.exit()" << endl;
        tout << "    eqloc = a.find('=')" << endl;
        tout << "    if (eqloc == -1):" << endl;
        tout << "        print 'Invalid argument (missing =): ', a" << endl;
        tout << "        use()" << endl;
        tout << "        sys.exit()" << endl;
        tout << "    else:" << endl;
        tout << "        var = a[0:eqloc]" << endl;
        tout << "        val = a[eqloc+1:]" << endl;
        tout << "        if var in options:" << endl;
        tout << "            options[var] = float(val)" << endl;
        tout << "        else:" << endl;
        tout << "            print 'Unknown argument: ', a" << endl;
        tout << "            use()" << endl;
        tout << "            sys.exit()" << endl;
        tout << endl;
        tout << "y_ = (";
        for (int i = 0; i < nv; ++i) {
            tout << "options['" << varname_list[i] << "']";
            if (nv == 1 | i < nv-1) {
                tout << ",";
            }
        }
        tout << ")" << endl;
        tout << "p_ = (";
        for (int i = 0; i < np; ++i) {
            tout << "options['" << parname_list[i] << "']";
            if (np == 1 | i < np-1) {
                tout << ",";
            }
        }
        tout << ")" << endl;
        tout << endl;
        tout << "# Create the GSL ODE solver" << endl;
        tout << "step    = odeiv.step_rk8pd(N_," + Name() + ".vectorfield," + Name() + ".jacobian,args=p_)" << endl;
        tout << "control = odeiv.control_y_new(step,options['abserr'],options['relerr'])" << endl;
        tout << "evolve  = odeiv.evolve(step,control,N_)" << endl;
        tout << endl;
        tout << "stoptime = options['stoptime']" << endl;
        tout << "# Initial step size is stoptime/500" << endl;
        tout << "h = stoptime/500.0" << endl;
        tout << "t = 0" << endl;
        tout << "# Print the initial condition\n";
        tout << "print t";
        for (int i = 0; i < nv; ++i) {
            tout << ",y_[" << i << "]";
        }
        if (options["func"] == "yes") {
            for (int i = 0; i < nf; ++i) {
                tout << "," << Name() << "." << funcname_list[i] << "(t,y_,args=p_)";
            }
        }
        tout << endl;
        tout << "# Call evolve.apply(...) until the solution reaches stoptime" << endl;
        tout << "while t < stoptime:" << endl;
        tout << "    t, h, y_ = evolve.apply(t,stoptime,h,y_)" << endl;
        // tout << "    y_ = y_[0]" << endl;
        tout << "    print t";
        for (int i = 0; i < nv; ++i) {
            tout << ",y_[" << i << "]";
        }
        if (options["func"] == "yes") {
            for (int i = 0; i < nf; ++i) {
                tout << "," << Name() << "." << funcname_list[i] << "(t,y_,args=p_)";
            }
        }
        tout << endl;
        tout.close();
    }
}
