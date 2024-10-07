
//
//  vf_pydstool.cpp
//
//  This file defines the VectorField::PyDSTool method.
//
//
//  Copyright (C) 2008-2014 Warren Weckesser
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
// PrintPyDSTool -- The PyDSTool Code Generator.
//

void VectorField::PrintPyDSTool(map<string, string> options)
{
    int nc, nv, np;

    nc = conname_list.nops();
    nv = varname_list.nops();
    np = parname_list.nops();

    string filename = Name() + ".py";
    ofstream fout;
    fout.open(filename.c_str());
    // fout << python << left;
    fout << csrc << left;

    //
    //  Print Python file header information.
    //
    fout << "#" << endl;
    fout << "# " << filename << endl;
    fout << "#" << endl;
    fout << "# PyDSTool Python file for the vector field named: " << Name() << endl;
    fout << "#" << endl;
    PrintVFGENComment(fout, "# ");
    fout << "#" << endl;
    fout << endl;
    fout << "import PyDSTool" << endl;

    if (HasPi) {
        fout << "from math import pi\n";
    }
    fout << endl;
    //
    // Create a function that returns an "args" object for the vector field.
    //
    fout << "#" << endl;
    fout << "# args()\n";
    fout << "#\n";
    fout << "# This function creates a PyDSTool 'args' object for the\n";
    fout << "# '" << Name() << "' vector field.\n";
    fout << "#\n";
    fout << endl;
    fout << "def args():" << endl;
    fout << "    \"\"\"\n";
    fout << "    This function creates a PyDSTool 'args' object for the\n";
    fout << "    '" << Name() << "' vector field.\n";
    fout << "    \"\"\"\n";
    fout << "    DSargs = PyDSTool.args()\n";
    fout << "    DSargs.name = '" << Name() << "'\n";
    if (HasPi) {
        fout << "    Pi = pi\n";
    }

    fout << "    DSargs.pars = {";
    for (int i = 0; i < np; ++i) {
        if (i > 0) {
            fout << ", ";
        }
        fout << "'" << parname_list[i] << "':" << pardefval_list[i];
    }
    fout << "}\n";
    /*
    fout << "    DSargs.reuseterms = {";
    for (int i = 0; i < nc; ++i) {

        fout << "'" << convalue_list[i] << "':'" << conname_list[i] << "'";
        if (i < nc-1 || na > 0) {
            fout << ", ";
        }
    }
    if (nc > 0) {
        fout << "\n
    }                  ";
    for (int i = 0; i < na; ++i) {
        fout << "'" << exprformula_list[i] << "':'" << exprname_list[i] << "'";
        if (i < na-1) {
            fout << ", ";
        }
    }
    fout << "}\n";
    */
    // The vector field.  I think I should be able to put the Constants and
    // Expressions in 'reuseterms', but my initial attempt to do that caused
    // errors when I ran the program.
    fout << "    DSargs.varspecs = {";
    for (int i = 0; i < nv; ++i) {
        if (i > 0) {
            fout << ", ";
        }
        ex f = iterated_subs(varvecfield_list[i], expreqn_list);
        for (int k = 0; k < nc; ++k) {
            f = f.subs(conname_list[k] == convalue_list[k]);
        }
        fout << "'" << varname_list[i] << "':'" << f << "'";
    }
    fout << "}\n";
    fout << "    DSargs.fnspecs = {'Jacobian': (['t'";
    for (int i = 0; i < nv; ++i) {
        fout << ", '" << varname_list[i] << "'";
    }
    fout << "],\n";
    for (int i = 0; i < nv; ++i) {
        ex f = iterated_subs(varvecfield_list[i], expreqn_list);
        for (int k = 0; k < nc; ++k) {
            f = f.subs(conname_list[k] == convalue_list[k]);
        }
        if (i == 0) {
            fout << "            \"\"\"[";
        }
        else {
            fout << "                ";
        }
        fout << "[";
        for (int j = 0; j < nv; ++j) {
            symbol v = ex_to<symbol>(varname_list[j]);
            ex df = f.diff(v);
            if (j > 0) {
                fout << ", ";
            }
            fout << df;
        }
        fout << "]";
        if (i < nv-1) {
            fout << ",\n";
        }
        else {
            fout << "]\"\"\")";
        }
    }
    fout << "}\n";
    fout << "    DSargs.ics = {";
    for (int i = 0; i < nv; ++i) {
        if (i > 0) {
            fout << ", ";
        }
        fout << "'" << varname_list[i] << "':" << vardefic_list[i];
    }
    fout << "}\n";
    fout << "    DSargs.tdomain = [0, 10]\n";
    fout << "    return DSargs\n";
    fout.close();

    if (options["demo"] == "yes") {
        // Create a script that uses Vode_ODEsystem to create and plot a solution.
        string tfilename = Name() + "_dst.py";
        ofstream tout;
        tout.open(tfilename.c_str());
        tout << csrc << left;
        tout << "#" << endl;
        tout << "# " << tfilename << endl;
        tout << "#\n" ;
        tout << "#" << endl;
        tout << "# This script uses PyDSTool to plot a solution to the\n";
        tout << "# differential equations defined in " << Name() << ".py\n";
        tout << "#" << endl;
        PrintVFGENComment(tout, "# ");
        tout << "#\n" ;
        tout << endl;
        // tout << "from math import *" << endl;
        // tout << "from matplotlib.font_manager import FontProperties\n";
        tout << "import sys\n";
        tout << "import matplotlib.pyplot as plt\n";
        tout << "import PyDSTool\n";
        tout << "import " << Name() << endl;
        tout << endl;
        tout << "# Compute the solution\n";
        tout << "ds = " << Name() << ".args()\n";
        tout << "ode = PyDSTool.Generator.Vode_ODEsystem(ds)\n";
        tout << "traj = ode.compute('traj')\n";
        tout << "sol = traj.sample(dt=0.05)\n";
        tout << endl;
        tout << "# Plot the solution\n";
        tout << "lw = 1.5\n";
        for (int i = 0; i < nv; ++i) {
            tout << "plt.plot(sol['t'], sol['" << varname_list[i] << "'], linewidth=lw)\n";
        }
        tout << "plt.xlabel('t')\n";
        tout << "plt.title('" << Name() << "')\n";
        tout << "plt.legend((";
        for (int i = 0; i < nv; ++i) {
            if (i > 0) {
                tout << ", ";
            }
            tout << "'" << varname_list[i] << "'";
        }
        // tout << "), prop=FontProperties(size=14))\n";
        tout << "))\n";
        tout << "plt.grid(True)\n";
        tout << "\n";
        tout << "if len(sys.argv) > 1:\n";
        tout << "    plt.savefig(sys.argv[1], dpi=72)\n";
        tout << "else:\n";
        tout << "    plt.show()\n";
        tout.close();
    }
}
