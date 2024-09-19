
//
//  vf_javamath.cpp
//
//  This file defines the VectorField::PrintGSL method.
//
//
//  Copyright (C) 2022 Warren Weckesser
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


//
// PrintJavaMath -- Java code generator for the Apache Commons Math ODE solver.
//

void VectorField::PrintJavaMath(map<string,string> options)
{
    int nc, nv, na;

    symbol t(IndependentVariable);
    nc = conname_list.nops();
    nv = varname_list.nops();
    na = exprname_list.nops();

    string packagename = Name() + "ODE";

    string filename = packagename + ".java";
    ofstream fout;
    fout.open(filename.c_str());
    fout << csrc << left;

    //
    // Print preamble stuff: comments, imports, etc.
    //
    fout << "//" << endl;
    fout << "//  " << filename << endl;
    fout << "//" << endl;
    fout << "//  System definition file for the vector field named: " << Name() << endl;
    fout << "//" << endl;
    PrintVFGENComment(fout,"//  ");
    fout << "//" << endl;
    fout << endl;

    fout << "package " << packagename << ";" << endl;
    fout << endl;
    fout << "import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;" << endl;
    fout << "import static java.lang.Math.*;" << endl;
    fout << endl;

    //
    // Print the class the defines the vector field.
    //
    fout << "public class " << packagename << " implements FirstOrderDifferentialEquations {" << endl;
    fout << endl;
    CDeclare(fout, "private double", parname_list);
    fout << endl;

    //
    // Print the constructors.
    //
    fout << "    public " << packagename << "(";
    PrintTransformedList(fout, "double $", parname_list);
    fout << ") {" << endl;
    for (unsigned i = 0; i < parname_list.nops(); ++i) {
        fout << "        this." << parname_list[i] << " = " << parname_list[i] << ";" << endl;
    }
    fout << "    }" << endl;
    fout << endl;
    fout << "    public " << packagename << "() {" << endl;
    AssignNameValueLists(fout, "        ", parname_list, "=", pardefval_list, ";");
    fout << "    }" << endl;
    fout << endl;

    //
    // Print the methods.
    //
    fout << "    public int getDimension() {" << endl;
    fout << "        return " << nv << ";" << endl;
    fout << "    }" << endl;
    fout << endl;

    fout << "    //" << endl;
    fout << "    //  The vector field." << endl;
    fout << "    //" << endl;
    fout << "    public void computeDerivatives(double t, double[] state, double[] deriv) {" << endl;
    if (HasPi) {
        fout << "        double Pi = PI;\n";
    }
    for (int i = 0; i < nc; ++i) {
        fout << "        double " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
    }
    GetFromVector(fout, "        double ", varname_list, "=", "state", "[]", 0, ";");
    fout << endl;
    for (int i = 0; i < na; ++i) {
        fout << "        double " << exprname_list[i] << " = " << exprformula_list[i] << ";" << endl;
    }
    if (na > 0) {
        fout << endl;
    }
    for (int i = 0; i < nv; ++i) {
        fout << "        deriv[" << i << "]" << " = " << varvecfield_list[i] << ";" << endl;
    }
    fout << "    }" << endl;
    fout << "}" << endl;
    fout << endl;

    fout.close();

    if (options["demo"] == "yes") {
        //
        //  Create a self-contained ODE solver for this vector field
        //  that allows the user to give the initial conditions,
        //  parameter values, and some solver control parameters
        //  on the command line.
        //

        string classname = Name() + "Solver";
        string tfilename = classname + ".java";
        ofstream tout;
        tout.open(tfilename.c_str());
        tout << csrc << left;
        tout << "import static java.lang.Math.*;" << endl;
        tout << "import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;" << endl;
        tout << "import org.apache.commons.math3.ode.FirstOrderIntegrator;" << endl;
        tout << "import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;" << endl;
        tout << "import " << packagename << "." << packagename << ";" << endl;
        tout << endl;
        tout << "public class " << classname << " {" << endl;
        tout << endl;
        tout << "    static void printHeader() {" << endl;
        tout << "        System.out.println(\"" << IndependentVariable << ", ";
        PrintList(tout, varname_list);
        tout << "\");" << endl;
        tout << "    }" << endl;
        tout << endl;
        tout << "    static void printState(double t, double[] state) {" << endl;
        tout << "        String s = String.format(\"%14.6f\", t);" << endl;
        tout << "        System.out.print(s);" << endl;
        tout << "        for (int i = 0; i < " << nv << "; ++i) {" << endl;
        tout << "            s = String.format(\", %20.15f\", state[i]);" << endl;
        tout << "            System.out.print(s);" << endl;
        tout << "        }" << endl;
        tout << "        System.out.println();" << endl;
        tout << "    }" << endl;
        tout << endl;
        tout << "    public static void main(String[] args) {" << endl;
        if (HasPi) {
            tout << "        double Pi = PI;" << endl;
        }
        tout << endl;
        tout << "        // Adjust these integrator control parameters as needed." << endl;
        tout << "        double minstep = 1e-13;" << endl;
        tout << "        double maxstep = 10.0;" << endl;
        tout << "        double abstol = 1.0e-12;" << endl;
        tout << "        double reltol = 1.0e-9;" << endl;
        tout << endl;
        tout << "        FirstOrderIntegrator dp853 = new DormandPrince853Integrator(minstep, maxstep, abstol, reltol);" << endl;
        tout << "        FirstOrderDifferentialEquations ode = new " << packagename << "();" << endl;
        tout << endl;
        tout << "        // Initial conditions" << endl;
        tout << "        double[] state = new double[] {";
        PrintList(tout, vardefic_list);
        tout << "};" << endl;
        tout << endl;
        tout << "        // Adjust the following time parameters as needed." << endl;
        tout << "        double t = 0.0;    // Start time." << endl;
        tout << "        double t1 = 10.0;  // End time." << endl;
        tout << "        double stepsize = 0.05;" << endl;
        tout << endl;
        tout << "        printHeader();" << endl;
        tout << "        printState(t, state);" << endl;
        tout << "        while (t < t1) {" << endl;
        tout << "            double ts = t + stepsize;" << endl;
        tout << "            if (Math.abs((ts - t1)/ts) < 1e-12) {" << endl;
        tout << "                ts = t1;" << endl;
        tout << "            }" << endl;
        tout << "            dp853.integrate(ode, t, state, ts, state);" << endl;
        tout << "            printState(ts, state);" << endl;
        tout << "            t = ts;" << endl;
        tout << "        }" << endl;
        tout << "    }" << endl;
        tout << "}" << endl;
        tout.close();

        string bfilename = "build-and-run-" + classname + ".sh";
        ofstream bout;
        bout.open(bfilename.c_str());
        bout << "# Simple bash shell script for systems where Apache Commons Math is" << endl;
        bout << "# in /usr/share/java.  Adjust as needed, or use a more sophisticated" << endl;
        bout << "# Java build tool or IDE." << endl;
        bout << "export CLASSPATH=.:/usr/share/java/commons-math3.jar" << endl;
        bout << "javac -d . " << filename << endl;
        bout << "javac " << tfilename << endl;
        bout << "java " << classname << endl;
        bout.close();
    }
}
