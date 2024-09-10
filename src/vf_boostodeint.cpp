
//
//  vf_boostodeint.cpp
//
//  This file defines the VectorField::PrintBoostOdeint method.
//
//
//  Copyright (C) 2019, 2024 Warren Weckesser
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


/*
 * Idea:
 * if method is default:
 *     if there are no parameters
 *         system is defined using a function
 *     else
 *         system is defined using a class with operator()
 * elif method is implicit
 *     if there are no parameters
 *         system is defined with two functions
 *     else
 *         system is defined with a class, with two member functions
 *         called vf and jac.  To create the system argument (a pair),
 *         first create an instance of the class, and then create the
 *         pair using lambdas. For example, something like this:
 *             auto lor = lorenz(3, 5, 7);
 *             auto vf = [&lor](...args...) {return lor.vf(...args...);};
 *             auto jac = [&lor](...args...) {return lor.jac(...args...);}
 *             auto sys = make_pair(vf, jac)
 *         Then pass sys to the integration function.
 */


//
// PrintBoostOdeint -- The Boost Odeint Code Generator.
//

void VectorField::PrintBoostOdeint(map<string,string> options)
{
    int np, nv, na, nf;

    np = parname_list.nops();
    nv = varname_list.nops();
    na = exprname_list.nops();
    nf = funcname_list.nops();

    string filename = Name()+"_vf.cpp";
    ofstream fout;
    fout.open(filename.c_str());
    fout << csrc << left;

    string pfilename = Name()+"_vf.h";
    ofstream pout;
    pout.open(pfilename.c_str());
    pout << csrc << left;

    //
    //  Print the function/method definitions to fout.
    //
    fout << "//" << endl;
    fout << "//  " << filename << endl;
    fout << "//" << endl;
    fout << "//  Implementation for the vector field named: " << Name() << endl;
    fout << "//" << endl;
    PrintVFGENComment(fout, "//  ");
    fout << "//" << endl;
    fout << endl;
    fout << "#include <cmath>" << endl;
    fout << "#include \"" << pfilename << "\"" << endl;
    fout << endl;
    fout << "void " << Name() << "_vf";
    fout << "::operator()";
    fout << "(const state_type &x_, state_type &dxdt_, const double t_)" << endl;
    fout << "{" << endl;
    if (HasPi) {
        fout << "    const double Pi = M_PI;\n";
    }
    AssignNameValueLists(fout, "    const double ", conname_list, "=", convalue_list, ";");
    CDeclare_double(fout, varname_list);
    CDeclare_double(fout, exprname_list);
    fout << endl;
    GetFromVector(fout, "    ", varname_list, "=", "x_", "[]", 0, ";");
    fout << endl;
    AssignNameValueLists(fout, "    ", exprname_list, "=", exprformula_list, ";");
    SetVectorFromNames(fout, "    ", "dxdt_", varvecfield_list, "[]", 0, ";");
    fout << "}" << endl;
    fout << endl;

    //
    // Print the Jacobian function.
    //
    fout << "//" << endl;
    fout << "//  The Jacobian." << endl;
    fout << "//" << endl;
    fout << endl;

    fout << "void " << Name() << "_vf::jac(";
    fout << "const state_type &x_, ";
    fout << "matrix_type &J_, ";
    fout << "const double &t_, ";
    fout << "state_type &dfdt_";
    fout << ")" << endl;
    fout << "{" << endl;
    if (HasPi) {
        fout << "    const double Pi = M_PI;\n";
    }
    AssignNameValueLists(fout, "    const double ", conname_list, "=", convalue_list, ";");

    CDeclare(fout, "double", varname_list);
    GetFromVector(fout, "    ", varname_list, "=", "x_", "[]", 0, ";");
    for (int i = 0; i < nv; ++i) {
        ex f = iterated_subs(varvecfield_list[i], expreqn_list);
        for (int j = 0; j < nv; ++j) {
            symbol v = ex_to<symbol>(varname_list[j]);
            ex df = f.diff(v);
            fout << "    J_(" << i << ", " << j << ") = " << f.diff(v) << ";" << endl;
        }
    }
    for (int i = 0; i < nv; ++i) {
        ex f = iterated_subs(varvecfield_list[i], expreqn_list);
        ex dfdt = f.diff(IndVar);
        fout << "    dfdt_(" << i << ") = " << dfdt << ";" << endl;
    }
    fout << "}" << endl;

    if ((options["func"] == "yes") && (nf > 0)) {
        //
        // Print the user-defined functions.
        //
        fout << endl;
        fout << "//" << endl;
        fout << "//  User-defined functions. " << endl;
        fout << "//" << endl;
        for (int n = 0; n < nf; ++n) {
            fout << endl;
            fout << "double " << Name() << "_vf::" << funcname_list[n] << "(const state_type &x_, double t_)" << endl;
            fout << "{" << endl;
            if (HasPi) {
                fout << "    const double Pi = M_PI;\n";
            }
            AssignNameValueLists(fout, "    const double ", conname_list, "=", convalue_list, ";");
            CDeclare(fout, "double", varname_list);
            CDeclare(fout, "double", exprname_list);
            fout << endl;
            GetFromVector(fout, "    ", varname_list, "=", "x_", "[]", 0, ";");
            fout << endl;
            for (int i = 0; i < na; ++i) {
                fout << "    " << exprname_list[i] << " = " << exprformula_list[i] << ";" << endl;
            }
            if (na > 0) {
                fout << endl;
            }
            fout << "    return " << funcformula_list[n] << ";" << endl;
            fout << "}" << endl;
        }
    }

    fout.close();

    //
    // Print the class declaration to pout.
    //
    pout << "//" << endl;
    pout << "//  " << pfilename << endl;
    pout << "//" << endl;
    pout << "//  Header file for the vector field " << Name() << endl;
    pout << "//" << endl;
    PrintVFGENComment(pout,"//  ");
    pout << "//" << endl;
    pout << endl;

    pout << "#ifndef " << Name() << "_VF_H" << endl;
    pout << "#define " << Name() << "_VF_H" << endl;
    pout << endl;
    
    pout << "#include <boost/numeric/odeint.hpp>" << endl;
    pout << endl;
    pout << "typedef boost::numeric::ublas::vector<double> state_type;" << endl;
    pout << "typedef boost::numeric::ublas::matrix<double> matrix_type;" << endl;
    pout << endl;

    //
    //  Print the vector field class.
    //
    pout << "//" << endl;
    pout << "//  The vector field." << endl;
    pout << "//" << endl;
    pout << endl;

    pout << "class " << Name() << "_vf" << endl;
    pout << "{" << endl;

    if (np > 0) {
        CDeclare_double(pout, parname_list);
        pout << endl;
    }

    pout << "public:" << endl;

    if (np > 0) {
        pout << "    " << Name() + "_vf" << "(";
        PrintTransformedList(pout, "double $_", parname_list);
        pout << ") : ";
        PrintTransformedList(pout, "$($_)", parname_list);
        pout << " {}" << endl;
        pout << endl;
    }

    pout << "    void operator()";
    pout << "(const state_type &x_, state_type &dxdt_, const double t_);" << endl;
    pout << "    void jac";
    pout << "(const state_type &x_, matrix_type &J_, const double &t_, state_type &dfdt_);" << endl;
    if (options["func"] == "yes") {
        for (int n = 0; n < nf; ++n) {
            pout << "    double " << funcname_list[n]
                 << "(const state_type &x_, double t_);" << endl;
        }
    }
    pout << "};" << endl;

    pout << endl;
    pout << "#endif" << endl;
    pout.close();

    if (options["demo"] == "yes") {
        string tfilename = Name()+"_demo.cpp";
        ofstream tout;
        tout.open(tfilename.c_str());
        tout << "#include <iostream>" << endl;
        tout << "#include <vector>" << endl;
        tout << "#include <cmath>" << endl;
        tout << "#include <boost/numeric/odeint.hpp>" << endl;
        if (options["system"] != "implicit") {
            tout << "#include <boost/numeric/odeint/stepper/bulirsch_stoer_dense_out.hpp>" << endl;
        }
        tout << endl;
        tout << "#include \"" << Name() << "_vf.h\"" << endl;
        tout << endl;
        tout << "using namespace std;" << endl;
        tout << "using namespace boost::numeric::odeint;" << endl;
        tout << endl;
        tout << "void writer(const state_type &y, const double t)" << endl;
        tout << "{" << endl;
        tout << "    cout << t;" << endl;
        tout << "    for (auto v : y) {" << endl;
        tout << "        cout << \" \" << v;" << endl;
        tout << "    }" << endl;
        tout << "    cout << endl;" << endl;
        tout << "}" << endl;
        tout << endl;
        tout << "int main(int argc, char *argv[])" << endl;
        tout << "{" << endl;

        if (HasPi) {
            tout << "    const double Pi = M_PI;\n";
        }
        AssignNameValueLists(tout, "    const double ", conname_list, "=", convalue_list, ";");

        tout << "    state_type y_(" << nv << ");" << endl;

        tout << "    double t0 = 0.0;" << endl;
        tout << "    double tfinal = 10.0;" << endl;
        tout << "    double dt = 0.01;" << endl;
        tout << endl;
        
        for (int i = 0; i < nv; ++i) {
            tout << "    y_[" << i << "] = " << vardefic_list[i] << ";" << endl;
        }
        tout << endl;

        if (options["system"] != "implicit") {
            tout << "    bulirsch_stoer_dense_out<state_type> stepper(1e-9, 1e-9, 1.0, 0.0);" << endl;
        }

        tout << "    auto " << Name() << " = " << Name() + "_vf(";
        PrintList(tout, pardefval_list);
        tout << ");" << endl;

        if (options["system"] == "implicit") {
            tout << "    integrate_const(make_dense_output<rosenbrock4<double>>(1e-9, 1e-9), ";
        }
        else {
            tout << "    integrate_const(stepper, ";
        }

        if (options["system"] == "implicit") {
            tout << endl;
            tout << "        make_pair(" << Name() << "," << endl;
            tout << "            [&" << Name() << "](const state_type &x_, matrix_type &J_, const double &t_, state_type &dfdt_)" << endl;
            tout << "            {" << Name() << ".jac(x_, J_, t_, dfdt_);}" << endl;
            tout << "        )";
        }
        else { // Not "implicit"
            tout << Name();
        }
        tout << ", y_, t0, tfinal, dt, writer);" << endl;

        tout << "}" << endl;
        tout.close();

        string makefilename = "Makefile." + Name() + "_demo";
        ofstream mout;
        mout.open(makefilename.c_str());
        mout << "# Makefile for " << Name() << "_demo" << endl;
        mout << endl;
        mout << "OBJFILES = " << Name() + "_demo.o" << " " << Name() + "_vf.o" << endl;
        mout << endl;
        mout << Name() + "_demo: $(OBJFILES)" << endl;
        mout << "\t$(CXX) $(OBJFILES) $(LDFLAGS) -o $@" << endl;
        mout << endl;
        mout << Name() + "_demo.o: " << tfilename << " " << Name() + "_vf.h" << endl;
        mout << "\t$(CXX) $(CXXFLAGS) -c " << tfilename << endl;
        mout << endl;
        mout << Name() + "_vf.o: " << filename << " " << Name() + "_vf.h" << endl;
        mout << "\t$(CXX) $(CXXFLAGS) -c " << filename << endl;
        mout.close();
    }
}
