
//
//  vf_boostodeint.cpp
//
//  This file defines the VectorField::PrintBoostOdeint method.
//
//
//  Copyright (C) 2019 Warren Weckesser
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
    int nc, np, nv, na, nf;

    symbol t(IndependentVariable);
    nc = conname_list.nops();
    nv = varname_list.nops();
    np = parname_list.nops();
    na = exprname_list.nops();
    nf = funcname_list.nops();

    if (options["system"] != "default") {
        cerr << "Sorry, only 'system=default' is implemented.\n";
        exit(-1);
    }
    string filename = Name()+"_vf.cpp";
    ofstream fout;
    fout.open(filename.c_str());
    fout << csrc << left;

    string pfilename = Name()+"_vf.h";
    ofstream pout;
    pout.open(pfilename.c_str());
    pout << csrc << left;

    //
    //  Print the class method definitions to fout.
    //
    fout << "//" << endl;
    fout << "//  " << filename << endl;
    fout << "//" << endl;
    fout << "//  Method definitions for the vector field named: " << Name() << endl;
    fout << "//" << endl;
    PrintVFGENComment(fout,"//  ");
    fout << "//" << endl;
    fout << endl;

    fout << "#include \"" << pfilename << "\"" << endl;
    fout << endl;
    fout << "void " << Name() << "_vf";
    if (np > 0) {
        fout << "::operator()";
    }
    fout << "(const state_type &x_, state_type &dxdt_, const double t_)" << endl;
    fout << "{" << endl;
    if (HasPi) {
        fout << "    const double Pi = M_PI;\n";
    }
    for (int i = 0; i < nc; ++i) {
        fout << "    const double " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
    }
    CDeclare_double(fout, varname_list);
    CDeclare_double(fout, exprname_list);
    fout << endl;
    GetFromVector(fout, "    ", varname_list, "=", "x_", "[]", 0, ";");
    fout << endl;
    for (int i = 0; i < na; ++i) {
        fout << "    " << exprname_list[i] << " = " << exprformula_list[i] << ";" << endl;
    }
    if (na > 0) {
        fout << endl;
    }
    for (int i = 0; i < nv; ++i) {
        fout << "    dxdt_[" << i << "]" << " = " << varvecfield_list[i] << ";" << endl;
    }
    fout << "}" << endl;
    fout << endl;

    //
    // Print the class declaration to pout.
    //
    pout << "//" << endl;
    pout << "//  " << pfilename << endl;
    pout << "//" << endl;
    pout << "//  Class definition for the vector field " << Name() << endl;
    pout << "//" << endl;
    PrintVFGENComment(pout,"//  ");
    pout << "//" << endl;
    pout << endl;

    pout << "#ifndef " << Name() << "_VF_H" << endl;
    pout << "#define " << Name() << "_VF_H" << endl;
    pout << endl;
    pout << "#include <vector>" << endl;
    pout << "#include <math.h>" << endl;
    pout << endl;
    pout << "typedef std::vector<double> state_type;" << endl;
    pout << endl;
    //
    //  Print the vector field class.
    //
    pout << "//" << endl;
    pout << "//  The vector field." << endl;
    pout << "//" << endl;
    pout << endl;
    if (np > 0) {
        pout << "class " << Name() << "_vf" << endl;
        pout << "{" << endl;
        CDeclare_double(pout, parname_list);
        pout << endl;
        pout << "public:" << endl;
        pout << "    " << Name() + "_vf" << "(";
        for (int i = 0; i < np; ++i) {
            if (i > 0) {
                pout << ", ";
            }
            pout << "double " << parname_list[i] << "_";
        }
        pout << ") : ";
        for (int i = 0; i < np; ++i) {
            if ( i > 0) {
                pout << ", ";
            }
            pout << parname_list[i] << "(" << parname_list[i] << "_" << ")";
        }
        pout << " {}" << endl;
        pout << endl;
    }
    if (np > 0) {
        pout << "    void operator()";
    }
    else {
        pout << "void " << Name() + "_vf";
    }
    pout << "(const state_type &x_, state_type &dxdt_, const double t_);" << endl;
    if (np > 0) {
        pout << "};" << endl;
    }
    pout << endl;
    pout << "#endif" << endl;

    fout.close();
    pout.close();

    if (options["demo"] == "yes") {
        string tfilename = Name()+"_demo.cpp";
        ofstream tout;
        tout.open(tfilename.c_str());
        tout << "#include <iostream>" << endl;
        tout << "#include <vector>" << endl;
        tout << "#include <boost/numeric/odeint.hpp>" << endl;
        tout << "#include <boost/numeric/odeint/stepper/bulirsch_stoer_dense_out.hpp>" << endl;
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
        tout << "    state_type y_{";
        for (int i = 0; i < nv; ++i) {
            tout << vardefic_list[i];
            if (i != nv-1) {
                tout << ", " ;
            }
        }
        tout << "};" << endl;
        tout << "    double t = 0.0;" << endl;
        tout << "    double dt = 0.01;" << endl;

        tout << endl;
        tout << "    bulirsch_stoer_dense_out<state_type> stepper(1e-9, 1e-9, 1.0, 0.0);" << endl;
        
        if (np > 0) {
            tout << "    auto " << Name() << " = " << Name() + "_vf(";
            for (int i = 0; i < np; ++i) {
                tout << pardefval_list[i];
                if (i != np-1) {
                    tout << ", ";
                }
            }
            tout << ");" << endl;
        }

        tout << "    integrate_const(stepper, ";
        if (np > 0) {
            tout << Name();
        }
        else {
            tout << Name() + "_vf";
        }
        tout << ", y_, t, 100.0, dt, writer);" << endl;
        tout << endl;
        tout << "}" << endl;
        tout.close(); 
    }
#ifdef NOTDEFINED
    if (options["demo"] == "yes") {
        //
        //  Create a self-contained ODE solver for this vector field
        //  that allows the user to give the initial conditions,
        //  parameter values, and some solver control parameters
        //  on the command line.
        //

        string tfilename = Name()+"_solve.c";
        ofstream tout;
        tout.open(tfilename.c_str());
        tout << csrc << left;
        tout << "/*\n";
        tout << " *  " << tfilename << endl;
        tout << " *\n" ;
        tout << " *" << endl;
        tout << " *  GSL ODE solver for the vector field named: " << Name() << endl;
        tout << " *" << endl;
        PrintVFGENComment(tout," *  ");
        tout << " *\n" ;
        tout << " *  To compile and run this program:\n" ;
        tout << " *      gcc -c " << Name() << "_gvf.c\n" ;
        tout << " *      gcc -c " << Name() << "_solve.c\n" ;
        tout << " *      gcc -o " << Name() << "_solve  " << Name() << "_solve.o " << Name() << "_gvf.o -lgsl -lgslcblas -lm\n" ;
        tout << " *  This creates an executable file called " << Name() << "_solve\n" ;
        tout << " */" << endl;
        tout << endl;
        tout << "#include <string.h>" << endl;
        tout << "#include <math.h>" << endl;
        tout << "#include <gsl/gsl_errno.h>" << endl;
        tout << "#include <gsl/gsl_matrix.h>" << endl;
        tout << "#include <gsl/gsl_odeiv.h>" << endl;
        tout << endl;
        tout << "#include \"" << pfilename << "\"\n";
        tout << endl;
        tout << endl << endl;
        tout << "int use(int argc, char *argv[], int nv, char *vname[], double y_[], int np, char *pname[], double p_[])\n" ;
        tout << "{\n" ;
        tout << "    int i;\n" ;
        tout << "    printf(\"use: %s [options]\\n\", argv[0]);\n" ;
        tout << "    printf(\"options:\\n\");\n" ;
        tout << "    printf(\"    -h    Print this help message.\\n\");\n" ;
        tout << "    for (i = 0; i < nv; ++i) {\n" ;
        tout << "        printf(\"    %s=<initial_condition>   Default value is %e\\n\", vname[i], y_[i]);\n";
        tout << "    }\n";
        tout << "    for (i = 0; i < np; ++i) {\n" ;
        tout << "        printf(\"    %s=<parameter_value>   Default value is %e\\n\", pname[i], p_[i]);\n";
        tout << "    }\n";
        tout << "    printf(\"    abserr=<absolute_error_tolerance>\\n\");\n" ;
        tout << "    printf(\"    relerr=<relative_error_tolerance>\\n\");\n" ;
        tout << "    printf(\"    stoptime=<stop_time>\\n\");\n" ;
        tout << "    return 0;\n";
        tout << "}\n" ;
        tout << endl << endl;
        tout << "int assign(char *str[], int ns, double v[], char *a)\n" ;
        tout << "{\n" ;
        tout << "    int i;\n" ;
        tout << "    char name[256];\n" ;
        tout << "    char *e;\n" ;
        tout << "\n" ;
        tout << "    e = strchr(a,'=');\n" ;
        tout << "    if (e == NULL) {\n" ;
        tout << "        return(-1);\n" ;
        tout << "    }\n";
        tout << "    *e = '\\0';\n" ;
        tout << "    strcpy(name, a);\n" ;
        tout << "    *e = '=';\n" ;
        tout << "    ++e;\n" ;
        tout << "    for (i = 0; i < ns; ++i) {\n" ;
        tout << "        if (strcmp(str[i], name)==0) {\n" ;
        tout << "            break;\n" ;
        tout << "        }\n";
        tout << "    }\n";
        tout << "    if (i == ns) {\n" ;
        tout << "        return -1;\n" ;
        tout << "    }\n";
        tout << "    v[i] = atof(e);\n" ;
        tout << "    return i;\n" ;
        tout << "}\n" ;
        tout << endl << endl;
        tout << "int main (int argc, char *argv[])\n" ;
        tout << "{\n";
        tout << "    int i;\n" ;
        if (HasPi) {
            tout << "    const double Pi = M_PI;\n";
        }
        for (int i = 0; i < nc; ++i) {
            tout << "    const double " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
        CDeclare_double(tout, parname_list);
        tout << "    const int P_ = " << np << ";\n" ;
        tout << "    double def_p_[" << np << "] = {" ;
        for (int i = 0; i < np; ++i) {
            tout << pardefval_list[i] ;
            if (i != np-1) {
                tout << ", " ;
            }
        }
        tout << "};\n" ;
        GetFromVector(tout, "    ", parname_list, "=", "def_p_", "[]", 0, ";");
        tout << "    double p_[" << np << "];\n" ;
        tout << "    const int N_ = " << nv << ";\n" ;
        tout << "    double def_y_[" << nv << "] = {";
        for (int i = 0; i < nv; ++i) {
            // tout << def_var_value.at(i) ;
            // tout << "0.0" ;
            tout << vardefic_list[i];
            if (i != nv-1) {
                tout << ", " ;
            }
        }
        tout << "};\n" ;
        tout << "    double y_[" << nv << "];\n" ;

        tout << "    double solver_param_[3] = {1.0e-6, 0.0, 10.0};\n" ;

        MakeCArrayOfStrings(tout,"varnames_", varname_list);
        if (np > 0) {
            MakeCArrayOfStrings(tout,"parnames_", parname_list);
        }
        else {
            tout << "    char *parnames_[] = {\"\"};\n";
        }
        tout << "    char *solver_param_names_[3] = {\"abserr\", \"relerr\", \"stoptime\"};\n" ;
        tout << endl;
        tout << "    for (i = 0; i < N_; ++i) {\n" ;
        tout << "        y_[i] = def_y_[i];\n" ;
        tout << "    }\n";
        tout << "    for (i = 0; i < P_; ++i) {\n" ;
        tout << "        p_[i] = def_p_[i];\n" ;
        tout << "    }\n";
        tout << "    for (i = 1; i < argc; ++i) {\n" ;
        tout << "        int j;\n" ;
        tout << "        if (strcmp(argv[i], \"-h\") == 0) {\n" ;
        tout << "            use(argc, argv, N_, varnames_, def_y_, P_, parnames_, def_p_);\n" ;
        tout << "            exit(0);\n" ;
        tout << "        }\n" ;
        tout << "        j = assign(varnames_, N_, y_, argv[i]);\n" ;
        tout << "        if (j == -1) {\n" ;
        tout << "            j = assign(parnames_, P_, p_, argv[i]);\n" ;
        tout << "            if (j == -1) {\n" ;
        tout << "                j = assign(solver_param_names_, 3, solver_param_, argv[i]);\n" ;
        tout << "                if (j == -1) {\n" ;
        tout << "                    fprintf(stderr, \"unknown argument: %s\\n\", argv[i]);\n" ;
        tout << "                    use(argc, argv, N_, varnames_, def_y_, P_, parnames_, def_p_); \n";
        tout << "                    exit(-1);\n" ;
        tout << "                }\n" ;
        tout << "            }\n" ;
        tout << "        }\n" ;
        tout << "    }\n" ;
        tout << endl;
        tout << "    const gsl_odeiv_step_type *T_  = gsl_odeiv_step_rk8pd;\n" ;
        tout << "    gsl_odeiv_step    *step_    = gsl_odeiv_step_alloc(T_, N_);\n" ;
        tout << "    gsl_odeiv_control *control_ = gsl_odeiv_control_y_new(solver_param_[0], solver_param_[1]);\n" ;
        tout << "    gsl_odeiv_evolve  *evolve_  = gsl_odeiv_evolve_alloc(N_);\n" ;
        tout << "    gsl_odeiv_system sys_ = {" << Name() << "_vf, " << Name() << "_jac, N_, &(p_[0])};\n";
        tout << endl;
        tout << "    double t_  = 0.0;\n" ;
        tout << "    double t1_ = solver_param_[2];\n" ;
        tout << "    double h_ = 1e-6;\n" ;
        tout << endl;
        tout << "    while (t_ < t1_) {\n" ;
        tout << "        int j_;\n" ;
        tout << "        int status_ = gsl_odeiv_evolve_apply(evolve_, control_, step_, &sys_, &t_, t1_, &h_, y_);\n" ;
        tout << "        if (status_ != GSL_SUCCESS) {\n" ;
        tout << "            fprintf(stderr, \"status=%d\\n\", status_);\n" ;
        tout << "            break;\n" ;
        tout << "        }\n" ;
        tout << "        printf(\"%.8e\", t_);\n" ;
        tout << "        for (j_ = 0; j_ < N_; ++j_) {\n" ;
        tout << "            printf(\" %.8e\", y_[j_]);\n" ;
        tout << "        }\n";
        if (options["func"] == "yes") {
            for (int i = 0; i < nf; ++i) {
                tout << "        printf(\" %.8e\", " << Name() << "_" << funcname_list[i] << "(t_, y_, p_));\n";
            }
        }
        tout << "        printf(\"\\n\");\n";
        tout << "    }\n" ;
        tout << endl;
        tout << "    gsl_odeiv_evolve_free(evolve_);\n" ;
        tout << "    gsl_odeiv_control_free(control_);\n" ;
        tout << "    gsl_odeiv_step_free(step_);\n" ;
        tout << "    return 0;\n" ;
        tout << "}\n" ;
        tout.close();
        //
        // Create a Makefile for the  GSL demo program
        //
        string mfilename = "Makefile-"+Name();
        ofstream mout;
        mout.open(mfilename.c_str());
        mout << "#\n";
        mout << "# " << mfilename << endl;
        mout << "#\n";
        mout << "# This is the Makefile for the " << Name() << "_solve program.\n";
        mout << "#\n";
        PrintVFGENComment(mout,"# ");
        mout << "#\n";
        mout << "# This Makefile is not guaranteed to work in all operating systems.\n";
        mout << "# You may have to edit this file to meet the conventions of your operating system.\n";
        mout << "#\n\n";
        mout << "all: " << Name() << "_solve" << endl;
        mout << endl;
        mout << Name() << "_solve: " << Name() << "_solve.o " << Name() << "_gvf.o" << endl;
        mout << "\t$(CC) $(LDFLAGS) -o " << Name() << "_solve ";
        mout <<       Name() << "_solve.o " << Name() << "_gvf.o -lgsl -lgslcblas -lm" << endl;
        mout << endl;
        mout << Name() << "_solve.o: " << Name() << "_solve.c " << Name() << "_gvf.h" << endl;
        mout << "\t$(CC) $(CPPFLAGS) -c " << Name() << "_solve.c" << endl;
        mout << endl;
        mout << Name() << "_gvf.o: " << Name() << "_gvf.c " << Name() << "_gvf.h" << endl;
        mout << "\t$(CC) $(CPPFLAGS) -c " << Name() << "_gvf.c" << endl;
        mout << endl;
        mout << "clean:\n";
        mout << "\trm -f " << Name() << "_solve " << Name() << "_solve.o " << Name() << "_gvf.o" << endl;
        mout << endl;
        mout.close();
    }
#endif
}
