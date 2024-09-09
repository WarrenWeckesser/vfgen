
//
//  vf_cvode7.cpp
//
//  This file defines the VectorField::PrintCVODE7 method.
//
//
//
//  Copyright (C) 2024 Warren Weckesser
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
#include <algorithm>
#include <ginac/ginac.h> 

#include "vf.h"
#include "vf_utils.h"

using namespace std;
using namespace GiNaC;


string to_upper(string s)
{
    string t = s;
    transform(t.begin(), t.end(), t.begin(), ::toupper);
    return t;
}

//
// PrintCVODE -- The CVODE Code Generator
//

void VectorField::PrintCVODE7(map<string,string> options)
{
    size_t nc, np, nv, na, nf;
    string include_guard_name;

    symbol t(IndependentVariable);
    nc = conname_list.nops();
    nv = varname_list.nops();
    np = parname_list.nops();
    na = exprname_list.nops();
    nf = funcname_list.nops();

    string filename = Name()+"_cv7.c";
    ofstream fout;
    fout.open(filename.c_str());
    fout << csrc << left;

    string pfilename = Name()+"_cv7.h";
    ofstream pout;
    pout.open(pfilename.c_str());
    pout << csrc << left;

    //
    //  Print C file header information.
    //
    fout << "/*" << endl;
    fout << " *  " << filename << endl;
    fout << " *" << endl;
    fout << " *  CVODE C file for the vector field named: " << Name() << endl;
    fout << " *" << endl;
    PrintVFGENComment(fout," *  ");
    fout << " */" << endl;
    fout << endl;

    pout << "/*" << endl;
    pout << " *  " << pfilename << endl;
    pout << " *" << endl;
    pout << " *  CVODE C prototype file for the functions defined in " << filename << endl;
    pout << " *" << endl;
    PrintVFGENComment(pout," *  ");
    pout << " */" << endl;
    pout << endl;
    include_guard_name = to_upper(Name());
    include_guard_name.append("_C7_H");
    pout << "#ifndef " << include_guard_name << endl;
    pout << "#define " << include_guard_name << endl;
    pout << endl;

    fout << "#include <math.h>" << endl;
    fout << endl;
    fout << "#include <cvode/cvode.h>" << endl;
    fout << "#include <nvector/nvector_serial.h>" << endl;
    fout << "#include <sunmatrix/sunmatrix_dense.h> /* dense SUNMatrix       */" << endl;
    fout << endl;
    //
    //  Print the vector field function.
    //
    fout << "/*" << endl;
    fout << " *  The vector field." << endl;
    fout << " */" << endl;
    fout << endl;
    string func_return_type = "int";
    fout << func_return_type << " " << Name() << "_vf(sunrealtype t, N_Vector y_, N_Vector f_, void *params)" << endl;
    pout << func_return_type << " " << Name() << "_vf(sunrealtype, N_Vector, N_Vector, void *);" << endl;
    fout << "{" << endl;
    if (IsAutonomous) {
        fout << "    (void) t;  /* Hack to avoid 'unused parameter' compiler warnings. */\n";
    }
    if (HasPi) {
        fout << "    const sunrealtype Pi = SUN_RCONST(M_PI);\n";
    }
    for (size_t i = 0; i < nc; ++i) {
        fout << "    const sunrealtype " << conname_list[i] << " = SUN_RCONST(" << convalue_list[i] << ");" << endl;
    }
    CDeclare(fout, "sunrealtype", varname_list);
    CDeclare(fout, "sunrealtype", parname_list);
    CDeclare(fout, "sunrealtype", exprname_list);
    fout << "    sunrealtype *p_;" << endl;
    fout << endl;
    fout << "    p_ = (sunrealtype *) params;" << endl;
    fout << endl;
    size_t width = max_expr_len(varname_list);
    for (size_t i = 0; i < nv; ++i) {
        fout << "    ";
        fout.width(width);
        fout << varname_list[i];
        fout.width(0);
        fout << " = NV_Ith_S(y_," << i << ");" << endl;
    }

    fout << endl;
    GetFromVector(fout, "    ", parname_list, "=", "p_", "[]", 0, ";");
    fout << endl;
    for (size_t i = 0; i < na; ++i) {
        fout << "    " << exprname_list[i] << " = " << exprformula_list[i] << ";" << endl;
    }
    if (na > 0) {
        fout << endl;
    }
    for (size_t i = 0; i < nv; ++i) {
        fout << "    NV_Ith_S(f_," << i << ") = " << varvecfield_list[i] << ";" << endl;
    }
    fout << "    return 0;\n";
    fout << "}" << endl;
    fout << endl;
    //
    // Print the Jacobian function.
    //
    fout << "/*" << endl;
    fout << " *  The Jacobian." << endl;
    fout << " */" << endl;
    fout << endl;
    fout << func_return_type << " " << Name() << "_jac(sunrealtype t, N_Vector y_, N_Vector fy_,";
    fout << " SUNMatrix jac_, void *params," << endl;
    fout << "                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)" << endl;

    pout << func_return_type << " " << Name() << "_jac(sunrealtype, N_Vector, N_Vector,";
    pout << " SUNMatrix, void *," << endl;
    pout << "                 N_Vector, N_Vector, N_Vector);" << endl;
    fout << "{" << endl;
    if (IsAutonomous) {
        fout << "    (void) t;     /* Hack to avoid 'unused parameter' compiler warnings. */\n";
    }
    fout << "    (void) fy_;   /* Hack to avoid 'unused parameter' compiler warnings. */\n";
    fout << "    (void) tmp1;  /* Hack to avoid 'unused parameter' compiler warnings. */\n";
    fout << "    (void) tmp2;  /* Hack to avoid 'unused parameter' compiler warnings. */\n";
    fout << "    (void) tmp3;  /* Hack to avoid 'unused parameter' compiler warnings. */\n";
    if (HasPi) {
        fout << "    const sunrealtype Pi = SUN_RCONST(M_PI);\n";
    }
    for (size_t i = 0; i < nc; ++i) {
        fout << "    const sunrealtype " << conname_list[i] << " = SUN_RCONST(" << convalue_list[i] << ");" << endl;
    }
    CDeclare(fout,"sunrealtype",varname_list);
    CDeclare(fout,"sunrealtype",parname_list);
    fout << "    sunrealtype *p_;" << endl;
    fout << endl;
    fout << "    p_ = (sunrealtype *) params;" << endl;
    fout << endl;
    for (size_t i = 0; i < nv; ++i) {
        fout << "    ";
        fout.width(width);
        fout << varname_list[i];
        fout.width(0);
        fout << " = NV_Ith_S(y_, " << i << ");" << endl;
    }
    fout << endl;
    GetFromVector(fout, "    ", parname_list, "=", "p_", "[]", 0, ";");
    fout << endl;
    for (size_t i = 0; i < nv; ++i) {
        ex f = iterated_subs(varvecfield_list[i],expreqn_list);
        for (size_t j = 0; j < nv; ++j) {
            symbol v = ex_to<symbol>(varname_list[j]);
            ex df = f.diff(v);
            // Skip zero elements.  CVODE initializes jac_ to zero before calling the Jacobian function.
            if (df != 0)
                fout << "    SM_ELEMENT_D(jac_, " << i << ", " << j << ") = " << f.diff(v) << ";" << endl;
        }
    }
    fout << "    return 0;\n";
    fout << "}" << endl;

    if ((options["func"] == "yes") && (nf > 0)) {
        //
        // Print the user-defined functions.
        // A single function is created that puts all the
        // user-defined function values in an array.  This
        // function is defined so that it can be used with
        // the CVODE rootfinding code.
        //
        fout << endl;
        fout << "/*" << endl;
        fout << " *  User-defined functions. " << endl;
        fout << " */" << endl;
        fout << endl;
        fout << func_return_type << " " << Name() << "_func(sunrealtype t, N_Vector y_, sunrealtype *func_, void *params)" << endl;
        pout << func_return_type << " " << Name() << "_func(sunrealtype, N_Vector, sunrealtype *, void *);" << endl;
        fout << "{" << endl;
        if (IsAutonomous) {
            fout << "    (void) t;  /* Hack to avoid 'unused parameter' compiler warnings. */\n";
        }
        if (HasPi) {
            fout << "    const sunrealtype Pi = SUN_RCONST(M_PI);\n";
        }
        for (size_t i = 0; i < nc; ++i) {
            fout << "    const sunrealtype " << conname_list[i] << " = SUN_RCONST(" << convalue_list[i] << ");" << endl;
        }
        CDeclare(fout,"sunrealtype",varname_list);
        CDeclare(fout,"sunrealtype",parname_list);
        CDeclare(fout,"sunrealtype",exprname_list);
        fout << "    sunrealtype *p_;" << endl;
        fout << endl;
        fout << "    p_ = (sunrealtype *) params;" << endl;
        fout << endl;
        for (size_t i = 0; i < nv; ++i) {
            fout << "    ";
            fout.width(width);
            fout << varname_list[i];
            fout.width(0);
            fout << " = NV_Ith_S(y_," << i << ");" << endl;
        }
        fout << endl;
        GetFromVector(fout, "    ", parname_list, "=", "p_", "[]", 0, ";");
        fout << endl;
        for (size_t i = 0; i < na; ++i) {
            fout << "    " << exprname_list[i] << " = " << exprformula_list[i] << ";" << endl;
        }
        if (na > 0) {
            fout << endl;
        }
        for (size_t n = 0; n < nf; ++n) {
            fout << "    /* " << funcname_list[n] << ":  */" << endl;
            fout << "    func_[" << n << "] = " << funcformula_list[n] << ";" << endl;
        }
        fout << "    return 0;\n";
        fout << "}" << endl;
    }
    fout.close();
    pout << endl;
    pout << "#endif /* " << include_guard_name << " */" << endl;
    pout.close();

    if (options["demo"] == "yes") {
        //
        //  Create a self-contained ODE solver for this vector field
        //  that allows the user to give the initial conditions,
        //  parameter values, and some solver control parameters
        //  on the command line.
        //

        string tfilename = Name()+"_cv7demo.c";
        ofstream tout;
        tout.open(tfilename.c_str());
        tout << csrc << left;
        tout << "/*\n";
        tout << " *  " << tfilename << endl;
        tout << " *\n" ;
        tout << " *" << endl;
        tout << " *  CVODE ODE solver for the vector field named: " << Name() << endl;
        tout << " *" << endl;
        PrintVFGENComment(tout," *  ");
        tout << " *\n";
        tout << " */" << endl;
        tout << endl;
        tout << "#include <stdlib.h>" << endl;
        tout << "#include <stdio.h>" << endl;
        tout << "#include <string.h>" << endl;
        tout << "#include <math.h>" << endl;
        tout << endl;

        tout << "/* Include headers for CVODE */" << endl;
        tout << "#include <cvode/cvode.h>" << endl;
        tout << "#include <nvector/nvector_serial.h>" << endl;
        tout << "#include <sunlinsol/sunlinsol_dense.h> /* dense SUNLinearSolver */" << endl;
        tout << "#include <sunmatrix/sunmatrix_dense.h> /* dense SUNMatrix       */" << endl;

        tout << endl;
        tout << "#include \"" << pfilename << "\"\n";
        tout << endl << endl;
        tout << "int check_int_status(int retval, char *funcname)\n";
        tout << "{\n";
        tout << "    if (retval != 0) {\n";
        tout << "        fprintf(stderr, \"SUNDIALS ERROR: %s() failed - returned %d\\n\", funcname, retval);\n";
        tout << "        return 1;\n";
        tout << "    }\n";
        tout << "    return 0;\n";
        tout << "}\n";
        tout << endl;
        tout << "int check_pointer(void *ptr, char *funcname)\n";
        tout << "{\n";
        tout << "    if (ptr == NULL) {\n";
        tout << "        fprintf(stderr, \"SUNDIALS ERROR: %s() failed - returned NULL\\n\", funcname);\n";
        tout << "        return 1;\n";
        tout << "    }\n";
        tout << "    return 0;\n";
        tout << "}\n";
        tout << endl;
        tout << "int use(const char *program_name, int nv, char *vname[], double y_[], int np, char *pname[], const double p_[])\n" ;
        tout << "{\n" ;
        tout << "    int i;\n" ;
        tout << "    printf(\"use: %s [options]\\n\", program_name);\n" ;
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
        tout << "    e = strchr(a, '=');\n" ;
        tout << "    if (e == NULL) {\n" ;
        tout << "        return -1;\n" ;
        tout << "    }\n";
        tout << "    *e = '\\0';\n" ;
        tout << "    strcpy(name, a);\n" ;
        tout << "    *e = '=';\n" ;
        tout << "    ++e;\n" ;
        tout << "    for (i = 0; i < ns; ++i) {\n" ;
        tout << "        if (strcmp(str[i], name) == 0) {\n" ;
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
        tout << "{\n" ;
        tout << "    SUNContext sunctx;\n";
        tout << "    int i, j;\n";
        tout << "    int retval;\n";
        tout << "    const int N_ = " << nv << ";\n" ;
        tout << "    const int P_ = " << np << ";\n" ;
        if (HasPi) {
            tout << "    const sunrealtype Pi = SUN_RCONST(M_PI);\n";
        }
        for (size_t i = 0; i < nc; ++i) {
            tout << "    const sunrealtype " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
        tout << "    const sunrealtype def_p_[" << np << "] = {\n" ;
        for (size_t i = 0; i < np; ++i) {
            tout << "        SUN_RCONST(" << pardefval_list[i] << ")" ;
            if (i != np-1) {
                tout << "," ;
            }
            tout << endl;
        }
        tout << "    };\n" ;
        tout << endl;
        tout << "    /* Create the SUNDIALS context */\n";
        tout << "    retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);\n";
        tout << "    if (check_int_status(retval, \"SUNContext_Create\")) {\n";
        tout << "        return 1;\n";
        tout << "    }\n";
        tout << endl;

        tout << "    sunrealtype def_y_[" << nv << "] = {";
        for (size_t i = 0; i < nv; ++i) {
            tout << "SUN_RCONST(" << vardefic_list[i] << ")";
            if (i != nv-1) {
                tout << ", " ;
            }
        }
        tout << "};\n" ;
        tout << "    sunrealtype y_[" << nv << "];\n" ;

        tout << "    sunrealtype p_[" << np << "];\n" ;
        tout << "    sunrealtype solver_param_[3] = {SUN_RCONST(1.0e-6), SUN_RCONST(0.0), SUN_RCONST(10.0)};\n" ;
        MakeCArrayOfStrings(tout,"varnames_",varname_list);
        if (np > 0) {
            MakeCArrayOfStrings(tout,"parnames_",parname_list);
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
        tout << "            use(argv[0], N_, varnames_, def_y_, P_, parnames_, def_p_);\n" ;
        tout << "            exit(0);\n" ;
        tout << "        }\n" ;
        tout << "        j = assign(varnames_, N_, y_, argv[i]);\n" ;
        tout << "        if (j == -1) {\n" ;
        tout << "            j = assign(parnames_, P_, p_, argv[i]);\n" ;
        tout << "            if (j == -1) {\n" ;
        tout << "                j = assign(solver_param_names_, 3, solver_param_, argv[i]);\n" ;
        tout << "                if (j == -1) {\n" ;
        tout << "                    fprintf(stderr, \"unknown argument: %s\\n\", argv[i]);\n" ;
        tout << "                    use(argv[0], N_, varnames_, def_y_, P_, parnames_, def_p_); \n";
        tout << "                    exit(-1);\n" ;
        tout << "                }\n" ;
        tout << "            }\n" ;
        tout << "        }\n" ;
        tout << "    }\n" ;
        tout << endl;
        tout << "    /* Initial conditions */\n";
        tout << "    N_Vector y0_ = N_VNew_Serial(N_, sunctx);\n";
        tout << "    if (check_pointer((void*)y0_, \"N_VNew_Serial\")) {\n";
        tout << "        return (1);\n";
        tout << "    }\n";
        tout << "    for (i = 0; i < N_; ++i) {\n";
        tout << "        NV_Ith_S(y0_, i) = y_[i];\n";
        tout << "    }\n";
        tout << endl;
        tout << "    /* Use CV_ADAMS for non-stiff problems, and CV_BDF for stiff problems:   */\n";
        tout << "    void *cvode_mem = CVodeCreate(CV_ADAMS, sunctx);\n";
        tout << "    if (check_pointer((void*)cvode_mem, \"CVodeCreate\")) {return 1;}\n";
        tout << endl;

        tout << "    sunrealtype t = SUN_RCONST(0.0);\n";

        tout << "    retval = CVodeInit(cvode_mem, " << Name() << "_vf, t, y0_);\n";
        tout << "    if (check_int_status(retval, \"CVodeInit\")) {\n";
        tout << "        return 1;\n";
        tout << "    }\n";
        tout << "    retval = CVodeSStolerances(cvode_mem, solver_param_[1], solver_param_[0]);\n";
        tout << "    if (check_int_status(retval, \"CVodeSStolerances\")) {return 1;}\n";
        tout << endl;
        tout << "    retval = CVodeSetUserData(cvode_mem, &(p_[0]));\n";
        tout << "    if (check_int_status(retval, \"CVodeSetUserData()\")) {return 1;}\n";
        tout << endl;
        tout << "    /* Create dense SUNMatrix for use in linear solves */\n";
        tout << "    SUNMatrix A = SUNDenseMatrix(N_, N_, sunctx);\n";
        tout << "    if (check_pointer((void*)A, \"SUNDenseMatrix()\")) {return 1;}\n";
        tout << endl;
        tout << "    /* Create dense SUNLinearSolver object for use by CVode */\n";
        tout << "    SUNLinearSolver LS = SUNLinSol_Dense(y0_, A, sunctx);\n";
        tout << "    if (check_pointer((void*)LS, \"SUNLinSol_Dense()\")) {return 1;}\n";
        tout << endl;
        tout << "    /* Attach the matrix and linear solver */\n";
        tout << "    retval = CVodeSetLinearSolver(cvode_mem, LS, A);\n";
        tout << "    if (check_int_status(retval, \"CVodeSetLinearSolver()\")) {return 1;}\n";
        tout << endl;
        tout << "    /* Set the Jacobian routine */\n";
        tout << "    retval = CVodeSetJacFn(cvode_mem, " << Name() << "_jac);\n";
        tout << "    if (check_int_status(retval, \"CVodeSetJacFn()\")) {return 1;}\n";
        tout << endl;
        tout << "    sunrealtype t1 = solver_param_[2];\n" ;
        tout << "    /* Print the initial condition  */\n";
        tout << "    printf(\"%.8e\", t);\n" ;
        tout << "    for (j = 0; j < N_; ++j) {\n" ;
        tout << "        printf(\" %.8e\", NV_Ith_S(y0_, j));\n";
        tout << "    }\n";
        if ((options["func"] == "yes") && (nf > 0)) {
            tout << "    sunrealtype funcval[" << nf << "];\n";
            tout << "    " << Name() << "_func(t, y0_, funcval, (void *) p_);\n";
            for (size_t i = 0; i < nf; ++i) {
                 tout << "    printf(\" %.8e\",funcval[" << i << "]);\n";
            }
        }
        tout << "    printf(\"\\n\");\n";
        tout << endl;
        tout << "    retval = CVodeSetStopTime(cvode_mem, t1);\n";
        tout << "    if (check_int_status(retval, \"CVodeSetStopTime()\")) {return 1;}\n";
        tout << endl;
        tout << "    while (t < t1) {\n";
        tout << "        /* Advance the solution. */\n";
        tout << "        retval = CVode(cvode_mem, t1, y0_, &t, CV_ONE_STEP);\n";
        tout << "        if (retval != CV_SUCCESS && retval != CV_TSTOP_RETURN) {\n";
        tout << "            fprintf(stderr, \"retval=%d\\n\", retval);\n" ;
        tout << "            retval = -1;\n";
        tout << "            break;\n" ;
        tout << "        }\n" ;
        tout << "        else {\n";
        tout << "            retval = 0;\n";
        tout << "        }\n";
        tout << "        /* Print the solution at the current time */\n";
        tout << "        printf(\"%.8e\", t);\n" ;
        tout << "        for (j = 0; j < N_; ++j) {\n" ;
        tout << "            printf(\" %.8e\", NV_Ith_S(y0_,j));\n" ;
        tout << "        }\n";

        if ((options["func"] == "yes") && (nf > 0)) {
            tout << "        " << Name() << "_func(t, y0_, funcval, (void *) p_);\n";
            for (size_t i = 0; i < nf; ++i) {
                 tout << "        printf(\" %.8e\", funcval[" << i << "]);\n";
            }
        }
        tout << "        printf(\"\\n\");\n";
        tout << "    }\n" ;
        tout << endl;
        tout << "    /* Free memory */\n";
        tout << "    N_VDestroy(y0_);\n";
        tout << "    CVodeFree(&cvode_mem);\n";
        tout << "    SUNLinSolFree(LS);\n";
        tout << "    SUNMatDestroy(A);\n";
        tout << "    SUNContext_Free(&sunctx);\n";
        tout << "    return retval;\n";
        tout << "}\n" ;
        tout.close();
        //
        // Create a Makefile for the CVODE demo program
        //
        string mfilename = "Makefile-"+Name()+"_cv7demo";
        ofstream mout;
        mout.open(mfilename.c_str());
        mout << "#\n";
        mout << "# " << mfilename << endl;
        mout << "#\n";
        mout << "# This is the Makefile for the " << Name() << "_cv7demo program.\n";
        mout << "# This file is configured for CVODE 7.x.\n";
        mout << "#\n";
        PrintVFGENComment(mout,"# ");
        mout << "#\n";
        mout << "# This Makefile is not guaranteed to work in all operating systems.\n";
        mout << "# You may have to edit this file to meet the conventions of your operating system.\n";
        mout << "#\n\n";
        mout << endl;
        mout << "# Set SUNDIALS_DIR externally or by editing this file,\n";
        mout << "# e.g. SUNDIALS_DIR=/usr/local\n";
        mout << "SUNDIALS_LIB_DIR=$(SUNDIALS_DIR)/lib\n";
        mout << "SUNDIALS_INC_DIR=$(SUNDIALS_DIR)/include\n";
        mout << "SUNDIALS_LIBS=-lsundials_cvode -lsundials_nvecserial -lsundials_sunmatrixdense -lsundials_core\n";
        mout << "SUNDIALS_INCS=-I$(SUNDIALS_INC_DIR)\n";
        mout << "LIBS=-lm\n";
        mout << endl;
        mout << "all: " << Name() << "_cv7demo" << endl;
        mout << endl;
        mout << Name() << "_cv7demo: " << Name() << "_cv7demo.o " << Name() << "_cv7.o" << endl;
        mout << "\t$(CC) $(LDFLAGS) -o " << Name() << "_cv7demo ";
        mout <<       Name() << "_cv7demo.o " << Name() << "_cv7.o -L$(SUNDIALS_LIB_DIR) $(SUNDIALS_LIBS) $(LIBS)" << endl;
        mout << endl;
        mout << Name() << "_cv7demo.o: " << Name() << "_cv7demo.c " << Name() << "_cv7.h" << endl;
        mout << "\t$(CC) $(CPPFLAGS) $(SUNDIALS_INCS) -c " << Name() << "_cv7demo.c" << endl;
        mout << endl;
        mout << Name() << "_cv7.o: " << Name() << "_cv7.c " << Name() << "_cv7.h" << endl;
        mout << "\t$(CC) $(CPPFLAGS) $(SUNDIALS_INCS) -c " << Name() << "_cv7.c" << endl;
        mout << endl;
        mout << "clean:\n";
        mout << "\trm -f " << Name() << "_cv7demo " << Name() << "_cv7demo.o " << Name() << "_cv7.o" << endl;
        mout << endl;
        mout.close();

        /*
         * Create a CMake file for the demo.
         */
        string cmfilename = "CMakeLists.txt";
        ofstream cmout;
        cmout.open(cmfilename.c_str());
        cmout << "#\n";
        cmout << "# " << cmfilename << endl;
        cmout << "#\n";
        cmout << "# This is the CMake file for the " << Name() << "_cv7demo program.\n";
        cmout << "#\n";
        cmout << "# If the SUNDIALS library is not installed in /usr or /usr/local,\n";
        cmout << "# set CMAKE_INCLUDE_PATH and CMAKE_LIBRARY_PATH to the appropriate\n";
        cmout << "# paths.\n";
        cmout << "#\n";
        PrintVFGENComment(cmout,"# ");
        cmout << "#\n";
        cmout << endl;
        cmout << "cmake_minimum_required(VERSION 3.22.1)\n";
        cmout << "project(" << Name() << "_cv7demo C)\n";
        cmout << "add_executable(" << Name() << "_cv7demo " << Name() << "_cv7demo.c " << Name() << "_cv7.c)";
        cmout << endl;
        cmout << "target_compile_options(" << Name() << "_cv7demo PUBLIC \"-Wall\" \"-Wextra\")\n";
        cmout << endl;
        cmout << "find_path(CVODE_INCLUDE_DIR cvode/cvode.h /usr/local/include /usr/include)\n";
        cmout << "find_path(NVECTOR_SERIAL_INCLUDE_DIR nvector/nvector_serial.h /usr/local/include /usr/include)\n";
        cmout << "find_path(SUNMATRIX_DENSE_INCLUDE_DIR sunmatrix/sunmatrix_dense.h /usr/local/include /usr/include)\n";
        cmout << "find_path(SUNLINSOL_DENSE_INCLUDE_DIR sunlinsol/sunlinsol_dense.h /usr/local/include /usr/include)\n";
        cmout << endl;
        cmout << "find_library(SUNDIALS_CVODE_LIBRARY sundials_cvode /usr/local/lib /usr/lib)\n";
        cmout << "find_library(SUNDIALS_NVECSERIAL_LIBRARY sundials_nvecserial /usr/local/lib /usr/lib)\n";
        cmout << "find_library(SUNDIALS_CORE_LIBRARY sundials_core /usr/local/lib /usr/lib)\n";
        cmout << "find_library(SUNDIALS_SUNMATRIXDENSE_LIBRARY sundials_sunmatrixdense /usr/local/lib /usr/lib)\n";
        cmout << endl;
        cmout << "include_directories(${CVODE_INCLUDE_DIR} ${NVECTOR_SERIAL_INCLUDE_DIR} ${SUNMATRIX_DENSE_INCLUDE_DIR} ${SUNLINSOL_DENSE_INCLUDE_DIR})\n";
        cmout << endl;
        cmout << "target_link_libraries(" << Name() << "_cv7demo ${SUNDIALS_CVODE_LIBRARY} ${SUNDIALS_NVECSERIAL_LIBRARY} ${SUNDIALS_CORE_LIBRARY} ${SUNDIALS_SUNMATRIXDENSE_LIBRARY} -lm)\n";
        cmout.close();
    }
}
