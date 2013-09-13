
//
//  vf_cvode.cpp
//
//  This file defines the VectorField::PrintCVODE method.
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
#include <string>
#include <ginac/ginac.h> 

#include "vf.h"
#include "vf_utils.h"

using namespace std;
using namespace GiNaC;

//
// PrintCVODE -- The CVODE Code Generator
//

void VectorField::PrintCVODE(map<string,string> options)
{

    int nc, np, nv, na, nf;

    if (options.count("version") > 0 && options["version"] != "2.3.0" && options["version"] != "2.4.0"
            && options["version"] != "2.5.0" && options["version"] != "2.6.0" && options["version"] != "2.7.0") {
        cerr << "vfgen CVODE command: unknown version specified: " << options["version"] << endl;
        cerr << "Versions of CVODE supported by VFGEN are 2.3.0, 2.4.0, 2.5.0 and 2.6.0. Default: version=2.6.0" << endl;
        exit(-1);
    } 
    if (options.count("version") == 0) {
        // Explicitly set the default value.
        options["version"] = "2.7.0";
    }

    symbol t(IndependentVariable);
    nc = conname_list.nops();
    nv = varname_list.nops();
    np = parname_list.nops();
    na = exprname_list.nops();
    nf = funcname_list.nops();

    string filename = Name()+"_cv.c";
    ofstream fout;
    fout.open(filename.c_str());
    fout << csrc << left;

    string pfilename = Name()+"_cv.h";
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

    fout << "#include <math.h>" << endl;
    fout << endl;
    if (options["version"] == "2.3.0") {
        fout << "/* Include headers for CVODE 2.3.0 */" << endl;
        fout << "#include <sundialstypes.h>" << endl;
        fout << "#include <cvode.h>" << endl;
        fout << "#include <cvdense.h>" << endl;
        fout << "#include <nvector_serial.h>" << endl;
        fout << "#include <dense.h>" << endl;
    }
    else if (options["version"] == "2.4.0") {
        fout << "/* Include headers for CVODE 2.4.0 */" << endl;
        fout << "#include <sundials/sundials_types.h>" << endl;
        fout << "#include <sundials/sundials_dense.h>" << endl;
        fout << "#include <cvode.h>" << endl;
        fout << "#include <cvode/cvode_dense.h>" << endl;
        fout << "#include <nvector_serial.h>" << endl;
    }
    else if (options["version"] == "2.5.0") {
        fout << "/* Include headers for CVODE 2.5.0 */" << endl;
        fout << "#include <sundials/sundials_types.h>" << endl;
        fout << "#include <sundials/sundials_dense.h>" << endl;
        fout << "#include <sundials/sundials_nvector.h>" << endl;
        fout << "#include <nvector/nvector_serial.h>" << endl;
        fout << "#include <cvode/cvode.h>" << endl;
        fout << "#include <cvode/cvode_dense.h>" << endl;
    }
    else {
        // CVODE 2.6.0 or 2.7.0
        fout << "/* Include headers for CVODE " << options["version"] << " */" << endl;
        fout << "#include <cvode/cvode.h>" << endl;
        fout << "#include <nvector/nvector_serial.h>" << endl;
        fout << "#include <cvode/cvode_dense.h>" << endl;
    }
    fout << endl;
    //
    //  Print the vector field function.
    //
    fout << "/*" << endl;
    fout << " *  The vector field." << endl;
    fout << " */" << endl;
    fout << endl;
    string func_return_type = "int";
    if (options["version"] == "2.3.0") {
        func_return_type = "void";
    }
    fout << func_return_type << " " << Name() << "_vf(realtype t, N_Vector y_, N_Vector f_, void *params)" << endl;
    pout << func_return_type << " " << Name() << "_vf(realtype, N_Vector, N_Vector, void *);" << endl;
    fout << "    {" << endl;
    if (HasPi) {
        fout << "    const realtype Pi = RCONST(M_PI);\n";
    }
    for (int i = 0; i < nc; ++i) {
        fout << "    const realtype " << conname_list[i] << " = RCONST(" << convalue_list[i] << ");" << endl;
    }
    CDeclare(fout,"realtype",varname_list);
    CDeclare(fout,"realtype",parname_list);
    CDeclare(fout,"realtype",exprname_list);
    fout << "    realtype *p_;" << endl;
    fout << endl;
    fout << "    p_ = (realtype *) params;" << endl;
    fout << endl;
    // GetFromVector(fout,"    ",varname_list,"y_","[]",0,";");
    for (int i = 0; i < nv; ++i) {
        fout << "    ";
        fout.width(10);
        fout << varname_list[i];
        fout.width(0);
        fout << " = NV_Ith_S(y_," << i << ");" << endl;
    }

    fout << endl;
    GetFromVector(fout,"    ",parname_list,"p_","[]",0,";");
    fout << endl;
    for (int i = 0; i < na; ++i) {
        fout << "    " << exprname_list[i] << " = " << exprformula_list[i] << ";" << endl;
    }
    if (na > 0) {
        fout << endl;
    }
    for (int i = 0; i < nv; ++i) {
        fout << "    NV_Ith_S(f_," << i << ") = " << varvecfield_list[i] << ";" << endl;
    }
    if (options["version"] != "2.3.0") {
        fout << "    return 0;\n";
    }
    fout << "    }" << endl;
    fout << endl;
    //
    // Print the Jacobian function.
    //
    fout << "/*" << endl;
    fout << " *  The Jacobian." << endl;
    fout << " */" << endl;
    fout << endl;
    if (options["version"] == "2.6.0" || options["version"] == "2.7.0") {
        string aux_type;
        if (options["version"] == "2.7.0") {
            aux_type = "long ";
        }
        else {
            aux_type = "";
        }
        fout << func_return_type << " " << Name() << "_jac(" << aux_type << "int N_, realtype t, N_Vector y_, N_Vector fy_,";
        fout << " DlsMat jac_, void *params," << endl;
        fout << "                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)" << endl;

        pout << func_return_type << " " << Name() << "_jac(" << aux_type << "int, realtype, N_Vector, N_Vector,";
        pout << " DlsMat, void *," << endl;
        pout << "                N_Vector, N_Vector, N_Vector);" << endl;
    }
    else {
        fout << func_return_type << " " << Name() << "_jac(long int N_, DenseMat jac_, realtype t," << endl;
        fout << "                N_Vector y_, N_Vector fy_, void *params," << endl;
        fout << "                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)" << endl;

        pout << func_return_type << " " << Name() << "_jac(long int, DenseMat, realtype," << endl;
        pout << "                N_Vector, N_Vector, void *," << endl;
        pout << "                N_Vector, N_Vector, N_Vector);" << endl;
    }
    fout << "    {" << endl;
    if (HasPi) {
        fout << "    const realtype Pi = RCONST(M_PI);\n";
    }
    for (int i = 0; i < nc; ++i) {
        fout << "    const realtype " << conname_list[i] << " = RCONST(" << convalue_list[i] << ");" << endl;
    }
    CDeclare(fout,"realtype",varname_list);
    CDeclare(fout,"realtype",parname_list);
    fout << "    realtype *p_;" << endl;
    fout << endl;
    fout << "    p_ = (realtype *) params;" << endl;
    fout << endl;
    // GetFromVector(fout,"    ",varname_list,"y_","[]",0,";");
    for (int i = 0; i < nv; ++i) {
        fout << "    ";
        fout.width(10);
        fout << varname_list[i];
        fout.width(0);
        fout << " = NV_Ith_S(y_," << i << ");" << endl;
    }
    fout << endl;
    GetFromVector(fout,"    ",parname_list,"p_","[]",0,";");
    fout << endl;
    for (int i = 0; i < nv; ++i) {
        ex f = iterated_subs(varvecfield_list[i],expreqn_list);
        for (int j = 0; j < nv; ++j) {
            symbol v = ex_to<symbol>(varname_list[j]);
            ex df = f.diff(v);
            // Skip zero elements.  CVODE initializes jac_ to zero before calling the Jacobian function.
            if (df != 0)
                fout << "    DENSE_ELEM(jac_, " << i << ", " << j << ") = " << f.diff(v) << ";" << endl;
        }
    }
    if (options["version"] != "2.3.0") {
        fout << "    return 0;\n";
    }
    fout << "    }" << endl;

    if (options["func"] == "yes" & nf > 0) {
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
        fout << func_return_type << " " << Name() << "_func(realtype t, N_Vector y_, realtype *func_, void *params)" << endl;
        pout << func_return_type << " " << Name() << "_func(realtype, N_Vector, realtype *, void *);" << endl;
        fout << "    {" << endl;
        if (HasPi) {
            fout << "    const realtype Pi = RCONST(M_PI);\n";
        }
        for (int i = 0; i < nc; ++i) {
            fout << "    const realtype " << conname_list[i] << " = RCONST(" << convalue_list[i] << ");" << endl;
        }
        CDeclare(fout,"realtype",varname_list);
        CDeclare(fout,"realtype",parname_list);
        CDeclare(fout,"realtype",exprname_list);
        fout << "    realtype *p_;" << endl;
        fout << endl;
        fout << "    p_ = (realtype *) params;" << endl;
        fout << endl;
        // GetFromVector(fout,"    ",varname_list,"y_","[]",0,";");
        for (int i = 0; i < nv; ++i) {
            fout << "    ";
            fout.width(10);
            fout << varname_list[i];
            fout.width(0);
            fout << " = NV_Ith_S(y_," << i << ");" << endl;
        }
        fout << endl;
        GetFromVector(fout,"    ",parname_list,"p_","[]",0,";");
        fout << endl;
        for (int i = 0; i < na; ++i) {
            fout << "    " << exprname_list[i] << " = " << exprformula_list[i] << ";" << endl;
        }
        if (na > 0) {
            fout << endl;
        }
        for (int n = 0; n < nf; ++n) {
            fout << "    /* " << funcname_list[n] << ":  */" << endl;
            fout << "    func_[" << n << "] = " << funcformula_list[n] << ";" << endl;
        }
        fout << "    return 0;\n";
        fout << "    }" << endl;
    }
    fout.close();
    pout.close();

    if (options["demo"] == "yes") {
        //
        //  Create a self-contained ODE solver for this vector field
        //  that allows the user to give the initial conditions,
        //  parameter values, and some solver control parameters
        //  on the command line.
        //

        string tfilename = Name()+"_cvdemo.c";
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
        tout << "#include <string.h>" << endl;
        tout << "#include <math.h>" << endl;
        tout << endl;
        if (options["version"] == "2.3.0") {
            tout << "/* Include headers for CVODE v2.3.0 */" << endl;
            tout << "#include <sundialstypes.h>" << endl;
            tout << "#include <cvode.h>" << endl;
            tout << "#include <cvdense.h>" << endl;
            tout << "#include <nvector_serial.h>" << endl;
            tout << "#include <dense.h>" << endl;
        }
        else if (options["version"] == "2.4.0") {
            tout << "/* Include headers for CVODE 2.4.0 */" << endl;
            tout << "#include <sundials/sundials_types.h>" << endl;
            tout << "#include <sundials/sundials_dense.h>" << endl;
            tout << "#include <cvode.h>" << endl;
            tout << "#include <cvode/cvode_dense.h>" << endl;
            tout << "#include <nvector_serial.h>" << endl;
        }
        else if (options["version"] == "2.5.0") {
            tout << "/* Include headers for CVODE 2.5.0 */" << endl;
            tout << "#include <sundials/sundials_types.h>" << endl;
            tout << "#include <sundials/sundials_dense.h>" << endl;
            tout << "#include <sundials/sundials_nvector.h>" << endl;
            tout << "#include <nvector/nvector_serial.h>" << endl;
            tout << "#include <cvode/cvode.h>" << endl;
            tout << "#include <cvode/cvode_dense.h>" << endl;
        }
        else {
            // CVODE >= 2.6.0
            tout << "/* Include headers for CVODE " << options["version"] << " */" << endl;
            tout << "#include <cvode/cvode.h>" << endl;
            tout << "#include <nvector/nvector_serial.h>" << endl;
            tout << "#include <cvode/cvode_dense.h>" << endl;
        }

        tout << endl;
        tout << "#include \"" << pfilename << "\"\n";
        tout << endl << endl;
        tout << "int use(int argc, char *argv[], int nv, char *vname[], double y_[], int np, char *pname[], const double p_[])\n" ;
        tout << "    {\n" ;
        tout << "    int i;\n" ;
        tout << "    printf(\"use: %s [options]\\n\",argv[0]);\n" ;
        tout << "    printf(\"options:\\n\");\n" ;
        tout << "    printf(\"    -h    Print this help message.\\n\");\n" ;
        tout << "    for (i = 0; i < nv; ++i)\n" ;
        tout << "        printf(\"    %s=<initial_condition>   Default value is %e\\n\",vname[i],y_[i]);\n";
        tout << "    for (i = 0; i < np; ++i)\n" ;
        tout << "        printf(\"    %s=<parameter_value>   Default value is %e\\n\",pname[i],p_[i]);\n";
        tout << "    printf(\"    abserr=<absolute_error_tolerance>\\n\");\n" ;
        tout << "    printf(\"    relerr=<relative_error_tolerance>\\n\");\n" ;
        tout << "    printf(\"    stoptime=<stop_time>\\n\");\n" ;
        tout << "    return 0;\n";
        tout << "    }\n" ;
        tout << endl << endl;
        tout << "int assign(char *str[], int ns, double v[], char *a)\n" ;
        tout << "    {\n" ;
        tout << "    int i;\n" ;
        tout << "    char name[256];\n" ;
        tout << "    char *e;\n" ;
        tout << "\n" ;
        tout << "    e = strchr(a,'=');\n" ;
        tout << "    if (e == NULL)\n" ;
        tout << "        return(-1);\n" ;
        tout << "    *e = '\\0';\n" ;
        tout << "    strcpy(name,a);\n" ;
        tout << "    *e = '=';\n" ;
        tout << "    ++e;\n" ;
        tout << "    for (i = 0; i < ns; ++i)\n" ;
        tout << "        if (strcmp(str[i],name)==0)\n" ;
        tout << "            break;\n" ;
        tout << "    if (i == ns)\n" ;
        tout << "        return(-1);\n" ;
        tout << "    v[i] = atof(e);\n" ;
        tout << "    return(i);\n" ;
        tout << "    }\n" ;
        tout << endl << endl;
        tout << "int main (int argc, char *argv[])\n" ;
        tout << "    {\n" ;
        tout << "    int i,j;\n";
        tout << "    int flag;\n";
        tout << "    const int N_ = " << nv << ";\n" ;
        tout << "    const int P_ = " << np << ";\n" ;
        if (HasPi) {
            tout << "    const realtype Pi = RCONST(M_PI);\n";
        }
        for (int i = 0; i < nc; ++i) {
            tout << "    const realtype " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
        tout << "    const realtype def_p_[" << np << "] = \n" ;
        tout << "        {\n";
        for (int i = 0; i < np; ++i) {
            tout << "        RCONST(" << pardefval_list[i] << ")" ;
            if (i != np-1) {
                tout << "," ;
            }
            tout << endl;
        }
        tout << "        };\n" ;
        // CDeclare(tout,"realtype",parname_list);
        GetFromVector(tout,"    const realtype ",parname_list,"def_p_","[]",0,";");
        tout << "    realtype def_y_[" << nv << "] = {";
        for (int i = 0; i < nv; ++i) {
            // tout << def_var_value.at(i) ;
            // tout << "RCONST(0.0)" ;
            tout << "RCONST(" << vardefic_list[i] << ")";
            if (i != nv-1) {
                tout << ", " ;
            }
        }
        tout << "};\n" ;
        tout << "    realtype y_[" << nv << "];\n" ;

        tout << "    realtype p_[" << np << "];\n" ;
        tout << "    realtype solver_param_[3] = {RCONST(1.0e-6),RCONST(0.0),RCONST(10.0)};\n" ;
        MakeCArrayOfStrings(tout,"varnames_",varname_list);
        if (np > 0) {
            MakeCArrayOfStrings(tout,"parnames_",parname_list);
        }
        else {
            tout << "    char *parnames_[] = {\"\"};\n";
        }
        tout << "    char *solver_param_names_[3] = {\"abserr\",\"relerr\",\"stoptime\"};\n" ;
        tout << endl;
        tout << "    for (i = 0; i < N_; ++i)\n" ;
        tout << "        y_[i] = def_y_[i];\n" ;
        tout << "    for (i = 0; i < P_; ++i)\n" ;
        tout << "        p_[i] = def_p_[i];\n" ;
        tout << "    for (i = 1; i < argc; ++i)\n" ;
        tout << "        {\n" ;
        tout << "        int j;\n" ;
        tout << "        if (strcmp(argv[i],\"-h\") == 0)\n" ;
        tout << "            {\n" ;
        tout << "            use(argc,argv,N_,varnames_,def_y_,P_,parnames_,def_p_);\n" ;
        tout << "            exit(0);\n" ;
        tout << "            }\n" ;
        tout << "        j = assign(varnames_,N_,y_,argv[i]);\n" ;
        tout << "        if (j == -1)\n" ;
        tout << "            {\n" ;
        tout << "            j = assign(parnames_,P_,p_,argv[i]);\n" ;
        tout << "            if (j == -1)\n" ;
        tout << "                {\n" ;
        tout << "                j = assign(solver_param_names_,3,solver_param_,argv[i]);\n" ;
        tout << "                if (j == -1)\n" ;
        tout << "                    {\n" ;
        tout << "                    fprintf(stderr,\"unknown argument: %s\\n\",argv[i]);\n" ;
        tout << "                    use(argc,argv,N_,varnames_,def_y_,P_,parnames_,def_p_); \n";
        tout << "                    exit(-1);\n" ;
        tout << "                    }\n" ;
        tout << "                }\n" ;
        tout << "            }\n" ;
        tout << "        }\n" ;
        tout << endl;
        tout << "    N_Vector y0_;\n";
        tout << "    y0_ = N_VNew_Serial(N_);\n";
        tout << "    for (i = 0; i < N_; ++i)\n";
        tout << "        NV_Ith_S(y0_,i) = y_[i];\n";
        tout << endl;

        tout << "    /* For non-stiff problems:   */\n";
        tout << "    void *cvode_mem = CVodeCreate(CV_ADAMS,CV_FUNCTIONAL);\n";
        tout << "    /* For stiff problems:       */\n";
        tout << "    /* void *cvode_mem = CVodeCreate(CV_BDF,CV_NEWTON); */\n";
        tout << endl;

        tout << "    realtype t = RCONST(0.0);\n";

        if (options["version"] == "2.6.0" || options["version"] == "2.7.0") {
            tout << "    flag = CVodeInit(cvode_mem," << Name() << "_vf, t, y0_);\n";
            tout << "    flag = CVodeSStolerances(cvode_mem, solver_param_[1], solver_param_[0]);\n";
            tout << "    flag = CVodeSetUserData(cvode_mem, &(p_[0]));\n";
        }
        else {
            tout << "    flag = CVodeMalloc(cvode_mem," << Name() << "_vf, t, y0_, CV_SS, solver_param_[1], &(solver_param_[0]));\n";
            tout << "    flag = CVodeSetFdata(cvode_mem, &(p_[0]));\n";
        }
        tout << "    flag = CVDense(cvode_mem, N_);\n";
        if (options["version"] == "2.6.0" || options["version"] == "2.7.0") {
            tout << "    flag = CVDlsSetDenseJacFn(cvode_mem, " << Name() << "_jac);\n";
        }
        else {
            tout << "    flag = CVDenseSetJacFn(cvode_mem, " << Name() << "_jac, &(p_[0]));\n";
        }
        tout << endl;
        tout << "    realtype t1 = solver_param_[2];\n" ;
        tout << "    /* Print the initial condition  */\n";
        tout << "    printf(\"%.8e\", t);\n" ;
        tout << "    for (j = 0; j < N_; ++j)\n" ;
        tout << "        printf(\" %.8e\", NV_Ith_S(y0_,j));\n" ;
        if (options["func"] == "yes" & nf > 0) {
            tout << "    realtype funcval[" << nf << "];\n";
            tout << "    " << Name() << "_func(t, y0_, funcval, (void *) p_);\n";
            for (int i = 0; i < nf; ++i) {
                 tout << "    printf(\" %.8e\",funcval[" << i << "]);\n";
            }
        }
        tout << "    printf(\"\\n\");\n";
        tout << "    flag = CVodeSetStopTime(cvode_mem, t1);\n";
        tout << "    while (t < t1)\n";
        tout << "        {\n" ;
        tout << "        /* Advance the solution */\n";
        if (options["version"] == "2.6.0" || options["version"] == "2.7.0") {
            tout << "        flag = CVode(cvode_mem, t1, y0_, &t, CV_ONE_STEP);\n";
        }
        else {
            tout << "        flag = CVode(cvode_mem, t1, y0_, &t, CV_ONE_STEP_TSTOP);\n";
        }        
        tout << "        if (flag != CV_SUCCESS && flag != CV_TSTOP_RETURN)\n";
        tout << "            {\n" ;
        tout << "            fprintf(stderr, \"flag=%d\\n\", flag);\n" ;
        tout << "            break;\n" ;
        tout << "            }\n" ;
        tout << "        /* Print the solution at the current time */\n";
        tout << "        printf(\"%.8e\", t);\n" ;
        tout << "        for (j = 0; j < N_; ++j)\n" ;
        tout << "            printf(\" %.8e\", NV_Ith_S(y0_,j));\n" ;

        if (options["func"] == "yes" & nf > 0) {
            tout << "        " << Name() << "_func(t, y0_, funcval, (void *) p_);\n";
            for (int i = 0; i < nf; ++i) {
                 tout << "        printf(\" %.8e\", funcval[" << i << "]);\n";
            }
        }
        tout << "        printf(\"\\n\");\n";
        tout << "        }\n" ;
        tout << endl;
        tout << "    N_VDestroy_Serial(y0_);\n";
        if (options["version"] == "2.3.0") {
            tout << "    CVodeFree(cvode_mem);\n";
        }
        else {
            tout << "    CVodeFree(&cvode_mem);\n";
        }
        tout << "    }\n" ;
        tout.close();
        //
        // Create a Makefile for the CVODE demo program
        //
        string mfilename = "Makefile-"+Name()+"_cvdemo";
        ofstream mout;
        mout.open(mfilename.c_str());
        mout << "#\n";
        mout << "# " << mfilename << endl;
        mout << "#\n";
        mout << "# This is the Makefile for the " << Name() << "_cvdemo program.\n";
        if (options["version"] == "2.3.0") {
            mout << "# This file is configured for CVODE 2.3.0.\n";
        }
        else if (options["version"] == "2.4.0") {
            mout << "# This file is configured for CVODE 2.4.0.\n";
        }
        else {
            // CVODE >= 2.5.0
            mout << "# This file is configured for CVODE " << options["version"] << ".\n";
        }
        mout << "#\n";
        PrintVFGENComment(mout,"# ");
        mout << "#\n";
        mout << "# This Makefile is not guaranteed to work in all operating systems.\n";
        mout << "# You may have to edit this file to meet the conventions of your operating system.\n";
        mout << "#\n\n";
        mout << endl;
        mout << "SUNDIALS_DIR=/usr/local\n";
        mout << "SUNDIALS_LIB_DIR=$(SUNDIALS_DIR)/lib\n";
        mout << "SUNDIALS_INC_DIR=$(SUNDIALS_DIR)/include\n";
        if (options["version"] == "2.3.0") {
            mout << "SUNDIALS_LIBS=-lsundials_cvode -lsundials_nvecserial -lsundials_shared\n";
        }
        else {
            mout << "SUNDIALS_LIBS=-lsundials_cvode -lsundials_nvecserial\n";
        }
        if (options["version"] == "2.4.0") {
            mout << "SUNDIALS_INCS=-I$(SUNDIALS_INC_DIR) -I$(SUNDIALS_INC_DIR)/cvode -I$(SUNDIALS_INC_DIR)/sundials\n";
        }
        else {
            // CVODE >= 2.5.0
            mout << "SUNDIALS_INCS=-I$(SUNDIALS_INC_DIR)\n";
        }
        mout << "LIBS=-lm\n";
        mout << endl;
        mout << "all: " << Name() << "_cvdemo" << endl;
        mout << endl;
        mout << Name() << "_cvdemo: " << Name() << "_cvdemo.o " << Name() << "_cv.o" << endl;
        mout << "\t$(CC) $(LDFLAGS) -o " << Name() << "_cvdemo ";
        mout <<       Name() << "_cvdemo.o " << Name() << "_cv.o -L$(SUNDIALS_LIB_DIR) $(SUNDIALS_LIBS) $(LIBS)" << endl;
        mout << endl;
        mout << Name() << "_cvdemo.o: " << Name() << "_cvdemo.c " << Name() << "_cv.h" << endl;
        mout << "\t$(CC) $(CPPFLAGS) $(SUNDIALS_INCS) -c " << Name() << "_cvdemo.c" << endl;
        mout << endl;
        mout << Name() << "_cv.o: " << Name() << "_cv.c " << Name() << "_cv.h" << endl;
        mout << "\t$(CC) $(CPPFLAGS) $(SUNDIALS_INCS) -c " << Name() << "_cv.c" << endl;
        mout << endl;
        mout << "clean:\n";
        mout << "\trm -f " << Name() << "_cvdemo " << Name() << "_cvdemo.o " << Name() << "_cv.o" << endl;
        mout << endl;
        mout.close();
    }
}
