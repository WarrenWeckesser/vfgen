
//
//  vf_gsl.cpp
//
//  This file defines the VectorField::PrintGSL method.
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


//
// PrintGSL -- The GSL Code Generator.
//

void VectorField::PrintGSL(map<string,string> options)
    {
    int nc, np, nv, na, nf;

    symbol t(IndependentVariable);
    nc = conname_list.nops();
    nv = varname_list.nops();
    np = parname_list.nops();
    na = exprname_list.nops();
    nf = funcname_list.nops();

    string filename = Name()+"_gvf.c";
    ofstream fout;
    fout.open(filename.c_str());
    fout << csrc << left;

    string pfilename = Name()+"_gvf.h";
    ofstream pout;
    pout.open(pfilename.c_str());
    pout << csrc << left;

    //
    //  Print C file header information.
    //
    fout << "/*" << endl;
    fout << " *  " << filename << endl;
    fout << " *" << endl;
    fout << " *  GSL C file for the vector field named: " << Name() << endl;
    fout << " *" << endl;
    PrintVFGENComment(fout," *  ");
    fout << " */" << endl;
    fout << endl;

    pout << "/*" << endl;
    pout << " *  " << pfilename << endl;
    pout << " *" << endl;
    pout << " *  GSL C prototype file for the functions defined in " << filename << endl;
    pout << " *" << endl;
    PrintVFGENComment(pout," *  ");
    pout << " */" << endl;
    pout << endl;

    fout << "#include <math.h>" << endl;
    fout << "#include <gsl/gsl_errno.h>" << endl;
    fout << "#include <gsl/gsl_matrix.h>" << endl;
    fout << "#include <gsl/gsl_odeiv.h>" << endl;
    fout << endl;
    //
    //  Print the vector field function.
    //
    fout << "/*" << endl;
    fout << " *  The vector field." << endl;
    fout << " */" << endl;
    fout << endl;
    fout << "int " << Name() << "_vf(double t, const double y_[], double f_[], void *params)" << endl;
    pout << "int " << Name() << "_vf(double, const double [], double [], void *);" << endl;
    fout << "    {" << endl;
    if (HasPi)
        {
        fout << "    const double Pi = M_PI;\n";
        }
    for (int i = 0; i < nc; ++i)
        {
        fout << "    const double " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    CDeclare_double(fout,varname_list);
    CDeclare_double(fout,parname_list);
    CDeclare_double(fout,exprname_list);
    fout << "    double *p_;" << endl;
    fout << endl;
    fout << "    p_ = (double *) params;" << endl;
    fout << endl;
    GetFromVector(fout,"    ",varname_list,"y_","[]",0,";");
    fout << endl;
    GetFromVector(fout,"    ",parname_list,"p_","[]",0,";");
    fout << endl;
    for (int i = 0; i < na; ++i)
        {
        fout << "    " << exprname_list[i] << " = " << exprformula_list[i] << ";" << endl;
        }
    if (na > 0)
        fout << endl;
    for (int i = 0; i < nv; ++i)
        {
        fout << "    f_[" << i << "]" << " = " << varvecfield_list[i] << ";" << endl;
        }
    fout << endl;
    fout << "    return GSL_SUCCESS;" << endl;
    fout << "    }" << endl;
    fout << endl;
    //
    // Print the Jacobian function.
    //
    fout << "/*" << endl;
    fout << " *  The Jacobian." << endl;
    fout << " */" << endl;
    fout << endl;
    fout << "int " << Name() << "_jac(double t, const double y_[], double *jac_, double *dfdt_, void *params)" << endl;
    pout << "int " << Name() << "_jac(double, const double [], double *, double *, void *);" << endl;
    fout << "    {" << endl;
    if (HasPi)
        {
        fout << "    const double Pi = M_PI;\n";
        }
    for (int i = 0; i < nc; ++i)
        {
        fout << "    const double " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    CDeclare_double(fout,varname_list);
    CDeclare_double(fout,parname_list);
    fout << "    double *p_;" << endl;
    fout << endl;
    fout << "    p_ = (double *) params;" << endl;
    fout << endl;
    GetFromVector(fout,"    ",varname_list,"y_","[]",0,";");
    fout << endl;
    GetFromVector(fout,"    ",parname_list,"p_","[]",0,";");
    fout << endl;
    fout << "    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(jac_," << nv << "," << nv << ");" << endl;
    fout << "    gsl_matrix *m_ = &dfdy_mat.matrix;" << endl;
    fout << endl;
    for (int i = 0; i < nv; ++i)
        {
        ex f = iterated_subs(varvecfield_list[i],expreqn_list);
        for (int j = 0; j < nv; ++j)
            {
            symbol v = ex_to<symbol>(varname_list[j]);
            fout << "    gsl_matrix_set(m_, " << i << ", " << j << ", " << f.diff(v) << ");" << endl;
            }
        }
    fout << endl;
    for (int i = 0; i < nv; ++i)
        {
        ex f = iterated_subs(varvecfield_list[i],expreqn_list);
        fout << "    dfdt_[" << i << "] = " << f.diff(t) << ";" << endl;
        }
    fout << endl;
    fout << "    return GSL_SUCCESS;" << endl;
    fout << "    }" << endl;
    //
    // Print the function that computes the Jacobian with respect
    // to the parameters.
    //
    fout << endl;
    fout << "/*" << endl;
    fout << " *  The Jacobian with respect to the parameters." << endl;
    fout << " */" << endl;
    fout << endl;
    fout << "int " << Name() << "_jacp(double t, const double y_[], double *jacp_, void *params)" << endl;
    pout << "int " << Name() << "_jacp(double, const double [], double *, void *);" << endl;
    fout << "    {" << endl;
    if (HasPi)
        {
        fout << "    const double Pi = M_PI;\n";
        }
    for (int i = 0; i < nc; ++i)
        {
        fout << "    const double " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    CDeclare_double(fout,varname_list);
    CDeclare_double(fout,parname_list);
    fout << "    double *p_;" << endl;
    fout << endl;
    fout << "    p_ = (double *) params;" << endl;
    fout << endl;
    GetFromVector(fout,"    ",varname_list,"y_","[]",0,";");
    fout << endl;
    GetFromVector(fout,"    ",parname_list,"p_","[]",0,";");
    fout << endl;
    fout << "    gsl_matrix_view dfdp_mat = gsl_matrix_view_array(jacp_," << nv << "," << np << ");" << endl;
    fout << "    gsl_matrix *m_ = &dfdp_mat.matrix;" << endl;
    fout << endl;
    for (int i = 0; i < nv; ++i)
        {
        ex f = iterated_subs(varvecfield_list[i],expreqn_list);
        for (int j = 0; j < np; ++j)
            {
            symbol p = ex_to<symbol>(parname_list[j]);
            fout << "    gsl_matrix_set(m_, " << i << ", " << j << ", " << f.diff(p) << ");" << endl;
            }
        }
    fout << endl;
    fout << "    return GSL_SUCCESS;" << endl;
    fout << "    }" << endl;

    if (options["func"] == "yes")
        {
        //
        // Print the user-defined functions.
        //
        for (int n = 0; n < nf; ++n)
            {
            fout << endl;
            fout << "/*" << endl;
            fout << " *  User function: " << funcname_list[n] << endl;
            fout << " */" << endl;
            fout << endl;
            fout << "double " << Name() << "_" << funcname_list[n] << "(double t, const double y_[], void *params)" << endl;
            pout << "double " << Name() << "_" << funcname_list[n] << "(double, const double [], void *);" << endl;
            fout << "    {" << endl;
            if (HasPi)
                {
                fout << "    const double Pi = M_PI;\n";
                }
            for (int i = 0; i < nc; ++i)
                {
                fout << "    const double " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
                }
            CDeclare_double(fout,varname_list);
            CDeclare_double(fout,parname_list);
            CDeclare_double(fout,exprname_list);
            fout << "    double *p_;" << endl;
            fout << endl;
            fout << "    p_ = (double *) params;" << endl;
            fout << endl;
            GetFromVector(fout,"    ",varname_list,"y_","[]",0,";");
            fout << endl;
            GetFromVector(fout,"    ",parname_list,"p_","[]",0,";");
            fout << endl;
            for (int i = 0; i < na; ++i)
                {
                fout << "    " << exprname_list[i] << " = " << exprformula_list[i] << ";" << endl;
                }
            if (na > 0)
                fout << endl;
            fout << "    return " << funcformula_list[n] << ";" << endl;
            fout << "    }" << endl;
            }
        }

    fout.close();
    pout.close();

    if (options["demo"] == "yes")
        {
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
        tout << "    int i;\n" ;
        if (HasPi)
            {
            tout << "    const double Pi = M_PI;\n";
            }
        for (int i = 0; i < nc; ++i)
            {
            tout << "    const double " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
            }
        CDeclare_double(tout,parname_list);
        tout << "    const int P_ = " << np << ";\n" ;
        tout << "    double def_p_[" << np << "] = {" ;
        for (int i = 0; i < np; ++i)
            {
            tout << pardefval_list[i] ;
            if (i != np-1)
                tout << ", " ;
            }
        tout << "};\n" ;
        GetFromVector(tout,"    ",parname_list,"def_p_","[]",0,";");
        tout << "    double p_[" << np << "];\n" ;
        tout << "    const int N_ = " << nv << ";\n" ;
        tout << "    double def_y_[" << nv << "] = {";
        for (int i = 0; i < nv; ++i)
            {
            // tout << def_var_value.at(i) ;
            // tout << "0.0" ;
            tout << vardefic_list[i];
            if (i != nv-1)
                tout << ", " ;
            }
        tout << "};\n" ;
        tout << "    double y_[" << nv << "];\n" ;

        tout << "    double solver_param_[3] = {1.0e-6,0.0,10.0};\n" ;

        MakeCArrayOfStrings(tout,"varnames_",varname_list);
        if (np > 0)
            MakeCArrayOfStrings(tout,"parnames_",parname_list);
        else
            tout << "    char *parnames_[] = {\"\"};\n";
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
        tout << "    const gsl_odeiv_step_type *T_  = gsl_odeiv_step_rk8pd;\n" ;
        tout << "    gsl_odeiv_step    *step_    = gsl_odeiv_step_alloc(T_,N_);\n" ;
        tout << "    gsl_odeiv_control *control_ = gsl_odeiv_control_y_new(solver_param_[0], solver_param_[1]);\n" ;
        tout << "    gsl_odeiv_evolve  *evolve_  = gsl_odeiv_evolve_alloc(N_);\n" ;
        tout << "    gsl_odeiv_system sys_ = {" << Name() << "_vf, " << Name() << "_jac, N_, &(p_[0])};\n";
        tout << endl;
        tout << "    double t_  = 0.0;\n" ;
        tout << "    double t1_ = solver_param_[2];\n" ;
        tout << "    double h_ = 1e-6;\n" ;
        tout << endl;
        tout << "    while (t_ < t1_)\n" ;
        tout << "        {\n" ;
        tout << "        int j_;\n" ;
        tout << "        int status_ = gsl_odeiv_evolve_apply(evolve_, control_, step_, &sys_, &t_, t1_, &h_, y_);\n" ;
        tout << "        if (status_ != GSL_SUCCESS)\n" ;
        tout << "            {\n" ;
        tout << "            fprintf(stderr,\"status=%d\\n\",status_);\n" ;
        tout << "            break;\n" ;
        tout << "            }\n" ;
        tout << "        printf(\"%.8e\",t_);\n" ;
        tout << "        for (j_ = 0; j_ < N_; ++j_)\n" ;
        tout << "            printf(\" %.8e\",y_[j_]);\n" ;
        if (options["func"] == "yes")
            for (int i = 0; i < nf; ++i)
                 tout << "        printf(\" %.8e\"," << Name() << "_" << funcname_list[i] << "(t_,y_,p_));\n";
        tout << "        printf(\"\\n\");\n";
        tout << "        }\n" ;
        tout << endl;
        tout << "    gsl_odeiv_evolve_free(evolve_);\n" ;
        tout << "    gsl_odeiv_control_free(control_);\n" ;
        tout << "    gsl_odeiv_step_free(step_);\n" ;
        tout << "    return 0;\n" ;
        tout << "    }\n" ;
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
    }
