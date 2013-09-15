//
//  vfgen.cpp -- Multi-format vector field file generator.
//
//  by Warren Weckesser
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
#include <vector>
#include <map>
#include <ginac/ginac.h>

#include "vf.h"

using namespace std;
using namespace GiNaC;


//
// This function overrides the default print method for the Zlags_ function.
// Ginac thinks Zlags_ is a function, but in fact, in the code generated in
// a couple of the VFGEN commands that handle delay equations, Zlags_ is a
// two-dimensional array.  So we want Zlags_(1,2) printed like that, not as
// Zlags_(1.0,2.0).
//

static void Zlags_print(const ex& arg1, const ex& arg2, const print_context& c) {
    c.s << "Zlags_(";
    if (is_a<numeric>(arg1)) {
        c.s << ex_to<numeric>(arg1).to_int();
    }
    else {
        arg1.print(c);
    }
    c.s << ",";
    if (is_a<numeric>(arg2)) {
        c.s << ex_to<numeric>(arg2).to_int();
    }
    else {
        arg2.print(c);
    }
    c.s << ")";
}

//
// Add a function called "delay" to the functions known by ginac.
// (The corresponding macro DECLARE_FUNCTION_2P(delay) is in
// ginac_declare_funcs.h.)
//

REGISTER_FUNCTION(delay, dummy())
//REGISTER_FUNCTION(Zlags_, dummy())
REGISTER_FUNCTION(Zlags_, print_func<print_csrc_float>(Zlags_print).
                          print_func<print_csrc_double>(Zlags_print).
                          print_func<print_python>(Zlags_print) )

#define NAMEWIDTH 9

const char *commands[] = {
        "adolc",
        "auto",
        "check",
        "cvode",
        "dde23",
        "ddebiftool",
        "dde_solver",
        "delay2ode",
        "dstool",
        "evf",
        "gsl",
        "help",
        "javascript",
        "latex",
        "lsoda",
        "matcont",
        "matlab",
        "octave",
        "pddecont",
        "pydstool",
        "pygsl",
        "r",
        "radau5",
        "scilab",
        "scipy",
        "taylor",
        "xml",
        "xpp",
        "end"};
        
map<string,vector<string> > command_options;

int checkcommand(const char *s)
{
    int i;

    i = 0;
    while (strcmp(commands[i],"end") != 0) {
        if (strcmp(s,commands[i]) == 0) {
            return i;
        }
        else {
            i = i + 1;
        }
    }
    return -1;
}
    
void printcommands(ostream &c)
{
    int i = 0;
    c << "    ";
    while (strcmp(commands[i],"end") != 0) {
        if (i > 0) {
            c << ", ";
        }
        if (i > 0 && i % 8 == 0) {
            c << "\n    ";
        }
        c << commands[i];
        i = i + 1;
    }
    c << endl;
}


void printuse()
{
    cerr << "VFGEN (Version:" << VERSION << ")" << endl;
    cerr << "Use: vfgen command  vector-field-file.vf" << endl;
    cerr << "or:  vfgen command:option=value,...,option=value vector-field-file.vf" << endl;
    cerr << "or:  vfgen help command\n";
    cerr << "where command is one of:\n";
    printcommands(cerr);
}

int help(char *command)
{
    int m;
    m = checkcommand(command);
    if (m < 0) {
        cout << "VFGEN error: \"" << command << "\" is not a valid command.\n";
        cout << "VFGEN commands are:\n";
        printcommands(cout);
        return -1;
    }
    if (strcmp(command,"adolc") == 0) {
        cout << "use: vfgen adolc vector_field_file.vf\n";
        cout << endl;
        cout << "This command creates a C++ function that can be used with the ADOL-C library.\n";
        cout << "Two files are created:\n";
        cout << "    [name]_adolc.cpp\n";
        cout << "    [name]_adolc.h\n";
        cout << "The file [name]_adolc.cpp contains the function\n";
        cout << "    [name]_vf(short int tag, double *y_, double *f_, double *params_)\n";
        cout << "The code in this function uses the ADOL-C data types adouble and adoublev,\n";
        cout << "and wraps the actual computation of the vector field inside trace_on(tag)\n";
        cout << "and trace_off(tag) statement.\n";
        cout << endl;
        cout << "Options: none.\n";
    }
    else if (strcmp(command,"auto") == 0) {
        cout << "use:  vfgen auto vector_field_file.vf\n";
        cout << endl;
        cout << "This command generates a C or FORTRAN file to be used with the AUTO\n";
        cout << "continuation software.  A single file is created, called [name]_avf.c\n";
        cout << "(or [name]_avf.f if the option lang=fortran is used) where [name] is the\n";
        cout << "Name attribute of the VectorField entity in the vector field file.\n";
        cout << endl;
        cout << "Options: (default is listed first)\n";
        cout << "lang=c|fortran\n";
        cout << "    Specify the language in which to write the AUTO equations file.\n";
        cout << "    For AUTO2000, use C (the default).  For AUTO07p, either C or FORTRAN\n";
        cout << "    may be used.\n";
    }
    else if (strcmp(command,"check") == 0) {
        cout << "use: vfgen check vector_field_file.vf\n";
        cout << endl;
        cout << "This command prints information from the vector field file.  It can be used to\n";
        cout << "to check for errors before using another command.\n";
        cout << endl;
        cout << "Options: none.\n";
    }
    else if (strcmp(command,"cvode") == 0) {
        cout << "use: vfgen cvode vector_field_file.vf\n";
        cout << "     vfgen cvode:option=value,...,option=value vector_field_file.vf\n";
        cout << endl;
        cout << "This command generates code to be used with the CVODE library (which is part\n";
        cout << "of the SUNDIALS suite).\n";
        cout << "With no options, two files are created:\n";
        cout << "    [name]_cv.c\n";
        cout << "    [name]_cv.h\n";
        cout << "where [name] is the Name attribute of the VectorField entity in the vector\n";
        cout << "field file.  The file [name]_cv.c will define the C functions [name]_vf() and\n";
        cout << "[name]_jac().\n";
        cout << endl;
        cout << "Options: (default is listed first)\n";
        cout << "version=2.7.0|2.6.0|2.5.0|2.4.0|2.3.0\n";
        cout << "    This option determines of version of CVODE for which code is generated.\n";
        cout << "func=no|yes\n";
        cout << "    By default, any user function defined in the vector field file is not\n";
        cout << "    defined in [name]_cv.c.  If the option func=yes is given, the function\n";
        cout << "    [name]_func() will be defined.  This function will fill in an array with\n";
        cout << "    values that are the user defined functions.\n";
        cout << "demo=no|yes\n";
        cout << "    If the option demo=yes is given, two more files are created:\n";
        cout << "        [name]_cvdemo.c\n";
        cout << "        Makefile-[name]_cvdemo\n";
        cout << "    where [name] is the Name attribute of the VectorField entity in the\n";
        cout << "    vector field file.  This program provides a simple command-line solver.\n";
        cout << "    It takes arguments in the form name=value, where name can be a state\n";
        cout << "    variable (to give an initial condition), a parameter (to give the value\n";
        cout << "    of a parameter), or one of abserr, relerr or stoptime, to control the ODE\n";
        cout << "    solver.  The output of the program consists of columns of numbers; the\n";
        cout << "    first column is the time, and the rest are the state variables.\n";
    }
    else if (strcmp(command,"dde23") == 0) {
        cout << "use: vfgen dde23 vector_field_file.vf\n";
        cout << "     vfgen dde23:option=value,...,option=value vector_field_file.vf\n";
        cout << endl;
        cout << "This command creates a MATLAB file that can be used with the DDE23 delay\n";
        cout << "equation solver.  The file is [name]_dde23.m, where [name] is the Name\n";
        cout << "attribute of the VectorField entity in the vector field file.\n";
        cout << endl;
        cout << "Options: (default is listed first)\n";
        cout << "parstyle=vector|list\n";
        cout << "    This option controls how parameters are passed to the function defined in\n";
        cout << "    the file [name]_dde23.m.  By default (parstyle=vector), the parameters are\n";
        cout << "    stored in a single array.  If the option parstyle=list is given, the\n";
        cout << "    parameters are all listed separately as arguments.\n";
        cout << "demo=no|yes\n";
        cout << "    If the option demo=yes is given, the file [name]_dde23_demo.m is created.\n";
        cout << "    This defines the MATLAB function [name]_dde23_demo(stoptime), which uses\n";
        cout << "    the DDE23 function with the default parameter values and initial conditions\n";
        cout << "    to plot a solution to the system.\n";
    }
    else if (strcmp(command,"ddebiftool") == 0) {
        cout << "use: vfgen ddebiftool vector_field_file.vf\n";
        cout << "     vfgen ddebiftool:path=/path/to/ddebiftool vector_field_file.vf\n";
        cout << endl;
        cout << "This command creates the files sys_init.m, sys_rhs.m, sys_deri.m and sys_tau.m,\n";
        cout << "to be used with the MATLAB software package DDE-BIFTOOL.\n";
        cout << "Note: Currently, VFGEN only generates DDE-BIFTOOL files for systems with\n";
        cout << "constant delays.\n";
        cout << endl;
        cout << "Options:\n";
        cout << "path=/path/to/ddebiftool\n";
        cout << "    The 'path' option allows the user to specify a directory that is added to\n";
        cout << "    the MATLAB search path in the file sys_init.m.\n";  
    }
    else if (strcmp(command,"dde_solver") == 0) {
        cout << "use: vfgen dde_solver vector_field_file.vf\n";
        cout << "     vfgen dde_solver:demo=yes vector_field_file.vf\n";
        cout << endl;
        cout << "This command creates a Fortran 90 file that can be used with the DDE solver\n";
        cout << "written by S. Thompson and L. F. Shampine.  The file created is [name].f90.\n";
        cout << endl;
        cout << "Options: (default is listed first)\n";
        cout << "demo=no|yes\n";
        cout << "    If the option demo=yes is given, the file [name]_demo.f90 is created.\n";
        cout << "    This program will use [name].f90 to generate a solution to the delay\n";
        cout << "    equations.\n";
    }
    else if (strcmp(command,"delay2ode") == 0) {
        cout << "use: vfgen delay2ode vector_field_file.vf\n";
        cout << "     vfgen delay2ode:option=value,...,option=value vector_field_file.vf\n";
        cout << endl;
        cout << "This command generates a new vector field file that is a finite dimensional\n";
        cout << "approximation to a delay equation.  The approximation is created by dividing\n";
        cout << "each delay into N smaller delays, and then approximating each fractional\n";
        cout << "delay with an Taylor series truncated at order p.\n";
        cout << endl;
        cout << "The new vector field file is written to 'standard output', so to use this\n";
        cout << "command, one usually redirects the output to a file. For example,\n";
        cout << "    vfgen delay2ode mydelay.vf > mydelay_2ode.vf\n";
        cout << endl;
        cout << "Options:\n";
        cout << "p=1|2|3\n";
        cout << "    The order of the Taylor approximation of the fractional delay step.\n";
        cout << "    The default is p=1.\n";
        cout << "N=integer\n";
        cout << "    The number of fractional delay steps into which each delay is split.\n";
        cout << "    The default is N=10.\n";
    }
    else if (strcmp(command,"dstool") == 0) {
        cout << "use: vfgen dstool vector_field_file.vf\n";
        cout << endl;
        cout << "This command creates a C definition file for the vector field to be used\n";
        cout << "with DSTOOL. The name of the file is [name]_def.c.\n";
        cout << endl;
        cout << "Options: none.\n";
    }
    else if (strcmp(command,"evf") == 0) {
        cout << "use: vfgen evf vector_field_file.vf\n";
        cout << endl;
        cout << "This command generates a new vector field file, in which the original vector\n";
        cout << "field has been extended with its variational equations.\n";
        cout << endl;
        cout << "The new vector field file is written to 'standard output', so to use this\n";
        cout << "command, one usually redirects the output to a file. For example,\n";
        cout << "    vfgen evf system.vf > system_evf.vf\n";
        cout << endl;
        cout << "Options:\n";
        cout << "par=p\n";
        cout << "    p must be the Name of a Parameter in the vector field file.\n";
        cout << "    If this option is given, the extension of the vector field is the variation\n";
        cout << "    with respect to the parameter.  That is, the vector field F(x) is extended\n";
        cout << "    with\n";
        cout << "        v' = (DF(x)/dx)v + DF(x)/dp\n";
        cout << "    This extended system allows for computing the sensitivity of a solution\n";
        cout << "    with respect to a parameter.\n";
        cout << "    If this option is not given, the vector field is simply extended with\n";
        cout << "        v' = (DF(x)/dx)v\n";
    }
    else if (strcmp(command,"gsl") == 0) {
        cout << "use: vfgen gsl vector_field_file.vf\n";
        cout << "     vfgen gsl:option=value,...,option=value vector_field_file.vf\n";
        cout << endl;
        cout << "This command generates a C file containing functions to be used with the\n";
        cout << "GNU Scientific Library (GSL) ODE suite.  The files that it creates are:\n";
        cout << "    [name]_gvf.c\n";
        cout << "    [name]_gvf.h\n";
        cout << "where [name] is the Name attribute of the VectorField entity in the vector\n";
        cout << "field file.  The file [name]_gvf.c will define the C functions [name]_vf(),\n";
        cout << "[name]_jac(), and [name]_jacp(). These compute the vector field, its Jacobian\n";
        cout << "matrix, and the derivatives of the vector field with respect to the parameters.\n";
        cout << endl;
        cout << "Options: (default is listed first)\n";
        cout << "func=no|yes\n";
        cout << "    If the option func=yes is given, the user functions defined in the\n";
        cout << "    vector field file are implemented as C functions in [name]_gvf.c.\n";
        cout << "    The name of the function will be [name]_[funcname](), where [name]\n";
        cout << "    is the Name attribute of the VectorField entity, and [funcname] is\n";
        cout << "    the Name attribute of the Function entity.\n";
        cout << "demo=no|yes\n";
        cout << "    If the option demo=yes is given, two more files are created:\n";
        cout << "        [name]_solver.c\n";
        cout << "        Makefile-[name]\n";
        cout << "    where [name] is the Name attribute of the VectorField entity in the\n";
        cout << "    vector field file.  This program provides a simple command-line solver.\n";
        cout << "    It takes arguments in the form name=value, where name can be a state\n";
        cout << "    variable (to give an initial condition), a parameter (to give the value\n";
        cout << "    of a parameter), or one of abserr, relerr or stoptime, to control the ODE\n";
        cout << "    solver.  The output of the program consists of columns of numbers; the\n";
        cout << "    first column is the time, and the rest are the state variables.\n";
    }
    else if (strcmp(command,"help") == 0) {
        cout << "use: vfgen help <command>\n";
        cout << endl;
        cout << "The help command will print a short description of the given <command>.\n";
        cout << "For example,\n";
        cout << "    vfgen help cvode\n";
        cout << "prints a description of the cvode command.\n";
    }
    else if (strcmp(command,"javascript") == 0) {
        cout << "use: vfgen javascript vector_field_file.vf\n";
        cout << "     vfgen javascript:option=value,... vector_field_file.vf\n";
        cout << endl;
        cout << "This command creates Javascript functions for computing a Taylor polynomial\n";
        cout << "approximation of the solution to the differential equations.\n";
        cout << "The file created is [name].js, where [name] is the Name attribute of the\n";
        cout << "VectorField entity in the vector field file.  The Javascript file contains\n";
        cout << "the functions [name]_vf and [name]_evaltaylor[r], where [r] is the order of\n";
        cout << "the Taylor polynomial. The first computes the vector field, and the second\n";
        cout << "computes the Taylor approximation. Also defined in the Javascript file are the\n";
        cout << "functions [name]_diff[k], for k=1, 2, ..., r-1. These functions compute the\n";
        cout << "k-linear differentials DkF(X)[V1,...,Vk]. These could be useful in programs\n";
        cout << "that analyze bifurcations or that compute invariant manifolds.\n";
        cout << endl;
        cout << "Options: (default is listed first)\n";
        cout << "demo=no|yes\n";
        cout << "    If the option demo=yes is given, two additional files are created:\n";
        cout << "    [name].html and [name]_solverdemo.js.  The code in [name]_solverdemo.js\n";
        cout << "    uses the functions defined in [name].js to implement a Taylor method solver\n";
        cout << "    with an adaptive step size.  The HTML file [name].html creates an interface\n";
        cout << "    to this solver.  The interfaces provides a form for defining the parameters\n";
        cout << "    and initial conditions, and allows some additional solver parameters to be\n";
        cout << "    changed.  It also has a graphical display of a phase-space projection of the\n";
        cout << "    solution.  The demo should work in any browser that has Javascript enabled\n";
        cout << "    and that supports the <canvas> tag in HTML.\n";
        cout << "order=[integer]\n";
        cout << "    The order of the Taylor polynomial. The default is order=5.\n";
    }
    else if (strcmp(command,"latex") == 0) {
        cout << "use: vfgen latex vector_field_file.vf\n";
        cout << endl;
        cout << "This command generate a LaTeX fragment for the vector field, in a file called\n";
        cout << "[name].tex, where [name] is the Name attribute of the VectorField entity in\n";
        cout << "the vector field file. The vector field will be defined inside nested\n";
        cout << "'equation' and 'split' environments. (The 'split' environment is provided\n";
        cout << "by the AMSMATH LaTeX package.)\n";
        cout << endl;
        cout << "Options: none.\n";
    }
    else if (strcmp(command,"lsoda") == 0) {
        cout << "use: vfgen lsoda vector_field_file.vf\n";
        cout << "     vfgen lsoda:option=value,option=value,... vector_field_file.vf\n";
        cout << endl;
        cout << "This command creates the Fortran file [name]_rhs.f, which defines the vector\n";
        cout << "field and its Jacobian in subroutines to be used with the Fortran 77 ODE solvers\n";
        cout << "LSODA, LSODAR or LSODE from the ODEPACK suite.  The subroutines defined in the\n";
        cout << "file are\n";
        cout << "    [name]_rhs    The vector field\n";
        cout << "    [name]_jac    The Jacobian of the vector field\n";
        cout << endl;
        cout << "Options: (default is listed first)\n";
        cout << "func=no|yes\n";
        cout << "    If the option func=yes is given, VFGEN creates a subroutine that computes\n";
        cout << "    the values of any Functions that were defined in the vector field file. The\n";
        cout << "    subroutine is designed to be used with the rootfinding capability of the\n";
        cout << "    LSODAR solver.\n";
        cout << "demo=no|yes\n";
        cout << "    If the option demo=yes is given, the file [name]_demo.f is also created.\n";
        cout << "    This program provides a sample driver for the LSODA subroutine. The initial\n";
        cout << "    conditions and parameter values are the default values defined in the vector\n";
        cout << "    field file.\n";
        cout << "parstyle=after|common\n";
        cout << "    There are two methods for passing parameters to the vector field and\n";
        cout << "    Jacobian functions:  include them in the array that holds the state\n";
        cout << "    variables, or put them in a named common block. The parstyle option\n";
        cout << "    indicates which method to use. If the option parstyle=common is given, the\n";
        cout << "    parameters are put in a common block with the name [name]_parameters, where\n";
        cout << "    [name] is the name of the vector field. If the option parstyle=after is\n";
        cout << "    given, the parameters are included in the state variable vector, beginning\n";
        cout << "    just after the last state variable.\n";
    }
    else if (strcmp(command,"matcont") == 0) {
        cout << "use: vfgen matcont vector_field_file.vf\n";
        cout << endl;
        cout << "This command generates a MATLAB file to be used with the MATCONT and\n";
        cout << "CL_MATCONT programs.  The name of the file is [name].m, where [name]\n";
        cout << "is the Name attribute of the VectorField entity in the vector field file.\n";
        cout << endl;
        cout << "Options: none.\n";
    }
    else if (strcmp(command,"matlab") == 0) {
        cout << "use: vfgen matlab vector_field_file.vf\n";
        cout << "     vfgen matlab:option=value,...,option=value vector_field_file.vf\n";
        cout << endl;
        cout << "This command generates MATLAB files to be used with the MATLAB ODE solvers.\n";
        cout << "The files created are:\n";
        cout << "    [name]_vf.m     The vector field\n";
        cout << "    [name]_jac.m    The Jacobian of the vector field\n";
        cout << "    [name]_jacp.m   Derivatives of the vector field with resperct to the\n";
        cout << "                    parameters.\n";
        cout << "    [name]_hess.m   The Hessians of the vector field\n";
        cout << "    [name]_der3.m   The third order derivatives of the vector field\n";
        cout << "where [name] is the Name attribute of the VectorField entity in the\n";
        cout << "vector field file. (Only the first two are actually used by MATLAB's\n";
        cout << "ODE solvers.)\n";
        cout << endl;
        cout << "Options: (default is listed first)\n";
        cout << "parstyle=vector|list\n";
        cout << "    This option controls how parameters are passed to the functions. By\n";
        cout << "    default (parstyle=vector), the parameters are in a single array.\n";
        cout << "    If the option parstyle=list is given, the parameters are all listed\n";
        cout << "    separately as arguments.\n";
        cout << "demo=no|yes\n";
        cout << "    If the option demo=yes is given, the file [name]_demo.m is created.\n";
        cout << "    This defines the MATLAB function [name]_demo, which creates a simple\n";
        cout << "    user interface with fields for the parameter values and initial\n";
        cout << "    conditions, and a button that will generate and plot a solution.\n";
        cout << "func=no|yes\n";
        cout << "    If the option func=yes is given, the file [name]_[funcname].m will be\n";
        cout << "    created for each user function.\n";
        cout << "evf=no|yes\n";
        cout << "    If the option evf=yes is given, the file [name]_evf.m is created. This\n";
        cout << "    file defines a function that computes the extended vector field consisting\n";
        cout << "    of the original vector field and the variational equations. This returns a\n";
        cout << "    2N dimensional column vector. (Note: The EVF command is probably a better\n";
        cout << "    alternative, because you can then use VFGEN to create the Jacobian of the\n";
        cout << "     extended vector field.)\n";
    }
    else if (strcmp(command,"octave") == 0) {
        cout << "use: vfgen octave vector_field_file.vf\n";
        cout << "     vfgen octave:option=value,...,option=value vector_field_file.vf\n";
        cout << endl;
        cout << "This command generates an OCTAVE file, [name].m, to be used with the OCTAVE\n";
        cout << "ODE solver, where [name] is the Name attribute of the VectorField defined in\n";
        cout << "the vector field file. The file contains functions for the vector field, the\n";
        cout << "Jacobian, and, if the option func=yes was given, the user-defined functions.\n";
        cout << endl;
        cout << "Options: (default is listed first)\n";
        cout << "demo=no|yes\n";
        cout << "    If the option demo=yes is given, the file [name]_demo.m is created.\n";
        cout << "    This OCTAVE script generates a plot of a solution to the differential\n";
        cout << "    equations.  The default parameter values and initial conditions are\n";
        cout << "    used.\n";
        cout << "func=no|yes\n";
        cout << "    If the option func=yes is given, the file [name].m will define a\n";
        cout << "    function for each user function given in the vector field file.\n";
    }
    else if (strcmp(command,"pddecont") == 0) {
        cout << "use: vfgen pddecont vector_field_file.vf\n";
        cout << endl;
        cout << "This command creates the file sys-[name].cpp to be used with the PDDE-CONT\n";
        cout << "software package.  PDDE-CONT does continuation and bifurcation computations\n";
        cout << "for delay differential equations.\n";
        cout << endl;
        cout << "Options: none.\n";
    }
	else if (strcmp(command,"pydstool") == 0) {
		cout << "use: vfgen pydstool vector_field_file.vf\n";
		cout << "     vfgen pydstool:demo=yes vector_field_file.vf\n";
		cout << endl;
		cout << "This command creates a Python file that can be used to define a differential\n";
	    cout << "equation for the software package PyDSTool.  The file [name].py defines the\n";
	    cout << "function args(), which creates an 'args' object for the vector field.\n";
		cout << endl;
		cout << "Options: (default is list first)\n";
		cout << "demo=no|yes\n";
		cout << "    If the option demo=yes is given, the file [name]_dst.py is also created.\n";
		cout << "    This file contains a script that uses the function defined in [name].py\n";
	    cout << "    to generate and plot a solution to the differential equation.\n";
	}
    else if (strcmp(command,"pygsl") == 0) {
        cout << "use: vfgen pygsl vector_field_file.vf\n";
        cout << "     vfgen pygsl:option=value,...,option=value vector_field_file.vf\n";
        cout << endl;
        cout << "This command creates Python files that can be used with the PyGSL Python\n";
        cout << "library. The files created are [name].py and, if the option demo=yes option\n";
        cout << "is given, [name]_demo.py, where [name] is the Name attribute of the VectorField\n";
        cout << "entity in the vector field file. The Python file [name].py will contain the\n";
        cout << "functions vectorfield(...), jacobian(...), and, if the func=yes option is\n";
        cout << "given, a function for each user-defined function.\n";
        cout << endl;
        cout << "Options: (default is listed first)\n";
        cout << "func=no|yes\n";
        cout << "    If the option func=yes is given, VFGEN also converts any user-defined\n";
        cout << "    functions in the vector field file into functions in the Python file.\n";
        cout << "    The names of the functions will be the same as those given in the vector\n";
        cout << "    field file.\n"; 
        cout << "demo=no|yes\n";
	    cout << "    If the option demo=yes is given, the file [name]_demo.py will contain a\n";
        cout << "    Python script for a command-line ODE solver for the vector field. The\n";
        cout << "    initial conditions and parameters can be specified on the command line.\n";
        cout << "    The program will print the solution data to the console.\n";
    }
    else if (strcmp(command,"r") == 0) {
        cout << "use: vfgen r vector_field_file.vf\n";
        cout << endl;
        cout << "This command creates the file [name].R, which defines the vector field\n";
        cout << "and its Jacobian in subroutines to be used with the R package deSolve.\n";
        cout << "XXX ignore the rest... \n";
        cout << "The subroutines defined in the file are\n";
        cout << "    [name]_rhs    The vector field\n";
        cout << "    [name]_jac    The Jacobian of the vector field\n";
        cout << "    [name]_sol    A subroutine that prints a point of the solution.\n";
        cout << endl;
        cout << "Options: (default is listed first)\n";
        cout << "demo=no|yes\n";
        cout << "    If the option demo=yes is given, the file [name]_dr5.f is also created.\n";
        cout << "    This program provides a driver that will generated a solution to equations\n";
        cout << "    by calling RADAU5.  The initial conditions and parameters of the solution\n";
        cout << "    are the default values given in the vector field file.\n";
    }
	else if (strcmp(command,"radau5") == 0) {
		cout << "use: vfgen radau5 vector_field_file.vf\n";
		cout << "     vfgen radau5:demo=yes vector_field_file.vf\n";
		cout << endl;
		cout << "This command creates the file [name]_rhs.f, which defines the vector field\n";
		cout << "and its Jacobian in subroutines to be used with the Fortran 77 ODE solver\n";
		cout << "RADAU5.  The subroutines defined in the file are\n";
		cout << "    [name]_rhs    The vector field\n";
		cout << "    [name]_jac    The Jacobian of the vector field\n";
		cout << "    [name]_sol    A subroutine that prints a point of the solution.\n";
		cout << endl;
	    cout << "Options: (default is listed first)\n";
		cout << "demo=no|yes\n";
		cout << "    If the option demo=yes is given, the file [name]_dr5.f is also created.\n";
		cout << "    This program provides a driver that will generated a solution to equations\n";
		cout << "    by calling RADAU5.  The initial conditions and parameters of the solution\n";
	    cout << "    are the default values given in the vector field file.\n";
	}
    else if (strcmp(command,"scilab") == 0) {
        cout << "use: vfgen scilab vector_field_file.vf\n";
        cout << "     vfgen scilab:option=value,...,option=value vector_field_file.vf\n";
        cout << endl;
        cout << "This command creates files to be used with Scilab. The main file created\n";
        cout << "is [name].sci, where [name] is the Name attribute of the VectorField\n";
        cout << "entity in the vector field file. The Scilab functions in this file are\n";
        cout << "    [name]_vf(t,x,p)\n";
        cout << "    [name]_jac(t,x,p)\n";
        cout << endl;
        cout << "Options: (default is listed first)\n";
        cout << "parstyle=vector|list\n";
        cout << "    The method by which parameters are passed to the Scilab functions\n";
        cout << "    generated by VFGEN is controlled by the parstyle option. When\n";
        cout << "    parstyle=list is given as\n";
        cout << "        $ vfgen scilab:parstyle=list vector_field_file.vf\n";
        cout << "    each parameter is listed explicitly as an argument of Scilab functions.\n";
        cout << "    That is, instead of [name]_vf(t,x,p), the vector field function will be\n";
        cout << "    [name]_vf(t,x,param1,param2,...,paramM).\n"; 
        cout << "func=no|yes\n";
        cout << "    If the option func=yes is given, a Scilab function will be created for\n";
        cout << "    each user-defined function.\n";
        cout << "demo=no|yes\n";
        cout << "    If the option demo=yes is given, a second file is created called\n";
        cout << "    [name]_demo.sce. This is a script that when run provides a simple GUI\n";
        cout << "    interface to the initial conditions, system parameters, and the ODE\n";
        cout << "    solver parameters. It will call the odeint function and plot the\n";
        cout << "    solution. Run the script in Scilab with the command\n";
        cout << "        -->exec [name]_demo.sce;\n"; 
    }
    else if (strcmp(command,"scipy") == 0) {
        cout << "use: vfgen scipy vector_field_file.vf\n";
        cout << "     vfgen scipy:option=value,...,option=value vector_field_file.vf\n";
        cout << endl;
        cout << "This command generates Python files to be used with the SciPy library.\n";
        cout << "The files created are [name].py and, if the option demo=yes option is\n";
        cout << "given, [name]_demo.py, where [name] is the Name attribute of the VectorField\n";
        cout << "entity in the vector field file. The Python file [name].py will contain the\n";
        cout << "functions vectorfield(...) and jacobian(...).\n";
        cout << endl;
        cout << "Options: (default is listed first)\n";
        cout << "func=no|yes\n";
        cout << "    If the option func=yes is given, VFGEN also converts any user-defined\n";
        cout << "    functions in the vector field file into functions in the Python file.\n";
        cout << "    The names of the functions will be the same as those given in the vector\n";
        cout << "    field file.\n"; 
        cout << "demo=no|yes\n";
	    cout << "    If the option demo=yes is given, the file [name]_demo.py will contain a\n";
        cout << "    Python script for a command-line ODE solver for the vector field. The\n";
        cout << "    initial conditions and parameters can be specified on the command line.\n";
        cout << "    The program will print the solution data to the console.\n";         
    }
    else if (strcmp(command,"taylor") == 0) {
        cout << "use: vfgen taylor vector_field_file.vf\n";
        cout << "     vfgen taylor:order=[integer] vector_field_file.vf\n";
        cout << endl;
        cout << "This command creates C functions for computing a Taylor polynomial\n";
        cout << "approximation of the solution to the differential equations.\n";
        cout << "The files created are [name]_taylor[r].c and [name]_taylor[r].h, where [name]\n";
        cout << "is the Name attribute of the VectorField entity in the vector field file, and\n";
        cout << "[r] is the order given with the order option.  The C file [name]_taylor[r].c\n";
        cout << "will contain the functions [name]_vf and [name]_taylor[r]. The first computes\n";
        cout << "the vector field, and the second computes the Taylor approximation. Also defined\n";
        cout << "in the C file are the functions [name]_diff[k], for k=1, 2, ..., r-1. These\n";
        cout << "functions compute the k-linear differentials DkF(X)[V1,...,Vk]. These could be\n";
        cout << "useful in programs that analyze bifurcations or that compute invariant\n";
        cout << "manifolds.\n";
        cout << endl;
        cout << "Options:\n";
        cout << "order=[integer]\n";
        cout << "    The order of the Taylor polynomial. The default is order=5.\n";
    }
    else if (strcmp(command,"xpp") == 0) {
        cout << "use: vfgen xpp vector_field_file.vf\n";
        cout << "     vfgen xpp:option=value,...,option=value vector_field_file.vf\n";
        cout << endl;
        cout << "This command creats a file to be used with XPP (aka XPP-AUT).\n";
        cout << "The name of the file created by this command is [name].ode, where\n";
        cout << "[name] is the Name attribute of the VectorField entity in the vector\n";
        cout << "field file.\n";
        cout << endl;
        cout << "Options:\n";
        cout << "extra=[text]\n";
        cout << "    The option extra=[text] allows the user to customize the ODE file with\n";
        cout << "    additional lines. The characters in text are added to a section of the\n";
        cout << "    ODE file. Any semi-colons in text are converted to newlines. For example,\n";
        cout << "    extra=\"@ total=100;@ maxstor=10000\" will add the lines\n";
        cout << "        @ total=100\n";
        cout << "        @ maxstor=10000\n";
        cout << "    to the ODE file. (A final semi-colon in text is optional.)\n";
        cout << "    Note: Currently, text must not contain any commas!\n";
    }
    else {
        cout << "Sorry, help for \"" << command << "\" is not available yet!" << endl;
    }
    return 0;
}

///////////////////////////////////////////////////////////
//  main
///////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    VectorField vf;
    string s, extrastr;

    if (argc != 3) {
        printuse();
        exit(-1);
    }

    //
    //  Allowed command options:
    //
    command_options["auto"].push_back("lang");
    command_options["cvode"].push_back("demo");
    command_options["cvode"].push_back("func");
    command_options["cvode"].push_back("version");
    command_options["dde23"].push_back("demo");
    command_options["dde23"].push_back("parstyle");
    command_options["ddebiftool"].push_back("path");
    command_options["dde_solver"].push_back("demo");
    command_options["delay2ode"].push_back("N");
    command_options["delay2ode"].push_back("p");
    command_options["evf"].push_back("par");
    command_options["gsl"].push_back("demo");
    command_options["gsl"].push_back("func");
    command_options["javascript"].push_back("order");
    command_options["javascript"].push_back("demo");
    command_options["lsoda"].push_back("demo");
    command_options["lsoda"].push_back("func");
    command_options["lsoda"].push_back("parstyle");
    command_options["matlab"].push_back("demo");
    command_options["matlab"].push_back("evf");
    command_options["matlab"].push_back("func");
    command_options["matlab"].push_back("parstyle");
    command_options["octave"].push_back("demo");
    command_options["octave"].push_back("func");
    command_options["pydstool"].push_back("demo");
    command_options["pygsl"].push_back("demo");
    command_options["pygsl"].push_back("func");
    command_options["r"].push_back("demo");
    command_options["r"].push_back("func");    
    command_options["radau5"].push_back("demo");
    command_options["scilab"].push_back("demo");
    command_options["scilab"].push_back("func");
    command_options["scilab"].push_back("parstyle");
    command_options["scipy"].push_back("demo");
    command_options["scipy"].push_back("func");
    command_options["taylor"].push_back("order");
    command_options["xpp"].push_back("extra");
    
    string commandstr(argv[1]);
    //
    // Check for the help command.
    //
    if (commandstr == "help") {
        exit(help(argv[2]));
    }

    //
    //  Check for any options appended to the command.
    //
    map<string,string> options;

    string::size_type loccolon = commandstr.find(":",0);
    if (loccolon != string::npos) {
        // extrastr holds all the options given after the :
        extrastr = commandstr.substr(loccolon+1,commandstr.length()-loccolon-1);
        commandstr.erase(loccolon,commandstr.length()-loccolon);
        // cerr << "Options \"" << extrastr << "\"" << endl;
        string::size_type pos = 0;
        do {
            string::size_type locsep = extrastr.find(",",pos);
            string current_option, option_name, option_value;
            if (locsep == string::npos) {
                locsep = extrastr.length();
            }
            current_option = extrastr.substr(pos,locsep-pos);
            // cerr << "current_option = \"" << current_option << "\"\n";
            string::size_type loceq = current_option.find("=",0);
            if (loceq == string::npos) {
                // No "=" given in the option.
                option_name = current_option;
                option_value = "";
            }
            else {
                if (loceq == 0) {
                    // The option was "=something"; not valid
                    printuse();
                    exit(-1);
                }
                option_name = current_option.substr(0,loceq);
                option_value = current_option.substr(loceq+1,current_option.length()-loceq);
            }
            options[option_name] = option_value;
            pos = locsep+1;
        } while (pos < extrastr.length());
    }
        
        
    int command = checkcommand(commandstr.c_str());
    if (command == -1) {
        cout << "vfgen: unknown command: " << commandstr << endl;
        printuse();
        exit(-1);
    }

    //
    //  Check that any options given are known options.
    //  (This doesn't check that the value of the option is valid; it just
    //  checks that the name of the option is one of the  allowed options for
    //  the given command.)
    //
    bool bad_opt = false;        
    map<string,string>::const_iterator opt;
    for (opt = options.begin(); opt != options.end(); ++opt) {
        string optstr = opt->first;
        // cerr << "Option: " << optstr;
        // if (opt->second != "")
        //     cerr << "=" << opt->second;
        // cerr << endl;
        bool validopt = false;
        vector<string>::const_iterator w;
        for (w = command_options[commandstr].begin(); w != command_options[commandstr].end(); ++w) {
            if (optstr == *w) {
                validopt = true;
                break;
            }
        }
        if (!validopt) {
            cerr << "Errror: \"" << optstr << "\" is not a valid option for the " << commandstr << " command.\n";
            bad_opt = true;
        }
    }
    if (bad_opt) {
        printuse();
        exit(-1);
    }

    //
    //  Read the vector field file.  This just puts the strings into the
    //  appropriate fields.  It doesn't do any symbolic processing.
    //
    vf.ReadXML(argv[2]);

    //
    //  Process the strings to create the GiNaC symbolic expressions in the object.
    //
    int pserr = vf.ProcessSymbols();
    if (pserr == -1) {
        exit(-1);
    }

    //
    // Call the appropriate output function based on the first
    // command line argument.
    //
    if (commandstr == "check") {
        vf.Print();
    }
    else if (commandstr == "xpp") {
        vf.PrintXPP(options);
    }
    else if (commandstr == "xml") {
        vf.PrintXML("xml");
    }
    else if (commandstr == "delay2ode") {
        if (vf.IsDelay == true) {
            vf.PrintDelay2ODE(options);
        }
        else {
            cerr << "This system is not a delay equation.\n";
        }
    }
    else if (commandstr == "dde23") {
        if (vf.IsDelay == true) {
            if (vf.HasNonconstantDelay) {
                cerr << "This system has nonconstant delays.  DDE23 can only be used with constant delays.\n";
            }
            else {
                vf.PrintDDE23(options);
            }
        }
        else {
            cerr << "This system is not a delay equation.\n";
        }
    }
    else if (commandstr == "ddebiftool") {
        if (vf.IsDelay == true) {
            if (vf.HasNonconstantDelay) {
                cerr << "This system has nonconstant delays. VFGEN does not (yet) generate DDE-BIFTOOL files for such systems.\n";
            }
            else {
                vf.PrintDDEBIFTOOL(options);
            }
        }
        else {
            cerr << "This system is not a delay equation.\n";
        }
    }
    else if (commandstr == "pddecont") {
        if (vf.HasNonconstantDelay) {
            cerr << "This system has nonconstant delays. PDDE-CONT is for systems with constant delays.\n";
        }
        else {
            vf.PrintPDDECONT(options);
        }
    }
    else if (commandstr == "dde_solver") {
        if (vf.IsDelay == true) {
            vf.PrintDDE_SOLVER(options);
        }
        else {
            cerr << "This system is not a delay equation.\n";
        }
    }
    else if (commandstr == "matlab") {
        if (vf.IsDelay == false) {
            vf.PrintMATLAB(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "octave") {
        if (vf.IsDelay == false) {
            vf.PrintOctave(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "scilab") {
        if (vf.IsDelay == false) {
            vf.PrintScilab(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "matcont") {
        if (vf.IsDelay == false) {
            vf.PrintMATCONT(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "dstool") {
        if (vf.IsDelay == false) {
            vf.PrintDSTool();
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "auto") {
        if (vf.IsDelay == false) {
            vf.PrintAUTO(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "gsl") {
        if (vf.IsDelay == false) {
            vf.PrintGSL(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "cvode") {
        if (vf.IsDelay == false) {
            vf.PrintCVODE(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "adolc") {
        if (vf.IsDelay == false) {
            vf.PrintADOLC(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "r") {
        if (vf.IsDelay == false) {
            vf.PrintR(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "radau5") {
        if (vf.IsDelay == false) {
            vf.PrintRadau5(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "lsoda") {
        if (vf.IsDelay == false) {
            vf.PrintLSODA(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "scipy") {
        if (vf.IsDelay == false) {
            vf.PrintSciPy(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "pydstool") {
        if (vf.IsDelay == false) {
            vf.PrintPyDSTool(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "taylor") {
        if (vf.IsDelay == false) {
            vf.PrintTaylor(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "latex") {
        vf.PrintLatex(options);
    }
    else if (commandstr == "pygsl") {
        if (vf.IsDelay == false) {
            vf.PrintPyGSL(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "evf") {
        if (vf.IsDelay == false) {
            vf.PrintEVF(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "javascript") {
        if (vf.IsDelay == false) {
            vf.PrintJavascript(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else {
        // This should not happen!!!
        cerr << "vfgen: Unknown command: " << argv[1] << endl;
        printuse();
        exit(-1);
    }
    return(0);
}  // end main()
