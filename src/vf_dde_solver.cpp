
//
//  vf_dde_solver.cpp
//
//  This file defines the VectorField::PrintDDE_SOLVER method.
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
#include <cassert>
#include <ginac/ginac.h> 

#include "vf.h"
#include "vf_utils.h"

using namespace std;
using namespace GiNaC;


//
// PrintDDE_SOLVER --
//

void VectorField::PrintDDE_SOLVER(map<string,string> options)
{
    int nc, np, nv, na;
    int version;

    nc = conname_list.nops();
    nv = varname_list.nops();
    np = parname_list.nops();
    na = exprname_list.nops();

    if ((options["version"] == "") || (options["version"] == "2")) {
        version = 2;
    }
    else if (options["version"] == "1") {
        version = 1;
    }
    else {
        cerr << "vfgen DDE_SOLVER command: unknown version specified: " << options["version"] << endl;
        cerr << "Versions of DDE_SOLVER supported by VFGEN are 1 and 2.  The default is 2." << endl;
        exit(-1);
    }
    //
    //  Create the vector field function.
    //
    string vf_filename = Name()+".f90";
    ofstream fout;
    fout.open(vf_filename.c_str());
    fout << left;
    fout << csrc;
    // Also override the csrc style for powers.
    // IMPORTANT: This means we can NOT subsequently print C/C++ code!
    set_print_func<power,print_csrc>(print_power_as_fortran);

    fout << "!" << endl;
    fout << "! " << vf_filename << endl;
    fout << "!" << endl;
    fout << "! Vector field function for: " << Name() << endl;
    fout << "! These functions are to be used with DDE_SOLVER_M.\n";
    fout << "!" << endl;
    PrintVFGENComment(fout,"! ");
    fout << "!" << endl;
    fout << "!" << endl;
    // Include a list of the lags in the comments.
    fout << "! The lags are: {";
    for (unsigned k = 0; k < Delays.size(); ++k) {
        fout << Delays[k];
        if (k < Delays.size()-1) {
            fout << ", "; 
        }
    }
    fout << "}" << endl;
    fout << endl;
    fout << "module define_" << Name() << "_ddes\n";
    fout << endl;
    fout << "    implicit none\n";
    fout << "    integer, parameter :: NEQN=" << nv << ", NLAGS=" << Delays.size() << endl;
    //
    // Constants. (Declared as global variables.)
    //
    // COMMENTED OUT--I'm not sure why I wanted these to be global...
    // if (conname_list.nops() > 0)
    //     {
    //     fout << "    ! Constants\n";
    //     Declare(&fout, "    ","DOUBLE PRECISION", conname_list,"");
    //     }
    //
    // Parameters. (Declared as global variables.)
    //
    if (parname_list.nops() > 0) {
        fout << "    ! Parameters\n";
        Declare(fout, "    ","double precision", parname_list,"");
    }
    fout << endl;
    fout << "contains\n";
    fout << endl;
    // DDE function definition starts here.
    fout << "    subroutine " << Name() << "_ddes(" << IndependentVariable << ", x_, Zlags_, vf_)\n";
    fout << "    ! Arguments\n";
    fout << "    double precision :: " << IndependentVariable << endl;
    fout << "    double precision, dimension(";
    if (version == 2) {
        fout << ":";
    }
    else {
        fout << "NEQN";
    }
    fout << ") :: x_, vf_\n";
    fout << "    double precision, dimension(";
    if (version == 2) {
        fout << ":,:";
    }
    else {
        fout << "NEQN, NLAGS";
    }
    fout << ") :: Zlags_\n";
    fout << "    intent(in) :: " << IndependentVariable << ", x_, Zlags_\n";
    fout << "    intent(out) :: vf_\n";
    fout << "    ! Local variables\n";
    if (nc > 0) {
        Declare(fout, "    ", "double precision", conname_list, "");
    }
    if (na > 0) {
        Declare(fout, "    ", "double precision", exprname_list, "");
    }
    Declare(fout, "    ", "double precision", varname_list, "");
    if (HasPi) {
        fout << "    double precision Pi\n";
    }
    if (HasPi) {
        fout << "    Pi = 3.1415926535897932385D0\n";
    }
    //
    // Constants...
    //
    if (nc > 0) {
        fout << "    ! Constants\n";
        AssignNameValueLists(fout, "    ", conname_list, "=", convalue_list, "");
    }
    fout << "    ! State variables\n";
    GetFromVector(fout, "    ", varname_list, "=", "x_", "()", 1, "");

    //
    // Expressions...
    //
    if (na > 0) {
        fout << "    ! Expressions\n";
    }
    for (int i = 0; i < na; ++i) {
        ex f = exprformula_list[i];
        if (f.has(delay(wild(1),wild(2)))) {
            ConvertDelaysToZlags(f, 1, 1);
        }
        ostringstream os;
        os << left << csrc;
        os << "    " << exprname_list[i] << " = " << f << endl;
        F90Write(fout, os.str());
    }
    //
    // StateVariables...
    //
    fout << "    ! The vector field\n";
    for (int i = 0; i < nv; ++i) {
        ex f = varvecfield_list[i];
        if (f.has(delay(wild(1),wild(2)))) {
            ConvertDelaysToZlags(f, 1, 1);
        }
        ostringstream os;
        os << left << csrc;
        os << "    vf_(" << (i+1) << ")" << " = " << f << endl;
        F90Write(fout, os.str());
    }
    fout << endl;
    fout << "    return\n";
    fout << "    end subroutine " << Name() << "_ddes\n";
    fout << endl;
    // DDE history function starts here
    fout << "    subroutine " << Name() << "_history(" << IndependentVariable << ",x_)\n";
    fout << "    double precision :: " << IndependentVariable << endl;
    fout << "    double precision, dimension(";
    if (version == 2) {
        fout << ":";
    }
    else {
        fout << "NEQN";
    }
    fout << ") :: x_\n";
    fout << "    intent(in) :: " << IndependentVariable << endl;
    fout << "    intent(out) :: x_\n";
    fout << endl;
    for (int i = 0; i < nv; ++i) {
        ostringstream os;
        os << left << csrc;
        os << "    x_(" << i+1 << ") = " << vardefhist_list[i] << endl;
        F90Write(fout, os.str());        
    }
    fout << endl;
    fout << "    return\n";
    fout << "    end subroutine " << Name() << "_history\n";
    fout << endl;
    if (HasNonconstantDelay) {
        // If there is a nonconstant delay, create the BETA subroutine
        fout << "    subroutine " << Name() << "_beta(" << IndependentVariable << ",x_,bval_)\n";
        fout << "    ! Arguments\n";
        fout << "    double precision :: " << IndependentVariable << endl;
        fout << "    double precision, dimension(";
        if (version == 2) {
            fout << ":";
        }
        else {
            fout << "NEQN";
        }
        fout << ") :: x_\n";
        fout << "    double precision, dimension(";
        if (version == 2) {
            fout << ":";
        }
        else {
            fout << "NLAGS";
        }
        fout << ") :: bval_\n";
        fout << "    intent(in) :: " << IndependentVariable << ", x_\n";
        fout << "    intent(out) :: bval_\n";
        fout << "    ! Local variables\n";
        if (nc > 0) {
            Declare(fout, "    ", "double precision", conname_list, "");
        }
        if (na > 0) {
            Declare(fout, "    ", "double precision", exprname_list, "");
        }
        Declare(fout, "    ", "double precision", varname_list, "");
        // Constants...
        if (nc > 0) {
            AssignNameValueLists(fout, "    ", conname_list, "=", convalue_list, "");
        }
        fout << endl;
        // State Variables...
        fout << "    ! State variables\n";
        GetFromVector(fout, "    ", varname_list, "=", "x_", "()", 1, "");
        // Expressions...
        if (na > 0) {
            fout << "    ! Expressions\n";
        }
        for (int i = 0; i < na; ++i) {
            ex f = exprformula_list[i];
            if (!f.has(delay(wild(1),wild(2)))) {
                ostringstream os;
                os << left << csrc;
                os << "    " << exprname_list[i] << " = " << f << endl;
                F90Write(fout, os.str());
            }
        }
        for (unsigned i = 0; i < Delays.size(); ++i) {
            ostringstream os;
            os << left << csrc;
            os << "    bval_(" << i+1 << ") = " << IndVar - Delays[i] << endl;
            F90Write(fout, os.str());
        }
        fout << "    end subroutine " << Name() << "_beta\n";
        fout << endl;
    }
    fout << "end module define_" << Name() << "_ddes\n";
    fout.close();

    if (options["demo"] == "yes") {
        //
        // Create the demo function.
        //
        string filename = Name() + "_demo.f90";
        ofstream fout;
        fout.open(filename.c_str());
        fout << left;
        fout << csrc;

        fout << "!" << endl;
        fout << "! " << filename << endl;
        fout << "!" << endl;
        fout << "! Fortran 90 program that will use DDE_SOLVER_M to solve the DDEs defined\n";
        fout << "! in the vector field: " << Name() << endl;
        fout << "!" << endl;
        PrintVFGENComment(fout,"! ");
        fout << "!" << endl;
        fout << "!" << endl;
        fout << endl;
        fout << "program " << Name() << "_demo\n";
        fout << endl;
        fout << "use define_" << Name() << "_ddes, only: &\n";
        fout << "            NEQN, NLAGS, ";
        if (parname_list.nops() > 0) {
            PrintList(fout, parname_list);
            fout << ", ";
        }
        fout << "&\n";
        fout << "            " << Name() << "_ddes, &\n";
        fout << "            " << Name() << "_history";
        if (HasNonconstantDelay) {
            fout << ", &\n";
            fout << "            " << Name() << "_beta\n";
        }
        else {
            fout << "\n";
        }
        fout << "use dde_solver_m, only: dde_sol, dde_opts, dde_set, dde_solver, release_arrays\n";
        fout << endl;
        fout << "implicit none\n";
        fout << endl;
        fout << "integer, dimension(2) :: NVAR = (/NEQN, NLAGS/)\n";
        fout << endl;
        fout << "type(dde_sol) :: sol\n";
        fout << "type(dde_opts) :: opts\n";
        fout << endl;
        if (!HasNonconstantDelay) {
            fout << "double precision, dimension(" << Delays.size() << ") :: lags\n";
        }
        if (np > 0) {
            fout << "double precision, dimension(" << np << ") :: p_\n";
        }
        fout << "double precision, dimension(2) :: tspan\n";
        fout << endl;
        fout << "integer :: i, j\n";
        fout << "character(7+6*NEQN) :: f\n";
        fout << "double precision :: relerr, abserr, stoptime\n";
        if (HasPi) {
            fout << "double precision :: Pi\n";
        }
        if (nc > 0) {
            Declare(fout, "", "double precision ::", conname_list, "");
        }
        fout << endl;

        if (HasPi) {
            fout << "Pi = 3.1415926535897932385D0\n";
        }

        if (nc > 0) {
            fout << "! Constants\n";
            AssignNameValueLists(fout, "", conname_list, "=", convalue_list, "");
            fout << endl;
        }

        if (np > 0) {
            fout << "! Set the parameters of the DDE\n";
            AssignNameValueLists(fout, "", parname_list, "=", pardefval_list, "");
        }

        fout << "! Set the solver parameters: relative error, abs. error, stop time\n";
        fout << "relerr = 1D-7\n";
        fout << "abserr = 1D-9\n";
        fout << "stoptime = 10.0\n";
        for (unsigned k = 0; k < parname_list.nops(); ++k) { 
            fout << "p_(" << k+1 << ") = " << parname_list[k] << "\n";
        }
        if (!HasNonconstantDelay) {
            // If there are only constant delays, put them in the array LAGS
            fout << "! Initialize the array of lags\n";
            for (unsigned k = 0; k < Delays.size(); ++k) {
                fout << "lags(" << k+1 << ") = " << Delays[k] << endl;
            }
        }
        fout << endl;
        fout << "tspan(1) = 0.0\n";
        fout << "tspan(2) = stoptime\n";
        fout << "opts = dde_set(re=relerr, ae=abserr)\n";
        fout << endl;
        string lags_arg;
        if (HasNonconstantDelay) {
            lags_arg = Name() + "_beta";
        }
        else {
            lags_arg = "lags";
        } 
        fout << "sol = dde_solver(NVAR, " << Name() << "_ddes, " << lags_arg <<
             ", " << Name() << "_history, tspan, options=opts)\n";
        fout << endl;
        fout << "f = \"(E17.8\"//REPEAT(\",E17.8\",NEQN)//\")\"\n";
        fout << "do I = 1, sol%npts\n";
        fout << "    write(*, fmt=f) sol%t(i), (sol%y(i,j), j=1,NEQN)\n";
        fout << "end do\n";
        fout << endl;
        fout << "call release_arrays(sol, opts)\n";
        fout << endl;
        fout << "end program " << Name() << "_demo\n";
        fout.close();
    }
}
