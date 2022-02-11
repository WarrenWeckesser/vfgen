
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
    fout << "MODULE DEFINE_" << Name() << "_DDEs\n";
    fout << endl;
    fout << "    IMPLICIT NONE\n";
    fout << "    INTEGER, PARAMETER :: NEQN=" << nv << ", NLAGS=" << Delays.size() << endl;
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
        Declare(fout, "    ","DOUBLE PRECISION", parname_list,"");
    }
    fout << endl;
    fout << "CONTAINS\n";
    fout << endl;
    // DDE function definition starts here.
    fout << "    SUBROUTINE " << Name() << "_ddes(" << IndependentVariable << ",x_,Zlags_,vf_)\n";
    fout << "    ! Arguments\n";
    fout << "    DOUBLE PRECISION :: " << IndependentVariable << endl;
    fout << "    DOUBLE PRECISION, DIMENSION(";
    if (version == 2) {
        fout << ":";
    }
    else {
        fout << "NEQN";
    }
    fout << ") :: x_, vf_\n";
    fout << "    DOUBLE PRECISION, DIMENSION(";
    if (version == 2) {
        fout << ":,:";
    }
    else {
        fout << "NEQN,NLAGS";
    }
    fout << ") :: Zlags_\n";
    fout << "    INTENT(IN) :: " << IndependentVariable << ", x_, Zlags_\n";
    fout << "    INTENT(OUT) :: vf_\n";
    fout << "    ! Local variables\n";
    if (nc > 0) {
        Declare(fout,"    ","DOUBLE PRECISION",conname_list,"");
    }
    if (na > 0) {
        Declare(fout,"    ","DOUBLE PRECISION",exprname_list,"");
    }
    Declare(fout,"    ","DOUBLE PRECISION",varname_list,"");
    if (HasPi) {
        fout << "    DOUBLE PRECISION Pi\n";
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
    fout << "    RETURN\n";
    fout << "    END SUBROUTINE " << Name() << "_ddes\n";
    fout << endl;
    // DDE history function starts here
    fout << "    SUBROUTINE " << Name() << "_history(" << IndependentVariable << ",x_)\n";
    fout << "    DOUBLE PRECISION :: " << IndependentVariable << endl;
    fout << "    DOUBLE PRECISION, DIMENSION(";
    if (version == 2) {
        fout << ":";
    }
    else {
        fout << "NEQN";
    }
    fout << ") :: x_\n";
    fout << "    INTENT(IN) :: " << IndependentVariable << endl;
    fout << "    INTENT(OUT) :: x_\n";
    fout << endl;
    for (int i = 0; i < nv; ++i) {
        ostringstream os;
        os << left << csrc;
        os << "    x_(" << i+1 << ") = " << vardefhist_list[i] << endl;
        F90Write(fout, os.str());        
    }
    fout << endl;
    fout << "    RETURN\n";
    fout << "    END SUBROUTINE " << Name() << "_history\n";
    fout << endl;
    if (HasNonconstantDelay) {
        // If there is a nonconstant delay, create the BETA subroutine
        fout << "    SUBROUTINE " << Name() << "_beta(" << IndependentVariable << ",x_,bval_)\n";
        fout << "    ! Arguments\n";
        fout << "    DOUBLE PRECISION :: " << IndependentVariable << endl;
        fout << "    DOUBLE PRECISION, DIMENSION(";
        if (version == 2) {
            fout << ":";
        }
        else {
            fout << "NEQN";
        }
        fout << ") :: x_\n";
        fout << "    DOUBLE PRECISION, DIMENSION(";
        if (version == 2) {
            fout << ":";
        }
        else {
            fout << "NLAGS";
        }
        fout << ") :: bval_\n";
        fout << "    INTENT(IN) :: " << IndependentVariable << ", x_\n";
        fout << "    INTENT(OUT) :: bval_\n";
        fout << "    ! Local variables\n";
        if (nc > 0) {
            Declare(fout,"    ","DOUBLE PRECISION",conname_list,"");
        }
        if (na > 0) {
            Declare(fout,"    ","DOUBLE PRECISION",exprname_list,"");
        }
        Declare(fout,"    ","DOUBLE PRECISION",varname_list,"");
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
        fout << "    END SUBROUTINE " << Name() << "_beta\n";
        fout << endl;
    }
    fout << "END MODULE DEFINE_" << Name() << "_DDEs\n";
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
        fout << "PROGRAM " << Name() << "_demo\n";
        fout << endl;
        fout << "USE DEFINE_" << Name() << "_DDEs\n";
        fout << "USE DDE_SOLVER_M\n";
        fout << endl;
        fout << "IMPLICIT NONE\n";
        fout << endl;
        fout << "INTEGER, DIMENSION(2) :: NVAR = (/NEQN,NLAGS/)\n";
        fout << endl;
        fout << "TYPE(DDE_SOL) :: SOL\n";
        fout << "TYPE(DDE_OPTS) :: OPTS\n";
        fout << endl;
        if (!HasNonconstantDelay) {
            fout << "DOUBLE PRECISION, DIMENSION(" << Delays.size() << ") :: LAGS\n";
        }
        if (np > 0) {
            fout << "DOUBLE PRECISION, DIMENSION(" << np << ") :: p_\n";
        }
        fout << "DOUBLE PRECISION, DIMENSION(2) :: TSPAN\n";
        fout << endl;
        fout << "INTEGER :: I,J\n";
        fout << "CHARACTER(7+6*NEQN) :: F\n";
        fout << "DOUBLE PRECISION :: relerr, abserr, stoptime\n";
        if (HasPi) {
            fout << "DOUBLE PRECISION :: Pi\n";
        }
        if (nc > 0) {
            Declare(fout,"","DOUBLE PRECISION ::",conname_list,"");
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
                fout << "LAGS(" << k+1 << ") = " << Delays[k] << endl;
            }
        }
        fout << endl;
        fout << "TSPAN(1) = 0.0\n";
        fout << "TSPAN(2) = stoptime\n";
        fout << "OPTS = DDE_SET(RE=relerr,AE=abserr)\n";
        fout << endl;
        string lags_arg;
        if (HasNonconstantDelay) {
            lags_arg = Name() + "_beta";
        }
        else {
            lags_arg = "LAGS";
        } 
        fout << "SOL = DDE_SOLVER(NVAR," << Name() << "_ddes," << lags_arg <<
             "," << Name() << "_history,TSPAN,OPTIONS=OPTS)\n";
        fout << endl;
        fout << "F = \"(E17.8\"//REPEAT(\",E17.8\",NEQN)//\")\"\n";
        fout << "DO I = 1, SOL%NPTS\n";
        fout << "    WRITE(*,FMT=F) SOL%T(I), (SOL%Y(I,J),J=1,NEQN)\n";
        fout << "END DO\n";
        fout << endl;
        fout << "CALL RELEASE_ARRAYS(SOL, OPTS)\n";
        fout << endl;
        fout << "END PROGRAM " << Name() << "_demo\n";
        fout.close();
    }
}
