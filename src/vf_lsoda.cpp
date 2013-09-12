
//
//  vf_lsoda.cpp
//
//  This file defines the VectorField::PrintLSODA method.
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
#include <sstream>
#include <string>
#include <cassert>
#include <ginac/ginac.h> 

#include "vf.h"
#include "vf_utils.h"

using namespace std;
using namespace GiNaC;


//
// PrintLSODA --
//

void VectorField::PrintLSODA(map<string,string> options)
    {
    int nc, np, nv, na, nf;

    nc = conname_list.nops();
    nv = varname_list.nops();
    np = parname_list.nops();
    na = exprname_list.nops();
    nf = funcname_list.nops();

    //
    //  Create the vector field function.
    //
    string vf_filename = Name()+"_rhs.f";
    ofstream fout;
    fout.open(vf_filename.c_str());
    fout << left;
    fout << csrc;
    // Also override the csrc style for powers.
    // IMPORTANT: This means we can NOT subsequently print C/C++ code!
    set_print_func<power,print_csrc>(print_power_as_fortran);

    fout << "c\n";
    fout << "c " << vf_filename << endl;
    fout << "c\n";
    fout << "c Vector field functions for the vector field '" << Name() << "'\n";
    fout << "c These functions are to be used with the Fortran ODE solver LSODA.\n";
    fout << "c\n";
    PrintVFGENComment(fout,"c ");
    fout << "c\n";
    F77Write(fout,"subroutine " + Name() + "_rhs(n_,t_,y_,f_)");
    fout << "      implicit none\n";
    fout << "      integer n_\n";
    fout << "      double precision t_, y_, f_\n";
    int ydim = nv;
    if (options["parstyle"] != "common")
        ydim = nv + np;
    fout << "      dimension y_(" << ydim << "), f_(" << nv << ")\n";
    if (options["parstyle"] == "common")
        {
        fout << "      double precision rpar_\n";
        fout << "      dimension rpar_(" << np << ")\n";
        fout << "      common /" << Name()+"_parameters/ rpar_\n";
        }
    if (nc > 0)
        F77Declare(fout,conname_list);
    if (np > 0)
        F77Declare(fout,parname_list);
    if (na > 0)
        F77Declare(fout,exprname_list);
    F77Declare(fout,varname_list);
    fout << endl;
    if (nc > 0)
        fout << "c     --- Constants ---\n";
    for (int i = 0; i < nc; ++i)
        {
        ostringstream os;
        os << left << csrc;
        os << conname_list[i] << " = " << convalue_list[i];
        F77Write(fout, os.str());
        }
    if (np > 0)
        {
        fout << "c     --- Parameters ---\n";
        if (options["parstyle"] == "common")
            GetFromVector(fout,"      ",parname_list,"rpar_","()",1,"");
        else
            GetFromVector(fout,"      ",parname_list,"y_","()",nv+1,"");
        }
    fout << "c     --- State variables ---\n";
    GetFromVector(fout,"      ",varname_list,"y_","()",1,"");
    if (na > 0)
        fout << "c     --- Expressions ---\n";
    for (int i = 0; i < na; ++i)
        {
        ostringstream os;
        os << left << csrc;
        os << exprname_list[i] << " = " << exprformula_list[i];
        F77Write(fout,os.str());
        }
    fout << "c     --- The vector field ---\n";
    for (int i = 0; i < nv; ++i)
        {
        ex f = varvecfield_list[i];
        ostringstream os;
        os << left << csrc;
        os << "f_(" << (i+1) << ")" << " = " << f;
        F77Write(fout,os.str());
        }
    fout << endl;
    fout << "      return\n";
    fout << "      end\n";
    fout << endl;
    fout << "c\n";
    fout << "c Jacobian of the vector field '" << Name() << "'\n";
    fout << "c\n";
    fout << "c The subroutine assumes that the full Jacobian is to be computed.\n";
    fout << "c ml_ and mu_ are ignored, and nrowpd_ is assumed to be n_.\n";
    fout << "c\n";
    fout << "      subroutine " << Name() << "_jac(n_,t_,y_,ml_,mu_,jac_,nrowpd_)\n";
    fout << "      implicit none\n";
    fout << "      integer n_, ml_, mu_, nrowpd_\n";
    fout << "      double precision t_, y_, jac_\n";
    fout << "      dimension y_(" << ydim << "), jac_(nrowpd_," << nv << ")\n";
    if (options["parstyle"] == "common")
        {
        fout << "      double precision rpar_\n";
        fout << "      dimension rpar_(" << np << ")\n";
        fout << "      common /" << Name()+"_parameters/ rpar_\n";
        }
    if (nc > 0)
        F77Declare(fout,conname_list);
    if (np > 0)
        F77Declare(fout,parname_list);
    F77Declare(fout,varname_list);
    fout << endl;
    if (nc > 0)
        fout << "c     --- Constants ---\n";
    for (int i = 0; i < nc; ++i)
        {
        ostringstream os;
        os << left << csrc;
        os << conname_list[i] << " = " << convalue_list[i];
        F77Write(fout,os.str());
        }
    if (np > 0)
        {
        fout << "c     --- Parameters ---\n";
        if (options["parstyle"] == "common")
            GetFromVector(fout,"      ",parname_list,"rpar_","()",1,"");
        else
            GetFromVector(fout,"      ",parname_list,"y_","()",nv+1,"");
        }
    fout << "c     --- State variables ---\n";
    GetFromVector(fout,"      ",varname_list,"y_","()",1,"");
    fout << "c     --- Jacobian ---\n";
    for (int i = 0; i < nv; ++i)
        {
        ex f = iterated_subs(varvecfield_list[i],expreqn_list);
        for (int j = 0; j < nv; ++j)
            {
            symbol v = ex_to<symbol>(varname_list[j]);
            ostringstream os;
            os << left << csrc;
            ex df = f.diff(v);
            if (df != 0)
                {
                os << "jac_(" << i+1 << ", " << j+1 << ") = " << df;
                F77Write(fout,os.str());
                }
            else
                fout << "c     jac_(" << i+1 << ", " << j+1 << ") = 0\n";
            }
        }
    fout << endl;
    fout << "      return\n";
    fout << "      end\n";
    fout << endl;
    
    if (options["func"] == "yes" & nf > 0)
        {
        //
        // Print the user-defined functions.
        // A single function is created that puts all the
        // user-defined function values in an array.  This
        // function is defined so that it can be used with
        // the LSODAR rootfinding capability.
        //
        fout << endl;
        fout << "c" << endl;
        fout << "c  User-defined functions. " << endl;
        fout << "c" << endl;
        fout << "      subroutine "<< Name() << "_func(neq_, t_, y_, nf_, func_)\n";
        fout << "      implicit none\n";
        fout << "      integer neq_, nf_\n";
        fout << "      double precision t_, y_, func_\n";
        fout << "      dimension y_(" << ydim << "), func_(" << nf << ")\n";
        if (options["parstyle"] == "common")
            {
            fout << "      double precision rpar_\n";
            fout << "      dimension rpar_(" << np << ")\n";
            fout << "      common /" << Name()+"_parameters/ rpar_\n";
            }
        if (nc > 0)
            F77Declare(fout,conname_list);
        if (np > 0)
            F77Declare(fout,parname_list);
        if (na > 0)
            F77Declare(fout,exprname_list);
        F77Declare(fout,varname_list);
        fout << endl;
        if (HasPi)
            {
            fout << "      double precision Pi\n";
            fout << "      Pi = ";
            PrintPi(fout);
            fout << "D0\n";
            }
        if (nc > 0)
            fout << "c     --- Constants ---\n";
        for (int i = 0; i < nc; ++i)
            {
            ostringstream os;
            os << left << csrc;
            os << conname_list[i] << " = " << convalue_list[i];
            F77Write(fout,os.str());
            }
        if (np > 0)
            {
            fout << "c     --- Parameters ---\n";
            if (options["parstyle"] == "common")
                GetFromVector(fout,"      ",parname_list,"rpar_","()",1,"");
            else
                GetFromVector(fout,"      ",parname_list,"y_","()",nv+1,"");
            }
        fout << "c     --- State variables ---\n";
        GetFromVector(fout,"      ",varname_list,"y_","()",1,"");
        if (na > 0)
            fout << "c     --- Expressions ---\n";
		for (int i = 0; i < na; ++i)
			{
			ostringstream os;
			os << left << csrc;
			os << exprname_list[i] << " = " << exprformula_list[i];
			F77Write(fout,os.str());
			}        

        fout << "c     --- Compute the functions ---\n";
        for (int n = 0; n < nf; ++n)
            {
            fout << "c     " << funcname_list[n] << ":" << endl;
            fout << "      func_(" << n+1 << ") = " << funcformula_list[n] << ";" << endl;
            }
        fout << endl;
        fout << "      return\n";
        fout << "      end\n";
        }
    
    
    fout.close();

    if (options["demo"] == "yes")
        {
        //
        // Create the demo function.
        //
        string filename = Name() + "_demo.f";
        ofstream fout;
        fout.open(filename.c_str());
        fout << left;
        fout << csrc;

        int lrn = 20 + 16*nv;
        int lrs = 22 + 9*nv + nv*nv;
        int lrw  = max(lrn,lrs);
        int liw  = 20 + nv; 
        fout << "c\n";
        fout << "c " << filename << endl;
        fout << "c" << endl;
        fout << "c Fortran 77 program that uses LSODA to solve the differential equations\n";
        fout << "c defined in the vector field '" << Name() << "'\n";
        fout << "c\n";
        PrintVFGENComment(fout,"c ");
        fout << "c\n";
        fout << "      program " << Name() << endl;
        fout << endl;
        fout << "      implicit none\n";
        fout << endl;
        fout << "      external " << Name() << "_rhs\n";
        fout << "      external " << Name() << "_jac\n";
        fout << endl;
        fout << "      double precision atol_, rtol_, y_, t_, tout_, tfinal_, rwork_\n";
        fout << "      integer iwork_\n";
        fout << "      dimension y_(" << ydim << "), rwork_(" << lrw << "), iwork_(" << liw << ")\n";
        if (options["parstyle"] == "common")
            {
            fout << "      double precision rpar_\n";
            fout << "      dimension rpar_(" << np << ")\n";
            fout << "      common /" << Name()+"_parameters/ rpar_\n";
            }
        fout << "      integer neq_, i_, j_, nsteps_\n";
        fout << "      integer itol_, iopt_, itask_, istate_, jt_, lrw_, liw_\n";
        if (nc > 0)
            F77Declare(fout,conname_list);
        if (np > 0)
            F77Declare(fout,parname_list);
        F77Declare(fout,varname_list);
        if (HasPi)
            fout << "      double precision Pi\n";
        if (HasPi)
            fout << "      Pi = 3.1415926535897932385D0\n";
        fout << "c     --- t range ---\n";
        fout << "      t_ = 0.0D0\n";
        fout << "      tfinal_  = 10.0D0\n";
        fout << "      nsteps_ = 100\n";
        if (nc > 0)
            fout << "c     --- Constants ---\n";
        for (int i = 0; i < nc; ++i)
            {
            ostringstream os;
            os << left << csrc;
            os << conname_list[i] << " = " << convalue_list[i];
            F77Write(fout,os.str());
            }
        if (np > 0)
            fout << "c     --- Parameters ---\n";
        for (int i = 0; i < np; ++i)
            {
            ostringstream os;
            os << left << csrc;
            os << parname_list[i] << " = " << pardefval_list[i];
            F77Write(fout,os.str()); 
            }
        for (int i = 0; i < np; ++i)
            if (options["parstyle"] == "common")
                fout << "      rpar_(" << i+1 << ") = " << parname_list[i] << endl;
            else
                fout << "      y_(" << nv+i+1 << ") = " << parname_list[i] << endl;
        fout << "c     --- Initial conditions ---\n";
        for (int i = 0; i < nv; ++i)
            {
            ostringstream os;
            os << left << csrc;
            os << varname_list[i] << " = " << vardefic_list[i];
            F77Write(fout,os.str());
            }
        for (int i = 0; i < nv; ++i)
            fout << "      y_(" << i+1 << ") = " << varname_list[i] << endl;
        fout << "c     --- Solver tolerances ---\n";
        fout << "      rtol_ = 1.0D-6\n";
        fout << "      atol_ = 1.0D-8\n";
        fout << "      itol_ = 1\n";
        fout << "c     --- Other LSODA parameters ---\n";
        fout << "      neq_ = " << nv << "\n";
        fout << "      itask_ = 1\n";
        fout << "      istate_ = 1\n";
        fout << "      iopt_ = 0\n";
        fout << "      lrw_ = " << lrw << endl;
        fout << "      liw_ = " << liw << endl;
        fout << "c\n";
        fout << "c     LSODE and LSODA take the same arguments, so either may\n";
        fout << "c     be used in the loop below.  jt_ must be set as follows:\n";
        fout << "c     jt_ =  1 for LSODA\n";
        fout << "c     jt_ = 10 for LSODE, non-stiff (Adams) method\n";
        fout << "c     jt_ = 21 for LSODE, stiff (BDF) method\n";
        fout << "c     See the documentation in the Fortran file for more details.\n";
        fout << "      jt_ = 1\n";
        fout << "c     --- Print the first point ---\n";
        fout << "      write (6,49) t_, (y_(j_), j_ = 1," << nv << ")\n";
        fout << "c     --- Call DLSODA in a loop to compute the solution ---\n";
        fout << "      do 40 i_ = 1,nsteps_\n";
        fout << "          tout_ = (i_*tfinal_)/nsteps_\n";
        fout << "          call DLSODA(" << Name() << "_rhs, neq_, y_, t_, tout_,\n";
        fout << "     &           itol_, rtol_, atol_, itask_, istate_, iopt_,\n";
        fout << "     &           rwork_, lrw_, iwork_, liw_,\n";
        fout << "     &           " << Name() << "_jac, jt_)\n";
        fout << "          if (istate_ .lt. 0) goto 80\n";
        fout << "40        write (6,49) t_, (y_(j_), j_ = 1," << nv << ")\n";
        fout << "49        format(1X,F10.6," << nv << "E18.10)\n";
        fout << "      stop\n";
        fout << "80    write (6,89) istate_\n";
        fout << "89    format(1X,\"Error: istate=\",I3)\n";
        fout << "      stop\n";
        fout << "      end\n";
        fout << endl;
        fout.close();
        }
    }
