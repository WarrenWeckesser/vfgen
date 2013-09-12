
//
//  vf_radau5.cpp
//
//  This file defines the VectorField::PrintRadau5 method.
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
// PrintRadau5 --
//

void VectorField::PrintRadau5(map<string,string> options)
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
    fout << "c These functions are to be used with the Fortran ODE solver RADAU5.\n";
    fout << "c\n";
    PrintVFGENComment(fout,"c ");
    fout << "c\n";
    // fout << "      subroutine " << Name() << "_rhs(n_,t_,y_,f_,rpar_,ipar_)\n";
    F77Write(fout,"subroutine " + Name() + "_rhs(n_,t_,y_,f_,rpar_,ipar_)");
    fout << "      implicit none\n";
    fout << "      integer n_, ipar_\n";
    fout << "      double precision t_, y_, f_, rpar_\n";
    fout << "      dimension y_(" << nv << "), f_(" << nv << "), rpar_(" << np << ")\n";
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
        fout << "c     --- Parameters ---\n";
    GetFromVector(fout,"      ",parname_list,"rpar_","()",1,"");
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
    fout << "      subroutine " << Name() << "_jac(n_,t_,y_,dfy_,ldfy_,\n";
    fout << "     &                             rpar_,ipar_)\n";
    fout << "      implicit none\n";
    fout << "      integer n_, ldfy_, ipar_\n";
    fout << "      double precision t_, y_, dfy_, rpar_\n";
    fout << "      dimension y_(" << nv << "), dfy_(ldfy_," << nv << "), rpar_(" << np << ")\n";
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
        fout << "c     --- Parameters ---\n";
    GetFromVector(fout,"      ",parname_list,"rpar_","()",1,"");
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
            os << "dfy_(" << i+1 << ", " << j+1 << ") = " << f.diff(v);
            F77Write(fout,os.str());
            }
        }
    fout << endl;
    fout << "      return\n";
    fout << "      end\n";
    fout << endl;
    fout << "c\n";
    fout << "c " << Name() << "_out -- output subroutine\n";
    fout << "c\n";
    fout << "      subroutine " << Name() << "_out(nr,told,t,y,cont,lrc,n,\n";
    fout << "     &                  rpar,ipar,irtrn)\n";
    fout << "      implicit none\n";
    fout << "      integer nr, lrc, n, ipar, irtrn\n";
    fout << "      double precision told, t, y, cont, rpar\n";
    fout << "      integer i\n";
    fout << "      dimension y(" << nv << "), rpar(" << np << ")\n";
    fout << endl;
    fout << "      write (6,99) t, (y(i), i = 1," << nv << ")\n";
    fout << "99    format(1x,f10.5," << nv << "E18.10)\n";
    fout << "      return\n";
    fout << "      end\n";
    fout.close();

    if (options["demo"] == "yes")
        {
        //
        // Create the demo function.
        //
        string filename = Name() + "_dr5.f";
        ofstream fout;
        fout.open(filename.c_str());
        fout << left;
        fout << csrc;

        fout << "c\n";
        fout << "c " << filename << endl;
        fout << "c" << endl;
        fout << "c Fortran 77 program that uses RADAU5 to solve the differential equations\n";
        fout << "c defined in the vector field '" << Name() << "'\n";
        fout << "c\n";
        PrintVFGENComment(fout,"c ");
        fout << "c\n";
        fout << "      program " << Name() << endl;
        fout << "      implicit none\n";
        fout << "      integer nd_, lwork_, liwork_\n";
        fout << "      parameter (nd_=" << nv << ", lwork_=" << 4*nv*nv + 12*nv + 20 << ", liwork_=" << 3*nv+20 << ")\n";
        fout << "      double precision y_, rpar_, work_, iwork_\n";
        fout << "      dimension y_(nd_), work_(lwork_), iwork_(liwork_)\n";
        fout << "      integer n_, ipar_, ijac_, mljac_, mujac_, imas_, mlmas_, mumas_\n";
        fout << "      integer iout_, itol_, idid_\n";
        fout << "      dimension rpar_(" << np << ")\n";
        fout << "      integer i_\n";
        fout << "      double precision t_, tstop_\n";
        fout << "      double precision atol_, rtol_, h_\n";
        if (nc > 0)
            F77Declare(fout,conname_list);
        if (np > 0)
            F77Declare(fout,parname_list);
        F77Declare(fout,varname_list);
        if (HasPi)
            fout << "      double precision Pi\n";
        fout << "      external " << Name() << "_rhs\n";
        fout << "      external " << Name() << "_jac\n";
        fout << "      external " << Name() << "_out\n";
        fout << endl;
        fout << "      n_ = " << nv << "\n";
        fout << "      ijac_ = 1\n";
        fout << "      mljac_ = n_\n";
        fout << "      imas_ = 0\n";
        fout << "      iout_ = 1\n";
        if (HasPi)
            fout << "      Pi = 3.1415926535897932385D0\n";
        fout << "c     --- t range ---\n";
        fout << "      t_ = 0.0D0\n";
        fout << "      tstop_  = 10.0D0\n";
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
            fout << "      rpar_(" << i+1 << ") = " << parname_list[i] << endl;
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
        fout << "      itol_ = 0\n";
        fout << "c     --- Initial step size ---\n";
        fout << "      h_ = 1.0D-5\n";
        fout << "c     --- Set default values ---\n";
        fout << "      do i_ = 1, 20\n";
        fout << "          iwork_(i_) = 0\n";
        fout << "          work_(i_) = 0.0D0\n";
        fout << "      end do\n";
        fout << "c     --- Call RADAU5 ---\n";
        fout << "      call radau5(n_," << Name() << "_rhs, t_, y_, tstop_, h_,\n";
        fout << "     &           rtol_, atol_, itol_,\n";
        fout << "     &           " << Name() << "_jac, ijac_, mljac_, mujac_,\n";
        fout << "     &           " << Name() << "_rhs, imas_, mlmas_, mumas_,\n";
        fout << "     &           " << Name() << "_out, iout_,\n";
        fout << "     &           work_, lwork_, iwork_, liwork_, rpar_, ipar_, idid_)\n";
        // fout << "      write (6,99) t_, (y_(i_), i_ = 1," << nv << ")\n";
        // fout << "99    format(1x,f10.5," << nv << "E18.10)\n";
        fout << "      stop\n";
        fout << "      end\n";
        fout << endl;
        fout.close();
        }
    }
