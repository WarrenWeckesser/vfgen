
//
//  vf_auto.cpp
//
//  This file defines the VectorField::PrintAUTO method.
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
#include <ginac/ginac.h> 

#include "vf.h"
#include "vf_utils.h"

using namespace std;
using namespace GiNaC;

//
// PrintAUTO -- The AUTO Code Generator.
//

void VectorField::PrintAUTO(map<string,string> options)
    {
    int nc, nv, np, na, nf;

    nc = conname_list.nops();
    nv = varname_list.nops();
    np = parname_list.nops();
    na = exprname_list.nops();
    nf = funcname_list.nops();

    if (options["lang"] == "fortran")
        {
        //
        // Generate FORTRAN code.
        //
        string filename = Name()+"_avf.f";
        ofstream fout;
        fout.open(filename.c_str());
        fout << csrc << left;
        // Also override the csrc style for powers.
        // IMPORTANT: This means we can NOT subsequently print C/C++ code!
        set_print_func<power,print_csrc>(print_power_as_fortran);

        fout << "c\n";
        fout << "c  AUTO FORTRAN file for the vector field named: " << Name() << endl;
        fout << "c\n";
        PrintVFGENComment(fout,"c  ");
        fout << "c\n";
        fout << "      subroutine func(n_,u_,icp_,par_,ijac_,F_,DFDU_,DFDP_)\n"; 
        fout << "      implicit none\n";
        fout << "      integer n_, icp_, ijac_\n";
        fout << "      double precision u_, par_, F_, DFDU_, DFDP_\n";
        fout << "      dimension u_(" << nv << "), F_(" << nv << ")\n";
        fout << "      dimension par_(" << np << ")\n";
        fout << "      dimension DFDU_(" << nv << "," << nv << "), DFDP_(" << nv << "," << np << ")\n";
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
            GetFromVector(fout,"      ",parname_list,"par_","()",1,"");
            }
        fout << "c     --- State variables ---\n";
        GetFromVector(fout,"      ",varname_list,"u_","()",1,"");
        if (na > 0)
            fout << "c     --- Expressions ---\n";
        for (int i = 0; i < na; ++i)
            {
            ostringstream os;
            os << left << csrc;
            os << exprname_list[i] << " = " << exprformula_list[i];
            F77Write(fout,os.str());
            }
        fout << endl;
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
        fout << "      if (ijac_ .eq. 0) return\n";
        fout << endl;
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
                os << "DFDU_(" << i+1 << ", " << j+1 << ") = " << df;
                F77Write(fout,os.str());
                }
            }
        fout << endl;
        fout << "      if (ijac_ .eq. 1) return\n";
        fout << endl;
        fout << "c     --- Derivatives with respect to the parameters ---\n";
        for (int i = 0; i < nv; ++i)
            {
            ex f = iterated_subs(varvecfield_list[i],expreqn_list);
            for (int j = 0; j < np; ++j)
                {
                symbol p = ex_to<symbol>(parname_list[j]);
                ostringstream os;
                os << left << csrc;
                ex df = f.diff(p);
                os << "DFDP_(" << i+1 << ", " << j+1 << ") = " << df;
                F77Write(fout,os.str());
                }
            }
        fout << endl;
        fout << "      return\n";
        fout << "      end\n";
        fout << "c\n";
        fout << "      subroutine stpnt(n_,u_,par_)\n";
        fout << "      implicit none\n";
        fout << "      integer n_;\n";
        fout << "      double precision u_, par_\n";
        fout << "      dimension u_(" << nv << "), par_(" << np << ")\n";

        if (nc > 0)
            F77Declare(fout,conname_list);
        if (np > 0)
            F77Declare(fout,parname_list);
        if (na > 0)
            F77Declare(fout,exprname_list);
        F77Declare(fout,varname_list);
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
            F77Write(fout, os.str());
            };
        /*
        if (na > 0)
            fout << "c     --- Expressions ---\n";
        for (int i = 0; i < na; ++i)
            {
            ostringstream os;
            os << left << csrc;
            os << exprname_list[i] << " = " << exprformula_list[i];
            F77Write(fout,os.str());
            }
        */
        fout << endl;
        fout << "c     --- Starting point -- Update these with correct values! ---\n";
        for (int i = 0; i < np; ++i)
            {
            ostringstream os;
            os << left << csrc;
            os << parname_list[i] << " = " << pardefval_list[i];
            F77Write(fout,os.str());
            }
        fout << "c\n";
        for (int i = 0; i < nv; ++i)
            {
            ostringstream os;
            os << left << csrc;
            os << varname_list[i] << " = " << vardefic_list[i];
            F77Write(fout,os.str());
            }
        fout << endl;
        for (int i = 0; i < np; ++i)
            fout << "      par_(" << i+1 << ") = " << parname_list[i] << endl;
        for (int i = 0; i < nv; ++i)
            fout << "      u_(" << i+1 << ") = " << varname_list[i] << endl;
        fout << "c\n";
        fout << "      return\n";
        fout << "      end\n";
        fout << "c" << endl;
        fout << "c  The remaining functions are just stubs." << endl;
        fout << "c  You will have to edit these by hand if you need these functions." << endl;
        fout << "c\n";
        fout << "      subroutine bcnd\n";
        fout << "      return\n";
        fout << "      end\n";
        fout << "c\n";
        fout << "      subroutine icnd\n";
        fout << "      return\n";
        fout << "      end\n";
        fout << "c\n";
        fout << "      subroutine fopt\n";
        fout << "      return\n";
        fout << "      end\n";
        fout << "c\n";
        fout << "      subroutine pvls\n";
        fout << "      return\n";
        fout << "      end\n";
        fout.close();
        }
    else
        {
        //
        // Generate C code.
        //
        string filename = Name()+"_avf.c";
        ofstream fout;
        fout.open(filename.c_str());
        fout << csrc << left;

        //
        //  Print the C file header information.
        //
        fout << "/*" << endl;
        fout << " *  " << filename << endl;
        fout << " *" << endl;
        fout << " *  AUTO C file for the vector field named: " << Name() << endl;
        fout << " *" << endl;
        PrintVFGENComment(fout," *  ");
        fout << " */" << endl;
        fout << endl;
        fout << "#include <math.h>\n";
        fout << "#include \"auto_f2c.h\"" << endl;
        fout << endl;
        //
        //  Print the AUTO function func(...).
        //
        fout << "/*" << endl;
        fout << " *  FUNC  Defines the vector field and its derivatives" << endl;
        fout << " */" << endl;
        fout << endl;
        fout << "int func(integer ndim_, const doublereal *u_, const integer *icp_," << endl;
        fout << "         const doublereal *par_, integer ijac_," << endl;
        fout << "         doublereal *f_, doublereal *dfdu_, doublereal *dfdp_)" << endl;
        fout << "    {" << endl;
        fout << "    integer dfdu__dim1, dfdp__dim1;" << endl;
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
        fout << endl;
        fout << "    dfdu__dim1 = ndim_;" << endl;
        fout << "    dfdp__dim1 = ndim_;" << endl;
        fout << endl;
        GetFromVector(fout,"    ",varname_list,"u_","[]",0,";");
        fout << endl;
        GetFromVector(fout,"    ",parname_list,"par_","[]",0,";");
        fout << endl;
        for (int i = 0; i < na; ++i)
            {
            fout << "    " << exprname_list[i] << " = " << exprformula_list[i] << ";" << endl;
            }
        fout << endl;
        fout << "    /*" << endl;
        fout << "     *  The vector field" << endl;
        fout << "     */" << endl;
        fout << endl;
        for (int i = 0; i < nv; ++i)
            {
            fout << "    f_[" << i << "]" << " = " << varvecfield_list[i] << ";" << endl;
            }
        fout << endl;
        fout << "    if (ijac_ == 0)" << endl;
        fout << "        return 0;" << endl;
        fout << endl;
        //
        // The Jacobian section...
        //
        fout << "    /*" << endl;
        fout << "     *  The Jacobian" << endl;
        fout << "     */" << endl;
        fout << endl;
        for (int i = 0; i < nv; ++i)
            {
            ex f = iterated_subs(varvecfield_list[i],expreqn_list);
            for (int j = 0; j < nv; ++j)
                {
                symbol v = ex_to<symbol>(varname_list[j]);

                fout << "    ARRAY2D(dfdu_," << i << "," << j << ") = " << f.diff(v) << ";" << endl;
                }
            }
        fout << endl;
        //
        // The Jacobian with respect to the parameters
        //
        fout << "    if (ijac_ == 1)" << endl;
        fout << "        return 0;" << endl;
        fout << endl;
        fout << "    /*" << endl;
        fout << "     *  The Jacobian with respect to the parameters." << endl;
        fout << "     */" << endl;
        fout << endl;
        for (int i = 0; i < nv; ++i)
            {
            ex f = iterated_subs(varvecfield_list[i],expreqn_list);
            for (int j = 0; j < np; ++j)
                {
                symbol p = ex_to<symbol>(parname_list[j]);
                fout << "    ARRAY2D(dfdp_," << i << "," << j << ") = " << f.diff(p) << ";" << endl;
                }
            }
        fout << endl;
        fout << "    return 0;" << endl;
        fout << "    }" << endl;
        fout << endl;
        fout << "/*" << endl;
        fout << " *  STPNT  Gives a starting point" << endl;
        fout << " */" << endl;
        fout << endl;
        fout << "int stpnt(integer ndim_, doublereal t_, doublereal *u_, doublereal *par_)" << endl;
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
        fout << endl;
        fout << "    /* Change the parameter values and the starting point to correct values! */" << endl;
        fout << endl;
        for (int i = 0; i < np; ++i)
            fout << "    " << parname_list[i] << " = " << pardefval_list[i] << ";\n" ;
        fout << endl;
        // Assign the DefaultInitialConditions to the variable. Not necessarily the
        // correct thing to do, but it might work for some purposes.
        for (int i = 0; i < nv; ++i)
            fout << "    " << varname_list[i] << " = " << vardefic_list[i] << ";\n" ;
        fout << endl;
        for (int i = 0; i < np; ++i)
            fout << "    par_[" << i << "] = " << parname_list[i] << ";" << endl;
        for (int i = 0; i < nv; ++i)
            fout << "    u_[" << i << "] = " << varname_list[i] << ";" << endl;
        fout << endl;
        fout << "    return 0;" << endl;
        fout << "    }" << endl;
        fout << endl;
        fout << "/*" << endl;
        fout << " *  The remaining functions are just stubs." << endl;
        fout << " *  You will have to edit these by hand if you need these functions." << endl;
        fout << " */" << endl;
        fout << endl;
        fout << "/*" << endl;
        fout << " *  BCND  Defines the boundary conditions" << endl;
        fout << " */" << endl;
        fout << endl;
        fout << "int bcnd(integer ndim_, const doublereal *par_, const integer *icp_," << endl;
        fout << "         integer nbc_, const doublereal *u0_, const doublereal *u1_, integer ijac_," << endl;
        fout << "         doublereal *fb_, doublereal *dbc_)" << endl;
        fout << "    {" << endl;
        fout << "    return 0;" << endl;
        fout << "    }" << endl;
        fout << endl; 
        fout << "/*" << endl;
        fout << " *  ICND  Defines the integral conditions" << endl;
        fout << " */" << endl;
        fout << endl;
        fout << "int icnd(integer ndim_, const doublereal *par_, const integer *icp_," << endl;
        fout << "         integer nint_, const doublereal *u_, const doublereal *uold_," << endl;
        fout << "         const doublereal *udot_, const doublereal *upold_, integer ijac_," << endl;
        fout << "         doublereal *fi_, doublereal *dint_)" << endl;
        fout << "    {" << endl;
        fout << "    return 0;" << endl;
        fout << "    }" << endl;
        fout << endl; 
        fout << "/*" << endl;
        fout << " *  FOPT" << endl;
        fout << " */" << endl;
        fout << endl;
        fout << "int fopt(integer ndim_, const doublereal *u_, const integer *icp_," << endl;
        fout << "         const doublereal *par_, integer ijac_," << endl;
        fout << "         doublereal *fs_, doublereal *dfdu_, doublereal *dfdp_)" << endl;
        fout << "    {" << endl;
        fout << "    return 0;" << endl;
        fout << "    }" << endl;
        fout << endl;
        fout << "/*" << endl;
        fout << " *  PVLS" << endl;
        fout << " */" << endl;
        fout << endl;
        fout << "int pvls(integer ndim_, const doublereal *u_, doublereal *par_)" << endl;
        fout << "    {" << endl;
        fout << "    return 0;" << endl;
        fout << "    }" << endl; 
        fout.close();
        }
    }
