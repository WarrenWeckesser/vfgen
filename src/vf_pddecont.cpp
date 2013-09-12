
//
//  vf_pddecont.cpp
//
//  This file defines the VectorField::PrintPDDECONT method.
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
// PDDEC_ConvertDelaysToZlags(ex& f)
//
// This function converts each subexpression of the form delay(delayexpr,del)
// in f to an expression in terms of Zlags_(i,j).
//

void VectorField::PDDEC_ConvertDelaysToZlags(ex& f)
    {
    exset dlist;
    f.find(delay(wild(1),wild(2)),dlist);
    for (exset::const_iterator iter = dlist.begin(); iter != dlist.end(); ++iter)
        {
        ex delayfunc = *iter;
        ex delayexpr = delayfunc.op(0);
        lst vars = FindVarsInEx(delayexpr);
        ex del = delayfunc.op(1);
        int dindex = FindDelay(del);
        assert(dindex != -1);
        for (lst::const_iterator viter = vars.begin(); viter != vars.end(); ++viter)
            {
            int vindex = FindVar(ex_to<symbol>(*viter));
            delayexpr = delayexpr.subs(*viter == Zlags_(vindex,dindex+1));
            }
        f = f.subs(delayfunc == delayexpr);
        }
    }

//
// PDDEC_ConvertStateToZlags(ex& f)
//
// This function converts all references to state variables--current time or
// delayed--to an expression in terms of Zlags_(i,j).
//

void VectorField::PDDEC_ConvertStateToZlags(ex& f)
    {
    exset dlist;
    f.find(delay(wild(1),wild(2)),dlist);
    for (exset::const_iterator iter = dlist.begin(); iter != dlist.end(); ++iter)
        {
        ex delayfunc = *iter;
        ex delayexpr = delayfunc.op(0);
        lst vars = FindVarsInEx(delayexpr);
        ex del = delayfunc.op(1);
        int dindex = FindDelay(del);
        assert(dindex != -1);
        for (lst::const_iterator viter = vars.begin(); viter != vars.end(); ++viter)
            {
            int vindex = FindVar(ex_to<symbol>(*viter));
            delayexpr = delayexpr.subs(*viter == Zlags_(vindex,dindex+1));
            }
        f = f.subs(delayfunc == delayexpr);
        }
    int nv = varname_list.nops();    
    for (int i = 0; i < nv; ++i)
        {
        f = f.subs(varname_list[i] == Zlags_(i,0));
        }
    }

void VectorField::PDDEC_PrintParDerivs(ofstream &dout, const vector<ex> &vf0)
    {
    // int nc = conname_list.nops();
    int nv = varname_list.nops();
    int np = parname_list.nops();
    // int na = exprname_list.nops();
    // int nf = funcname_list.nops();


    for (int j = 0; j < np; ++j)
        {
        if (j == 0)
            {    
            dout << "    if (j == 1)\n";
            dout << "        {\n";
            }
        else
            {
            dout << "    else if (j == " << j+1 << ")\n";
            dout << "        {\n";
            }
        dout << "        // Derivative wrt " << parname_list[j] << endl;
        for (int i = 0; i < nv; ++i)
            {
            symbol p = ex_to<symbol>(parname_list[j]);
            ex df = vf0[i].diff(p);
            dout << "        jac_(" << i << ",0)" << " = " << df << ";" << endl;
            }
        dout << "        }\n";
        }
    }
    

void VectorField::PDDEC_PrintJacobians(ofstream &dout, const vector<ex> &vf0)
    {
    // int nc = conname_list.nops();
    int nv = varname_list.nops();
    // int np = parname_list.nops();
    // int na = exprname_list.nops();
    // int nf = funcname_list.nops();

    int nd = Delays.size();
    for (int k = 0; k < nd+1; ++k)
        {
        if (k == 0)
            {
            dout << "    if (k == 0)\n";
            dout << "        {\n";
            dout << "        // Derivatives wrt the state variables\n";
            }
        else
            {
            dout << "    else if (k == " << k << ")\n";
            dout << "        {\n";
            dout << "        // Derivatives wrt state variables with delay " << Delays[k-1] << endl;
            }
        for (int i = 0; i < nv; ++i)
            {
            ex f = vf0[i];
            for (int j = 0; j < nv; ++j)
                {
                symbol vtmp_("vtmp_");
                ex fj = f.subs(Zlags_(j,k) == vtmp_);
                ex df = fj.diff(vtmp_);
                df = df.subs(vtmp_ == Zlags_(j,k));
                dout << "        jac_(" << i << "," << j << ")" << " = " << df << ";" << endl;
                }
            }
        dout << "        }\n";
        }
    }

void VectorField::PDDEC_PrintXandParJacobians(ofstream &dout, const vector<ex> &vf0)
    {
    // int nc = conname_list.nops();
    int nv = varname_list.nops();
    int np = parname_list.nops();
    // int na = exprname_list.nops();
    // int nf = funcname_list.nops();

    int nd = Delays.size();
    for (int k = 0; k < nd+1; ++k)
        {
        if (k == 0)
            {
            dout << "    if (k == 0)\n";
            dout << "        {\n";
            dout << "        // Derivatives wrt the state variables\n";
            }
        else
            {
            dout << "        }\n";
            dout << "    else if (k == " << k << ")\n";
            dout << "        {\n";
            dout << "        // Derivatives wrt state variables with delay " << Delays[k-1] << endl;
            }
        for (int m = 0; m < np; ++m)
            {
            if (m == 0)
                {
                dout << "        if (j == " << m+1 << ")\n";
                dout << "            {\n";
                }
            else
                {
                dout << "        else if (j == " << m+1 << ")\n";
                dout << "            {\n";
                }
            dout << "            // Derivative wrt " << parname_list[m] << "\n";
            for (int i = 0; i < nv; ++i)
                {
                ex f = vf0[i];
                for (int j = 0; j < nv; ++j)
                    {
                    symbol vtmp_("vtmp_");
                    ex fj = f.subs(Zlags_(j,k) == vtmp_);
                    ex df = fj.diff(vtmp_);
                    df = df.subs(vtmp_ == Zlags_(j,k));
                    symbol p = ex_to<symbol>(parname_list[m]);
                    df = df.diff(p);
                    // if (df != 0)
                    dout << "            jac_(" << i << "," << j << ")" << " = " << df << ";" << endl;
                    }
                }
            dout << "            }\n";
            }
        }
    dout << "        }\n";
    }

ex pddec_second_deriv(const ex &f, int lag1, int var1, int lag2, int var2)
    {
    symbol Z_("Z_");
    ex fj = f.subs(Zlags_(var1,lag1) == Z_);
    ex df = fj.diff(Z_);
    df = df.subs(Z_ == Zlags_(var1,lag1));
    df = df.subs(Zlags_(var2,lag2) == Z_);
    df = df.diff(Z_);
    df = df.subs(Z_ == Zlags_(var2,lag2));
    return df;
    }

void VectorField::PDDEC_PrintHessiansTimesV(ofstream &dout, const vector<ex> &vf0)
    {
    // int nc = conname_list.nops();
    int nv = varname_list.nops();
    // int np = parname_list.nops();
    // int na = exprname_list.nops();
    // int nf = funcname_list.nops();

    int nd = Delays.size();
    for (int k1 = 0; k1 < nd+1; ++k1)
        {
        if (k1 == 0)
            {
            dout << "    if (k1 == 0)\n";
            dout << "        {\n";
            dout << "        // Derivatives wrt the state variables\n";
            }
        else
            {
            dout << "        }\n";
            dout << "    else if (k1 == " << k1 << ")\n";
            dout << "        {\n";
            dout << "        // Derivatives wrt state variables with delay " << Delays[k1-1] << endl;
            }

        for (int k2 = 0; k2 < nd+1; ++k2)
            {
            if (k2 == 0)
                {
                dout << "        if (k2 == 0)\n";
                dout << "            {\n";
                dout << "            // Derivatives wrt the state variables\n";
                }
            else
                {
                dout << "        else if (k2 == " << k2 << ")\n";
                dout << "            {\n";
                dout << "            // Derivatives wrt state variables with delay " << Delays[k2-1] << endl;
                }

            for (int i = 0; i < nv; ++i)
                {
                ex f = vf0[i];
                for (int j = 0; j < nv; ++j)
                    {
                    ostringstream os;
                    os << csrc;
                    for (int h = 0; h < nv; ++h)
                        {
                        ex d2f = pddec_second_deriv(f,k1,j,k2,h);
                        if (d2f != 0)
                            {
                            if (os.str() != "")
                                os << " + ";
                            os << "(" << d2f << ")*v_(" << h << ",m)";
                            }
                        }
                    if (os.str() != "")
                        dout << "            jac_(" << i << "," << j << ")" << " = " << os.str() << ";\n";
                    else
                        dout << "            jac_(" << i << "," << j << ") = 0.0;\n";
                    }

                }
            dout << "            }\n";
            }
        }
    dout << "        }\n";
    }
    
//
// PrintPDDECONT -- The PDDE-CONT code generator.
//

void VectorField::PrintPDDECONT(map<string,string> options)
    {
    int nc = conname_list.nops();
    int nv = varname_list.nops();
    int np = parname_list.nops();
    int na = exprname_list.nops();
    // int nf = funcname_list.nops();

    // cerr << "Delays:\n";
    // for (vector<ex>::iterator di = Delays.begin(); di != Delays.end(); ++di)
    //     cerr << "   " << *di << endl;

    //
    //  Create the system definition file.
    //
    string sys_filename = "sys-" + Name() + ".cpp";
    ofstream sys_out;
    sys_out.open(sys_filename.c_str());
    sys_out << csrc;
    sys_out << "//" << endl;
    sys_out << "// " << sys_filename << endl;
    sys_out << "//" << endl;
    sys_out << "// PDDE-CONT System Definition file for the VFGEN vector field: " << Name() << endl;
    sys_out << "//" << endl;
    PrintVFGENComment(sys_out,"// ");
    sys_out << "//" << endl;
    sys_out << endl;
    // sys_out << "#include <cstdlib>\n";
    sys_out << "#include <cmath>\n";
    sys_out << "#include \"pddesys.h\"\n";
    sys_out << endl;
    sys_out << "extern \"C\" {\n";
    sys_out << endl;
    sys_out << "int sys_ndim()  { return " << nv << "; }  // System dimension\n";
    sys_out << "int sys_npar()  { return " << 1+np << "; }  // Number of parameters, plus one (for the period)\n";
    sys_out << "int sys_ntau()  { return " << 1+Delays.size() << "; }  // Number of delays, plus one\n";
    sys_out << "int sys_nderi() { return 2; }  // Order of derivatives computed here\n";
    sys_out << endl;

    //
    //  Create the vector field function sys_rhs.m
    //
    sys_out << "//" << endl;
    sys_out << "// sys_rhs(...) computes the vector field.\n";
    sys_out << "//" << endl;
    // Include a list of the lags in the comments.
    sys_out << "// The lags are: {";
    for (unsigned k = 0; k < Delays.size(); ++k)
        {
        sys_out << Delays[k];
        if (k < Delays.size()-1)
            sys_out << ", "; 
        }
    sys_out << "}" << endl;
    sys_out << "//\n";
    sys_out << "// If X(t) is the state vector at time t, then\n";
    sys_out << "//    Zlags_ = [ X(t) ";
    for (unsigned k = 0; k < Delays.size(); ++k)
        {
        sys_out << "X(t-" << Delays[k] << ")";
        if (k < Delays.size()-1)
            sys_out << " "; 
        }
    sys_out << " ]\n";
    // Function definition starts here.
    sys_out << "//" << endl;
    sys_out << "void sys_rhs(Vector& out, double t, const Matrix& Zlags_, const Vector& par_)\n";
    sys_out << "    {\n";
    if (HasPi)
        {
        sys_out << "    const double Pi = M_PI;\n";
        }
    //
    // Constants...
    //
    for (int i = 0; i < nc; ++i)
        {
        sys_out << "    const double " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    // sys_out << endl;
    sys_out << "    // State variables\n";
    GetFromVector2(sys_out,"    const double ",varname_list,"Zlags_","(",",0)",0,";");
    //
    // Parameters...
    //
    if (parname_list.nops() > 0)
        {
        sys_out << "    // Parameters (par_(0) is the period)\n";
        GetFromVector(sys_out,"    const double ",parname_list,"par_","()",1,";");
        // sys_out << endl;
        }
    //
    // The following code assumes that the delays are single parameters,
    // and not mathematical expressions.
    //
    // Expressions...
    //
    if (na > 0)
        sys_out << "    // Expressions\n";
    for (int i = 0; i < na; ++i)
        {
        ex f = exprformula_list[i];
        if (f.has(delay(wild(1),wild(2))))
            PDDEC_ConvertDelaysToZlags(f);
        sys_out << "    const double " << exprname_list[i] << " = " << f << ";" << endl;
        }
    //
    // StateVariables...
    //
    // sys_out << "    vf_ = zeros(" << nv << ",1);" << endl;
    sys_out << "    // Compute the vector field\n";
    for (int i = 0; i < nv; ++i)
        {
        ex f = varvecfield_list[i];
        if (f.has(delay(wild(1),wild(2))))
            PDDEC_ConvertDelaysToZlags(f);
        sys_out << "    out(" << i << ")" << " = " << f << ";" << endl;
        }
    sys_out << endl;
    sys_out << "    return;\n";
    sys_out << "    }\n";

    //
    //  Create the derivative functions to be used in sys_deri(...)
    //
    sys_out << endl;
    sys_out << "//\n";
    sys_out << "// " << Name() << "_jacx\n";
    sys_out << "//\n";
    sys_out << "// This function computes the Jacobian of the vector field\n";
    sys_out << "// with respect to the k-th delayed vector.\n";
    sys_out << "//\n";
    sys_out << endl;
    sys_out << "static inline void " << Name() << "_jacx(Matrix& jac_, double t, int k, const Matrix& Zlags_, const Vector& par_)\n";
    sys_out << "    {\n";
    if (HasPi)
        {
        sys_out << "    const double Pi = M_PI;\n";
        }
    //
    // Constants...
    //
    for (int i = 0; i < nc; ++i)
        {
        sys_out << "    const double " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    // sys_out << endl;

    //
    // Parameters...
    //
    if (parname_list.nops() > 0)
        {
        sys_out << "    // Parameters\n";
        GetFromVector(sys_out,"    const double ",parname_list,"par_","()",1,";");
        sys_out << endl;
        }

    //
    // After the following loop, vf0 will hold a version of the vector field
    // in which all occurrences of delay(expr,lag) will have been converted
    // to use Zlags_(i,j).  The i^th non-delayed state variable will be converted
    // to Zlags_(i,1).
    //
    vector<ex> vf0;
    for (int i = 0; i < nv; ++i)
        {
        // Get the i^th formula, and substitute all expressions.
        ex f = varvecfield_list[i].subs(expreqn_list);
        // Convert the state variables and delay expressions to Zlags_
        PDDEC_ConvertStateToZlags(f);
        vf0.push_back(f);
        }
    PDDEC_PrintJacobians(sys_out,vf0);
    sys_out << "    }\n";
    sys_out << endl;
    
    sys_out << endl;
    sys_out << "//\n";
    sys_out << "// " << Name() << "_jacp\n";
    sys_out << "//\n";
    sys_out << "// This function computes the Jacobian of the vector field\n";
    sys_out << "// with respect to the j-th parameter.\n";
    sys_out << "//\n";
    sys_out << endl;
    sys_out << "static inline void " << Name() << "_jacp(Matrix& jac_, double t, int j, const Matrix& Zlags_, const Vector& par_)\n";
    sys_out << "    {\n";
    if (HasPi)
        {
        sys_out << "    const double Pi = M_PI;\n";
        }
    //
    // Constants...
    //
    for (int i = 0; i < nc; ++i)
        {
        sys_out << "    const double " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    sys_out << endl;

    //
    // Parameters...
    //
    if (parname_list.nops() > 0)
        {
        sys_out << "    // Parameters\n";
        GetFromVector(sys_out,"    const double ",parname_list,"par_","()",1,";");
        sys_out << endl;
        }
    PDDEC_PrintParDerivs(sys_out,vf0);
    sys_out << "    }\n";
    sys_out << endl;

    sys_out << endl;
    sys_out << "//\n";
    sys_out << "// " << Name() << "_jacxp\n";
    sys_out << "//\n";
    sys_out << "// This function computes the Jacobian of the vector field\n";
    sys_out << "// with respect to the k-th delayed state vector and the j-th parameter.\n";
    sys_out << "//\n";
    sys_out << endl;
    sys_out << "static inline void " << Name() << "_jacxp(Matrix& jac_, double t, int k, int j, const Matrix& Zlags_, const Vector& par_)\n";
    sys_out << "    {\n";
    if (HasPi)
        {
        sys_out << "    const double Pi = M_PI;\n";
        }
    //
    // Constants...
    //
    for (int i = 0; i < nc; ++i)
        {
        sys_out << "    const double " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    sys_out << endl;

    //
    // Parameters...
    //
    if (parname_list.nops() > 0)
        {
        sys_out << "    // Parameters\n";
        GetFromVector(sys_out,"    const double ",parname_list,"par_","()",1,";");
        sys_out << endl;
        }
    PDDEC_PrintXandParJacobians(sys_out,vf0);
    sys_out << "    }\n";
    sys_out << endl;


    sys_out << endl;
    sys_out << "//\n";
    sys_out << "// " << Name() << "_hess_times_v\n";
    sys_out << "//\n";
    sys_out << "// This function computes the Hessian of the vector field\n";
    sys_out << "// with respect to the k1-th and k2-th delayed state vectors,\n";
    sys_out << "// then multiplies this by the m-th column of v to obtain a matrix.\n";
    sys_out << "//\n";
    sys_out << endl;
    sys_out << "static inline void " << Name() << "_hess_times_v(Matrix& jac_, double t, int k1, int k2, int m, const Matrix& v_, const Matrix& Zlags_, const Vector& par_)\n";
    sys_out << "    {\n";
    if (HasPi)
        {
        sys_out << "    const double Pi = M_PI;\n";
        }
    //
    // Constants...
    //
    for (int i = 0; i < nc; ++i)
        {
        sys_out << "    const double " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    sys_out << endl;

    //
    // Parameters...
    //
    if (parname_list.nops() > 0)
        {
        sys_out << "    // Parameters\n";
        GetFromVector(sys_out,"    const double ",parname_list,"par_","()",1,";");
        sys_out << endl;
        }
    PDDEC_PrintHessiansTimesV(sys_out,vf0);
    sys_out << "    }\n";
    sys_out << endl;


    //
    //  Create the derivatives function sys_deri(...)
    //

    sys_out << "//" << endl;
    sys_out << "// sys_deri(...)\n";
    sys_out << "//" << endl;
    // Include a list of the lags in the comments.
    sys_out << "// The lags are: {";
    for (unsigned k = 0; k < Delays.size(); ++k)
        {
        sys_out << Delays[k];
        if (k < Delays.size()-1)
            sys_out << ", "; 
        }
    sys_out << "}" << endl;
    sys_out << "//\n";
    sys_out << "// If X(t) is the state vector at time t, then\n";
    sys_out << "//    Zlags_ = [ X(t) ";
    for (unsigned k = 0; k < Delays.size(); ++k)
        {
        sys_out << "X(t-" << Delays[k] << ")";
        if (k < Delays.size()-1)
            sys_out << " "; 
        }
    sys_out << " ]\n";
    sys_out << "//\n";
    sys_out << "// The state vector:\n";
    GetFromVector2(sys_out,"// ",varname_list,"Zlags_","(",",0)",0,";");
    sys_out << "//\n";
    // Function definition starts here.
    sys_out << "//" << endl;
    sys_out << "void sys_deri(Matrix& jac_, double t, const Matrix& Zlags_, const Vector& par_, int nx_, const int* vx_, int np_, const int* vp_, const Matrix& v_)\n";
    sys_out << "    {\n";
    sys_out << "    if (nx_ == 1 & np_ == 0)\n";
    sys_out << "        " << Name() << "_jacx(jac_, t, vx_[0], Zlags_, par_);\n"; 
    sys_out << "    else if (nx_ == 0 & np_ == 1)\n";
    sys_out << "        " << Name() << "_jacp(jac_, t, vp_[0], Zlags_, par_);\n"; 
    sys_out << "    else if (nx_ == 1 & np_ == 1)\n";
    sys_out << "        " << Name() << "_jacxp(jac_,t,vx_[0],vp_[0],Zlags_,par_);\n";
    sys_out << "    else if (nx_ == 2 & np_ == 0)\n";
    sys_out << "        " << Name() << "_hess_times_v(jac_,t,vx_[0],vx_[1],vx_[0],v_,Zlags_,par_);\n";
    sys_out << "    else\n";
    sys_out << "        {\n";
    sys_out << "        std::cerr << \"sys_deri: Requested derivative has not been implemented.\\n\";\n";
    sys_out << "        exit(-1);\n";
    sys_out << "        }\n";
    sys_out << "    return;\n";
    sys_out << "    }\n";

    //
    // Create sys_tau()
    //
    sys_out << endl;
    sys_out << "//" << endl;
    sys_out << "// sys_tau(...) computes the vector of delays.\n";
    sys_out << "//" << endl;
    sys_out << endl;
    sys_out << "void sys_tau(Vector& out, double t, const Vector& par_)\n";
    sys_out << "    {\n";
    if (HasPi)
        {
        sys_out << "    const double Pi = M_PI;\n";
        }
    //
    // Constants...
    //
    for (int i = 0; i < nc; ++i)
        {
        sys_out << "    const double " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    sys_out << "    // par_(0) is the period.\n";
    GetFromVector(sys_out,"    const double ",parname_list,"par_","()",1,";");
    sys_out << endl;
    sys_out << "    out(0) = 0.0;\n";
    int j = 1;
    for (vector<ex>::iterator p = Delays.begin(); p != Delays.end(); ++p)
        {
        if (!is_a<symbol>(*p))
            cerr << "Error: the delay expression " << *p << " is not a single symbol.\n";
        else
            {
            int k = 0;
            lst::const_iterator q;
            for (q = parname_list.begin(); q != parname_list.end(); ++q)
                if (*q == *p)
                    break;
                else
                    ++k;
            if (q == parname_list.end())
                cerr << "Error: the delay expression " << *p << " is not a parameter.\n";
            else
                sys_out << "    out(" << j << ") = " <<  *p << ";\n";
            }
        ++j;
        }
    sys_out << "    }\n";
    sys_out << endl;
    //
    // Create sys_dtau()
    //
    sys_out << "//\n";
    sys_out << "// sys_dtau(...) computes the derivatives of the delays with respect to the parameters.\n";
    sys_out << "//\n";
    sys_out << "// The delays are: ";
    for (unsigned k = 0; k < Delays.size(); ++k)
        {
        sys_out << Delays[k];
        if (k < Delays.size()-1)
            sys_out << ", "; 
        }
    sys_out << endl;
    sys_out << "//\n";
    sys_out << endl;
    sys_out << "void sys_dtau(Vector& out, double t, const Vector& par_, int p_)\n";
    sys_out << "    {\n";
    if (HasPi)
        {
        sys_out << "    const double Pi = M_PI;\n";
        }
    //
    // Constants...
    //
    for (int i = 0; i < nc; ++i)
        {
        sys_out << "    const double " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }

    sys_out << "    // par_(0) is the period.\n";
    GetFromVector(sys_out,"    const double ",parname_list,"par_","()",1,";");
    sys_out << endl;
    sys_out << "    out(0) = 0.0;\n";
    sys_out << "    if (p_ == 0)\n";
    sys_out << "        {\n";
    for (unsigned k = 0; k < Delays.size(); ++k)
        {
        sys_out << "        out(" << k+1 << ") = 0.0;\n";
        }
    sys_out << "        }\n";
    for (int j = 0; j < np; ++j)
        {
        sys_out << "    else if (p_ == " << j+1 << ")\n";
        sys_out << "        {\n";
        sys_out << "        // Derivative wrt " << parname_list[j] << endl;
        for (unsigned k = 0; k < Delays.size(); ++k)
            {
            symbol p = ex_to<symbol>(parname_list[j]);
            ex df = Delays[k].diff(p);
            sys_out << "        out(" << k+1 << ") = " << df << ";\n";
            }
        sys_out << "        }\n";
        }

    sys_out << "    }\n";
    sys_out << endl;
    
    sys_out << "void sys_stpar(Vector& par_)\n";
    sys_out << "    {\n";
    if (HasPi)
        {
        sys_out << "    const double Pi = M_PI;\n";
        }
    sys_out << "    // VFGEN used the DefaultValues of the Parameters.\n";
    sys_out << "    // Change the following values to match your known solution.\n";
    sys_out << "    par_(0) = 1.0;\n";
    for (int j = 0; j < np; ++j)
        {
        sys_out << "    par_(" << j+1 << ") = " << pardefval_list[j] << ";\n";
        }
    sys_out << "    }\n";
    sys_out << endl;
    sys_out << "void sys_stsol(Vector& out, double t)\n";
    sys_out << "    {\n";
    if (HasPi)
        {
        sys_out << "    const double Pi = M_PI;\n";
        }
    sys_out << "    // VFGEN used the DefaultInitialConditions of the StateVariables.\n";
    sys_out << "    // Change the following values to implement your known solution.\n";
    for (int j = 0; j < nv; ++j)
        {
        sys_out << "    out(" << j << ") = " << vardefic_list[j] << ";\n";
        }
    sys_out << "    }\n";
    sys_out << endl;
    sys_out << "}  // extern \"C\"\n";
    sys_out.close();

    return;
    }
