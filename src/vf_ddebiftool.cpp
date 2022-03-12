
//
//  vf_ddebiftool.cpp
//
//  This file defines the VectorField::PrintDDEBIFTOOL method.
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



void VectorField::DDEBT_PrintParDerivs(ofstream &dout, const vector<ex> &vf0)
{
    // int nc = conname_list.nops();
    int nv = varname_list.nops();
    int np = parname_list.nops();
    // int na = exprname_list.nops();
    // int nf = funcname_list.nops();

    for (int j = 0; j < np; ++j) {
        if (j == 0) {
            dout << "        if np_ == 1\n";
        }
        else {
            dout << "        elseif np_ == " << j+1 << endl;
        }
        dout << "            % Derivative wrt " << parname_list[j] << endl;
        for (int i = 0; i < nv; ++i) {
            symbol p = ex_to<symbol>(parname_list[j]);
            ex df = vf0[i].diff(p);
            if (df != 0) {
                dout << "            jac_(" << (i+1) << ")" << " = " << df << ";" << endl;
            }
        }
    }
    dout << "        end\n";
}


void VectorField::DDEBT_PrintJacobians(ofstream &dout, const vector<ex> &vf0)
{
    // int nc = conname_list.nops();
    int nv = varname_list.nops();
    // int np = parname_list.nops();
    // int na = exprname_list.nops();
    // int nf = funcname_list.nops();

    int nd = Delays.size();
    for (int k = 0; k < nd+1; ++k) {
        if (k == 0) {
            dout << "        if nx_ == 0\n";
            dout << "            % Derivatives wrt the state variables\n";
        }
        else {
            dout << "        elseif nx_ == " << k << endl;
            dout << "            % Derivatives wrt state variables with delay " << Delays[k-1] << endl;
        }
        for (int i = 0; i < nv; ++i) {
            ex f = vf0[i];
            for (int j = 0; j < nv; ++j) {
                symbol vtmp_("vtmp_");
                ex fj = f.subs(Zlags_(j+1,k+1) == vtmp_);
                ex df = fj.diff(vtmp_);
                df = df.subs(vtmp_ == Zlags_(j+1,k+1));
                if (df != 0) {
                    dout << "            jac_(" << (i+1) << "," << (j+1) << ")" << " = " << df << ";" << endl;
                }
            }
        }
    }
    dout << "        end\n";
}


void VectorField::DDEBT_PrintXandParJacobians(ofstream &dout, const vector<ex> &vf0)
{
    // int nc = conname_list.nops();
    int nv = varname_list.nops();
    int np = parname_list.nops();
    // int na = exprname_list.nops();
    // int nf = funcname_list.nops();

    int nd = Delays.size();
    for (int k = 0; k < nd+1; ++k) {
        if (k == 0) {
            dout << "        if nx_ == 0\n";
            dout << "            % Derivatives wrt the state variables\n";
        }
        else {
            dout << "        elseif nx_ == " << k << endl;
            dout << "            % Derivatives wrt state variables with delay " << Delays[k-1] << endl;
        }
        for (int m = 0; m < np; ++m) {
            if (m == 0) {
                dout << "            if np_ == " << m+1 << "\n";
            }
            else {
                dout << "            elseif np_ == " << m+1 << "\n";
            }
            dout << "                % Derivative wrt " << parname_list[m] << "\n";
            for (int i = 0; i < nv; ++i) {
                ex f = vf0[i];
                for (int j = 0; j < nv; ++j) {
                    symbol vtmp_("vtmp_");
                    ex fj = f.subs(Zlags_(j+1,k+1) == vtmp_);
                    ex df = fj.diff(vtmp_);
                    df = df.subs(vtmp_ == Zlags_(j+1,k+1));
                    symbol p = ex_to<symbol>(parname_list[m]);
                    df = df.diff(p);
                    if (df != 0) {
                        dout << "                jac_(" << (i+1) << "," << (j+1) << ")" << " = " << df << ";" << endl;
                    }
                }
            }
        }
        dout << "            end    % if np_ == ...\n";
    }
    dout << "        end    % if nx_ == ...\n";
}


ex second_deriv(const ex &f, int lag1, int var1, int lag2, int var2)
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


void VectorField::DDEBT_PrintHessiansTimesV(ofstream &dout, const vector<ex> &vf0)
{
    // int nc = conname_list.nops();
    int nv = varname_list.nops();
    // int np = parname_list.nops();
    // int na = exprname_list.nops();
    // int nf = funcname_list.nops();

    int nd = Delays.size();
    for (int k1 = 0; k1 < nd+1; ++k1) {
        if (k1 == 0) {
            dout << "        if nx_(1) == 0\n";
            dout << "            % Derivatives wrt the state variables\n";
        }
        else {
            dout << "        elseif nx_(1) == " << k1 << endl;
            dout << "            % Derivatives wrt state variables with delay " << Delays[k1-1] << endl;
        }

        for (int k2 = 0; k2 < nd+1; ++k2) {
            if (k2 == 0) {
                dout << "            if nx_(2) == 0\n";
                dout << "                % Derivatives wrt the state variables\n";
            }
            else {
                dout << "            elseif nx_(2) == " << k2 << endl;
                dout << "                % Derivatives wrt state variables with delay " << Delays[k2-1] << endl;
            }

            for (int i = 0; i < nv; ++i) {
                ex f = vf0[i];
                for (int j = 0; j < nv; ++j) {
                    ostringstream os;
                    for (int h = 0; h < nv; ++h) {
                        ex d2f = second_deriv(f,k1+1,j+1,k2+1,h+1);
                        if (d2f != 0) {
                            if (os.str() != "") {
                                os << " + ";
                            }
                            os << "(" << d2f << ")*v_(" << h+1 << ")";
                        }
                    }
                    if (os.str() != "") {
                        dout << "                jac_(" << (i+1) << "," << (j+1) << ")" << " = " << os.str() << ";\n";
                    }
                }
            }
        }
        dout << "            end % if nx_(2) == ...\n";
    }
    dout << "        end  % if nx_(1) == ...\n";
}


//
// PrintDDEBIFTOOL -- The DDEBIFTOOL Matlab function code generator.
//

void VectorField::PrintDDEBIFTOOL(map<string,string> options)
{
    int nc = conname_list.nops();
    int nv = varname_list.nops();
    // int np = parname_list.nops();
    int na = exprname_list.nops();
    // int nf = funcname_list.nops();

    // cerr << "Delays:\n";
    // for (vector<ex>::iterator di = Delays.begin(); di != Delays.end(); ++di)
    //     cerr << "   " << *di << endl;

    //
    //  Create the initialization file sys_init.m
    //
    // string si_filename = Name() + "_sys_init.m";
    string si_filename = "sys_init.m";
    ofstream si_out;
    si_out.open(si_filename.c_str());
    si_out << "%" << endl;
    si_out << "% " << si_filename << endl;
    si_out << "%" << endl;
    si_out << "% DDEBIFTOOL MATLAB sys_init() function for: " << Name() << endl;
    si_out << "%" << endl;
    PrintVFGENComment(si_out,"% ");
    si_out << "%" << endl;
    si_out << "%" << endl;
    si_out << "function [name,dim] = sys_init()\n";
    si_out << "    name = '" << Name() << "';\n";
    si_out << "    dim = " << nv << ";\n";
    if (options["path"] != "") {
        si_out << "    path(path,'" << options["path"] << "');\n";
    }
    else {
        si_out << "    % Add the DDE-BIFTOOL path here, if necessary.\n";
        si_out << "    % path(path,'" << options["path"] << "');\n";
    }
    si_out << "    return;\n";
    si_out.close();

    //
    //  Create the vector field function sys_rhs.m
    //
    // string vf_filename = Name()+"_sys_rhs.m";
    string vf_filename = "sys_rhs.m";
    ofstream fout;
    fout.open(vf_filename.c_str());
    fout << left;

    fout << "%" << endl;
    fout << "% " << vf_filename << endl;
    fout << "%" << endl;
    fout << "% DDEBIFTOOL MATLAB vector field function for: " << Name() << endl;
    fout << "%" << endl;
    PrintVFGENComment(fout,"% ");
    fout << "%" << endl;
    fout << "%" << endl;
    // Include a list of the lags in the comments.
    fout << "% The lags are: {";
    for (unsigned k = 0; k < Delays.size(); ++k) {
        fout << Delays[k];
        if (k < Delays.size()-1)
            fout << ", "; 
    }
    fout << "}" << endl;
    fout << "%\n";
    fout << "% If X(t) is the state vector at time t, then\n";
    fout << "%    Zlags_ = [ X(t) ";
    for (unsigned k = 0; k < Delays.size(); ++k) {
        fout << "X(t-" << Delays[k] << ")";
        if (k < Delays.size()-1) {
            fout << " "; 
        }
    }
    fout << " ]\n";
    // Function definition starts here.
    fout << "%" << endl;
    fout << "function vf_ = sys_rhs(Zlags_, par_)";

    if (HasPi) {
        fout << "    Pi = pi;\n";
    }
    //
    // Constants...
    //
    for (int i = 0; i < nc; ++i) {
        fout << "    " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
    }
    fout << endl;
    fout << "    % State variables\n";
    GetFromVector2(fout, "    ", varname_list, "=", "Zlags_", "(", ",1)", 1, ";");
    //
    // Parameters...
    //
    if (parname_list.nops() > 0) {
        fout << "    % Parameters\n";
        GetFromVector(fout, "    ", parname_list, "=", "par_", "()", 1, ";");
        fout << endl;
    }
    //
    // The following code assumes that the delays are single parameters,
    // and not mathematical expressions.
    //
    // Expressions...
    //
    for (int i = 0; i < na; ++i) {
        ex f = exprformula_list[i];
        if (f.has(delay(wild(1),wild(2)))) {
            ConvertDelaysToZlags(f, 1, 2);
        }
        fout << "    " << exprname_list[i] << " = " << f << ";" << endl;
    }
    //
    // StateVariables...
    //
    fout << "    vf_ = zeros(" << nv << ",1);" << endl;
    for (int i = 0; i < nv; ++i) {
        ex f = varvecfield_list[i];
        if (f.has(delay(wild(1),wild(2)))) {
            ConvertDelaysToZlags(f, 1, 2);
        }
        fout << "    vf_(" << (i+1) << ")" << " = " << f << ";" << endl;
    }
    fout << endl;
    fout << "    return\n";
    fout.close();

    //
    //  Create the derivatives file sys_deri.m
    //
    // string d_filename = Name()+"_sys_deri.m";
    string d_filename = "sys_deri.m";
    ofstream dout;
    dout.open(d_filename.c_str());
    dout << left;

    dout << "%" << endl;
    dout << "% " << d_filename << endl;
    dout << "%" << endl;
    dout << "% DDE-BIFTOOL MATLAB vector field function for: " << Name() << endl;
    dout << "%" << endl;
    PrintVFGENComment(dout,"% ");
    dout << "%" << endl;
    dout << "%" << endl;
    // Include a list of the lags in the comments.
    dout << "% The lags are: {";
    for (unsigned k = 0; k < Delays.size(); ++k) {
        dout << Delays[k];
        if (k < Delays.size()-1) {
            dout << ", "; 
        }
    }
    dout << "}" << endl;
    dout << "%\n";
    dout << "% If X(t) is the state vector at time t, then\n";
    dout << "%    Zlags_ = [ X(t) ";
    for (unsigned k = 0; k < Delays.size(); ++k) {
        dout << "X(t-" << Delays[k] << ")";
        if (k < Delays.size()-1) {
            dout << " "; 
        }
    }
    dout << " ]\n";
    dout << "%\n";
    dout << "% The state vector:\n";
    GetFromVector2(dout, "% ", varname_list, "=", "Zlags_", "(", ",1)", 1, ";");
    dout << "%\n";
    dout << "% (In the comments below, \"wrt\" means \"with respect to\".)\n";
    // Function definition starts here.
    dout << "%" << endl;
    dout << "function jac_ = sys_deri(Zlags_, par_, nx_, np_, v_)";

    if (HasPi) {
        dout << "    Pi = pi;\n";
    }
    //
    // Constants...
    //
    for (int i = 0; i < nc; ++i) {
        dout << "    " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
    }
    dout << endl;

    //
    // Parameters...
    //
    if (parname_list.nops() > 0) {
        dout << "    % Parameters\n";
        GetFromVector(dout, "    ", parname_list, "=", "par_", "()", 1, ";");
        dout << endl;
    }
    //
    // We do not print the Expressions, since the vector field formulas will
    // have all their Expressions replaced by their formulas.
    //

    //
    // After the following loop, vf0 will hold a version of the vector field
    // in which all occurrences of delay(expr,lag) will have been converted
    // to use Zlags_(i,j).  The i^th non-delayed state variable will be converted
    // to Zlags_(i,1).
    //
    vector<ex> vf0;
    for (int i = 0; i < nv; ++i) {
        // Get the i^th formula, and substitute all expressions.
        ex f = varvecfield_list[i].subs(expreqn_list);
        // Convert the state variables and delay expressions to Zlags_
        ConvertStateToZlags(f, 1);
        vf0.push_back(f);
    }

    dout << "    if length(nx_) == 1 && length(np_) == 0 && isempty(v_)\n";
    dout << "        jac_ = zeros(" << nv << "," << nv << ");\n";
    DDEBT_PrintJacobians(dout, vf0);
    dout << "    elseif length(nx_) == 0 && length(np_) == 1 && isempty(v_)\n";
    dout << "        jac_ = zeros(" << nv << ",1);\n";
    DDEBT_PrintParDerivs(dout, vf0);
    dout << "    elseif length(nx_) == 1 && length(np_) == 1 && isempty(v_)\n";
    dout << "        % mixed state variable and parameter derivatives\n";
    dout << "        jac_ = zeros(" << nv << "," << nv << ");\n";
    DDEBT_PrintXandParJacobians(dout,vf0);
    dout << "    elseif length(nx_) == 2 && length(np_) == 0 && ~isempty(v_)\n";
    dout << "        jac_ = zeros(" << nv << "," << nv << ");\n";
    DDEBT_PrintHessiansTimesV(dout,vf0);
    dout << "    else\n";
    dout << "        error('sys_deri: Requested derivative has not been implemented.')\n";
    dout << "    end\n";
    dout << "    return\n";

    dout.close();

    //
    // Create sys_tau.m
    //
    // (The code currently assumes that all the delays are constants.)
    //
    // string t_filename = Name()+"_sys_tau.m";
    string t_filename = "sys_tau.m";
    ofstream tout;
    tout.open(t_filename.c_str());
    tout << left;

    tout << "%" << endl;
    tout << "% " << t_filename << endl;
    tout << "%" << endl;
    tout << "% DDE-BIFTOOL MATLAB sys_tau function for: " << Name() << endl;
    tout << "%" << endl;
    PrintVFGENComment(tout,"% ");
    tout << "%" << endl;
    tout << endl;
    tout << "function tau = sys_tau()\n";
    tout << endl;
    tout << "tau = [ ";
    for (vector<ex>::iterator p = Delays.begin(); p != Delays.end(); ++p) {
        if (!is_a<symbol>(*p)) {
            cerr << "Error: the delay expression " << *p << " is not a single symbol.\n";
        }
        else {
            int k = 0;
            lst::const_iterator q;
            for (q = parname_list.begin(); q != parname_list.end(); ++q) {
                if (*q == *p) {
                    break;
                }
                else {
                    ++k;
                }
            }
            if (q == parname_list.end()) {
                cerr << "Error: the delay expression " << *p << " is not a parameter.\n";
            }
            else {
                tout << k+1 << " ";
            }
        }
    }
    tout << "];\n";
    tout <<"return\n";
    tout.close();

    return;
}
