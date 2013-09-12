
//
//  vf_delay2ode.cpp
//
//  This file defines the VectorField::PrintDelay2ODE method.
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
// Delay2ODE_ConvertExprToDefHist(ex& f)
//
// Replace each occurrence of a StateVariable in f with its DefaultHistory.
//

void VectorField::Delay2ODE_ConvertExprToDefHist(ex& f)
    {
    for (unsigned k = 0; k < varname_list.nops(); ++k)
        {
        if (f.has(varname_list[k]))
            f = f.subs(varname_list[k] == vardefhist_list[k]);
        }
    }

//
// Delay2ODE_ConvertAndExtend(ex& f, int N, int p)
//
// This function is only used in PrintDelay2ODE.
//

void VectorField::Delay2ODE_ConvertAndExtend(ex& f, int N, int p)
    {
    int k = 0;
    exset dlist;
    f.find(delay(wild(1),wild(2)),dlist);
    // dlist is now a ginac lst of expressions of the form delay(delayexpr,del)
    for (exset::const_iterator diter = dlist.begin(); diter != dlist.end(); ++diter)
        {
        // diter points to a ginac expression of the form delay(delayexpr,del).
        ex delayfunc = *diter;
        ostringstream symbol_name;
        symbol_name << "lag" << ++k;
        string sym = symbol_name.str();
        ex lag = delayfunc.op(1);
        ex hist = delayfunc.op(0);
        Delay2ODE_ConvertExprToDefHist(hist);
        string delayed_var_name;
        ostringstream os_N_over_delta;
        os_N_over_delta << "(" << N << "/(" << lag << "))";
        string N_over_delta = os_N_over_delta.str();
        for (int k = 0; k < N; ++k)
            {
            if (p == 1)
                {
                string vstr1, prev_vstr1;
                ostringstream os;
                os << k+1;
                vstr1 = sym + "_1_" + os.str();
                os.str("");
                os << k;
                prev_vstr1 = sym + "_1_" + os.str();
                //
                StateVariable *var = new StateVariable(vstr1);
                ostringstream os_varformula;
                os_varformula << N_over_delta << "*(";
                if (k == 0)
                    os_varformula << delayfunc.op(0);
                else
                    os_varformula << prev_vstr1;
                os_varformula << " - " << vstr1 << ")";
                var->Formula(os_varformula.str());
                ostringstream os_vardefic;
                os_vardefic << hist.subs(IndVar==-(k+1)*lag/N);
                var->DefaultInitialCondition(os_vardefic.str());
                AddStateVariable(var);
                if (k == N-1)
                    delayed_var_name = var->Name();
                }
            else if (p == 2)
                {
                string vstr1,vstr2, prev_vstr1;
                ostringstream os;
                os << k+1;
                vstr1 = sym + "_1_" + os.str();
                vstr2 = sym + "_2_" + os.str();
                os.str("");
                os << k;
                prev_vstr1 = sym + "_1_" + os.str();
                //
                StateVariable *var = new StateVariable(vstr1);
                var->Formula(vstr2);
                ostringstream os_vardefic;
                os_vardefic << hist.subs(IndVar==-(k+1)*lag/N);
                var->DefaultInitialCondition(os_vardefic.str());
                AddStateVariable(var);
                if (k == N-1)
                    delayed_var_name = var->Name();
                //
                var = new StateVariable(vstr2);
                ostringstream os_varformula;
                os_varformula << "2*" << N_over_delta << "*";
                os_varformula << "( -" << vstr2 << " + " << N_over_delta << "*";
                os_varformula << "(";
                if (k == 0)
                    os_varformula << delayfunc.op(0);
                else
                    os_varformula << prev_vstr1;
                os_varformula << " - " << vstr1 << "))";
                var->Formula(os_varformula.str());
                ex icderiv = hist.diff(IndVar);
                ostringstream os_vardefic_deriv;
                os_vardefic_deriv << icderiv.subs(IndVar==-(k+1)*lag/N);
                var->DefaultInitialCondition(os_vardefic_deriv.str());
                AddStateVariable(var);
                }
            else // p == 3
                {
                string vstr1,vstr2,vstr3, prev_vstr1;
                ostringstream os;
                os << k+1;
                vstr1 = sym + "_1_" + os.str();
                vstr2 = sym + "_2_" + os.str();
                vstr3 = sym + "_3_" + os.str();
                os.str("");
                os << k;
                prev_vstr1 = sym + "_1_" + os.str();
                //
                StateVariable *var = new StateVariable(vstr1);
                var->Formula(vstr2);
                ostringstream os_vardefic;
                os_vardefic << hist.subs(IndVar==-(k+1)*lag/N);
                var->DefaultInitialCondition(os_vardefic.str());
                AddStateVariable(var);
                if (k == N-1)
                    delayed_var_name = var->Name();
                //
                var = new StateVariable(vstr2);
                var->Formula(vstr3);
                ex icderiv = hist.diff(IndVar);
                ostringstream os_vardefic_deriv;
                os_vardefic_deriv << icderiv.subs(IndVar==-(k+1)*lag/N);
                var->DefaultInitialCondition(os_vardefic_deriv.str());
                AddStateVariable(var);
                //
                var = new StateVariable(vstr3);
                ostringstream os_varformula;
                os_varformula << N_over_delta << "*";
                os_varformula << "(-3*" << vstr3;
                os_varformula << " - 6*" << N_over_delta << "*";
                os_varformula << "(" << vstr2 << " - " << N_over_delta << "*";
                os_varformula << "(";
                if (k == 0)
                    os_varformula << delayfunc.op(0);
                else
                    os_varformula << prev_vstr1;
                os_varformula << " - " << vstr1 << ")))";
                var->Formula(os_varformula.str());
                ex icderiv2 = icderiv.diff(IndVar);
                ostringstream os_vardefic_deriv2;
                os_vardefic_deriv2 << icderiv2.subs(IndVar==-(k+1)*lag/N);
                var->DefaultInitialCondition(os_vardefic_deriv2.str());
                AddStateVariable(var);
                }                    
            } // end k loop
        symbol s(delayed_var_name);
        f = f.subs(delayfunc == s);
        } // end dlist loop
    }


//
// PrintDelay2ODE -- The Delay2ODE Code Generator.
//

void VectorField::PrintDelay2ODE(map<string,string> options)
    {
    //
    // Note that this function MODIFIES THE OBJECT!!!
    // It then calls PrintXML to output the extended vector field to the
    // standard output in XML format.
    //

    int ne = exprname_list.nops();
    int nv = varname_list.nops();
    int N, p;

    if (options.find("N") == options.end())
        N = 10;
    else
        {
        N = string_to_int(options["N"]);
        if (N < 1)
            {
            cerr << "N must be at least 1.\n";
            return;
            }
        }
    if (options.find("p") == options.end())
        p = 1;
    else
        {
        p = string_to_int(options["p"]);
        if (p < 1 || p > 3)
            {
            cerr << "Only p=1, p=2 or p=3 are supported.\n";
            return;
            }
        }

    for (int i = 0; i < ne; ++i)
        {
        ex f = exprformula_list[i];
        if (f.has(delay(wild(1),wild(2))))
            {
            Delay2ODE_ConvertAndExtend(f,N,p);
            ostringstream os_newformula;
            os_newformula << f;
            Expressions[i]->Formula(os_newformula.str());
            }
        }

    for (int k = 0; k < nv; ++k)
        {
        ex f = varvecfield_list[k];
        if (f.has(delay(wild(1),wild(2))))
            {
            Delay2ODE_ConvertAndExtend(f,N,p);
            ostringstream os_newformula;
            os_newformula << f;
            StateVariables[k]->Formula(os_newformula.str());
            }
        }

    for (unsigned k = 0; k < vardefic_list.nops(); ++k)
        {
        if (has(vardefic_list[k],IndVar))
            {
            vardefic_list[k] = vardefic_list[k].subs(IndVar==0);
            ostringstream os;
            os << vardefic_list[k];
            StateVariables[k]->DefaultInitialCondition(os.str());
            }
        }
    Name(Name()+"_2ode");
    PrintXML("delay2ode");
    }
