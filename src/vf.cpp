
//
// vf.cpp
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

#include <vector>
#include <fstream>
#include <string>
#include <ginac/ginac.h> 

#include "vf.h"
#include "vf_utils.h"

using namespace std;
using namespace GiNaC;

//
// Symbol Methods
//

Symbol::Symbol()
    {
    }

Symbol::Symbol(string n)
    {
    name = n;
    }

Symbol::Symbol(string n, string descr)
    {
    name = n;
    description = descr;
    }

void Symbol::Name(string n)
    {
    name = n;
    }

string Symbol::Name(void)
    {
    return name;
    }

void Symbol::Description(string descr)
    {
    description = descr;
    }

string Symbol::Description(void)
    {
    return description;
    }


void Symbol::Latex(string l)
    {
    latex = l;
    }

string Symbol::Latex(void)
    {
    return latex;
    }

//
// FormulaSymbol Methods
//

FormulaSymbol::FormulaSymbol() : Symbol()
    {
    }

FormulaSymbol::FormulaSymbol(string name) : Symbol(name)
    {
    }

FormulaSymbol::FormulaSymbol(string name, string descr) : Symbol(name,descr)
    {
    }

void FormulaSymbol::Formula(string f)
    {
    formula = f;
    }

string FormulaSymbol::Formula(void)
    {
    return formula;
    }


//
// Constant Methods
//

Constant::Constant(string name) : Symbol(name)
    {
    value = "";
    }

Constant::Constant(string name, string descr) : Symbol(name,descr)
    {
    value = "";
    }

void Constant::Value(string val)
    {
    value = val;
    }

string Constant::Value(void)
    {
    return value;
    }


//
// Parameter Methods
//

Parameter::Parameter(string name) : Symbol(name)
    {
    defaultvalue = "";
    }

Parameter::Parameter(string name, string descr) : Symbol(name,descr)
    {
    defaultvalue = "";
    }

void Parameter::DefaultValue(string val)
    {
    defaultvalue = val;
    }

string Parameter::DefaultValue(void)
    {
    return defaultvalue;
    }

//
// Expression Methods
//
Expression::Expression(string name) : FormulaSymbol(name)
    {
    }

Expression::Expression(string name, string descr) : FormulaSymbol(name,descr)
    {
    }


//
// StateVariable Methods
//
StateVariable::StateVariable(string name) : FormulaSymbol(name)
    {
    periodicfrom = "";
    periodicto   = "";
    default_ic   = "";
    }

StateVariable::StateVariable(string name, string descr) : FormulaSymbol(name, descr)
    {
    periodicfrom = "";
    periodicto   = "";
    default_ic   = "";
    }


void StateVariable::PeriodicFrom(string pfrom)
    {
    periodicfrom = pfrom;
    }

string StateVariable::PeriodicFrom(void)
    {
    return periodicfrom;
    }


void StateVariable::PeriodicTo(string pto)
    {
    periodicto = pto;
    }

string StateVariable::PeriodicTo(void)
    {
    return periodicto;
    }

bool StateVariable::IsPeriodic(void)
    {
    return periodicfrom != "";
    }

void StateVariable::DefaultInitialCondition(string ic)
    {
    default_ic = ic;
    }

string StateVariable::DefaultInitialCondition(void)
    {
    return default_ic;
    }

void StateVariable::DefaultHistory(string hist)
    {
    default_history = hist;
    }

string StateVariable::DefaultHistory()
    {
    return default_history;
    }

//
// Function Methods
//
Function::Function(string name) : FormulaSymbol(name)
    {
    }

Function::Function(string name, string descr) : FormulaSymbol(name, descr)
    {
    }


//
// VectorField Methods
//

VectorField::VectorField(void)
    {
    IndependentVariable = "t";
    IsAutonomous = true;
    }

VectorField::VectorField(string name, string descr) : Symbol(name,descr)
    {
    IndependentVariable = "t";
    IsAutonomous = true;
    }

VectorField::VectorField(string name, string descr, string indvar) : Symbol(name,descr)
    {
    IndependentVariable = indvar;
    IsAutonomous = true;
    }

// TO DO: Create the destructor; it must delete the memory used
//        by the vector<>s in VectorField.


void VectorField::AddConstant(Constant *p)
    {
    Constants.push_back(p);
    }

void VectorField::AddParameter(Parameter *p)
    {
    Parameters.push_back(p);
    }

void VectorField::AddExpression(Expression *e)
    {
    Expressions.push_back(e);
    }

void VectorField::AddStateVariable(StateVariable *sv)
    {
    StateVariables.push_back(sv);
    }

int VectorField::FindVar(const symbol &var)
    {
    bool found = false;
    unsigned k;
    for (k = 0; k < varname_list.nops(); ++k)
        if (varname_list[k] == var)
            {
            // cerr << "Found the variable " << var << " at varname_list[" << k << "]\n";
            found = true;
            break;
            }
    int vindex;
    if (found == false)
        vindex = -1;
    else
        vindex = k;
    return vindex;
    }

GiNaC::lst VectorField::FindVarsInEx(const GiNaC::ex &e)
    {
    unsigned k;
    GiNaC::lst vlist;
    for (k = 0; k < varname_list.nops(); ++k)
        {
        if (e.has(varname_list[k]))
            vlist.append(varname_list[k]);
        }
    return vlist;
    }

void VectorField::AddFunction(Function *f)
    {
    Functions.push_back(f);
    }


ex VectorField::SubsAllExpressions(const ex& e)
    {
    int na = exprname_list.nops();
    ex s = e;
    for (int k = na-1; k >= 0; --k)
        {
        s = s.subs( ex_to<symbol>(exprname_list[k]) == exprformula_list[k] );
        }
    return s;
    }

int VectorField::FindDelay(ex &del)
    {
    bool found = false;
    unsigned k;
    for (k = 0; k < Delays.size(); ++k)
        {
        // cerr << "FindDelay: Delays[" << k << "] = " << Delays[k] << endl;
        if (Delays[k] == del)
            {
            // cerr << "FindDelay: Found the delay " << del << " at Delays[" << k << "]\n";
            found = true;
            break;
            }
         }
    int dindex;
    if (found == false)
        dindex = -1;
    else
        dindex = k;
    return dindex;
    }


int VectorField::AddDelay(ex &del)
    {
    if (FindDelay(del) != -1)
        return 1;
    Delays.push_back(del);

    ex s = SubsAllExpressions(del);
    ex delayedtime = IndVar - s;
    bool nonconstant = false;
    // cout << "AddDelay: del=" << del << "  s=" << s << "  delayedtime=" << delayedtime << endl;
    if (s.has(IndVar))
        nonconstant = true;
    int nv = varname_list.nops();
    for (int i = 0; i < nv; ++i)
        {
        if (delayedtime.has(varname_list[i]))
            {
            nonconstant = true;
            break;
            }
        }
    // if (nonconstant)
    //     cout << "AddDelay: " << del << " = " << s << " is a nonconstant delay.\n";
    return 0;
    }

//
// CheckForDelay(const ex& f)
//
// Look for expressions of the form delay(delayexpr,del) in f.
// For each of the these, call AddDelay.
// Also set the flags IsDelay and HasNonconstantDelay as appropriate.
//
// This function is called by ProcessSymbols for each Expression and
// StateVariable Formula.
//

void VectorField::CheckForDelay(const ex& f)
    {
    exset occurrences;
    if (f.find(delay(wild(1),wild(2)),occurrences))
        {
        IsDelay = true;
        for (exset::const_iterator iter = occurrences.begin(); iter != occurrences.end(); ++iter)
            {
            ex del = iter->op(1);
            AddDelay(del);
            if (del.has(IndVar))
                HasNonconstantDelay = true;  // time-dependent delay
            for (lst::const_iterator viter = varname_list.begin(); viter != varname_list.end(); ++viter)
                {
                if (del.has(*viter))
                    HasNonconstantDelay = true; // state-dependent delay
                }
            }
        }
    }


int VectorField::ProcessSymbols(void)
    {
    int rval = 0;
    
    IsDelay = false;
    HasNonconstantDelay = false;
    HasPi   = false;

    // Process the constants
    vector<Constant *>::iterator c;
    for (c = Constants.begin(); c != Constants.end(); ++c)
        {
        symbol con( (*c)->Latex() != ""
               ? symbol( (*c)->Name(),(*c)->Latex() )
               : symbol( (*c)->Name() ) );
        conname_list.append(con);
        allsymbols.append(con);
        string val = (*c)->Value();
        try
            {
            ex e(val,allsymbols);
            convalue_list.append(e);
            if (has(e,Pi))
                HasPi = true;
            }
        catch (exception &p)
            {
            // cerr << "VectorField:ProcessSymbols: exception while processing constants\n";
            cerr << "The Value \"" << (*c)->Value() << "\" for the Constant " << (*c)->Name() << " has an error: " << p.what() << endl;
            rval = -1;
            }
        }

    // Process the parameters
    vector<Parameter *>::iterator p;
    for (p = Parameters.begin(); p != Parameters.end(); ++p)
        {
        symbol par( (*p)->Latex() != ""
               ? symbol( (*p)->Name(),(*p)->Latex() )
               : symbol( (*p)->Name() ) );
        // symbol par((*p)->Name());
        parname_list.append(par);
        allsymbols.append(par);
        string defval = (*p)->DefaultValue();
        if (defval == "")
            defval = "0";
        try
            {
            ex e(defval,allsymbols);
            pardefval_list.append(e);
            if (has(e,Pi))
                HasPi = true;
            }
        catch (exception &ep)
            {
            // cerr << "VectorField:ProcessSymbols: exception while processing parameters\n";
            cerr << "The DefaultValue \"" << (*p)->DefaultValue() << "\" for the Parameter " << (*p)->Name() << " has an error: " << ep.what() << endl;
            rval = -1;
            }
        }

    // At this point, allsymbols is a list of the ginac
    // symbols of the constants and parameters.
    // Now add IndVar to the list, because the expressions
    // for the default initial conditions, vector field
    // formulas, expressions, and functions can all be
    // functions of the independent variable.

    IndVar = symbol(IndependentVariable);
    allsymbols.append(IndVar);

    // Process the state variable names and default ICs
    //   (but not the formulas for the vector field)
    int nv = 0;
    vector<StateVariable *>::iterator sv;
    for (sv = StateVariables.begin(); sv != StateVariables.end(); ++sv, ++nv)
        {
        symbol var( (*sv)->Latex() != ""
               ? symbol( (*sv)->Name(),(*sv)->Latex() )
               : symbol( (*sv)->Name() ) );
        varname_list.append(var);
        string defic = (*sv)->DefaultInitialCondition();
        if (defic == "")
            defic = "0";
        try
            {
            ex e(defic,allsymbols);
            vardefic_list.append(e);
            if (has(e,Pi))
                HasPi = true;
            }
        catch (exception &p)
            {
            // cerr << "VectorField::ProcessSymbols: exception while processing DefaultInitialCondition " << (*sv)->DefaultInitialCondition() << endl;
            cerr << "The DefaultInitialCondition \"" << (*sv)->DefaultInitialCondition() << "\" for the StateVariable " << (*sv)->Name() << " has an error: " << p.what() << endl;
            rval = -1;
            }
        string defhist = (*sv)->DefaultHistory();
        if (defhist == "")
            defhist = "0";  // The default DefaultHistory.  Shouldn't this be DefaultInitialCondition?
                            // (And only be 0 if DefaultInitialCondition was also not given.)
        try
            {
            ex e(defhist,allsymbols);
            vardefhist_list.append(e);
            if (has(e,Pi))
                HasPi = true;
            }
        catch (exception &p)
            {
            // cerr << "VectorField::ProcessSymbols: exception while processing DefaultHistory " << (*sv)->DefaultHistory() << endl;
            cerr << "The DefaultHistory \"" << (*sv)->DefaultHistory() << "\" for the StateVariable " << (*sv)->Name() << " has an error: " << p.what() << endl;
            rval = -1;
            }

        }
    //
    // Now add the state variable symbols to allsymbols.
    // We didn't do this in the above loop because we don't want
    // to allow the default initial condition of one state variable
    // to be a function of another state variable. (In other words,
    // the default initial conditions can be functions of only the
    // constants, parameters, and the independent variable [for
    // delay equations].)
    //
    for (int j = 0; j < nv; ++j)
        {
        allsymbols.append(varname_list.op(j));
        }

    // Process the expressions
    vector<Expression *>::iterator e;
    for (e = Expressions.begin(); e != Expressions.end(); ++e)
        {
        symbol auxe( (*e)->Latex() != ""
               ? symbol( (*e)->Name(),(*e)->Latex() )
               : symbol( (*e)->Name() ) );
        exprname_list.append(auxe);
        allsymbols.append(auxe);
        try
            {
            ex f((*e)->Formula(),allsymbols);
            exprformula_list.append(f);
            expreqn_list.append( auxe==f );
            CheckForDelay(f);
            if (has(f,Pi))
                HasPi = true;
            }
        catch (exception &p)
            {
            // cerr << "VectorField:ProcessSymbols: exception while processing expressions\n";
            cerr << "The Formula \"" << (*e)->Formula() << "\" for the Expression " << (*e)->Name() << " has an error: " << p.what() << endl;
            rval = -1;
            }
        }
    /*
    cerr << "*** Delays = { ";
    vector<ex>::iterator pos;
    for (pos = Delays.begin(); pos != Delays.end(); ++pos)
        {
        cerr << *pos << ' ';
        }
    cerr << "}" << endl;
    */

    //
    // At this point, allsymbols holds the ginac symbols for
    // the independent variable, the constants, the parameters,
    // the state variables, and the auxiliary expressions.
    // No more symbols will be added to allsymbols.
    //

    // Process the vector field formulas
    for (sv = StateVariables.begin(); sv != StateVariables.end(); ++sv)
        {
        try
            {
            ex f((*sv)->Formula(),allsymbols);
            varvecfield_list.append(f);
            CheckForDelay(f);
            if (has(f,Pi))
                HasPi = true;
            if (has(f,IndVar))
                IsAutonomous = false;
            }
        catch (exception &p)
            {
            // cerr << "VectorField:ProcessSymbols: exception while processing vector field formulas\n";
            cerr << "The Formula \"" << (*sv)->Formula() << "\" for the StateVariable " << (*sv)->Name() << " has an error: " << p.what() << endl;
            rval = -1;
            }
        }

    // Functions
    vector<Function *>::iterator f;
    for (f = Functions.begin(); f != Functions.end(); ++f)
        {
        symbol funcname((*f)->Name());
        funcname_list.append(funcname);
        try
            {
            ex funcexpr((*f)->Formula(),allsymbols);
            funcformula_list.append(funcexpr);
            if (has(funcexpr,Pi))
                HasPi = true;
            }
        catch (exception &p)
            {
            // cerr << "VectorField:ProcessSymbols: exception while processing function formulas\n";
            cerr << "The Formula \"" << (*f)->Formula() << "\" for the  Function " << (*f)->Name() << " has an error: " << p.what() << endl;
            rval = -1;
            }
        }
    return rval;
    }

void VectorField::Print(void)
    {
    cout << "Name:       " << Name();
    cout << "   Independent Variable: " << IndependentVariable << endl;
    cout << "Constants:  " << conname_list;
    cout << "  Values: " << convalue_list << endl;

    cout << "Parameters: " << parname_list;
    cout << "  Default values: " << pardefval_list << endl;

    cout << "Variables:  " << varname_list << endl;
    cout << "  DefaultICs:     " << vardefic_list << endl;
    cout << "  DefaultHistory: " << vardefhist_list << endl;

    cout << "Expressions: " << endl;

    for (unsigned i = 0; i < exprname_list.nops(); ++i)
        {
        cout << "   " << exprname_list[i] << "=" << exprformula_list[i] << endl;
        }

    cout << "Expressions (equation list): " ;
    cout << expreqn_list << endl;

    cout << "Vector field: " << endl;
    for (unsigned i = 0; i < varvecfield_list.nops(); ++i)
        {
        cout << "   " << varvecfield_list[i] << endl;
        }
    cout << "Functions: " << endl;
    for (unsigned i = 0; i < funcname_list.nops(); ++i)
        {
        cout << "   " << funcname_list[i] << "=" << funcformula_list[i] << endl;
        }
    string tf[2] = {"false","true"};
    cout << "IsDelay: " << tf[IsDelay] << endl;
    cout << "HasNonconstantDelay: " << tf[HasNonconstantDelay] << endl;
    cout << "IsAutonomous: " << tf[IsAutonomous] << endl;
    cout << "HasPi: " << tf[HasPi] << endl;
    }

