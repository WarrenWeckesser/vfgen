
//
// vf.h
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


#ifndef VF_H_INCLUDED_
#define VF_H_INCLUDED_

#include <string>
#include <vector>
#include <map>
#include <ginac/ginac.h>

#include "ginac_declare_funcs.h"

//
// Symbol
//     FormulaSymbol
//         Expression
//         StateVariable
//         Function
//     Constant
//     Parameter
//     VectorField
//


//
// Symbol is the base class for the objects used to define a vector field.
// The data associated with a Symbol is the name of the symbol, a
// description, and a latex representation.
//

class Symbol
{
    private:

    std::string name;
    std::string description;
    std::string latex;

    public:

    // Constructors
    Symbol();
    Symbol(std::string name);
    Symbol(std::string name, std::string descr);

    // Get/Set methods
    void Name(std::string n);
    std::string Name(void);
    void Description(std::string descr);
    std::string Description(void);
    void Latex(std::string l);
    std::string Latex(void);
};


class FormulaSymbol : public Symbol
{
    private:

    std::string formula;

    public:

    // Constructors
    FormulaSymbol();
    FormulaSymbol(std::string name);
    FormulaSymbol(std::string name, std::string descr);

    // Get/Set methods
    void Formula(std::string f);
    std::string Formula(void);
};


class Constant : public Symbol
{
    private:

    std::string value;

    public:

    // Constructors
    Constant(std::string name);
    Constant(std::string name, std::string descr);

    // Get/Set methods
    void Value(std::string val);
    std::string Value(void);
};


class Parameter : public Symbol
{
    private:

    std::string defaultvalue;

    public:

    // Constructors
    Parameter(std::string name);
    Parameter(std::string name, std::string descr);

    // Get/Set methods
    void DefaultValue(std::string val);
    std::string DefaultValue(void);
};


class Expression : public FormulaSymbol
{
    public:

    // Constructors
    Expression(std::string name);
    Expression(std::string name, std::string descr);

};


class StateVariable : public FormulaSymbol
{
    private:

    std::string periodicfrom;
    std::string periodicto;
    std::string default_ic;
    std::string default_history;

    public:

    // Constructors
    StateVariable(std::string name);
    StateVariable(std::string name, std::string descr);

    // Get/Set methods
    void PeriodicFrom(std::string pfrom);
    std::string PeriodicFrom(void);
    void PeriodicTo(std::string pto);
    std::string PeriodicTo(void);
    bool IsPeriodic();
    void DefaultInitialCondition(std::string ic);
    std::string DefaultInitialCondition(void);
    void DefaultHistory(std::string hist);
    std::string DefaultHistory(void);
};


class Function : public FormulaSymbol
{
    public:

    // Constructors
    Function(std::string name);
    Function(std::string name, std::string descr);
};


class VectorField : public Symbol
{
    private:

    // There is no implementation of the copy constructor.
    VectorField(const VectorField& vf);
    // There is no implementation of the assignment operator.
    VectorField operator=(const VectorField& vf);

    protected:

    GiNaC::lst conname_list;
    GiNaC::lst convalue_list;
    GiNaC::lst parname_list;
    GiNaC::lst pardefval_list;
    GiNaC::lst exprname_list;
    GiNaC::lst exprformula_list;
    GiNaC::lst expreqn_list;
    GiNaC::lst varname_list;
    GiNaC::lst varvecfield_list;
    GiNaC::lst vardefic_list;
    GiNaC::lst vardefhist_list;
    GiNaC::lst funcname_list;
    GiNaC::lst funcformula_list;

    GiNaC::lst allsymbols;
    GiNaC::symbol IndVar;

    // Everything is public, for now.
    public:

    std::string IndependentVariable;
    std::vector<Constant *>      Constants;
    std::vector<Parameter *>     Parameters;
    std::vector<Expression *>    Expressions;
    std::vector<StateVariable *> StateVariables;
    std::vector<Function *>      Functions;
    bool IsDelay;
    bool HasNonconstantDelay;
    bool HasPi;
    bool IsAutonomous;

    std::vector<GiNaC::ex> Delays;

    // Constructors
    VectorField();
    VectorField(std::string name, std::string descr);
    VectorField(std::string name, std::string descr, std::string indvar);

    void AddConstant(Constant *c);
    void AddParameter(Parameter *p);
    void AddExpression(Expression *e);
    void AddStateVariable(StateVariable *sv);
    int  FindVar(const GiNaC::symbol&);
    GiNaC::lst FindVarsInEx(const GiNaC::ex &e);
    GiNaC::ex SubsAllExpressions(const GiNaC::ex &e);

    void AddFunction(Function *f);

    int FindDelay(GiNaC::ex&);
    int AddDelay(GiNaC::ex&);

    void ConvertDelaysToZlags(GiNaC::ex& f, int i_offset, int j_offset);
    void ConvertStateToZlags(GiNaC::ex& f, int offset);
    void convert_delay_to_lagvalue(GiNaC::ex& f, GiNaC::lst& lags);

    void PrintXML(std::string cmdstr);
    int  ReadXML(std::string xmlfilename);

    void CheckForDelay(const GiNaC::ex& f);
    int  ProcessSymbols(void);

    void Print(void);
    void PrintGSL(std::map<std::string,std::string> options);
    void PrintBoostOdeint(std::map<std::string,std::string> options);
    void PrintCVODE(std::map<std::string,std::string> options);
    void PrintCVODE7(std::map<std::string,std::string> options);
    void PrintScilab(std::map<std::string,std::string> options);
    void PrintADOLC(/* std::map<std::string,std::string> options */);
    void PrintAUTO(std::map<std::string,std::string> options);
    void PrintDSTool(void);
    void PrintEVF(std::map<std::string,std::string> options);
    void PrintJavaMath(std::map<std::string,std::string> options);
    void PrintJavascript(std::map<std::string,std::string> options);

    void PrintJuliaFuncStart(std::ofstream &fout);
    void PrintJulia(std::map<std::string,std::string> options);

    void PrintMATCONT(/* std::map<std::string,std::string> options */);
    void PrintMATLAB(std::map<std::string,std::string> options);
    void PrintOctave(std::map<std::string,std::string> options);
    void PrintRode(std::map<std::string,std::string> options);
    void PrintRdede(std::map<std::string,std::string> options);
    void PrintRadau5(std::map<std::string,std::string> options);
    void PrintLSODA(std::map<std::string,std::string> options);
    void PrintSciPy(std::map<std::string,std::string> options);
    void PrintPyDSTool(std::map<std::string,std::string> options);
    void PrintTaylor(std::map<std::string,std::string> options);
    void PrintPyGSL(std::map<std::string,std::string> options);
    void PrintLatex(/* std::map<std::string,std::string> options */);
    void PrintXPP(std::map<std::string,std::string> options);

    void Delay2ODE_ConvertExprToDefHist(GiNaC::ex& f);
    void Delay2ODE_ConvertAndExtend(GiNaC::ex& f, int N, int p);
    void PrintDelay2ODE(std::map<std::string,std::string> options);

    void PrintDDE23(std::map<std::string,std::string> options);

    void PrintDDE_SOLVER(std::map<std::string,std::string> options);

    void DDEBT_PrintParDerivs(std::ofstream &dout, const std::vector<GiNaC::ex> &e);
    void DDEBT_PrintJacobians(std::ofstream &dout, const std::vector<GiNaC::ex> &e);
    void DDEBT_PrintXandParJacobians(std::ofstream &dout, const std::vector<GiNaC::ex> &e);
    void DDEBT_PrintHessiansTimesV(std::ofstream &dout, const std::vector<GiNaC::ex> &e);
    void PrintDDEBIFTOOL(std::map<std::string,std::string> options);

    void PDDEC_PrintParDerivs(std::ofstream &dout, const std::vector<GiNaC::ex> &e);
    void PDDEC_PrintJacobians(std::ofstream &dout, const std::vector<GiNaC::ex> &e);
    void PDDEC_PrintXandParJacobians(std::ofstream &dout, const std::vector<GiNaC::ex> &e);
    void PDDEC_PrintHessiansTimesV(std::ofstream &dout, const std::vector<GiNaC::ex> &e);
    void PrintPDDECONT(/* std::map<std::string,std::string> options */);
};

#endif
