
//
//  vf_julia.cpp
//
//  This file defines the VectorField::Julia method.
//
//
//  Copyright (C) 2024 Warren Weckesser
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


#include <cassert>
#include <fstream>
#include <string>
#include <ginac/ginac.h> 

#include "vf.h"
#include "vf_utils.h"

using namespace std;
using namespace GiNaC;

static void
generate_solver_struct(ofstream& out)
{
    out << "struct Solver\n";
    out << "    abstol::Float64\n";
    out << "    reltol::Float64\n";
    out << "    stoptime::Float64\n";
    out << "end\n";
}

static void
generate_solver_options(ofstream& out)
{
    out << "function get_solver_options()\n";
    out << "    abstol = 1e-12\n";
    out << "    reltol = 1e-8\n";
    out << "    stoptime = 10.0\n";
    out << "    for arg in ARGS\n";
    out << "        fields = split(arg, \"=\")\n";
    out << "        if length(fields) == 1\n";
    out << "            println(\"invalid argument '\", arg, \"'; missing '='\")\n";
    out << "            exit()\n";
    out << "        end\n";
    out << "        if length(fields) > 2\n";
    out << "            println(\"invalid argument '\", arg, \"'; too many '=' characters\")\n";
    out << "            exit()\n";
    out << "        end\n";
    out << "        name, value = fields\n";
    out << "        if name == \"abstol\"\n";
    out << "            abstol = parse(Float64, value)\n";
    out << "        elseif name == \"reltol\"\n";
    out << "            reltol = parse(Float64, value)\n";
    out << "        elseif name == \"stoptime\"\n";
    out << "            stoptime = parse(Float64, value)\n";
    out << "        else\n";
    out << "            println(\"Bad argument '\", arg, \"'; only 'abstol', 'reltol' and 'stoptime' are handled.\")\n";
    out << "            exit()\n";
    out << "        end\n";
    out << "    end\n";
    out << "    Solver(abstol, reltol, stoptime)\n";
    out << "end\n";
    out << endl;
}

static void
generate_julia_function_def_hack(ofstream &out)
{
    out << endl;
    out << "# In case these are used in any expressions..." << endl;
    out << "atan2(y, x) = atan(y, x)" << endl;
    out << "pow(x, p) = x^p" << endl;
    out << endl;
}

static void
generate_lag_assignment(ofstream& fout, const lst& lag, symbol& indvar)
{
    fout << "    " << lag.op(0) << " = h_(p_, " << indvar - lag.op(3) <<  "; idxs = " <<  ex_to<numeric>(lag.op(1)).to_int() << ")\n";
}

static void
generate_common_plot_code(ofstream& out, string name,
                          string IndependentVariable,
                          lst &varname_list)
{
    size_t nv = varname_list.nops();
    if (nv > 1) {
        out << "ylabels = ";
        ListOfStrings(out, varname_list, "[]", " ");
        out << endl;
        out << "plot(sol, ";
        out <<      "xaxis=\"" << IndependentVariable << "\", ";
        out <<      "ylabel=ylabels, ";
        out <<      "layout=(" << nv << ", 1), ";
        out <<      "legend=false)" << endl;
    }
    else {
        out << "plot(sol, ";
        out <<      "xaxis=\"" << IndependentVariable << "\", ";
        out <<      "ylabel=\"" << varname_list[0] << "\", ";
        out <<      "legend=false)" << endl;
    }

    out << "savefig(\"" << name << "_plot.svg\")" << endl;
}

void VectorField::PrintJuliaFuncStart(ofstream &fout)
{
    size_t np = parname_list.nops();

    if (HasPi) {
        fout << "    Pi = pi\n";
    }
    AssignNameValueLists(fout, "    ", conname_list, "=", convalue_list, "");
    GetFromVector(fout, "    ", varname_list, "=", "u_", "[]", 1, "");
    fout << endl;
    if (np > 0) {
        GetFromVector(fout, "    ", parname_list, "=", "p_", "[]", 1, "");
        fout << endl;
    }

    if (!IsDelay) {
        AssignNameValueLists(fout, "    ", exprname_list, "=", exprformula_list, "");
    }
    else {
        vector<lst> references;
        size_t nv = varname_list.nops();
        size_t na = exprname_list.nops();

        assert((na + nv) > 0);

        ex tuple_expr = to_nested_tuple(exprformula_list, varvecfield_list);
        lst tmp;
        convert_delay_to_lagvalue(tuple_expr, tmp);
        lst new_list = to_list(tuple_expr);
        exprformula_list = lst();
        for (size_t k = 0; k < na; ++k) {
            exprformula_list.append(new_list.op(k));
        }
        varvecfield_list = lst();
        for (size_t k = na; k < na + nv; ++k) {
            varvecfield_list.append(new_list.op(k));
        }

        int k = 0;
        for (lst::const_iterator iter = new_list.begin(); iter != new_list.end(); ++iter) {
            lst ref;
            ex e = *iter;
            for (lst::const_iterator p = tmp.begin(); p != tmp.end(); ++p) {
                lst lag_data = ex_to<lst>(*p);
                if (e.has(lag_data.op(0))) {
                    ref.append(lag_data);
                }
            }
            ++k;
            references.push_back(ref);
        }

        ex generated = 0;
        for (size_t i = 0; i < na; ++i) {
            for (lst::const_iterator iter = references[i].begin(); iter != references[i].end(); ++iter) {
                lst e = ex_to<lst>(*iter);
                if (!generated.has(e.op(0))) {
                    generate_lag_assignment(fout, e, IndVar);
                    generated = generated + e.op(0);
                }
            }

            fout << "    " << exprname_list[i] << " = " << exprformula_list[i] << endl;
        }

        for (size_t i = na; i < na + nv; ++i) {
            for (lst::const_iterator iter = references[i].begin(); iter != references[i].end(); ++iter) {
                lst e = ex_to<lst>(*iter);
                if (!generated.has(e.op(0))) {
                    generate_lag_assignment(fout, e, IndVar);
                    generated = generated + e.op(0);
                }
            }
        }
    }
}

//
// PrintJulia -- The Julia Code Generator.
//

void VectorField::PrintJulia(map<string,string> options)
{
    size_t nv, np, nf;   // currently na is not used.

    nv = varname_list.nops();
    np = parname_list.nops();
    // na = exprname_list.nops();
    nf = funcname_list.nops();

    string filename = Name() + ".jl";
    ofstream fout;
    fout.open(filename.c_str());
    fout << csrc << left;

    //
    //  Print Python file header information.
    //
    fout << "#" << endl;
    fout << "# " << filename << endl;
    fout << "#" << endl;
    fout << "# Function for the vector field named: " << Name() << endl;
    fout << "#\n";
    PrintVFGENComment(fout, "# ");
    fout << "#\n";

    generate_julia_function_def_hack(fout);

    //
    //  Print the vector field function.
    //
    fout << endl;
    fout << "function " << Name() << "!(du_, u_, p_, " << IndependentVariable << ")" << endl;
    PrintJuliaFuncStart(fout);
    fout << endl;
    for (size_t i = 0; i < nv; ++i) {
        fout << "    du_[" << i+1 << "] = " << varvecfield_list[i] << endl;
    }
    fout << "end" << endl;

    if (options["func"] == "yes" && nf > 0) {
        // Generate functions...
        fout << endl;
        fout << "# Functions..." << endl;
        fout << endl;
        for (size_t n = 0; n < nf; ++n) {
            fout << "function " << Name() << "_" << funcname_list[n] << "(u_, p_, t)" << endl;
            PrintJuliaFuncStart(fout);
            fout << endl;
            fout << "    " << funcformula_list[n] << endl;
            fout << "end" << endl;
            fout << endl;
        }
    }
    fout.close();

    if (options["demo"] == "yes") {
        //
        //  Create a self-contained ODE solver for this vector field.
        //

        string tfilename = Name() + "_demo.jl";
        ofstream tout;
        tout.open(tfilename.c_str());
        tout << csrc << left;
        tout << "#" << endl;
        tout << "# " << tfilename << endl;
        tout << "#" << endl;
        tout << "# Demo program for the vector field named: " << Name() << endl;
        tout << "#" << endl;
        PrintVFGENComment(tout, "# ");
        tout << "#\n" ;
        tout << endl;
        tout << "using DifferentialEquations" << endl;
        tout << "using Plots" << endl;
        tout << endl;
        tout << "include(\"" << Name() << ".jl\")" << endl;
        tout << endl << endl;
        generate_solver_struct(tout);
        tout << endl;
        generate_solver_options(tout);
        tout << endl;
        if (HasPi) {
            tout << "Pi = pi\n";
        }
        AssignNameValueLists(tout, "", conname_list, "=", convalue_list, "");
        tout << "u0 = [";
        PrintList(tout, vardefic_list);
        tout << "]\n" ;
        if (np > 0) {
            tout << "p = [" ;
            PrintList(tout, pardefval_list);
            tout << "]\n" ;
        }
        tout << endl;
        tout << "solver_options = get_solver_options()\n";
        tout << endl;
        tout << "tspan = (0.0, solver_options.stoptime)" << endl;
        tout << endl;

        tout << endl;
        tout << "prob = ODEProblem(" << Name() << "!, u0, tspan";
        if (np > 0) {
            tout << ", p";
        }
        tout << ")" << endl; 
        tout << "sol = solve(prob, abstol=solver_options.abstol, reltol=solver_options.reltol)" << endl;
        tout << endl;
        generate_common_plot_code(tout, Name(), IndependentVariable, varname_list);
        tout.close();
    }
}


void VectorField::PrintJuliaDelay(map<string,string> options)
{
    size_t nv, np, nf;

    nv = varname_list.nops();
    np = parname_list.nops();
    //na = exprname_list.nops();
    nf = funcname_list.nops();

    string filename = Name() + ".jl";
    ofstream fout;
    fout.open(filename.c_str());
    fout << csrc << left;

    //
    //  Print Python file header information.
    //
    fout << "#" << endl;
    fout << "# " << filename << endl;
    fout << "#" << endl;
    fout << "# Function for the vector field named: " << Name() << endl;
    fout << "#\n";
    PrintVFGENComment(fout, "# ");
    fout << "#\n";

    generate_julia_function_def_hack(fout);

    //
    //  Print the vector field function.
    //
    fout << endl;
    fout << "function " << Name() << "!(du_, u_, h_, p_, " << IndependentVariable << ")" << endl;
    PrintJuliaFuncStart(fout);

    fout << endl;
    for (size_t i = 0; i < nv; ++i) {
        fout << "    du_[" << i+1 << "] = " << varvecfield_list[i] << endl;
    }
    fout << "end" << endl;

    if (options["func"] == "yes" && nf > 0) {
        // Generate functions...
        fout << endl;
        fout << "# Functions..." << endl;
        fout << endl;
        for (size_t n = 0; n < nf; ++n) {
            fout << "function " << Name() << "_" << funcname_list[n] << "(u_, p_, t)" << endl;
            PrintJuliaFuncStart(fout);
            fout << endl;
            fout << "    " << funcformula_list[n] << endl;
            fout << "end" << endl;
            fout << endl;
        }
    }
    fout.close();

    if (options["demo"] == "yes") {
        //
        //  Create a self-contained ODE solver for this vector field.
        //

        string tfilename = Name() + "_demo.jl";
        ofstream tout;
        tout.open(tfilename.c_str());
        tout << csrc << left;
        tout << "#" << endl;
        tout << "# " << tfilename << endl;
        tout << "#" << endl;
        tout << "# Demo program for the vector field named: " << Name() << endl;
        tout << "#" << endl;
        PrintVFGENComment(tout, "# ");
        tout << "#\n" ;
        tout << endl;
        tout << "using DifferentialEquations" << endl;
        tout << "using Plots" << endl;
        tout << endl;
        tout << "include(\"" << Name() << ".jl\")" << endl;
        tout << endl << endl;
        generate_solver_struct(tout);
        tout << endl;
        generate_solver_options(tout);
        tout << endl;
        tout << "function history(p_, " << IndependentVariable << "; idxs = nothing)\n";
        AssignNameValueLists(tout, "    ", conname_list, "=", convalue_list, "");
        if (np > 0) {
            GetFromVector(tout, "    ", parname_list, "=", "p_", "[]", 1, "");
            tout << endl;
        }
        tout << "    if typeof(idxs) <: Number\n";
        tout << "        ";
        for (size_t k = 0; k < nv - 1; ++k) {
            tout << "if idxs == " << k + 1 << endl;
            tout << "            " << vardefhist_list[k] << endl;
            tout << "        else";
            if (k == nv - 2) {
                tout << endl;
                tout << "            ";
            }
        }
        tout << vardefhist_list[nv - 1] << endl;
        if (nv > 1) {
            tout << "        end\n";
        }
        tout << "    else\n";
        tout << "        [";
        for (size_t k = 0; k < nv; ++k) {
            if (k > 0) {
                tout << "         ";
            }
            tout << vardefhist_list[k];
            if (k < nv - 1) {
                tout << ",\n";
            }
        }
        tout << "]\n";
        tout << "    end\n";
        tout << "end\n\n";

        vector<ex> constant_delays, dependent_delays;
        for (size_t k = 0; k < Delays.size(); ++k) {
            ex f = iterated_subs(Delays[k], expreqn_list);
            bool dependent = false;
            if (has(f, IndVar)) {
                dependent = true;
            }
            for (size_t j = 0; j < Delays.size(); ++j) {
                if (has(f, varname_list[j])) {
                    dependent = true;
                }
            }
            if (dependent) {
                dependent_delays.push_back(f);
            }
            else {
                constant_delays.push_back(f);
            }
        }

        for (size_t k = 0; k < dependent_delays.size(); ++k) {
            tout << "function dependent_lag" << k + 1 << "(u_, p_, " << IndependentVariable << ")\n";
            AssignNameValueLists(tout, "    ", conname_list, "=", convalue_list, "");
            GetFromVector(tout, "    ", varname_list, "=", "u_", "[]", 1, "");
            tout << endl;
            if (np > 0) {
                GetFromVector(tout, "    ", parname_list, "=", "p_", "[]", 1, "");
                tout << endl;
            }
            tout << "    " << dependent_delays[k] << endl;
            tout << "end" << endl;
            tout << endl;
        }

        if (HasPi) {
            tout << "Pi = pi\n";
        }
        AssignNameValueLists(tout, "", conname_list, "=", convalue_list, "");
        AssignNameValueLists(tout, "", parname_list, "=", pardefval_list, "");
        tout << endl;
        tout << "u0 = [";
        PrintList(tout, vardefic_list);
        tout << "]\n" ;
        if (np > 0) {
            tout << "p = [" ;
            PrintList(tout, parname_list);
            tout << "]\n" ;
        }
        tout << endl;
        tout << "solver_options = get_solver_options()\n";
        tout << "tspan = (0.0, solver_options.stoptime)" << endl;
        tout << endl;
        tout << "prob = DDEProblem(" << Name() << "!, u0, history, tspan";
        if (np > 0) {
            tout << ", p;\n";
        }
        else {
            tout << ";\n";
        }

        if (constant_delays.size() > 0) {
            tout << "                  constant_lags = [";
            for (size_t k = 0; k < constant_delays.size(); ++k) {
                if (k > 0) {
                    tout << ", ";
                }
                tout << constant_delays[k];
            }
            tout << "]";
        }
        if (dependent_delays.size() == 0) {
            tout << ")" << endl;
        }
        else {
            if (constant_delays.size() > 0) {
                tout << ",\n";
            }
            tout << "                  dependent_lags = [";
            for (size_t k = 0; k < dependent_delays.size(); ++k) {
                if (k > 0) {
                    tout << ", ";
                }
                tout << "dependent_lag" << k + 1;
            }
            tout << "])" << endl;
        }
        tout << "alg = MethodOfSteps(Tsit5())" << endl;
        tout << "sol = solve(prob, alg, abstol=solver_options.abstol, reltol=solver_options.reltol)" << endl;
        tout << endl;
        generate_common_plot_code(tout, Name(), IndependentVariable, varname_list);
        tout.close();
    }
}
