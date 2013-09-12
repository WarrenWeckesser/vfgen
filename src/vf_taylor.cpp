

//
//  vf_taylor.cpp
//
//  This file defines the VectorField::PrintTaylor method.
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
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
// #include <ctime>
#include <ginac/ginac.h> 

#include "vf.h"
#include "vf_utils.h"
#include "MyVec.h"

using namespace std;
using namespace GiNaC;


void generate_deriv(string lang, ofstream &fout, ofstream &pout, string name, int r, lst vf, lst expreqn, lst vars, lst params);


long int factorial(long int);


//
// PrintTaylor -- this is the main function that prints the Taylor C function.
//

void VectorField::PrintTaylor(map<string,string> options)
    {
    int nc, np, nv, na, nf;

    symbol t(IndependentVariable);
    nc = conname_list.nops();
    nv = varname_list.nops();
    np = parname_list.nops();
    na = exprname_list.nops();
    nf = funcname_list.nops();

    int Order;
    if (options.find("order") == options.end())
        {
        Order = 5;
        cerr << "The order option was not specified; the default order=" << Order << " will be used.\n";
        options["order"] = "5";
        }
    else
        {
        Order = string_to_int(options["order"]);
        if (Order < 1)
            {
            cerr << "Error: Bad value for the order parameter: " << Order << "; this must be a positive integer.\n";
            exit(-1);  // Should handle this in a better way?
            }
        }

    string filename = Name()+"_taylor" + options["order"] + ".c";
    ofstream fout;
    fout.open(filename.c_str());
    fout << csrc << left;

    string pfilename = Name()+"_taylor" + options["order"] + ".h";
    ofstream pout;
    pout.open(pfilename.c_str());
    pout << csrc << left;

    //
    //  Print C file header information.
    //
    fout << "/*" << endl;
    fout << " *  " << filename << endl;
    fout << " *" << endl;
    fout << " *  C file with functions for computing the Taylor series approximate solution" << endl;
    fout << " *  for the vector field named: " << Name() << endl;
    fout << " *" << endl;
    PrintVFGENComment(fout," *  ");
    fout << " */" << endl;
    fout << endl;

    pout << "/*" << endl;
    pout << " *  " << pfilename << endl;
    pout << " *" << endl;
    pout << " *  C prototype file for the functions defined in " << filename << endl;
    pout << " *" << endl;
    PrintVFGENComment(pout," *  ");
    pout << " */" << endl;
    pout << endl;

    fout << "#include <math.h>" << endl;
    fout << "#include \"" << Name() << "_taylor" << Order << ".h\"\n";
    fout << endl;
    //
    //  Print the vector field function.
    //
    fout << "/*" << endl;
    fout << " *  The vector field." << endl;
    fout << " */" << endl;
    fout << endl;
    fout << "void " << Name() << "_vf(double t, const double y_[], double f_[], double params[])" << endl;
    pout << "void " << Name() << "_vf(double, const double [], double [], double []);" << endl;
    fout << "    {" << endl;
    for (int i = 0; i < nc; ++i)
        {
        fout << "    const double " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    CDeclare_double(fout,varname_list);
    CDeclare_double(fout,parname_list);
    CDeclare_double(fout,exprname_list);
    fout << "    double *p_;" << endl;
    fout << endl;
    fout << "    p_ = (double *) params;" << endl;
    fout << endl;
    GetFromVector(fout,"    ",varname_list,"y_","[]",0,";");
    fout << endl;
    GetFromVector(fout,"    ",parname_list,"p_","[]",0,";");
    fout << endl;
    for (int i = 0; i < na; ++i)
        {
        fout << "    " << exprname_list[i] << " = " << exprformula_list[i] << ";" << endl;
        }
    if (na > 0)
        fout << endl;
    for (int i = 0; i < nv; ++i)
        {
        fout << "    f_[" << i << "]" << " = " << varvecfield_list[i] << ";" << endl;
        }
    fout << endl;
    fout << "    }" << endl;
    fout << endl;


    //
    // Call generate_deriv for each integer from 1 to Order-1.
    //

    for (int i = 1; i < Order; ++i)
        generate_deriv("c",fout, pout, Name(), i, varvecfield_list, expreqn_list,
                                       varname_list, parname_list);

    //
    // Create the Taylor series function.
    // First we need to compute some coefficients.
    //

    vectormap *table = new vectormap[Order-1];

    vector<int> *initial = new vector<int>;
    initial->push_back(1);  // initial = [1]
    table[0][initial] = 1.0;

    // cout << "Generating table of coefficients\n";
    // clock_t tgen = clock();
    for (int k = 1; k < Order-1; ++k)
        {
        for (vectormap::iterator v = table[k-1].begin(); v != table[k-1].end(); ++v)
            {
            vector<int> *a = v->first;
            double coeff = v->second;

            vector<int> *a1 = new vector<int>;
            CopyMyVec(a,a1);
            ++(*a1)[0];
            table[k][a1] = table[k][a1] + coeff;

            for (size_t i = 0; i < a->size(); ++i)
                {
                int m = (*a)[i];
                if (m > 0)
                    {
                    vector<int> *a2 = new vector<int>;
                    CopyMyVec(a,a2);
                    --(*a2)[i];
                    if (i == a->size()-1)
                        a2->push_back(1);
                    else
                        ++(*a2)[i+1];               
                    table[k][a2] = table[k][a2] + m*coeff;
                    }
                }
            }
        }
    // tgen = clock()-tgen;
    // double tgen_secs = ((double) tgen)/CLOCKS_PER_SEC;
    // cerr << "Coeff gen time: " << tgen_secs << " seconds\n";

    fout << endl;
    fout.precision(0);
    fout.setf(ios_base::fixed);
    fout.setf(ios_base::showpoint);
    fout << "/*\n";
    fout << " *  " << Name() << "_derivs" << Order << endl;
    fout << " *\n";
    fout << " *  Compute the coefficients in the Taylor polynomial at X.\n";
    fout << " *  These are just the derivatives; they have not been scaled\n";
    fout << " *  by the appropriate factorial.\n";
    fout << " *\n";
    fout << " */\n";
    fout << endl;
    fout << "void " << Name() << "_derivs" << Order << "(double Xderiv[" << Order << "][" << nv << "], double X[], double params[])\n";
    pout << "void " << Name() << "_derivs" << Order << "(double Xderiv[" << Order << "][" << nv << "], double X[], double params[]);\n";
    fout << "    {\n";
    fout << "    int i;\n";
    fout << "    double s;\n";
    // fout << "    double Xderiv[" << Order << "][" << nv << "];\n";
    fout << "    double Q[" << nv << "];\n";
    fout << endl;
    fout << "    " << Name() << "_vf(0.0,X,Xderiv[0],params);\n";
    for (int k  = 0; k < Order-1; ++k)
        {
        fout << endl;
        fout << "    for (i = 0; i < " << nv << "; ++i)\n";
        fout << "        Xderiv[" << k+1 << "][i] = 0.0;\n";
        for (vectormap::iterator v = table[k].begin(); v != table[k].end(); ++v)
            {
            vector<int> *a = v->first;
            double coeff = v->second;
            fout << "    /*    [";
            PrintMyVec(fout,a);
            fout << "]  coeff = " << coeff << "  */\n";
            // int r = a->size();
            int r = SumVec(a);
            fout << "    " << Name() << "_diff" << r << "(Q,X,params";
            int i = 0;
            for (vector<int>::iterator iter = a->begin(); iter != a->end(); ++iter)
                {
                for (int j = 0; j < *iter; ++j)
                    fout << ",Xderiv[" << i << "]";
                ++i;
                }
            fout << ");\n";
            fout << "    for (i = 0; i < " << nv << "; ++i)\n";
            fout << "        Xderiv[" << k+1 << "][i] += ";
            if (coeff > 1)
                fout << coeff << "*";
            fout << "Q[i];\n";
            }
        }
    fout << "    }\n";
    fout << endl;
    delete [] table;
    fout << "/*\n";
    fout << " *  " << Name() << "_evaltaylor" << Order << endl;
    fout << " *\n";
    fout << " *  Use the Taylor method to approximate X(t+h) given X(t).\n";
    fout << " *  This function uses a Taylor polynomial of order " << Order << ".\n";
    fout << " *  " << Name() << "_derivs" << Order << "(...) must be called first to fill in the array Xderiv.\n";
    fout << " */\n";
    fout << endl;
    fout << "void " << Name() << "_evaltaylor" << Order << "(double Xnew[], double h, double X[], double Xderiv[" << Order << "][" << nv << "])\n";
    pout << "void " << Name() << "_evaltaylor" << Order << "(double Xnew[], double h, double X[], double Xderiv[" << Order << "][" << nv << "]);\n";
    fout << "    {\n";
    fout << "    int i;\n";
    fout << "    double s;\n";
    fout << "    /*  Xnew = X  */\n";
    fout << "    for (i = 0; i < " << nv << "; ++i)\n";
    fout << "        Xnew[i] = X[i];\n";
    for (int k = 1; k <= Order; ++k)
        {
        if (k == 1)
            fout << "    s = h;\n";
        else
            fout << "    s = s * h;\n";
        fout << "    /* Add order " << k << " term to Xnew */\n";
        fout << "    for (i = 0; i < " << nv << "; ++i)\n";
        fout << "        Xnew[i] += (1.0/" << factorial(k) << ")*Xderiv[" << k-1 << "][i]*s;" << endl;
        }
    fout << "    }\n";

    fout.close();
    pout.close();
    }
