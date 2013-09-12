
//
// vf_differentials.cpp
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
#include <iostream>
#include <sstream>
#include <map>
#include <algorithm>
#include <ginac/ginac.h>

using namespace std;
using namespace GiNaC;

#include "vf_utils.h"


unsigned long factorial(unsigned n)
    {
    if (n == 0 || n == 1)
        return 1;
    unsigned long f = 2;
    for (unsigned k = 3; k <= n; ++k)
        f = k*f;
    return f;
    }

//
// ipow -- Integer power function.
//
int ipow(int base, int p)
    {
    if (p < 0)
        {
        cerr << "Error: ipow given negative power.\n";
        exit(-1);
        }
    if (base < 0)
        {
        cerr << "Error: ipow given negative base.\n";
        exit(-1);
        }
    if (base == 0 && p == 0)
        {
        cerr << "Error: ipow given base 0 and power 0.\n";
        exit(-1);
        }
    if (base == 0 && p > 0)
        return 0;
    if (base > 0 && p == 0)
        return 1;
    int m = base;
    for (int k = 1; k < p; ++k)
        m = m * base;
    return m;
    }

//
// next_partition
//
// Input:
//     r is an integer, r >= 1.
//     n is an integer, n >= 1.
//     part is an integer array of length n.  Each element of
//     part must be nonnegative, and the sum of the elements
//     must be r.
//
// Algorithm:
// Given an n-partition of r, [part[0],...,part[n-1]], do this:
//     last = part[n-1]
//     part[n-1] = 0
//     find the right-most nonzero entry part[k]
//     subtract 1 from the right-most nonzero entry, and set
//     part[k+1] = 1+last
//  If this process starts with [r,0,0,...,0], it will generate
//  all the partitions [r-1,1,0,...,0], [r-1,0,1,0,...,0], ...,
//  [r-2,2,0,...,0],...,[1,1,...,1],...,[0,0,...,1,r-1],[0,0,...,0,r].
//
//  Return:
//  The function returns 0 if part = [0,0,...,0,r].
//  Otherwise, it replaces part with the next partition, and returns 1.
//     

int next_partition(int r, int n, int part[])
    {
    // Check args
    if (r < 1)
        {
        cerr << "Error: next_partition given r < 1.\n";
        exit(-1);
        }
    if (r < 1)
        {
        cerr << "Error: next_partition given n < 1.\n";
        exit(-1);
        }
    // Be safe--check that part actually contains a valid partition.
    int sum = 0;
    for (int k = 0; k < n; ++k)
        {
        if (part[k] < 0)
            {
            cerr << "Error: bad partition given to next_partition: an entry is negative.\n";
            exit(-1);
            }
        else
            sum = sum + part[k];
        }
    if (sum != r)
        {
        cerr << "Error: bad partition given to next_partition: sum is not r.\n";
        exit(-1);
        }
    // OK, the partition is valid.  Find the next one.
    int last = part[n-1];
    if (last == r)
        return 0;
    part[n-1] = 0;
    int k = n-2;    // OK, since if n=1, this point would not be reached.
    while (part[k] == 0)
        k = k - 1;
    part[k] = part[k]-1;
    part[k+1] = 1+last;
    return 1;
    }

void print_partition(int n, int part[])
    {
    for (int k = 0; k < n; ++k)
        {
        cout.width(3);
        cout << part[k];
        }
    }

/*
long unsigned num_comb_partition(int r, int n, int part[])
    {
    long unsigned m = factorial(r);
    for (int k = 0; k < n; ++k)
        m = m / factorial(part[k]);
    return m;
    }
*/

void var_coeffs(int m, int p[], int r, int n)
    {
    for (int k = r; k > 0; --k)
        {
        int rem = m % n;
        p[k-1] = rem;
        m = m/n;
        }
    }

string p_to_deriv_string(lst vars, int p[])
    {
    int n = vars.nops();
    ostringstream dname;
    for (int j = 0; j < n; ++j)
        for (int k = 0; k < p[j]; ++k)
            dname << "_" << vars[j];
    return dname.str();
    }

string q_to_deriv_string(lst vars, int r, int q[])
    {
    ostringstream dname;
    for (int j = 0; j < r; ++j)
            dname << "_" << vars[q[j]];
    return dname.str();
    }

//
// generate_deriv --
//
// lang    -- Either "c" or "javascript"
// fout    -- Source file
// pout    -- C header file (ignored if lang="javascript")
// r       -- order of the derivative to generate
// vf      -- list containing the vector field formulas
// vars    -- vector field variables
// params  -- parameters used in the vector field
//

void generate_deriv(string lang, ofstream &fout, ofstream &pout, string name, int r,
                         lst vf, lst expreqn_list,
                         lst vars, lst params)
    {
    if (r < 1)
        {
        cerr << "Error: generate_deriv called with order=" << r << ". order must be greater than 0.\n";
        exit(-1);
        }
    //
    // Print the function header
    //
    fout << endl;
    fout << "/*\n";
    if (r == 1)
        fout << " *     deriv = Df(x)[v1]\n";
    else if (r == 2)
        fout << " *     deriv = D^2f(x)[v1,v2]\n";
    else if (r == 3)
        fout << " *     deriv = D^3f(x)[v1,v2,v3]\n";
    else
        fout << " *     deriv = D^" << r << "f(x)[v1,...,v" << r << "]\n";
    fout << " */\n";
    if (lang == "c")
        {
        fout << "void " << name << "_diff" << r << "(double deriv[],double x_[], double p_[],";
        pout << "void " << name << "_diff" << r << "(double deriv[], double x_[], double p_[],";
        }
    else
        {
        fout << "function " << name << "_diff" << r << "(x_, p_,";
        }
    for (int k = 0; k < r; ++k)
        {
        if (lang == "c")
            {
            fout << "double v" << k+1 << "_[]";
            pout << "double v" << k+1 << "_[]";
            }
        else
            {
            fout << "v" << k+1 << "_";
            }
        if (k < r-1)
            {
            fout << ",";
            if (lang == "c")
                pout << ",";
            }
        }
    fout << ")\n";
    if (lang == "c")
        pout << ");\n";
    fout << "    {\n";
    if (lang == "c")
        {
        CDeclare(fout,"double",vars);
        CDeclare(fout,"double",params);
        fout << endl;
        }
    const char* pre = "    ";
    if (lang == "javascript")
        pre = "    var ";
    GetFromVector(fout,pre,vars,"x_","[]",0,";");
    fout << endl;
    GetFromVector(fout,pre,params,"p_","[]",0,";");
    fout << endl;
/*
    int na = exprnames.nops();
    for (int i = 0; i < na; ++i)
        {
        fout << "    " << exprnames[i] << " = " << exprformulas[i] << ";" << endl;
        }
    if (na > 0)
        fout << endl;
*/

    if (lang == "javascript")
        fout << "    var deriv = [];\n";
    int n = vars.nops();
    int *p = new int[n];
    for (int i = 0; i < n; ++i)
        {
        if (lang == "c")
            fout << "    {\n";
        // ex f = vf[i];
        ex f = iterated_subs(vf[i],expreqn_list);
        map<string,ex> derivatives;
        // Initialize p to [r,0,0,...,0]
        p[0] = r;
        for (int k = 1; k < n; ++k)
            p[k] = 0;
        fout << "    /*\n";
        fout << "     *  Partial derivatives of vf[" << i << "]. \n";
        fout << "     *  Any derivative not listed here is zero.\n";
        fout << "     */\n";
        do 
            {
            // Compute the partial derivative

            ex df = f;
            for (int j = 0; j < n; ++j)
                df = df.diff(ex_to<symbol>(vars[j]),p[j]);
            if (df != 0)
                {
                string dname = p_to_deriv_string(vars,p);
                derivatives[dname] = df;
                if (lang == "c")
                    fout << "    double vf";
                else
                    fout << "    var vf";
                fout << dname;
                fout << " = ";
                fout << df << ";\n";
                }
            }
        while (next_partition(r,n,p));

        int A = ipow(n,r);
        int *q = new int[r];
        int *sq = new int[r];
        fout << "    deriv[" << i << "] = 0.0;\n";
        for (int m = 0; m < A; ++m)
            {
            var_coeffs(m,q,r,n);
            var_coeffs(m,sq,r,n);
            sort(sq,sq+r);
            string dname = q_to_deriv_string(vars,r,sq);
            map<string,ex>::iterator lookupdname = derivatives.find(dname);
            if (lookupdname != derivatives.end())
                {
                fout << "    deriv[" << i << "] += vf";
                fout << q_to_deriv_string(vars,r,sq);
                fout << " * ";
                for (int k = 0; k < r; ++k)
                    {
                    fout << "v" << k+1 << "_[" << q[k] << "]";
                    if (k == r-1)
                        fout << ";\n";
                    else
                        fout << "*";
                    }
                }
            }
        if (lang == "c")
            fout << "    }\n";
        delete [] q;
        delete [] sq;
        }
    delete [] p;
    fout << endl;
    if (lang == "c")
        fout << "    return;\n";
    else
        fout << "    return deriv;\n";
    fout << "    }\n";
    }
