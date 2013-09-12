

//
//  vf_javascript.cpp
//
//  This file defines the VectorField::PrintJavascript method.
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
#include <ginac/ginac.h> 
#include <cmath>

#include "vf.h"
#include "vf_utils.h"
#include "MyVec.h"

using namespace std;
using namespace GiNaC;


void generate_deriv(string lang, ofstream &fout, ofstream &pout, string name, int r, lst vf, lst expreqn, lst vars, lst params);


long int factorial(long int);


//
// PrintJavascript -- this is the main function that prints the Javascript code.
//

void VectorField::PrintJavascript(map<string,string> options)
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

    string filename = Name() + ".js";
    ofstream fout;
    fout.open(filename.c_str());
    fout << csrc << left;

    //
    //  Print header information.
    //
    fout << "/*" << endl;
    fout << " *  " << filename << endl;
    fout << " *" << endl;
    PrintMultilineComment(fout,
        "Javascript code with functions for computing the Taylor series approximate solutions\nfor the vector field named: "+Name()+"\n",
        ""," *  ");
    PrintVFGENComment(fout," *  ");
    fout << " */" << endl;
    fout << endl;

    fout << "/*\n";
    fout << " *  Make several Math functions available without the Math prefix.\n";
    fout << " */\n";
    fout << "acos  = Math.acos\n";
    fout << "asin  = Math.asin\n";
    fout << "atan  = Math.atan\n";
    fout << "atan2 = Math.atan2\n";
    fout << "cos   = Math.cos\n";
    fout << "exp   = Math.exp\n";
    fout << "log   = Math.log\n";
    fout << "pow   = Math.pow\n";
    fout << "sin   = Math.sin\n";
    fout << "sqrt  = Math.sqrt\n";
    fout << "tan   = Math.tan\n";
    fout << endl;

    //
    //  Print the vector field function.
    //
    fout << "/*" << endl;
    fout << " *  The vector field." << endl;
    fout << " */" << endl;
    fout << endl;
    fout << "function " << Name() << "_vf(t, y_, params)" << endl;
    fout << "    {" << endl;
    for (int i = 0; i < nc; ++i)
        {
        fout << "    var " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    GetFromVector(fout,"    var ",varname_list,"y_","[]",0,";");
    fout << endl;
    GetFromVector(fout,"    var ",parname_list,"params","[]",0,";");
    fout << endl;
    for (int i = 0; i < na; ++i)
        {
        fout << "    " << exprname_list[i] << " = " << exprformula_list[i] << ";" << endl;
        }
    if (na > 0)
        fout << endl;
    fout << "    var f_ = [];\n";
    for (int i = 0; i < nv; ++i)
        {
        fout << "    f_[" << i << "]" << " = " << varvecfield_list[i] << ";" << endl;
        }
    fout << endl;
    fout << "    return f_;\n";
    fout << "    }" << endl;
    fout << endl;


    //
    // Call generate_deriv for each integer from 1 to Order-1.
    //

    for (int i = 1; i < Order; ++i)
        generate_deriv("javascript",fout, fout, Name(), i, varvecfield_list, expreqn_list,
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
    fout << "function " << Name() << "_derivs" << Order << "(X, params)\n";
    fout << "    {\n";
    fout << "    var i;\n";
    fout << "    var s;\n";
    fout << "    var Q = [];\n";
    fout << endl;
    fout << "    var Xderiv = [];\n";
    fout << "    Xderiv[0] = " << Name() << "_vf(0.0,X,params);\n";
    for (int k  = 0; k < Order-1; ++k)
        {
        fout << endl;
        fout << "    Xderiv[" << k+1 << "] = []\n";
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
            fout << "    Q = " << Name() << "_diff" << r << "(X,params";
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
    fout << endl;
    fout << "    return Xderiv;\n";
    fout << "    }\n";
    fout << endl;
    delete [] table;
    fout << "/*\n";
    fout << " *  " << Name() << "_evaltaylor" << Order << endl;
    fout << " *\n";
    fout << " *  Use the Taylor method to approximate X(t+h) given X(t).\n";
    fout << " *  This function uses a Taylor polynomial of order " << Order << ".\n";
    fout << " *  Xderiv must be created by calling " << Name() << "_derivs" << Order << "(...).\n";
    fout << " */\n";
    fout << endl;
    fout << "function " << Name() << "_evaltaylor" << Order << "(h, X, Xderiv)\n";
    fout << "    {\n";
    fout << "    var i;\n";
    fout << "    var s;\n";
    fout << "    var Xnew = X\n";
    // fout << "    for (i = 0; i < " << nv << "; ++i)\n";
    // fout << "        Xnew[i] = X[i];\n";
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
    fout << "    return Xnew;\n";
    fout << "    }\n";

    fout.close();
    
    if (options["demo"] == "yes")
        {
        //
        // Create the HTML file.
        //
        string htmlfilename = Name()+".html";
        ofstream hout;
        hout.open(htmlfilename.c_str()); // FIXME!  Check for failure.
        hout << "<html>\n";
        hout << "<head>\n";
        hout << "<script src=\"" << Name() << ".js\"></script>\n";
        hout << "<script src=\"" << Name() << "_solverdemo.js\"></script>\n";
        hout << "</head>\n";
        hout << "<body>\n";
        hout << "<h3>javascript solver (Taylor method, order " << Order << ") for the <i>" << Name() << "</i> vector field</h3>\n";
        hout << "<div style=\"float: left; width: 502px;  margin-right: 8px;\">\n";
        hout << "<div style=\"background: #F0F0E0; border: 1px solid red;\">\n";
        hout << "<canvas id=\"canvas\" width=\"500\" height=\"500\"></canvas>\n";
        hout << "</div>\n";
        hout << "Generated by <a href=\"http://www.warrenweckesser.net/vfgen\">VFGEN</a> (version " << VERSION << ")\n";
        hout << "</div>\n";
        hout << "<div>\n";

        hout << "<table border=\"1\" >\n";
        hout << "<tr>\n";
        hout << "<td><i>Parameter</i></td><td><i>Value</i></td>\n";
        hout << "</tr>\n";
        for (int i = 0; i < np; ++i)
            {
            ex val = pardefval_list[i];
            for (int j = i-1; j >= 0; j--)
                val = val.subs(parname_list[j]==pardefval_list[j]);
            for (int j = nc-1; j >= 0; j--)
                val = val.subs(conname_list[j]==convalue_list[j]);
            val = val.subs(Pi==M_PI);

            hout << "<tr>\n";
            hout << "<td>" << parname_list[i] << "</td>\n";
            hout << "<td><input id=\"" << parname_list[i] << "_value\"  type=\"text\" size=\"9\" value=\"" << val << "\" /></td>\n";
            hout << "</tr>\n";
            }
        hout << "</table>\n";
        hout << "<br />\n";
        hout << "<table border=\"1\" >\n";
        hout << "<tr>\n";
        hout << "<td>time</td><td>&nbsp;0</td><td><input id=\"t\" type=\"text\" size=\"9\" value=\"0\" readonly /></td>\n";
        hout << "</tr>\n";
        for (int i = 0; i < nv; ++i)
            {
            ex ic = vardefic_list[i];
            for (int j = np-1; j >= 0; j--)
                ic = ic.subs(parname_list[j]==pardefval_list[j]);
            for (int j = nc-1; j >= 0; j--)
                ic = ic.subs(conname_list[j]==convalue_list[j]);
            ic = ic.subs(Pi==M_PI);
            hout << "<tr>\n";
            hout << "<td>" << varname_list[i] << "</td>\n";
            hout << "<td><input id=\"" << varname_list[i] << "0\"  type=\"text\" size=\"9\" value=\"" << ic << "\" /></td>\n";
            hout << "<td><input id=\"" << varname_list[i] << "_t\" type=\"text\" size=\"9\" value=\"" << ic << "\" readonly /></td>\n";
            hout << "</tr>\n";
            }
        hout << "</table>\n";
        hout << "<br />\n";
        hout << "<div>\n";
        hout << "Stop time:\n";
        hout << "<input id=\"stoptime\" type=\"text\" size=\"9\" value=\"100.0\" />\n";
        hout << "<br />\n";
        hout << "Solver step tolerance:\n";
        hout << "<input id=\"tolerance\" type=\"text\" size=\"9\" value=\"1.0e-6\" />\n";
        hout << "<br />\n";
        hout << "Maximum time step:\n";
        hout << "<input id=\"hmax\" type=\"text\" size=\"9\" value=\"0.5\" />\n";
        hout << "<br /><br/>\n";
        hout << "Display X = ";
        for (int i = 0; i < nv; ++i)
            {
            double scale;
            if (i == 0)
                scale = 10.0;
            else
                scale = 0.0;
            hout << "<input id=\"" << varname_list[i] << "_xscale\" type=\"text\" size=\"4\" value=\"" << scale << "\" />*" << varname_list[i] << " + ";
            }
        hout << "<input id=\"xoffset\" type=\"text\" size=\"4\" value=\"0\" />\n";
        hout << "<br />\n";
        hout << "Display Y = ";
        for (int i = 0; i < nv; ++i)
            {
            double scale;
            if (i == 1)
                scale = 10.0;
            else
                scale = 0.0;
            hout << "<input id=\"" << varname_list[i] << "_yscale\" type=\"text\" size=\"4\" value=\"" << scale << "\" />*" << varname_list[i] << " + ";
            }
        hout << "<input id=\"yoffset\" type=\"text\" size=\"4\" value=\"0\" />\n";
        hout << "<br />\n";
        hout << "Display X and Y ranges are (-250,250).\n";
        hout << "</div>\n";
        hout << "<br />\n";
        hout << "Number of steps per display update:\n";
        hout << "<input id=\"numsteps\" type=\"text\" size=\"5\" value=\"100\" />\n";
        hout << "<br />\n";
        hout << "Color:\n";
        hout << "<input  id=\"black\"      type=\"radio\" name=\"palette\" value=\"black\" checked>black</input> &nbsp;\n";
        hout << "<input  id=\"blue+green\" type=\"radio\" name=\"palette\" value=\"blue+green\">blue+green</input> &nbsp;\n";
        hout << "<input  id=\"black+red\"  type=\"radio\" name=\"palette\" value=\"black+red\">black+red</input>\n";
        hout << "<br /><br />\n";
        hout << "<input type=\"button\" value=\"Restart\"  onclick=\"init()\"  />\n";
        hout << "<input type=\"button\" value=\"Continue\" onclick=\"advance()\"  />\n";
        hout << "<input type=\"button\" value=\"Clear\"    onclick=\"clear_canvas()\" />\n";
        hout << "<br />\n";
        hout << "<input type=\"button\" value=\"Start Animation\"    onclick=\"start_animation()\" />\n";
        hout << "<input type=\"button\" value=\"Stop Animation\"    onclick=\"stop_animation()\" />\n";
        hout << "</div>\n";
        hout << "</body>\n";
        hout << "</html>\n";
        hout.close();
        //
        // Create the Javascript ODE solver demonstration code.
        //
        string jsdemo_filename = Name() + "_solverdemo.js";
        ofstream dout;
        dout.open(jsdemo_filename.c_str());  // FIXME! Check for failure.
        dout << "/*\n";
        dout << " *  " << jsdemo_filename << endl;
        dout << " *" << endl;
        PrintVFGENComment(dout," *  ");
        dout << " */" << endl;
        dout << endl;
        // dout << "var c = 0\n";
        // dout << "var dc = 0.08\n";
        dout << endl;
        // dout << "var params = [";
        // for (int i = 0; i < np; ++i)
        //     dout << pardefval_list[i] << ((i == np-1) ? "]\n" : ",") ;
        dout << "var state = [";
        for (int i = 0; i < nv; ++i)
            {
            ex ic = vardefic_list[i];
            for (int j = np-1; j >= 0; j--)
                ic = ic.subs(parname_list[j]==pardefval_list[j]);
            for (int j = nc-1; j >= 0; j--)
                ic = ic.subs(conname_list[j]==convalue_list[j]);
            ic = ic.subs(Pi==M_PI);
            dout << ic << ((i == nv-1) ? "]\n" : ",") ;
            }
        dout << "var t = 0.0\n";
        dout << "var intervalID = \"none\"\n";
        dout << endl;
        dout << "function clear_canvas()\n";
        dout << "    {\n";
        dout << "    var ctx = document.getElementById('canvas').getContext('2d')\n";
        dout << "    ctx.clearRect(0,0,500,500)\n";
        dout << "    }\n";
        dout << endl;
        dout << "function init()\n";
        dout << "    {\n";
        dout << "    clear_canvas()\n";
        dout << "    t = 0.0\n";
        for (int i = 0; i < nv; ++i)
            dout << "    " << varname_list[i] << "0 = parseFloat(document.getElementById('"
                 << varname_list[i] << "0').value)\n";
        for (int i = 0; i < nv; ++i)
            dout << "    state[" << i << "] = " << varname_list[i] << "0\n";
        for (int i = 0; i < nv; ++i)
            dout << "    document.getElementById('" << varname_list[i] << "_t').value = document.getElementById('" << varname_list[i] << "0').value\n";
        dout << "    document.getElementById('t').value = t\n";
        dout << "    }\n";
        dout << endl;
        dout << "function start_animation()\n";
        dout << "    {\n";
        dout << "    if (intervalID == \"none\")\n";
        dout << "        intervalID = setInterval(advance,0)\n";
        dout << "    }\n";
        dout << endl;
        dout << "function stop_animation()\n";
        dout << "    {\n";
        dout << "    if (intervalID != \"none\")\n";
        dout << "        {\n";
        dout << "        clearInterval(intervalID)\n";
        dout << "        intervalID = \"none\"\n";
        dout << "        }\n";
        dout << "    }\n";
        dout << endl;
        dout << "function pw6(x)\n";
        dout << "    {\n";
        dout << "    var x1 = x-6*Math.floor(x/6.0)\n";
        dout << "    if (x1 < 1.0)\n";
        dout << "        return x1\n";
        dout << "    if (x1 < 3.0)\n";
        dout << "        return 1.0\n";
        dout << "    if (x1 < 4.0)\n";
        dout << "        return 1.0-(x1-3.0)\n";
        dout << "    return 0.0\n";
        dout << "    }\n";
        dout << endl;
        dout << "function displayX()\n";
        dout << "    {\n";
        dout << "    return ";
        for (int i = 0; i < nv; ++i)
            dout << "document.getElementById('" << varname_list[i] << "_xscale').value*state[" << i << "] + ";
        dout << "1.0*document.getElementById('xoffset').value\n";
        dout << "    }\n";
        dout << endl;
        dout << "function displayY()\n";
        dout << "    {\n";
        dout << "    return ";
        for (int i = 0; i < nv; ++i)
            dout << "- document.getElementById('" << varname_list[i] << "_yscale').value*state[" << i << "] ";
        dout << "- document.getElementById('yoffset').value\n";        
        dout << "    }\n";
        dout << endl;
        dout << "function advance()\n";
        dout << "    {\n";
        dout << "    var params = [";
        for (int i = 0; i < np; ++i)
            {
            dout << "parseFloat(document.getElementById('" << parname_list[i] << "_value').value)" << ((i == np-1) ? "]\n" : ",");
            }
        dout << "    var numsteps = parseInt(document.getElementById('numsteps').value)\n";
        dout << "    var stoptime = parseFloat(document.getElementById('stoptime').value)\n";
        dout << "    var tol = parseFloat(document.getElementById('tolerance').value)\n";
        dout << "    var hmax = parseFloat(document.getElementById('hmax').value)\n";
        dout << "    var palette_black      = document.getElementById('black').checked\n";
        dout << "    var palette_blue_green = document.getElementById('blue+green').checked\n";
        dout << "    var palette_black_red  = document.getElementById('black+red').checked\n";
        dout << "    var ctx = document.getElementById('canvas').getContext('2d')\n";
        dout << "    ctx.save()\n";
        dout << "    ctx.translate(250,250)\n";
        dout << "    ctx.beginPath()\n";
        dout << "    if (palette_black)\n";
        dout << "        {\n";
        dout << "        var red = 0.0\n";
        dout << "        var green = 0.0\n";
        dout << "        var blue = 0.0\n";
        dout << "        }\n";
        dout << "    else if (palette_blue_green)\n";
        dout << "        {\n";
        dout << "        var red = 0.0\n";
        dout << "        var green = 255*pw6(t)\n";
        dout << "        var blue = 255*pw6(t+7.0)\n";
        dout << "        }\n";
        dout << "    else if (palette_black_red)\n";
        dout << "        {\n";
        dout << "        var red = 255*pw6(t)\n";
        dout << "        var green = 0.0\n";
        dout << "        var blue = 0.0\n";
        dout << "        }\n";
        dout << "    ctx.strokeStyle = \"rgb(\"+parseInt(red)+\",\"+parseInt(green)+\",\"+parseInt(blue)+\")\"\n";
        dout << "    ctx.moveTo(displayX(),displayY())\n";
        dout << "    for (var i = 0; i < numsteps; ++i)\n";
        dout << "        {\n";
        dout << "        var Xderiv = " << Name() << "_derivs" << Order << "(state,params)\n";
        dout << "        var m = 0.0\n";
        dout << "        for (var j = 0; j < " << nv << "; ++j)\n";
        dout << "            {\n";
        dout << "            var a = Xderiv[" << Order-1 << "][j]\n";
        dout << "            m += a*a\n";
        dout << "            }\n";
        dout << "        m = Math.sqrt(m)/" << 1.0*factorial(Order) << endl;
        dout << "        var stepsize = Math.min(hmax,Math.pow(tol/m," << 1.0/(1.0*Order) << "))\n";
        dout << "        if (t + stepsize > stoptime)\n";
        dout << "            stepsize = stoptime - t\n";
        dout << "        state = " << Name() << "_evaltaylor" << Order << "(stepsize,state,Xderiv)\n";
        dout << "        t += stepsize\n";

        dout << "        ctx.lineTo(displayX(),displayY())\n";
        dout << "        if (t >= stoptime)\n";
        dout << "            {\n";
        dout << "            stop_animation()\n";
        dout << "            break\n";
        dout << "            }\n";
        dout << "        }\n";
        dout << "    ctx.stroke()\n";
        dout << "    ctx.restore()\n";
        for (int i = 0; i < nv; ++i)
            dout << "    document.getElementById('" << varname_list[i] << "_t').value = state[" << i << "]\n";
        dout << "    document.getElementById('t').value = t\n";
        dout << "    }\n";

        dout.close();
        }
    }

