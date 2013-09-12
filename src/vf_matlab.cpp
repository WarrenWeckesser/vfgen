
//
//  vf_matlab.cpp
//
//  This file defines the VectorField::PrintMATLAB method.
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
#include <string>
#include <ginac/ginac.h> 

#include "vf.h"
#include "vf_utils.h"

using namespace std;
using namespace GiNaC;

//
// PrintMATLAB -- The MATLAB Code Generator.
//

void VectorField::PrintMATLAB(map<string,string> options)
    {
    int nc, np, nv, na, nf;

    nc = conname_list.nops();
    nv = varname_list.nops();
    np = parname_list.nops();
    na = exprname_list.nops();
    nf = funcname_list.nops();

    //
    //  Create the vector field function.
    //
    string vf_filename = Name()+"_vf.m";
    ofstream fout;
    fout.open(vf_filename.c_str());
    fout << left;

    fout << "%" << endl;
    fout << "% " << vf_filename << endl;
    fout << "%" << endl;
    fout << "% MATLAB vector field function for: " << Name() << endl;
    fout << "%" << endl;
    PrintVFGENComment(fout,"% ");
    fout << "%" << endl;
    fout << "%" << endl;
    fout << "function vf_ = " << Name() << "_vf(" << IndependentVariable << ",x_";
    if (np > 0)
        {
        fout << ",";
        if (options["parstyle"] == "list")
            PrintNameList(fout,parname_list);
        else
            fout << "p_";
        }
    fout << ")" << endl;
    if (HasPi)
        {
        fout << "    Pi = pi;\n";
        }
    for (int i = 0; i < nc; ++i)
        {
        fout << "    " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    GetFromVector(fout,"    ",varname_list,"x_","()",1,";");
    if (options["parstyle"] != "list")
        GetFromVector(fout,"    ",parname_list,"p_","()",1,";");
    for (int i = 0; i < na; ++i)
        {
        fout << "    " << exprname_list[i] << " = " << exprformula_list[i] << ";" << endl;
        }
    fout << "    vf_ = zeros(" << nv << ",1);" << endl;
    for (int i = 0; i < nv; ++i)
        {
        fout << "    vf_(" << (i+1) << ")" << " = " << varvecfield_list[i] << ";" << endl;
     
   }
    fout << endl;
    fout.close();
    //
    //  Create the Jacobian function.
    //
    string jac_filename = Name()+"_jac.m";
    fout.open(jac_filename.c_str());
    fout << left;

    fout << "%" << endl;
    fout << "% " << jac_filename << endl;
    fout << "%" << endl;
    fout << "% This MATLAB function computes the Jacobian of the vector field" << endl;
    fout << "% defined in " << vf_filename << "." << endl;
    fout << "%" << endl;
    PrintVFGENComment(fout,"% ");
    fout << "%" << endl;
    fout << "%" << endl;
    fout << "function jac_ = " << Name() << "_jac(" << IndependentVariable << ",x_";
    if (np > 0)
        {
        fout << ",";
        if (options["parstyle"] == "list")
            PrintNameList(fout,parname_list);
        else
            fout << "p_";
        }
    fout << ")" << endl;
    if (HasPi)
        {
        fout << "    Pi = pi;\n";
        }
    for (int i = 0; i < nc; ++i)
        {
        fout << "    " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    GetFromVector(fout,"    ",varname_list,"x_","()",1,";");
    if (options["parstyle"] != "list")
        GetFromVector(fout,"    ",parname_list,"p_","()",1,";");
    fout << "    jac_ = zeros(" << nv << "," << nv << ");" << endl;
    for (int i = 0; i < nv; ++i)
        {
        ex f = iterated_subs(varvecfield_list[i],expreqn_list);
        for (int j = 0; j < nv; ++j)
            {
            symbol v = ex_to<symbol>(varname_list[j]);
            ex df = f.diff(v);
            if (df != 0)
                fout << "    jac_(" << (i+1) << "," << (j+1) << ")" << " = " << df << ";" << endl;
            }
        }
    fout << endl;
    fout.close();

    //
    //  Create the function that computes the Jacobian with respect
    //  to the parameters.
    //
    string jacp_filename = Name()+"_jacp.m";
    fout.open(jacp_filename.c_str());
    fout << left;

    fout << "%" << endl;
    fout << "% " << jacp_filename << endl;
    fout << "%" << endl;
    fout << "% This MATLAB function computes the Jacobian with respect to the parameters" << endl;
    fout << "% of the vector field defined in " << vf_filename << "." << endl;
    fout << "%" << endl;
    PrintVFGENComment(fout,"% ");
    fout << "%" << endl;
    fout << "%" << endl;
    fout << "function jacp_ = " << Name() << "_jacp(" << IndependentVariable << ",x_";
    if (np > 0)
        {
        fout << ",";
        if (options["parstyle"] == "list")
            PrintNameList(fout,parname_list);
        else
            fout << "p_";
        }
    fout << ")" << endl;
    if (HasPi)
        {
        fout << "    Pi = pi;\n";
        }
    for (int i = 0; i < nc; ++i)
        {
        fout << "    " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    GetFromVector(fout,"    ",varname_list,"x_","()",1,";");
    if (options["parstyle"] != "list")
        GetFromVector(fout,"    ",parname_list,"p_","()",1,";");
    fout << "    jacp_ = zeros(" << nv << "," << np << ");" << endl;
    for (int i = 0; i < nv; ++i)
        {
        ex f = iterated_subs(varvecfield_list[i],expreqn_list);
        for (int j = 0; j < np; ++j)
            {
            symbol p = ex_to<symbol>(parname_list[j]);
            ex df = f.diff(p);
            if (df != 0)
                fout << "    jacp_(" << (i+1) << "," << (j+1) << ")" << " = " << df << ";" << endl;
            }
        }
    fout << endl;
    fout.close();

    //
    //  Create the Hessian function.
    //
    string hess_filename = Name()+"_hess.m";
    fout.open(hess_filename.c_str());
    fout << left;

    fout << "%" << endl;
    fout << "% " << hess_filename << endl;
    fout << "%" << endl;
    fout << "% This MATLAB function computes the Hessians of the vector field" << endl;
    fout << "% defined in " << vf_filename << "." << endl;
    fout << "%" << endl;
    fout << "% hess_(n,i,j) is the second partial derivative of the n-th component" << endl;
    fout << "% of the vector field, taken with respect to the i-th and j-th variables." << endl;
    fout << "%" << endl;
    PrintVFGENComment(fout,"% ");
    fout << "%" << endl;
    fout << "%" << endl;
    fout << "function hess_ = " << Name() << "_hess(" << IndependentVariable << ",x_";
    if (np > 0)
        {
        fout << ",";
        if (options["parstyle"] == "list")
            PrintNameList(fout,parname_list);
        else
            fout << "p_";
        }
    fout << ")" << endl;
    if (HasPi)
        {
        fout << "    Pi = pi;\n";
        }
    for (int i = 0; i < nc; ++i)
        {
        fout << "    " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    GetFromVector(fout,"    ",varname_list,"x_","()",1,";");
    if (options["parstyle"] != "list")
        GetFromVector(fout,"    ",parname_list,"p_","()",1,";");
    fout << endl;
    fout << "    hess_ = zeros(" << nv << "," << nv << "," << nv << ");" << endl;
    for (int n = 0; n < nv; ++n)
        {
        fout << endl;
        // Get the n-th component of the vector field.
        ex f = iterated_subs(varvecfield_list[n],expreqn_list);
        for (int i = 0; i < nv; ++i)
            {
            // Get the i-th variable
            symbol x_i = ex_to<symbol>(varname_list[i]);
            // Differentiate f once
            ex df_i = f.diff(x_i);
            for (int j = i; j < nv; ++j)
                {
                // Get the j-th variable
                symbol x_j = ex_to<symbol>(varname_list[j]);
                // Differentiate again
                ex df_ij = df_i.diff(x_j);
                if (df_ij != 0)
                    {
                    fout << "    hess_(" << (n+1) << "," << (i+1) << "," << (j+1) << ")" << " = " << df_ij << ";" << endl;
                    if (j > i)
                        fout << "    hess_(" << (n+1) << "," << (j+1) << "," << (i+1) << ") = hess_(" << (n+1) << "," << (i+1) << "," << (j+1) << ");" << endl;
                    }
                }
            }
        }
    fout << endl;
    fout.close();

    //
    //  Create the third derivatives function function.
    //
    string d3_filename = Name()+"_der3.m";
    fout.open(d3_filename.c_str());
    fout << left;

    fout << "%" << endl;
    fout << "% " << d3_filename << endl;
    fout << "%" << endl;
    fout << "% This MATLAB function computes the third derivatives of the vector field" << endl;
    fout << "% defined in " << vf_filename << "." << endl;
    fout << "%" << endl;
    fout << "% der3_(n,i,j,k) is the third partial derivative of the n-th component" << endl;
    fout << "% of the vector field, taken with respect to the i-th, j-th and k-th variables." << endl;
    fout << "%" << endl;
    PrintVFGENComment(fout,"% ");
    fout << "%" << endl;
    fout << "%" << endl;
    fout << "function der3_ = " << Name() << "_der3(" << IndependentVariable << ",x_";
    if (np > 0)
        {
        fout << ",";
        if (options["parstyle"] == "list")
            PrintNameList(fout,parname_list);
        else
            fout << "p_";
        }
    fout << ")" << endl;
    if (HasPi)
        {
        fout << "    Pi = pi;\n";
        }
    for (int i = 0; i < nc; ++i)
        {
        fout << "    " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
        }
    GetFromVector(fout,"    ",varname_list,"x_","()",1,";");
    if (options["parstyle"] != "list")
        GetFromVector(fout,"    ",parname_list,"p_","()",1,";");
    fout << endl;
    fout << "    der3_ = zeros(" << nv << "," << nv << "," << nv << "," << nv << ");" << endl;
    for (int n = 0; n < nv; ++n)
        {
        fout << endl;
        // Get the n-th component of the vector field.
        ex f = iterated_subs(varvecfield_list[n],expreqn_list);
        for (int i = 0; i < nv; ++i)
            {
            // Get the i-th variable
            symbol x_i = ex_to<symbol>(varname_list[i]);
            // Differentiate f once
            ex df_i = f.diff(x_i);
            for (int j = i; j < nv; ++j)
                {
                // Get the j-th variable
                symbol x_j = ex_to<symbol>(varname_list[j]);
                // Differentiate again
                ex df_ij = df_i.diff(x_j); 
                for (int k = j; k < nv; ++k)
                    {
                    // Get the k-th variable
                    symbol x_k = ex_to<symbol>(varname_list[k]);
                    // Differentiate
                    ex df_ijk = df_ij.diff(x_k);
                    if (df_ijk != 0)
                        {
                        fout << "    der3_(" << (n+1) << "," << (i+1) << "," << (j+1) << "," << (k+1) << ") = " << df_ijk << ";" << endl;
                        if (j < k)
                            fout << "    der3_(" << (n+1) << "," << (i+1) << "," << (k+1) << "," << (j+1) << ") = der3_(" << (n+1) << "," << (i+1) << "," << (j+1) << "," << (k+1) << ");" << endl;
                        if (i < k)
                            fout << "    der3_(" << (n+1) << "," << (k+1) << "," << (j+1) << "," << (i+1) << ") = der3_(" << (n+1) << "," << (i+1) << "," << (j+1) << "," << (k+1) << ");" << endl;
                        if (i < j)
                            fout << "    der3_(" << (n+1) << "," << (j+1) << "," << (i+1) << "," << (k+1) << ") = der3_(" << (n+1) << "," << (i+1) << "," << (j+1) << "," << (k+1) << ");" << endl;
                        }
                    }
                }
            }
        }
    fout << endl;
    fout.close();

    if (np > 0)
        {
		//
		//  Create the function that computes the derivatives with respect to one
		//  variable and one parameter.
		//
		string hessp_filename = Name()+"_hessp.m";
		fout.open(hessp_filename.c_str());
		fout << left;
	
		fout << "%" << endl;
		fout << "% " << hessp_filename << endl;
		fout << "%" << endl;
		fout << "% This MATLAB function computes the derivatives of the vector field" << endl;
		fout << "% defined in " << vf_filename << "." << endl;
		fout << "%" << endl;
		fout << "% hessp_(n,i,j) is the second partial derivative of the n-th component" << endl;
		fout << "% of the vector field, taken with respect to the i-th variable" << endl;
		fout << "% and the j-th parameter." << endl;
		fout << "%" << endl;
		PrintVFGENComment(fout,"% ");
		fout << "%" << endl;
		fout << "%" << endl;
		fout << "function hessp_ = " << Name() << "_hessp(" << IndependentVariable << ",x_,";
		if (options["parstyle"] == "list")
			PrintNameList(fout,parname_list);
		else
			fout << "p_";
		fout << ")" << endl;
		if (HasPi)
			{
			fout << "    Pi = pi;\n";
			}
		for (int i = 0; i < nc; ++i)
			{
			fout << "    " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
			}
		GetFromVector(fout,"    ",varname_list,"x_","()",1,";");
		if (options["parstyle"] != "list")
			GetFromVector(fout,"    ",parname_list,"p_","()",1,";");
		fout << endl;
		fout << "    hessp_ = zeros(" << nv << "," << nv << "," << np << ");" << endl;
		for (int n = 0; n < nv; ++n)
			{
			fout << endl;
			// Get the n-th component of the vector field.
			ex f = iterated_subs(varvecfield_list[n],expreqn_list);
			for (int i = 0; i < nv; ++i)
				{
				// Get the i-th variable
				symbol x_i = ex_to<symbol>(varname_list[i]);
				// Differentiate f once
				ex df_i = f.diff(x_i);
				for (int j = 0; j < np; ++j)
					{
					// Get the j-th parameter
					symbol p_j = ex_to<symbol>(parname_list[j]);
					// Differentiate again
					ex df_ij = df_i.diff(p_j);
					if (df_ij != 0) 
						fout << "    hessp_(" << (n+1) << "," << (i+1) << "," << (j+1) << ")" << " = " << df_ij << ";" << endl;
					}
				}
			}
		fout << endl;
		fout.close();
        }


    if (options["func"] == "yes")
        {
        //
        //  Create the user-defined functions.
        //
        for (int n = 0; n < nf; ++n)
            {
            symbol fn = ex_to<symbol>(funcname_list[n]);
            string funcname = fn.get_name();
            string filename = Name() + "_" + funcname+".m";
            ofstream fout;
            fout.open(filename.c_str());
            fout << left;

            fout << "%" << endl;
            fout << "% " << filename << endl;
            fout << "%" << endl;
            fout << "% MATLAB user function for the vector field: " << Name() << endl;
            fout << "%" << endl;
            PrintVFGENComment(fout,"% ");
            fout << "%" << endl;
            fout << "%" << endl;
            fout << "function r_ = " << Name() << "_" << funcname << "(" << IndependentVariable << ",x_";
            if (np > 0)
                {
                fout << ",";
                if (options["parstyle"] == "list")
                    PrintNameList(fout,parname_list);
                else
                    fout << "p_";
                }
            fout << ")" << endl;
            if (HasPi)
                {
                fout << "    Pi = pi;\n";
                }
            for (int i = 0; i < nc; ++i)
                {
                fout << "    " << conname_list[i] << " = " << convalue_list[i] << ";" << endl;
                }
            GetFromVector(fout,"    ",varname_list,"x_","()",1,";");
            if (options["parstyle"] != "list")
                GetFromVector(fout,"    ",parname_list,"p_","()",1,";");
            for (int i = 0; i < na; ++i)
                {
                fout << "    " << exprname_list[i] << " = " << exprformula_list[i] << ";" << endl;
                }
            fout << "    r_ = " << funcformula_list[n] << ";" << endl;
            fout.close();
            }
        }

    if (options["evf"] == "yes")
        {
        //
        //  Create the extended vector field.  This has the vector field
        //  along with variational equations.
        //
        string evf_filename = Name()+"_evf.m";
        fout.open(evf_filename.c_str());
        fout << left;
 
        fout << "%" << endl;
        fout << "% " << evf_filename << endl;
        fout << "%" << endl;
        fout << "% MATLAB vector field and first variation for: " << Name() << endl;
        fout << "%" << endl;
        PrintVFGENComment(fout,"% ");
        fout << "%" << endl;
        fout << "%" << endl;
        fout << "function evf_ = " << Name() << "_evf(" << IndependentVariable << ",w_";
        if (np > 0)
            {
            fout << ",";
            if (options["parstyle"] == "list")
                PrintNameList(fout,parname_list);
            else
                fout << "p_";
            }
        fout << ")" << endl;
        fout << "    x_ = w_(1:" << nv << ");" << endl;
        fout << "    v_ = w_(" << (nv+1) << ":" << (2*nv) << ");" << endl;
        fout << "    evf_ = zeros(" << (2*nv) << ",1);" << endl;
        fout << "    evf_(1:" << nv << ") = " << Name() << "_vf(t,x_";
        if (np > 0)
            {
            fout << ",";
            if (options["parstyle"] == "list")
                PrintNameList(fout,parname_list);
            else
                fout << "p_";
            }
        fout << ");" << endl;
        fout << "    evf_(" << (nv+1) << ":" << (2*nv) << ") = " << Name() << "_jac(t,x_";
        if (np > 0)
            {
            if (options["parstyle"] == "list")
                PrintNameList(fout,parname_list);
            else
                fout << "p_";
            }
        fout << ") * v_;" << endl;
        fout.close();
        }

    if (options["demo"] == "yes")
        {
        //
        // Create the demo function.
        //
        int ht = 18;
        string filename = Name() + "_demo.m";
        ofstream fout;
        fout.open(filename.c_str());
        fout << left;

        fout << "%" << endl;
        fout << "% " << filename << endl;
        fout << "%" << endl;
        fout << "% MATLAB demo function for the vector field: " << Name() << endl;
        fout << "%" << endl;
        PrintVFGENComment(fout,"% ");
        fout << "%" << endl;
        fout << "%" << endl;
        fout << endl;
        int yinc = -(ht+2);
        int y = (nv+np+4)*(-yinc);
        int labelwidth = 100;
        fout << "function " << Name() << "_demo\n";
        fout << "    figure;\n";
        fout << "    clf\n";
        fout << "    set(gcf,'Position',[0 0 250 " << y << "]);\n";
        fout << "    v = [];\n";
        // Variables -- initial conditions
        symbol pi("pi");
        for (unsigned n = 0; n < varname_list.nops(); ++n)
            {
            y = y + yinc;
            fout << "    uicontrol('Style','text','Position',[10 " << y << " " << labelwidth << " " << ht << "],'String','" << varname_list[n] << "(0)');\n";
            fout << "    ui = uicontrol('Style','edit','Position',[" << labelwidth+15 << " " << y << " 100 " << ht << "],'String','" << vardefic_list[n].subs(Pi==pi) << "');\n";
            fout << "    v = [v; ui];\n";
            }
        // Parameters
        for (unsigned n = 0; n < parname_list.nops(); ++n)
            {
            y = y + yinc;
            fout << "    uicontrol('Style','text','Position',[10 " << y << " " << labelwidth << " " << ht << "],'String','" << parname_list[n] << "');\n";
            fout << "    ui = uicontrol('Style','edit','Position',[" << labelwidth+15 << " " << y << " 100 " << ht << "],'String','" << pardefval_list[n].subs(Pi==pi) << "');\n";
            fout << "    v = [v; ui];\n";
            }
        y = y + yinc;
        // Stop Time
        fout << "    uicontrol('Style','text','Position',[10 " << y << " " << labelwidth << " " << ht << "],'String','Stop Time');\n";
        fout << "    ui = uicontrol('Style','edit','Position',[" << labelwidth+15 << " " << y << " 100 " << ht << "],'String','10');\n";
        fout << "    v = [v; ui];\n";
        y = y + yinc;
        // Separate Axes checkbox
        fout << "    ui = uicontrol('Style','checkbox','Position',[10 " << y << " 200 " << ht << "],'String','Separate Axes');\n";
        fout << "    v = [v; ui];\n";
        y = y + yinc;
        // Go button
        fout << "    uicontrol('Style','pushbutton','Position',[10 " << y << " 40 " << ht << "],'String','Go','Callback',@go_cb,'UserData',v);\n";

        fout << "end\n";
        // Write the callback function for the "Go" button.
        fout << endl;
        fout << "function go_cb(arg1,arg2)\n";
        fout << "    v = get(arg1,'UserData');\n";
        fout << "    ic = zeros(size(v,1)-2,1);\n";
        fout << "    for k = 1:size(v,1)-2,\n";
        // fout << "        ic(k) = str2double(get(v(k),'String'));\n";
        fout << "        ic(k) = eval(get(v(k),'String'));\n";
        fout << "        if (isnan(ic(k)))\n";
        fout << "            ic(k) = 0.0;\n";
        fout << "        end;\n";
        fout << "    end;\n";
        // fout << "    disp(ic);\n";
        fout << "    stoptime = str2double(get(v(end-1),'String'));\n";
        fout << "    if (isnan(stoptime))\n";
        fout << "        stoptime = 0.0;\n";
        fout << "    end;\n";
        fout << "    abstol = 1e-9;\n";
        fout << "    reltol = 1e-7;\n";
        fout << "    opts = odeset('AbsTol',abstol,'RelTol',reltol,'Jacobian',@" << Name() << "_jac" << ");\n";
        fout << "    % Change ode45 to ode15s for stiff differential equations.\n";
        fout << "    [t,z_] = ode45(@" << Name() << "_vf,[0 stoptime],ic(1:" << varname_list.nops() << "),opts";
        if (np > 0)
            {
            if (options["parstyle"] != "list")
                fout << ",ic(" << varname_list.nops()+1 << ":end)";
            else
                {
                for (int n = 0; n < np; ++n)
                    fout << ",ic(" << varname_list.nops()+1+n << ")";
                }
            }
        fout << ");\n";
        fout << "    figure;\n";
        fout << "    clf;\n";
        fout << "    a = get(v(end),'Value');\n";
        fout << "    if (a == 0),\n";
        fout << "        plot(t,z_);\n";
        fout << "        xlabel('t');\n";
        fout << "        legend(";
        for (unsigned n = 0; n < varname_list.nops(); ++n)
            {
            if (n > 0)
                fout << ",";
            fout << "'" << varname_list[n] << "'";
            }
        fout << ");\n";
        fout << "    grid on\n";
        fout << "    else\n";
        for (unsigned n = 0; n < varname_list.nops(); ++n)
            {
            fout << "        subplot(" << varname_list.nops() << ",1," << n+1 << ");\n";
            fout << "        plot(t,z_(:," << n+1 << "))\n";
            // fout << "        legend('" << varname_list[n] << "')\n";
            fout << "        xlabel('t');\n";
            fout << "        ylabel('" << varname_list[n] << "')\n";
            fout << "        grid on\n";
            }
        fout << "    end\n";
        fout << "end\n";
        fout.close();
        }
    }

