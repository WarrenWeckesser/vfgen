//
//  vfgen.cpp -- Multi-format vector field file generator.
//
//  by Warren Weckesser
//
//
//  Copyright (C) 2008-2014 Warren Weckesser
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
#include <vector>
#include <map>
#include <ginac/ginac.h>

#include "vf.h"
#include "vf_help_text.h"

using namespace std;
using namespace GiNaC;

//
// This function overrides the default print method for the Zlags_ function.
// Ginac thinks Zlags_ is a function, but in fact, in the code generated in
// a couple of the VFGEN commands that handle delay equations, Zlags_ is a
// two-dimensional array.  So we want Zlags_(1,2) printed like that, not as
// Zlags_(1.0,2.0).
//

static void Zlags_print(const ex& arg1, const ex& arg2, const print_context& c) {
    c.s << "Zlags_(";
    if (is_a<numeric>(arg1)) {
        c.s << ex_to<numeric>(arg1).to_int();
    }
    else {
        arg1.print(c);
    }
    c.s << ",";
    if (is_a<numeric>(arg2)) {
        c.s << ex_to<numeric>(arg2).to_int();
    }
    else {
        arg2.print(c);
    }
    c.s << ")";
}

//
// Add a function called "delay" to the functions known by ginac.
// (The corresponding macro DECLARE_FUNCTION_2P(delay) is in
// ginac_declare_funcs.h.)
//

REGISTER_FUNCTION(delay, dummy())
//REGISTER_FUNCTION(Zlags_, dummy())
REGISTER_FUNCTION(Zlags_, print_func<print_csrc_float>(Zlags_print).
                          print_func<print_csrc_double>(Zlags_print).
                          print_func<print_python>(Zlags_print) )
REGISTER_FUNCTION(lagvalue, dummy())


#define NAMEWIDTH 9

const char *commands[] = {
        "adolc",
        "auto",
        "check",
        "cvode",
        "dde23",
        "ddebiftool",
        "dde_solver",
        "delay2ode",
        "dstool",
        "evf",
        "gsl",
        "help",
        "javascript",
        "latex",
        "lsoda",
        "matcont",
        "matlab",
        "octave",
        "pddecont",
        "pydstool",
        "pygsl",
        "r",
        "radau5",
        "scilab",
        "scipy",
        "taylor",
        "xml",
        "xpp",
        "end"};
        
map<string,vector<string> > command_options;

int checkcommand(const char *s)
{
    int i;

    i = 0;
    while (strcmp(commands[i],"end") != 0) {
        if (strcmp(s,commands[i]) == 0) {
            return i;
        }
        else {
            i = i + 1;
        }
    }
    return -1;
}
    
void printcommands(ostream &c)
{
    int i = 0;
    c << "    ";
    while (strcmp(commands[i],"end") != 0) {
        if (i > 0) {
            c << ", ";
        }
        if (i > 0 && i % 8 == 0) {
            c << "\n    ";
        }
        c << commands[i];
        i = i + 1;
    }
    c << endl;
}


void printuse()
{
    cerr << "VFGEN (Version:" << VERSION << ")" << endl;
    cerr << "Use: vfgen command  vector-field-file.vf" << endl;
    cerr << "or:  vfgen command:option=value,...,option=value vector-field-file.vf" << endl;
    cerr << "or:  vfgen help command\n";
    cerr << "where command is one of:\n";
    printcommands(cerr);
}

int help(char *command)
{
    int m;
    m = checkcommand(command);
    if (m < 0) {
        cout << "VFGEN error: \"" << command << "\" is not a valid command.\n";
        cout << "VFGEN commands are:\n";
        printcommands(cout);
        return -1;
    }
    if (strcmp(command, "adolc") == 0) {
        cout << help_adolc;
    }
    else if (strcmp(command, "auto") == 0) {
        cout << help_auto;
    }
    else if (strcmp(command, "check") == 0) {
        cout << help_check;
    }
    else if (strcmp(command, "cvode") == 0) {
        cout << help_cvode;
    }
    else if (strcmp(command, "dde23") == 0) {
        cout << help_dde23;
    }
    else if (strcmp(command, "ddebiftool") == 0) {
        cout << help_ddebiftool;
    }
    else if (strcmp(command, "dde_solver") == 0) {
        cout << help_dde_solver;
    }
    else if (strcmp(command, "delay2ode") == 0) {
        cout << help_delay2ode;
    }
    else if (strcmp(command, "dstool") == 0) {
        cout << help_dstool;
    }
    else if (strcmp(command, "evf") == 0) {
        cout << help_evf;
    }
    else if (strcmp(command, "gsl") == 0) {
        cout << help_gsl;
    }
    else if (strcmp(command, "help") == 0) {
        cout << help_help;
    }
    else if (strcmp(command, "javascript") == 0) {
        cout << help_javascript;
    }
    else if (strcmp(command, "latex") == 0) {
        cout << help_latex;
    }
    else if (strcmp(command, "lsoda") == 0) {
        cout << help_lsoda;
    }
    else if (strcmp(command, "matcont") == 0) {
        cout << help_matcont;
    }
    else if (strcmp(command, "matlab") == 0) {
        cout << help_matlab;
    }
    else if (strcmp(command, "octave") == 0) {
        cout << help_octave;
    }
    else if (strcmp(command, "pddecont") == 0) {
        cout << help_pddecont;
    }
    else if (strcmp(command, "pydstool") == 0) {
        cout << help_pydstool;
    }
    else if (strcmp(command, "pygsl") == 0) {
        cout << help_pygsl;
    }
    else if (strcmp(command, "r") == 0) {
        cout << help_r;
    }
    else if (strcmp(command, "radau5") == 0) {
        cout << help_radau5;
    }
    else if (strcmp(command, "scilab") == 0) {
        cout << help_scilab;
    }
    else if (strcmp(command,"scipy") == 0) {
        cout << help_scipy;
    }
    else if (strcmp(command,"taylor") == 0) {
        cout << help_taylor;
    }
    else if (strcmp(command,"xpp") == 0) {
        cout << help_xpp;
    }
    else {
        cout << "Sorry, help for \"" << command << "\" is not available yet!" << endl;
    }
    return 0;
}

///////////////////////////////////////////////////////////
//  main
///////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    VectorField vf;
    string s, extrastr;

    if (argc != 3) {
        printuse();
        exit(-1);
    }

    //
    //  Allowed command options:
    //
    command_options["auto"].push_back("lang");
    command_options["cvode"].push_back("demo");
    command_options["cvode"].push_back("func");
    command_options["cvode"].push_back("version");
    command_options["dde23"].push_back("demo");
    command_options["dde23"].push_back("parstyle");
    command_options["ddebiftool"].push_back("path");
    command_options["dde_solver"].push_back("demo");
    command_options["delay2ode"].push_back("N");
    command_options["delay2ode"].push_back("p");
    command_options["evf"].push_back("par");
    command_options["gsl"].push_back("demo");
    command_options["gsl"].push_back("func");
    command_options["javascript"].push_back("order");
    command_options["javascript"].push_back("demo");
    command_options["lsoda"].push_back("demo");
    command_options["lsoda"].push_back("func");
    command_options["lsoda"].push_back("parstyle");
    command_options["matlab"].push_back("demo");
    command_options["matlab"].push_back("evf");
    command_options["matlab"].push_back("func");
    command_options["matlab"].push_back("parstyle");
    command_options["octave"].push_back("demo");
    command_options["octave"].push_back("func");
    command_options["octave"].push_back("parstyle");
    command_options["pydstool"].push_back("demo");
    command_options["pygsl"].push_back("demo");
    command_options["pygsl"].push_back("func");
    command_options["r"].push_back("demo");
    command_options["r"].push_back("func");    
    command_options["radau5"].push_back("demo");
    command_options["scilab"].push_back("demo");
    command_options["scilab"].push_back("func");
    command_options["scilab"].push_back("parstyle");
    command_options["scipy"].push_back("demo");
    command_options["scipy"].push_back("func");
    command_options["taylor"].push_back("order");
    command_options["xpp"].push_back("extra");
    
    string commandstr(argv[1]);
    //
    // Check for the help command.
    //
    if (commandstr == "help") {
        exit(help(argv[2]));
    }

    //
    //  Check for any options appended to the command.
    //
    map<string,string> options;

    string::size_type loccolon = commandstr.find(":",0);
    if (loccolon != string::npos) {
        // extrastr holds all the options given after the :
        extrastr = commandstr.substr(loccolon+1,commandstr.length()-loccolon-1);
        commandstr.erase(loccolon,commandstr.length()-loccolon);
        // cerr << "Options \"" << extrastr << "\"" << endl;
        string::size_type pos = 0;
        do {
            string::size_type locsep = extrastr.find(",",pos);
            string current_option, option_name, option_value;
            if (locsep == string::npos) {
                locsep = extrastr.length();
            }
            current_option = extrastr.substr(pos,locsep-pos);
            // cerr << "current_option = \"" << current_option << "\"\n";
            string::size_type loceq = current_option.find("=",0);
            if (loceq == string::npos) {
                // No "=" given in the option.
                option_name = current_option;
                option_value = "";
            }
            else {
                if (loceq == 0) {
                    // The option was "=something"; not valid
                    printuse();
                    exit(-1);
                }
                option_name = current_option.substr(0,loceq);
                option_value = current_option.substr(loceq+1,current_option.length()-loceq);
            }
            options[option_name] = option_value;
            pos = locsep+1;
        } while (pos < extrastr.length());
    }
        
        
    int command = checkcommand(commandstr.c_str());
    if (command == -1) {
        cout << "vfgen: unknown command: " << commandstr << endl;
        printuse();
        exit(-1);
    }

    //
    //  Check that any options given are known options.
    //  (This doesn't check that the value of the option is valid; it just
    //  checks that the name of the option is one of the  allowed options for
    //  the given command.)
    //
    bool bad_opt = false;        
    map<string,string>::const_iterator opt;
    for (opt = options.begin(); opt != options.end(); ++opt) {
        string optstr = opt->first;
        // cerr << "Option: " << optstr;
        // if (opt->second != "")
        //     cerr << "=" << opt->second;
        // cerr << endl;
        bool validopt = false;
        vector<string>::const_iterator w;
        for (w = command_options[commandstr].begin(); w != command_options[commandstr].end(); ++w) {
            if (optstr == *w) {
                validopt = true;
                break;
            }
        }
        if (!validopt) {
            cerr << "Errror: \"" << optstr << "\" is not a valid option for the " << commandstr << " command.\n";
            bad_opt = true;
        }
    }
    if (bad_opt) {
        printuse();
        exit(-1);
    }

    //
    //  Read the vector field file.  This just puts the strings into the
    //  appropriate fields.  It doesn't do any symbolic processing.
    //
    vf.ReadXML(argv[2]);

    //
    //  Process the strings to create the GiNaC symbolic expressions in the object.
    //
    int pserr = vf.ProcessSymbols();
    if (pserr == -1) {
        exit(-1);
    }

    //
    // Call the appropriate output function based on the first
    // command line argument.
    //
    if (commandstr == "check") {
        vf.Print();
    }
    else if (commandstr == "xpp") {
        vf.PrintXPP(options);
    }
    else if (commandstr == "xml") {
        vf.PrintXML("xml");
    }
    else if (commandstr == "delay2ode") {
        if (vf.IsDelay == true) {
            vf.PrintDelay2ODE(options);
        }
        else {
            cerr << "This system is not a delay equation.\n";
        }
    }
    else if (commandstr == "dde23") {
        if (vf.IsDelay == true) {
            if (vf.HasNonconstantDelay) {
                cerr << "This system has nonconstant delays.  DDE23 can only be used with constant delays.\n";
            }
            else {
                vf.PrintDDE23(options);
            }
        }
        else {
            cerr << "This system is not a delay equation.\n";
        }
    }
    else if (commandstr == "ddebiftool") {
        if (vf.IsDelay == true) {
            if (vf.HasNonconstantDelay) {
                cerr << "This system has nonconstant delays. VFGEN does not (yet) generate DDE-BIFTOOL files for such systems.\n";
            }
            else {
                vf.PrintDDEBIFTOOL(options);
            }
        }
        else {
            cerr << "This system is not a delay equation.\n";
        }
    }
    else if (commandstr == "pddecont") {
        if (vf.HasNonconstantDelay) {
            cerr << "This system has nonconstant delays. PDDE-CONT is for systems with constant delays.\n";
        }
        else {
            vf.PrintPDDECONT(options);
        }
    }
    else if (commandstr == "dde_solver") {
        if (vf.IsDelay == true) {
            vf.PrintDDE_SOLVER(options);
        }
        else {
            cerr << "This system is not a delay equation.\n";
        }
    }
    else if (commandstr == "matlab") {
        if (vf.IsDelay == false) {
            vf.PrintMATLAB(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "octave") {
        if (vf.IsDelay == false) {
            vf.PrintOctave(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "scilab") {
        if (vf.IsDelay == false) {
            vf.PrintScilab(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "matcont") {
        if (vf.IsDelay == false) {
            vf.PrintMATCONT(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "dstool") {
        if (vf.IsDelay == false) {
            vf.PrintDSTool();
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "auto") {
        if (vf.IsDelay == false) {
            vf.PrintAUTO(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "gsl") {
        if (vf.IsDelay == false) {
            vf.PrintGSL(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "cvode") {
        if (vf.IsDelay == false) {
            vf.PrintCVODE(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "adolc") {
        if (vf.IsDelay == false) {
            vf.PrintADOLC(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "r") {
        if (vf.IsDelay) {
            vf.PrintRdede(options);
        }
        else {
            vf.PrintRode(options);
        }
    }
    else if (commandstr == "radau5") {
        if (vf.IsDelay == false) {
            vf.PrintRadau5(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "lsoda") {
        if (vf.IsDelay == false) {
            vf.PrintLSODA(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "scipy") {
        if (vf.IsDelay == false) {
            vf.PrintSciPy(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "pydstool") {
        if (vf.IsDelay == false) {
            vf.PrintPyDSTool(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "taylor") {
        if (vf.IsDelay == false) {
            vf.PrintTaylor(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "latex") {
        vf.PrintLatex(options);
    }
    else if (commandstr == "pygsl") {
        if (vf.IsDelay == false) {
            vf.PrintPyGSL(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "evf") {
        if (vf.IsDelay == false) {
            vf.PrintEVF(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else if (commandstr == "javascript") {
        if (vf.IsDelay == false) {
            vf.PrintJavascript(options);
        }
        else {
            cerr << "Delay equations can not be handled by the " << commandstr << " command.\n";
        }
    }
    else {
        // This should not happen!!!
        cerr << "vfgen: Unknown command: " << argv[1] << endl;
        printuse();
        exit(-1);
    }
    return(0);
}  // end main()
