
//
//  vf_xml.cpp
//
//  This file defines the VectorField::PrintXML and VectorField::ReadXML methods.
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
#include <mxml.h>
#include <ginac/ginac.h> 

#include "vf.h"
#include "vf_utils.h"

using namespace std;
using namespace GiNaC;

//
// PrintXML
//

void VectorField::PrintXML(string cmdstr)
{
    cout << "<?xml version=\"1.0\" ?>" << endl;
    cout << "<!--\n";
    cout << "  This file was generated by the program VFGEN (Version: " << VERSION << ") with the " << cmdstr << " command.\n";
    cout << "  " << DateTimeMsg() << endl;
    cout << "-->\n";
    cout << "<VectorField";
    cout << " Name=\"" << Name() << "\"";
    if (Description() != "") {
        cout << " Description=\"" << Description() << "\"";
    }
    cout << " IndependentVariable=\"" << IndependentVariable << "\"";
    cout << ">" << endl;

    vector<Constant *>::iterator c;
    for (c = Constants.begin(); c != Constants.end(); ++c) {
        cout << "<Constant";
        cout << " Name=\"" << (*c)->Name() << "\"";
        if ((*c)->Description() != "") {
            cout << " Description=\"" << (*c)->Description() << "\"";
        }
        cout << " Value=\"" << (*c)->Value() << "\"";
        if ((*c)->Latex() != "") {
            cout << " Latex=\"" << (*c)->Latex() << "\"";
        }
        cout << "/>" << endl;
    }

    vector<Parameter *>::iterator p;
    for (p = Parameters.begin(); p != Parameters.end(); ++p) {
        cout << "<Parameter";
        cout << " Name=\"" << (*p)->Name() << "\"";
        if ((*p)->Description() != "") {
            cout << " Description=\"" << (*p)->Description() << "\"";
        }
        cout << " DefaultValue=\"" << (*p)->DefaultValue() << "\"";
        if ((*p)->Latex() != "") {
            cout << " Latex=\"" << (*p)->Latex() << "\"";
        }
        cout << "/>" << endl;
    }

    vector<Expression *>::iterator e;
    for (e = Expressions.begin(); e != Expressions.end(); ++e) {
        cout << "<Expression";
        cout << " Name=\"" << (*e)->Name() << "\"";
        if ((*e)->Description() != "") {
            cout << " Description=\"" << (*e)->Description() << "\"";
        }
        cout << " Formula=\"" << (*e)->Formula() << "\"";
        if ((*e)->Latex() != "") {
            cout << " Latex=\"" << (*e)->Latex() << "\"";
        }
        cout << "/>" << endl;
    }

    vector<StateVariable *>::iterator sv;
    for (sv = StateVariables.begin(); sv != StateVariables.end(); ++sv) {
        cout << "<StateVariable";
        cout << " Name=\"" << (*sv)->Name() << "\"";
        if ((*sv)->Description() != "") {
            cout << " Description=\"" << (*sv)->Description() << "\"";
        }
        cout << " Formula=\"" << (*sv)->Formula() << "\"";
        if ((*sv)->PeriodicFrom() != "") {
            cout << " PeriodFrom=\"" << (*sv)->PeriodicFrom() << "\"";
        }
        if ((*sv)->PeriodicTo() != "") {
            cout << " PeriodTo=\"" << (*sv)->PeriodicTo() << "\"";
        }
        if ((*sv)->DefaultInitialCondition() != "") {
            cout << " DefaultInitialCondition=\"" << (*sv)->DefaultInitialCondition() << "\"";
        }
        if ((*sv)->DefaultHistory() != "") {
            cout << " DefaultHistory=\"" << (*sv)->DefaultHistory() << "\"";
        }
        if ((*sv)->Latex() != "") {
            cout << " Latex=\"" << (*sv)->Latex() << "\"";
        }
        cout << "/>" << endl;
    }

    vector<Function *>::iterator f;
    for (f = Functions.begin(); f != Functions.end(); ++f) {
        cout << "<Function";
        cout << " Name=\"" << (*f)->Name() << "\"";
        if ((*f)->Description() != "") {
            cout << " Description=\"" << (*f)->Description() << "\"";
        }
        cout << " Formula=\"" << (*f)->Formula() << "\"";
        cout << "/>" << endl;
    }

    cout << "</VectorField>" << endl;

}

void check_bad_name(string value, string attr, string element)
{
    if (value == "I") {
        cerr << "Error: bad symbol for '" << attr << "' of '" << element << "'." << endl;
        cerr << "The symbol 'I' is not allowed.  It conflicts with a predefined" << endl;
        cerr << "constant of the symbolic processor used by VFGEN, and it will" << endl;
        cerr << "conflict with a predefined constant when used in C code that also " << endl;
        cerr << "includes <complex.h>." << endl;
        exit(-1);
    }
}


//
// ReadXML
//

int VectorField::ReadXML(string xmlfilename)
{
    FILE *xmlfile;
    mxml_node_t *tree;
    mxml_node_t *node;
    bool bad_attr;

    xmlfile = fopen(xmlfilename.c_str(), "r");
    if (xmlfile == nullptr) {
       // Failed to open the file.
       cerr << "Error: Unable to open " << xmlfilename << "\n";
       exit(-1);
    }
    tree = mxmlLoadFile(nullptr, xmlfile, MXML_NO_CALLBACK);
    fclose(xmlfile);
    if (tree == nullptr) {
        cerr << "Error: Unable to load the vector field from the file " << xmlfilename << ".\n";
        cerr << "There may be an error in the XML definition of the vector field.\n";
        mxmlDelete(tree);
        exit(-1);
    }

    node = mxmlFindElement(tree, tree, "VectorField", nullptr, nullptr, MXML_DESCEND);
    if (node == nullptr) {
        cerr << "Error: No VectorField element found in XML defintion.\n";
        mxmlDelete(tree);
        exit(-1);
    }
    else {
        bad_attr = false;
        for (int i = 0; i < mxmlElementGetAttrCount(node); ++i) {
            const char *name;
            mxmlElementGetAttrByIndex(node, i, &name);
            string attr(name);
            if (attr != "Name" && attr != "IndependentVariable" && attr != "Description") {
                cerr << "Error: The VectorField element has an unknown attribute: " << attr << endl;
                bad_attr = true;
            }
        }
        if (bad_attr) {
            exit(-1);
        }
        const char *attr;
        attr = mxmlElementGetAttr(node, "Name");
        if (attr == nullptr) {
            cerr << "Error: The VectorField element has no Name attribute.\n";
            mxmlDelete(tree);
            exit(-1);
        }
        else {
            if (!isValidName(attr)) {
                cerr << "Error: The VectorField Name \"" << attr << "\" is not valid.\n";
                mxmlDelete(tree);
                exit(-1);
            }
            string s(attr);
            Name(s);
        }
        attr = mxmlElementGetAttr(node, "Description");
        if (attr != nullptr) {
            string s(attr);
            Description(s);
        }
        attr = mxmlElementGetAttr(node, "IndependentVariable");
        if (attr == nullptr) {
            IndependentVariable = "t";
        }
        else {
            if (!isValidName(attr)) {
                cerr << "Error: The VectorField IndependentVariable \"" << attr << "\" is not valid.\n";
                mxmlDelete(tree);
                exit(-1);
            }
            string s(attr);
            check_bad_name(s, "IndependentVariable", "VectorField");
            IndependentVariable = s;
        }
    }

    //
    // Get the constants
    //
    for (node = mxmlFindElement(tree, tree, "Constant", nullptr, nullptr, MXML_DESCEND);
         node != nullptr;
         node = mxmlFindElement(node, tree, "Constant", nullptr, nullptr, MXML_DESCEND)) {
        bad_attr = false;
        for (int i = 0; i < mxmlElementGetAttrCount(node); ++i) {
            const char *name;
            mxmlElementGetAttrByIndex(node, i, &name);
            string attr(name);
            if (attr != "Name" && attr != "Value" && attr != "Description" && attr != "Latex") {
                cerr << "Error: A Constant element has an unknown attribute: " << attr << endl;
                bad_attr = true;
            }
        }
        if (bad_attr) {
            cerr << "Valid Constant attributes are: Name, Value, Description, Latex.\n";
            exit(-1);
        } 
        const char *attr;
        attr = mxmlElementGetAttr(node, "Name");
        if (attr == nullptr) {
            cerr << "Error: A Constant element has no Name attribute.\n";
            mxmlDelete(tree);
            exit(-1);
        }
        else {
            if (!isValidName(attr)) {
                cerr << "Error: The Constant Name \"" << attr << "\" is not valid.\n";
                mxmlDelete(tree);
                exit(-1);
            }
            string name(attr);
            check_bad_name(name, "Name", "Constant");
            Constant *c = new Constant(name);
            AddConstant(c);
            attr = mxmlElementGetAttr(node, "Description");
            if (attr != nullptr) {
                string descr(attr);
                c->Description(descr);
            }
            attr = mxmlElementGetAttr(node, "Value");
            if (attr == nullptr) {
                cerr << "Error: The Constant element with Name=\"" << c->Name() << "\" has no Value attribute.\n";
                mxmlDelete(tree);
                exit(-1);
            }
            else {
                string val(attr);
                c->Value(val);
            }
            attr = mxmlElementGetAttr(node, "Latex");
            if (attr != nullptr) {
                string latex(attr);
                c->Latex(latex);
            }
        }
    }

    //
    // Get the parameters
    //
    for (node = mxmlFindElement(tree, tree, "Parameter", nullptr, nullptr, MXML_DESCEND);
         node != nullptr;
         node = mxmlFindElement(node, tree, "Parameter", nullptr, nullptr, MXML_DESCEND)) {
        bad_attr = false;
        for (int i = 0; i < mxmlElementGetAttrCount(node); ++i) {
            const char *name;
            mxmlElementGetAttrByIndex(node, i, &name);
            string attr(name);
            // string attr = node->value.element.attrs[i].name;
            if (attr != "Name" && attr != "DefaultValue" && attr != "Description" && attr != "Latex") {
                cerr << "Error: A Parameter element has an unknown attribute: " << attr << endl;
                bad_attr = true;
            }
        }
        if (bad_attr) {
            cerr << "Valid Parameter attributes are: Name, DefaultValue, Description, Latex.\n";
            exit(-1);
        }
        const char *attr;
        attr = mxmlElementGetAttr(node, "Name");
        if (attr == nullptr) {
            cerr << "Error: A Parameter element has no Name attribute.\n";
            mxmlDelete(tree);
            exit(-1);
        }
        else {
            if (!isValidName(attr)) {
                cerr << "Error: The Parameter Name \"" << attr << "\" is not valid.\n";
                mxmlDelete(tree);
                exit(-1);
            }
            string name(attr);
            check_bad_name(name, "Name", "Parameter");
            Parameter *p = new Parameter(name);
            AddParameter(p);
            attr = mxmlElementGetAttr(node, "Description");
            if (attr != nullptr) {
                string descr(attr);
                p->Description(descr);
            }
            attr = mxmlElementGetAttr(node, "DefaultValue");
            if (attr != nullptr) {
                string defval(attr);
                p->DefaultValue(defval);
            }
            attr = mxmlElementGetAttr(node, "Latex");
            if (attr != nullptr) {
                string latex(attr);
                p->Latex(latex);
            }
        }
    }

    //
    // Get the auxiliary expressions
    //
    for (node = mxmlFindElement(tree, tree, "Expression", nullptr, nullptr, MXML_DESCEND);
         node != nullptr;
         node = mxmlFindElement(node, tree, "Expression", nullptr, nullptr, MXML_DESCEND)) {
        bad_attr = false;
        for (int i = 0; i < mxmlElementGetAttrCount(node); ++i) {
            const char *name;
            mxmlElementGetAttrByIndex(node, i, &name);
            string attr(name);
            if (attr != "Name" && attr != "Formula" && attr != "Description" && attr != "Latex") {
                cerr << "Error: An Expression element has an unknown attribute: " << attr << endl;
                bad_attr = true;
            }
        }
        if (bad_attr) {
            cerr << "Valid Expression attributes are: Name, Formula, Description, Latex.\n";
            exit(-1);
        }
        const char *attr;
        attr = mxmlElementGetAttr(node, "Name");
        if (attr == nullptr) {
            cerr << "Error: An Expression element has no Name attribute.\n";
            mxmlDelete(tree);
            exit(-1);
        }
        else {
            if (!isValidName(attr)) {
                cerr << "Error: The Expression Name \"" << attr << "\" is not valid.\n";
                mxmlDelete(tree);
                exit(-1);
            }
            string name(attr);
            check_bad_name(name, "Name", "Expression");
            Expression *e = new Expression(name);
            AddExpression(e);
            attr = mxmlElementGetAttr(node, "Description");
            if (attr != nullptr) {
                string descr(attr);
                e->Description(descr);
            }
            attr = mxmlElementGetAttr(node, "Formula");
            if (attr == nullptr) {
                cerr << "Error: The Expression with Name=\"" << e->Name() << "\" has no Formula attribute.\n";
                mxmlDelete(tree);
                exit(-1);
            }
            else {
                string f(attr);
                e->Formula(f);
            }
            attr = mxmlElementGetAttr(node, "Latex");
            if (attr != nullptr) {
                string latex(attr);
                e->Latex(latex);
            }
        }
    }

    //
    // Get the state variables
    //
    for (node = mxmlFindElement(tree, tree, "StateVariable", nullptr, nullptr, MXML_DESCEND);
         node != nullptr;
         node = mxmlFindElement(node, tree, "StateVariable", nullptr, nullptr, MXML_DESCEND)) {
        bad_attr = false;
        for (int i = 0; i < mxmlElementGetAttrCount(node); ++i) {
            const char *name;
            mxmlElementGetAttrByIndex(node, i, &name);
            string attr(name);
            if (attr != "Name" && attr != "DefaultInitialCondition" && attr != "Description"
                  && attr != "Formula" && attr != "PeriodFrom" && attr != "PeriodTo"
                  && attr != "DefaultHistory" && attr != "Latex") {
                cerr << "Error: A StateVariable element has an unknown attribute: " << attr << endl;
                bad_attr = true;
            }
        }
        if (bad_attr) {
            cerr << "Valid StateVariable attributes are: Name, Formula, Description, DefaultInitialCondition, DefaultHistory, PeriodFrom, PeriodTo, Latex.\n";
            exit(-1);
        }
        const char *attr;
        attr = mxmlElementGetAttr(node, "Name");
        if (attr == nullptr) {
            cerr << "Error: A StateVariable element has no Name attribute.\n";
            mxmlDelete(tree);
            exit(-1);
        }
        else {
            if (!isValidName(attr)) {
                cerr << "Error: The StateVariable Name \"" << attr << "\" is not valid.\n";
                mxmlDelete(tree);
                exit(-1);
            }
            string name(attr);
            check_bad_name(name, "Name", "StateVariable");
            StateVariable *sv = new StateVariable(name);
            AddStateVariable(sv);
            attr = mxmlElementGetAttr(node, "Description");
            if (attr != nullptr) {
                string descr(attr);
                sv->Description(descr);
            }
            attr = mxmlElementGetAttr(node, "Formula");
            if (attr == nullptr) {
                cerr << "Error: The StateVariable with Name=\"" << sv->Name() << "\" has no Formula attribute.\n";
                mxmlDelete(tree);
                exit(-1);
            }
            else {
                string f(attr);
                sv->Formula(f);
            }
            attr = mxmlElementGetAttr(node, "PeriodFrom");
            if (attr != nullptr) {
                string pfrom(attr);
                sv->PeriodicFrom(pfrom);
            }
            attr = mxmlElementGetAttr(node, "PeriodTo");
            if (attr != nullptr) {
                string pto(attr);
                sv->PeriodicTo(pto);
            }
            if (sv->PeriodicFrom() != "" && sv->PeriodicTo() == "") {
                cerr << "Error: The StateVariable with Name=\"" << sv->Name() << "\" has a PeriodicFrom attribute but no PeriodicTo attribute.\n";
                mxmlDelete(tree);
                exit(-1);
            }
            if (sv->PeriodicFrom() == "" && sv->PeriodicTo() != "") {
                cerr << "Error: The StateVariable with Name=\"" << sv->Name() << "\" has a PeriodTo attribute but no PeriodicFrom attribute.\n";
                mxmlDelete(tree);
                exit(-1);
            }
            attr = mxmlElementGetAttr(node, "DefaultInitialCondition");
            if (attr != nullptr) {
                string ic(attr);
                sv->DefaultInitialCondition(ic);
            }
            attr = mxmlElementGetAttr(node, "DefaultHistory");
            if (attr != nullptr) {
                string hist(attr);
                sv->DefaultHistory(hist);
            }
            attr = mxmlElementGetAttr(node, "Latex");
            if (attr != nullptr) {
                string latex(attr);
                sv->Latex(latex);
            }
        }
    }

    //
    // Get the functions
    //
    for (node = mxmlFindElement(tree, tree, "Function", nullptr, nullptr, MXML_DESCEND);
         node != nullptr;
         node = mxmlFindElement(node, tree, "Function", nullptr, nullptr, MXML_DESCEND)) {
        bad_attr = false;
        for (int i = 0; i < mxmlElementGetAttrCount(node); ++i) {
            const char *name;
            mxmlElementGetAttrByIndex(node, i, &name);
            string attr(name);
            if (attr != "Name" && attr != "Formula" && attr != "Description") {
                cerr << "Error: A Function element has an unknown attribute: " << attr << endl;
                bad_attr = true;
            }
        }
        if (bad_attr) {
            cerr << "Valid Function attributes are: Name, Formula, Description.\n";
            exit(-1);
        }
        const char *attr;
        attr = mxmlElementGetAttr(node, "Name");
        if (attr == nullptr) {
            cerr << "Error: A Function element has no Name attribute.\n";
            mxmlDelete(tree);
            exit(-1);
        }
        else {
            if (!isValidName(attr)) {
                cerr << "Error: The Function Name \"" << attr << "\" is not valid.\n";
                mxmlDelete(tree);
                exit(-1);
            }
            string name(attr);
            check_bad_name(name, "Name", "Function");
            Function *func = new Function(name);
            AddFunction(func);
            attr = mxmlElementGetAttr(node, "Description");
            if (attr != nullptr) {
                string descr(attr);
                func->Description(descr);
            }
            attr = mxmlElementGetAttr(node, "Formula");
            if (attr == nullptr) {
                cerr << "Error: The Function element with Name=\"" << func->Name() << "\" has no Formula attibute.\n";
                mxmlDelete(tree);
                exit(-1);
            }
            else {
                string f(attr);
                func->Formula(f);
            }
        }
    }

    mxmlDelete(tree);
    return 0;
}
