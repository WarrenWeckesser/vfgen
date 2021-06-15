
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

//
// ReadXML
//

int VectorField::ReadXML(string xmlfilename)
{
    FILE *xmlfile;
    mxml_node_t *tree;
    mxml_node_t *node;
    bool bad_attr;

    xmlfile = fopen(xmlfilename.c_str(),"r");
    if (xmlfile == NULL) {
       // Failed to open the file.
       cerr << "Error: Unable to open " << xmlfilename << "\n";
       exit(-1);
    }
    tree = mxmlLoadFile(NULL,xmlfile,MXML_NO_CALLBACK);
    fclose(xmlfile);
    if (tree == NULL) {
        cerr << "Error: Unable to load the vector field from the file " << xmlfilename << ".\n";
        cerr << "There may be an error in the XML definition of the vector field.\n";
        mxmlDelete(tree);
        exit(-1);
    }

    node = mxmlFindElement(tree,tree,"VectorField",NULL,NULL,MXML_DESCEND);
    if (node == NULL) {
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
        attr = mxmlElementGetAttr(node,"Name");
        if (attr == NULL) {
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
        attr = mxmlElementGetAttr(node,"Description");
        if (attr != NULL) {
            string s(attr);
            Description(s);
        }
        attr = mxmlElementGetAttr(node,"IndependentVariable");
        if (attr == NULL) {
            IndependentVariable = "t";
        }
        else {
            if (!isValidName(attr)) {
                cerr << "Error: The VectorField IndependentVariable \"" << attr << "\" is not valid.\n";
                mxmlDelete(tree);
                exit(-1);
            }
            string s(attr);
            IndependentVariable = s;
        }
    }

    //
    // Get the constants
    //
    for (node = mxmlFindElement(tree,tree,"Constant",NULL,NULL,MXML_DESCEND);
         node != NULL;
         node = mxmlFindElement(node,tree,"Constant",NULL,NULL,MXML_DESCEND)) {
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
        attr = mxmlElementGetAttr(node,"Name");
        if (attr == NULL) {
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
            Constant *c = new Constant(name);
            AddConstant(c);
            attr = mxmlElementGetAttr(node,"Description");
            if (attr != NULL) {
                string descr(attr);
                c->Description(descr);
            }
            attr = mxmlElementGetAttr(node,"Value");
            if (attr == NULL) {
                cerr << "Error: The Constant element with Name=\"" << c->Name() << "\" has no Value attribute.\n";
                mxmlDelete(tree);
                exit(-1);
            }
            else {
                string val(attr);
                c->Value(val);
            }
            attr = mxmlElementGetAttr(node,"Latex");
            if (attr != NULL) {
                string latex(attr);
                c->Latex(latex);
            }
        }
    }

    //
    // Get the parameters
    //
    for (node = mxmlFindElement(tree,tree,"Parameter",NULL,NULL,MXML_DESCEND);
         node != NULL;
         node = mxmlFindElement(node,tree,"Parameter",NULL,NULL,MXML_DESCEND)) {
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
        attr = mxmlElementGetAttr(node,"Name");
        if (attr == NULL) {
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
            Parameter *p = new Parameter(name);
            AddParameter(p);
            attr = mxmlElementGetAttr(node,"Description");
            if (attr != NULL) {
                string descr(attr);
                p->Description(descr);
            }
            attr = mxmlElementGetAttr(node,"DefaultValue");
            if (attr != NULL) {
                string defval(attr);
                p->DefaultValue(defval);
            }
            attr = mxmlElementGetAttr(node,"Latex");
            if (attr != NULL) {
                string latex(attr);
                p->Latex(latex);
            }
        }
    }

    //
    // Get the auxiliary expressions
    //
    for (node = mxmlFindElement(tree,tree,"Expression",NULL,NULL,MXML_DESCEND);
         node != NULL;
         node = mxmlFindElement(node,tree,"Expression",NULL,NULL,MXML_DESCEND)) {
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
        attr = mxmlElementGetAttr(node,"Name");
        if (attr == NULL) {
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
            Expression *e = new Expression(name);
            AddExpression(e);
            attr = mxmlElementGetAttr(node,"Description");
            if (attr != NULL) {
                string descr(attr);
                e->Description(descr);
            }
            attr = mxmlElementGetAttr(node,"Formula");
            if (attr == NULL) {
                cerr << "Error: The Expression with Name=\"" << e->Name() << "\" has no Formula attribute.\n";
                mxmlDelete(tree);
                exit(-1);
            }
            else {
                string f(attr);
                e->Formula(f);
            }
            attr = mxmlElementGetAttr(node,"Latex");
            if (attr != NULL) {
                string latex(attr);
                e->Latex(latex);
            }
        }
    }

    //
    // Get the state variables
    //
    for (node = mxmlFindElement(tree,tree,"StateVariable",NULL,NULL,MXML_DESCEND);
         node != NULL;
         node = mxmlFindElement(node,tree,"StateVariable",NULL,NULL,MXML_DESCEND)) {
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
        attr = mxmlElementGetAttr(node,"Name");
        if (attr == NULL) {
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
            StateVariable *sv = new StateVariable(name);
            AddStateVariable(sv);
            attr = mxmlElementGetAttr(node,"Description");
            if (attr != NULL) {
                string descr(attr);
                sv->Description(descr);
            }
            attr = mxmlElementGetAttr(node,"Formula");
            if (attr == NULL) {
                cerr << "Error: The StateVariable with Name=\"" << sv->Name() << "\" has no Formula attribute.\n";
                mxmlDelete(tree);
                exit(-1);
            }
            else {
                string f(attr);
                sv->Formula(f);
            }
            attr = mxmlElementGetAttr(node,"PeriodFrom");
            if (attr != NULL) {
                string pfrom(attr);
                sv->PeriodicFrom(pfrom);
            }
            attr = mxmlElementGetAttr(node,"PeriodTo");
            if (attr != NULL) {
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
            attr = mxmlElementGetAttr(node,"DefaultInitialCondition");
            if (attr != NULL) {
                string ic(attr);
                sv->DefaultInitialCondition(ic);
            }
            attr = mxmlElementGetAttr(node,"DefaultHistory");
            if (attr != NULL) {
                string hist(attr);
                sv->DefaultHistory(hist);
            }
            attr = mxmlElementGetAttr(node,"Latex");
            if (attr != NULL) {
                string latex(attr);
                sv->Latex(latex);
            }
        }
    }

    //
    // Get the functions
    //
    for (node = mxmlFindElement(tree,tree,"Function",NULL,NULL,MXML_DESCEND);
         node != NULL;
         node = mxmlFindElement(node,tree,"Function",NULL,NULL,MXML_DESCEND)) {
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
        attr = mxmlElementGetAttr(node,"Name");
        if (attr == NULL) {
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
            Function *func = new Function(name);
            AddFunction(func);
            attr = mxmlElementGetAttr(node,"Description");
            if (attr != NULL) {
                string descr(attr);
                func->Description(descr);
            }
            attr = mxmlElementGetAttr(node,"Formula");
            if (attr == NULL) {
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
