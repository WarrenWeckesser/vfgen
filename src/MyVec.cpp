//
// MyVec.cpp
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
#include "MyVec.h"

using namespace std;

bool MyVectorCmpClass::operator()(const vector<int> *a1, const vector<int> *a2) const
    {
    if (a1->size() > a2->size())
        return true;
    if (a1->size() < a2->size())
        return false;
    for (size_t k = 0; k < a1->size(); ++k)
        {
        if ((*a1)[k] < (*a2)[k])
            return true;
        if ((*a1)[k] > (*a2)[k])
            return false;
        }
    return false;
    }



typedef map< vector<int> *, double, MyVectorCmpClass> vectormap;

void PrintMyVec(ofstream &out, vector<int> *a)
    {
    for (size_t j = 0; j < a->size(); ++j)
        out << " " << (*a)[j];
    }

//
// Copy from a to b
//

void CopyMyVec(vector<int> *a, vector<int> *b)
    {
    for (size_t k = 0; k < a->size(); ++k)
        b->push_back((*a)[k]);
    }

int SumVec(vector<int> *a)
    {
    int sum = 0;
    for (vector<int>::iterator ai = a->begin(); ai != a->end(); ++ai)
        sum = sum + *ai;
    return sum;
    }

void MyVecExtend1(vector<int> *a)
    {
    a->push_back(1);
    }
