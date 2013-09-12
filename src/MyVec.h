//
// MyVec.h
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

#include <vector>
#include <map>

struct MyVectorCmpClass
    {
    bool operator()(const std::vector<int> *a1, const std::vector<int> *a2) const;
    };


typedef std::map< std::vector<int> *, double, MyVectorCmpClass> vectormap;

void PrintMyVec(std::ofstream &out, std::vector<int> *a);
void CopyMyVec(std::vector<int> *a, std::vector<int> *b);
int SumVec(std::vector<int> *a);
void MyVecExtend1(std::vector<int> *a);

