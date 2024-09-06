//
//  codegen_utils.h
//
//
//  Prototypes for functions in codegen_utils.cpp
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

#ifndef CODEGEN_UTILS
#define CODEGEN_UTILS 1

char *DateTimeMsg();
void PrintVFGENComment(std::ofstream &fout, const char *prefix);
void Declare(std::ofstream &fout, std::string prefixstr, std::string typestr, GiNaC::lst name, std::string termstr);
void CDeclare(std::ofstream &fout, std::string typestr, GiNaC::lst names);
void CDeclare_double(std::ofstream &fout, GiNaC::lst names);
void MakeCArrayOfStrings(std::ofstream &fout, const char *var, GiNaC::lst names);
void MakePythonListOfStrings(std::ofstream &fout, const char *var, GiNaC::lst names, const char *pre);
size_t max_expr_len(GiNaC::lst names);
void GetFromVector(std::ofstream &fout, const char *skip, GiNaC::lst names, const char *assignop,
                    const char *vector, const char *braces, int istart, const char *term);
void GetFromVector2(std::ofstream &fout, const char *skip, GiNaC::lst names, const char *assignop,
                    const char *vector, const char *bropen, const char *brclose, int istart, const char *term);
void SetVectorFromNames(std::ofstream &fout, const char *skip, const char *vector, GiNaC::lst names,
                      const char *braces, int istart, const char *term);
void AssignNameValueLists(std::ofstream &fout, const char *skip,
                          GiNaC::lst names, const char *assignop, GiNaC::lst values,
                          const char *term);
void PrintList(std::ofstream &fout, GiNaC::lst names);
void PrintTransformedList(std::ofstream &fout, const std::string, GiNaC::lst names);
void PrintNameList(std::ofstream &fout, GiNaC::lst names);
void PrintPi(std::ofstream &fout);
void PrintMultilineComment(std::ofstream &fout, const std::string &comment,
        const std::string &pre, const std::string &cmark);

void print_power_as_fortran(const GiNaC::power& p, const GiNaC::print_csrc& c, unsigned level);
std::string fix_exp_notation(std::string &s);
void F77Write(std::ofstream &fout, std::string s);
void F77Declare(std::ofstream &fout, GiNaC::lst names);
void F90Write(std::ofstream &fout, std::string s);

#endif
