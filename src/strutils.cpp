//
// File:    strutils.cpp
// Author:  Warren Weckesser
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

#include <string>
#include <sstream>
#include <stdexcept>

using namespace std;


string TrimSpaces(const string str)
    {
    string::size_type first,last;
    string newstr;

    first = str.find_first_not_of(" ");
    if (first == string::npos)
        {
        newstr = string("");
        }
    else
        {
        last  = str.find_last_not_of(" ");
        newstr = str.substr(first,last-first+1);
        }
    return newstr;
    }


bool isValidName(const string s)
    {
    const string lower("abcdefghijklmnopqrstuvwxyz");
    const string upper("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
    const string letters = lower + upper + "_";
    const string digits("0123456789");
    const string validchars = letters + digits;

    if (s.empty())
        {
        // cerr << "s is empty\n";
        return false;
        }
    if (letters.find(s[0]) == string::npos)
        {
        // cerr << "char '" << s[0] << "' is not a valid first character\n";
        return false;
        }
    for (unsigned k = 1; k < s.length(); ++k)
        if (validchars.find(s[k]) == string::npos)
            {
            // cerr << "char '" << s[k] << "' is not valid\n";
            return false;
            }
    return true;  
    }


bool is_int(const string s)
    {
    for (unsigned i=0; i < s.length(); ++i)
        {
        if (s[i] < '0' || s[i] > '9')
            return false;
        }
    return true;
    }


class BadConversion : public std::runtime_error {
    public:
    BadConversion(const std::string& s) : std::runtime_error(s)
        {
        }
};
 
int string_to_int(const std::string& s)
   {
   std::istringstream i(s);
   int x;
   if (!(i >> x))
       throw BadConversion("string_to_int(\"" + s + "\")");
   return x;
   }
