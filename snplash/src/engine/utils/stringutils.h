//      stringutils.h
//      
//      Copyright 2010 Richard T. Guy <guyrt7@wfu.edu>
//      
//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//      
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//      
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//      MA 02110-1301, USA.

/*
 * String utilities designed to return a string object using different formatting
 * patterns. 
 */

#ifndef STRNUTILS_H
#define STRNUTILS_H

#define DB_V_STRNUTILS 0

#include <string>
#include <sstream>
#include <math.h>
#include <iostream>

namespace strnutils{

std::string spaced_string(std::string s, int numspaces, int buffer=0);
std::string spaced_number(int number, int numspaces, int buffer=0);
std::string spaced_number(double number, int numspaces, int precision, int buffer=0);
int int_log(int i);
int int_log(double i);
int small_log(double i);

/** test the strnutils methods. **/
void test();

}

#endif
