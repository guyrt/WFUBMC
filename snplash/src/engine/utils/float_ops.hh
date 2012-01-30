//      float_ops.hh
//      
//      Copyright 2010 Richard T. Guy <guy@cs.toronto.edu>
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


/**
 * @file float_ops.hh
 * 
 * Contains equals operator for floating point numbers.
 * 
 * 
 */ 
#ifndef FLOAT_OPS_H
#define FLOAT_OPS_H

#define EPS_EQUALS 1e-10

bool equal(double i, double d);

#endif
