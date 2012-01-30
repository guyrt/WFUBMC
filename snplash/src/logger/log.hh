//      log.hh
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
 * Pattern taken from yoLinux.com.
 * 
 * Singleton logger, meaning that there will only be one logger.
 * 
 * Called using Logger::Instance()->writeLine("foo");
 */

#ifndef LOGGER_H
#define LOGGER_H

#include <string>
#include <fstream> // for file output.

using namespace std;

class Logger{
public:
   static Logger* Instance();
   bool init(string logFile);
   void writeLine(string line);
   void close();

private:
   Logger(){};  // Private so that it can  not be called
   Logger(Logger const&){};             // copy constructor is private
   Logger& operator=(Logger const& ){return *this;}  // assignment operator is private
   static Logger* m_pInstance;
   
   ofstream outstream;
};

#endif
