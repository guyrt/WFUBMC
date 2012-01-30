#include "log.hh"

// Global static pointer used to ensure a single instance of the class.
Logger* Logger::m_pInstance = NULL;  

/** This function is called to create an instance of the class. 
    Calling the constructor publicly is not allowed. The constructor 
    is private and is only called by this Instance function.
*/
  
Logger* Logger::Instance()
{
   if (!m_pInstance)   // Only allow one instance of class to be generated.
      m_pInstance = new Logger;

   return m_pInstance;
}

bool Logger::init(string _logFile)
{
	outstream.open(_logFile.c_str());
	return outstream.is_open();
}

void Logger::writeLine(string s){
	
	#pragma omp critical
	outstream << s;
	
}

void Logger::close(){
	outstream.close();
}
