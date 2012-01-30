#include "snpinfo_out.hh"

SnpInfo::SnpInfo(){
	index = -1;
	name = "";
	position = -1;
	diff = 0;
	chr = "";
	majAllele = '0';
	minAllele = '0';
	refAllele = '0';
}

void SnpInfo::print(std::ostream& out, int maxMapSize, bool printMap) // output
{
	if (printMap){
		out << strnutils::spaced_string(chr, 3,0);
		out << strnutils::spaced_string(name, maxMapSize, 1);
		out << strnutils::spaced_number(position, 10, 1);
		out << strnutils::spaced_number(diff, 10, 1);
	}
    else{
		out << strnutils::spaced_number(index, 10);
	}
}
