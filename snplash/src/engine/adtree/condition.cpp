/*
 *      condition.cpp
 *
 *      Copyright 2009 Richard T. Guy <guyrt7@wfu.edu>
 *
 *      This program is free software; you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation; either version 2 of the License, or
 *      (at your option) any later version.
 *
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *
 *      You should have received a copy of the GNU General Public License
 *      along with this program; if not, write to the Free Software
 *      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *      MA 02110-1301, USA.
 */
#include "condition.h"

Condition::Condition(){
	this->genotype_conditional = GE;
	this->genotype_reference = 1;
	this->attribute_index = 1;
}

/**
 * Actually evaluate a condition given the data in vector
 *
 * Returns: the boolean value for this Condition on this data.
 */
bool Condition::evaluate(vector<short> *vec){

	if(attribute_index == -1){return true;} // A fake value.

	short a = vec->at(this->attribute_index);
	
	// Missing is always false.
	if(a == 0){return false;}

	// Treat each 1 2 as 2 1.
	if(a == 3){a = 2;}

	switch (this->genotype_conditional)
	{
		case GE :  // >=
			return a >= genotype_reference;
		break;
		case GT : // >
			return a > genotype_reference;
		break;
		case LE : // <=
			return a <= genotype_reference;
		break;
		case LT : // <
			return a < genotype_reference;
		break;
		case EQ : // ==
			return a == genotype_reference;
		break;
		case NE : // !=
			return a != genotype_reference;
		break;
	}
	return false;
}

/**
 * Return a formatted string.
 *
 * TODO: consider changing to to_string
 *
 * @return string representation of the condition.
 */
string Condition::print(DataAccess *data){

	stringstream ss;
	string cond;
	switch(genotype_conditional){
		case GE :  // >=
			cond = ">=";
		break;
		case GT : // >
			cond = ">";
		break;
		case LE : // <=
			cond = "<=";
		break;
		case LT : // <
			cond = "<";
		break;
		case EQ : // ==
			cond = "==";
		break;
		case NE : // !=
			cond = "!=";
		break;
	}
	if (attribute_index < 0){
		ss << "Base: " ;
	}else{
		ss << data->snp_name(attribute_index) << " " << cond << " " << genotype_reference << " " ;
	}
	return ss.str();
}

string Condition::print_reverse(DataAccess *data){

	stringstream ss;
	string cond;
	switch(genotype_conditional){
		case GE :  // >=
			cond = "<";
		break;
		case GT : // >
			cond = "<=";
		break;
		case LE : // <=
			cond = ">";
		break;
		case LT : // <
			cond = ">=";
		break;
		case EQ : // ==
			cond = "!=";
		break;
		case NE : // !=
			cond = "==";
		break;
	}
	if (attribute_index < 0){
		ss << "Base: " ;
	}else{
		ss << data->snp_name(attribute_index) << " " << cond << " " << genotype_reference << " " ;
	}
	return ss.str();
}

string Condition::hash(){
	stringstream ss;
	string cond;
	switch(genotype_conditional){
		case GE :  // >=
			cond = ">=";
		break;
		case GT : // >
			cond = ">";
		break;
		case LE : // <=
			cond = "<=";
		break;
		case LT : // <
			cond = "<";
		break;
		case EQ : // ==
			cond = "==";
		break;
		case NE : // !=
			cond = "!=";
		break;
	}
	ss << attribute_index << cond << genotype_reference;
	return ss.str();
}

/**
 * Switch the sign of the conditional to its opposite
 */
void Condition::inverse(){
	switch (genotype_conditional){
		case Condition::GE :
			genotype_conditional = Condition::LT;
		break;
		case Condition::GT :
			genotype_conditional = Condition::LE;
		break;
		case Condition::LE :
			genotype_conditional = Condition::GT;
		break;
		case Condition::LT :
			genotype_conditional = Condition::GE;
		break;
		case Condition::EQ :
			genotype_conditional = Condition::NE;
		break;
		case Condition::NE :
			genotype_conditional = Condition::EQ;
		break;
	}
}

/**
 * Check for equality of all parts of the condition.
 *
 * @return int 1 if equal, 0 otherwise
 */
int Condition::operator==(const Condition& right){

	if(attribute_index == right.attribute_index && genotype_reference == right.genotype_reference
	&& genotype_conditional == right.genotype_conditional){

			return 1;

	}
	else{
		return 0;
	}

}
