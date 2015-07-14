/*
 * Individual.cpp
 *
 *  Created on: Jul 7, 2015
 *      Author: Navs
 */

#include "Individual.h"
#include <stdexcept>

Individual::Individual() {
	// TODO Auto-generated constructor stub
	rawFitness = 0.0;
	valid = false;
}

Individual::~Individual() {
	// TODO Auto-generated destructor stub
}

bool Individual::isValid() {
	return valid;
}

bool Individual::getValid(){
	return isValid();
}

void Individual::setValid(bool valid){
	this->valid=valid;
	if (valid == false) this->rawFitness=0.0;
}

double Individual::getFitness(){
	return rawFitness;
}

void Individual::setFitness(double rawFitness){
if(rawFitness < 0)
	throw std::invalid_argument("received negative value");
	this->rawFitness = rawFitness;
}
