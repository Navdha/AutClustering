/*
 * Individual.cpp
 *
 *  Created on: Jul 7, 2015
 *      Author: Navs
 */

#include "Individual.h"
#include <stdexcept>

Individual::Individual(int kmax, int dim) {
	// TODO Auto-generated constructor stub
	rawFitness = 0.0;
	valid = false;
	k = kmax;
	clusCenter = new double* [kmax];
	clusters = new vector<int>*[kmax];
	for (int count = 0; count < k; count++)
	{
	    clusCenter[count] = new double[dim];
	    clusters[count] = new vector<int>;
	}
	threshold =  new double[kmax];
	active = new bool[kmax];
	active_ctr = 0;
}

Individual::~Individual() {
	// TODO Auto-generated destructor stub
	for(int i = 0; i < k; ++i) {
	    delete [] clusCenter[i];
	}
	delete [] clusCenter;
	delete [] active;
	delete [] threshold;
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

bool Individual::operator<=(const Individual& right)
{//overloaded operator <=
  return (rawFitness <= right.rawFitness)? true : false;
}//operator<=

