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
	cout << "Individual class constructor called." << endl;
	rawFitness = 0.0;
	valid = false;
	k = kmax;
	clusCenter = new double* [k];
	//clusters = new vector<int>*[k];
	for (int count = 0; count < k; count++)
	{
	    clusCenter[count] = new double[dim];
	   // clusters[count] = new vector<int>;
	}
	threshold =  new double[k];
	active = new bool[k];
	active_ctr = 0;
}

Individual::~Individual() {
	// TODO Auto-generated destructor stub
	//cout << "Destructor called" << endl;
	for(int i = 0; i < k; ++i) {
	    delete [] clusCenter[i];
	   // delete clusters[i];
	}
	delete [] clusCenter;
	delete [] active;
	delete [] threshold;
//	delete [] clusters;

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

/*bool Individual::operator<=(const Individual& right)
{//overloaded operator <=
  return (rawFitness <= right.rawFitness)? true : false;
}//operator<=*/

