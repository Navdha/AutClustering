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
//	cout << "Individual class constructor called." << endl;
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
	dimension = dim;
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

ostream& operator<<(ostream& o, const Individual& org)
{//overloaded insertion operator
  o << "#active clusters: " << org.active_ctr << endl;
  o << "fitness: " << org.rawFitness << endl;
  o << "validity: " << org.valid << endl;

  for (int i = 0; i < org.k; i++){
    if (org.active[i]){
      o << "cluster " << i << " centroid: (" << org.clusCenter[i][0];
      for (int j = 1; j < org.dimension; j++){
	o << ", " << org.clusCenter[i][j];
      }
      o << ")" << endl;
    }//end if
  }//end for i
  return o;
}//operator<<

Individual::Individual(const Individual& obj){
	rawFitness = obj.rawFitness;
	valid = obj.valid;
	k = obj.k;
	dimension = obj.dimension;
	clusCenter = new double* [k];
	for (int count = 0; count < k; count++)
	{
	    clusCenter[count] = new double[dimension];
	    for(int i = 0; i < dimension; i++) {
	      clusCenter[count][i] = obj.clusCenter[count][i];
	    }
	}
	threshold =  new double[k];
	active = new bool[k];
	for(int  i = 0 ; i < k; i ++){
	  threshold[i] = obj.threshold[i];
	  active[i] = obj.active[i];
	}
	active_ctr = obj.active_ctr;
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

