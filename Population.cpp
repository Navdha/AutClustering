/*
 * Population.cpp
 *
 *  Created on: Jul 8, 2015
 *      Author: Navs
 */

#include "Population.h"
#include <cassert>
#include <iostream>


Population::Population(int maxNumClusters, int numFeatures, int popFactorScale){
  //	cout << "Population class constructor called" << endl;
  popSize = popFactorScale * numFeatures;
  org = new Individual*[popSize];
  bestOrgIndex = -1;
  worstOrgIndex = -1;
  assert(org != NULL);
}



ostream& operator<<(ostream& o, const Population& pop)
{// overloaded insertion operator for Population
  o << "Population content " << endl;
  o << "\tBest Individiual index: " << pop.bestOrgIndex << endl;
  
  
  for (int i = 0; i < pop.popSize; i++){
    o << "Individual " << i << endl;
    o << *(pop.org[i]) << endl;
  }
  return o;
}//operator<<


Population::~Population() {
  // TODO Auto-generated destructor stub
  cout << "Pop destructor" << endl;
  for (int i = 0; i < popSize; i++)
    delete org[i];

  delete [] org;
}


