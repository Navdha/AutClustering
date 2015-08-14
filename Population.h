/*
 * Population.h
 *
 *  Created on: Jul 8, 2015
 *      Author: Navs
 */

#ifndef POPULATION_H_
#define POPULATION_H_

#include "Individual.h"


class Population {
 public:
  Population(int maxNumClusters, int numFeatures, int popFactorScale);
  ~Population();

  int popSize;
  Individual** org;
  int bestOrgIndex;
  int worstOrgIndex;
  friend ostream& operator<<(ostream& o, const Population& pop);
	
};

#endif /* POPULATION_H_ */
