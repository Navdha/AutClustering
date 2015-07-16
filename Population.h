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
	Population(int kmax, int dim);
	virtual ~Population();
	void setIndividial(Individual* individual, int index);
	//Individual getIndividual(int index);

	int size;
	Individual** chromosome;//change name
	//int bestChromosomeIndex;
	//int worstChromosomeIndex;
};

#endif /* POPULATION_H_ */
