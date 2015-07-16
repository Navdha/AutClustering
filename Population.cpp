/*
 * Population.cpp
 *
 *  Created on: Jul 8, 2015
 *      Author: Navs
 */

#include "Population.h"
#include <cassert>
#include <iostream>


Population::Population(int kmax, int dim){
	size = 10*dim;
	chromosome = new Individual*[size];
	//bestChromosomeIndex = -1;
	//worstChromosomeIndex = -1;
	assert(chromosome != NULL);
	  /*for (int i = 0; i < size; i++){
		  chromosome[i] = new Individual(kmax, dim);
	    assert(chromosome[i] != NULL);
	  }*/
}
Population::~Population() {
	// TODO Auto-generated destructor stub
	for (int i = 0; i < size; i++)
	    delete chromosome[i];

	  delete [] chromosome;
}

void Population:: setIndividial(Individual* individual, int index){
	delete chromosome[index];
	chromosome[index] = individual;
}

/*Individual Population::getIndividual(int index){
	return &chromosome[index];
}*/

