/*
 * DEMain.h
 *
 *  Created on: Jul 13, 2015
 *      Author: Navdha Sah
 */

#ifndef DEMAIN_H_
#define DEMAIN_H_

#define stRand1Bin 1

#include "Individual.h"
#include "Population.h"
#include "Item.h"


class DEMain {
public:
	DEMain(int kmax, int dim, int gen, int** placeholder, const Item** items, int itemSize);
	virtual ~DEMain();
	void setup(double min[],double max[],int deStrategy,double diffScale,double crossoverProb);
	double calcFitness(Individual* org, int index);
	double dist(double* x, double* y);
	double* avgDist(Individual* org);
	void selectSamples(int org, int *s1, int *s2, int *s3);
	Individual* crossover(int org, int generation);
	void run();


	Population* p;

	int strategy;
	double scale;
	double probability;
	int generations;
	int pSize;
	int kmax;
	int dim;
	int numItems;
	int** tracker;
	const Item** attr;

private:
	void Rand1Bin(int candIndex);
};

#endif /* DEMAIN_H_ */
