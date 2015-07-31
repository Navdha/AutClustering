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
#include <utility>
#include <vector>
#include <cstdlib>

class Dist_IC {
public:
		double distance;
		int itemIndex;
		int clustIndex;
	};

//typedef double (*compfn)(const void*, const void*);
class DEMain {
public:
	DEMain(int kmax, int dim, int gen, int** placeholder, Item** items, int itemSize, bool calcDB);
	virtual ~DEMain();
	void setup(double min[],double max[]);
	double calcFitness(Individual* org, int index, bool isInitial, int genNum);
	double dist(double* x, double* y);
	double* avgDist(Individual* org);
	void selectSamples(int org, int *s1, int *s2, int *s3);
	Individual* crossover(int org, int generation, double min[], double max[]);
	void run(double min[], double max[]);
	void report(int index);
	void reshuffle(Individual* org, int size, int index,  bool isInitial);
	void calcDistBtwnItems();

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
	Item** attr;
	vector<int>** clusters;
	Dist_IC* knn;
	int* offspring_arr;
	double** distItem;
	bool isDB;
private:
	void Rand1Bin(int candIndex);
};

#endif /* DEMAIN_H_ */
