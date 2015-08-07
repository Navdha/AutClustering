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
#include "Parameters.h"
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
	DEMain(int dim, int** placeholder, Item** items, int itemSize, int validityIndex, Parameters param);
	~DEMain();
	void setup(double min[],double max[]);
	double calcFitness(Individual* org, int index, bool isInitial, int genNum);
	double dist(double* x, double* y);
	double* avgDist(Individual* org);
	void selectSamples(int org, int *s1, int *s2, int *s3);
	Individual* crossover(int org, int generation, double min[], double max[]);
	void run(double min[], double max[]);
	void report(int index, int worstInd);
	void reshuffle(Individual* org, int size, int index,  bool isInitial);
	void calcDistBtwnItems();
	double calcDBIndex(Individual* org);
	double calcCSIndex(Individual* org);
	double calcPBIndex(Individual* org);
	double calcSD();

	Population* p;
	int strategy;
	double scale;
	double probability;
	long generations;
	int pSize;
	int kmax;
	int kmin;
	int dim;
	int numItems;
	double thresholdVal;
	int** tracker;
	Item** attr;
	vector<int>** clusters;
	Dist_IC* knn;
	int* offspring_arr;
	double** distItem;
	int indexForFit;
	bool* ItemUsed;
	bool* ClusFull;
	double* avgArr;
	double** newClustCenters;
	double* sumArr;
	int* ItemCounter;
	bool * new_pop;
private:
	void Rand1Bin(int candIndex);
};

#endif /* DEMAIN_H_ */
