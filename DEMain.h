/*
 * DEMain.h
 *
 *  Created on: Aug 12, 2015
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

class DistItemCluster {
 public:
  double distance;
  int itemIndex;
  int clustIndex;
};

namespace std {

class DEMain {
public:
  DEMain(int dim, int** tracebackArr, Item** items, int itemSize, int validityIndex, Parameters param);
  ~DEMain();
  void setup(double min[],double max[]);
  double calcFitness(Individual* org, int popIndex, bool isInitial, int genNum, double min[], double max[]);
  double dist(double* x, double* y);
  double* avgDist(Individual* org);
  void selectSamples(int org, int &s1, int &s2, int &s3);
  Individual* crossover(int orgIndex, double genNum, double min[], double max[]);
  void run(double min[], double max[], string filename);
  void report(int orgIndex, int worstOrgInd, string filename);
  void reshuffle(Individual* org, int numTriplesArray, int orgIndex,  bool isInitial);
  void calcDistBtwnItems();
  double calcDBIndex(Individual* org);
  double calcCSIndex(Individual* org);
  double calcPBIndex(Individual* org);
  double calcSD();

  Population* popObject;
  int strategy;
  double scale;
  double crossoverProbability;
  double numGenerations;
  int popSize;
  int maxNumClusters;
  int minNumClusters;
  int numFeatures;
  int numItems;
  double activationThreshold;
  int** trackerArray;
  Item** itemsArray;
  vector<int>** clusters;
  DistItemCluster* nearestNeighborTriples;
  int* offspringArray;
  double** distItem;
  int indexForFitness;
  bool* ItemUsed;
  bool* ClusFull;
  double* avgArr;
  double** newClustCenters;
  double* sumArr;
  int* ItemCounter;
  int popScaleFactor;
  bool * isReplaceOrg;
 private:
  void Rand1Bin(int candIndex);
};

} /* namespace std */

#endif /* DEMAIN_H_ */
