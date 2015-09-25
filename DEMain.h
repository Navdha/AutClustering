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

class ClusterSort{
public:
double size;
int clusIndex;
};

namespace std {

class DEMain {
public:
  DEMain(int dim, Item** items, int itemSize, int validityIndex, Parameters param);
  ~DEMain();
  void permuteBaseArray();
  void setup(double min[],double max[]);
  void computeClustering(Individual* org);
  double calcFitness(Individual* org);
  double dist(double* x, double* y);
  double* avgDist(Individual* org);
  void selectSamples(int org, int &s1, int &s2, int &s3);
  Individual* crossover(int orgIndex, double genNum, double min[], double max[]);
  void run(double min[], double max[], string filename);
  void report(int orgIndex, int worstOrgInd, string filename);
  Individual* replacement(Individual* offspring, double min[], double max[]);
  void centroidAddition(Individual* org, ClusterSort* objClus);
  void cleanIndividual(Individual* org, double min[], double max[]);
  void calcDistBtwnItems(double min[], double max[]);
  double calcDBIndex(Individual* org);
  double calcCSIndex(Individual* org);
  double calcPBIndex(Individual* org);
  double calcSD();
  void perturbPop(int cycle, double min[], double max[]);
  void restart();
  void initializePopCycle(Individual* temp, double min[], double max[]);
  double MI(int popInd1, int popInd2, bool isFinal);
  double randIndex(int popInd1,  bool isARI);

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
  int numRepeat;
  double activationThreshold;
  Item** itemsArray;
 // vector<int>** clusters;
  DistItemCluster* nearestNeighborTriples;
  double** distItem;
  int indexForFitness;
  bool* ItemUsed;
  bool* ClusFull;
  double* avgArr;
  double* scalingArr;
  double** newClustCenters;
  double* sumArr;
  int* ItemCounter;
  int popScaleFactor;
  bool * isReplaceOrg;
  int numClasses;
 private:
  void Rand1Bin(int candIndex);
};

} /* namespace std */

#endif /* DEMAIN_H_ */
