/*
 * DEMain.cpp
 *
 *  Created on: Aug 12, 2015
 *      Author: Navdha Sah
 */

#include "DEMain.h"
#include "Heap.h"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <cassert>
#include <random>     // needed for the mersenne twister random number generator
#include <functional> // needed for bind
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <set>
using namespace std;

random_device seed;
mt19937 mersenne(seed());  //Mersenne Twister random number generator
int uniformInRange(int from, int to) { // generate a random uniformly in the range [from, to]
  uniform_int_distribution<int> uniformFT(from, to);
  return uniformFT(mersenne);
}
double uniformInRange(double from, double to) { // generate a random uniformly in the range [from, to]
  uniform_real_distribution<double> uniformFT(from, to);
  return uniformFT(mersenne);
}
//uniform01() generates a random double uniformly in the range [0,1)
auto uniform01 = bind(uniform_real_distribution<double>(0, 1), mt19937(seed()));


/*
 * holds the max and min value obtained for the fitness validity index
 * among all members of the population
 * over all generations
 */
double minDB = std::numeric_limits<double>::max();
double maxDB= std::numeric_limits<double>::min();
double minCS = std::numeric_limits<double>::max();
double maxCS= std::numeric_limits<double>::min();
double minPB = std::numeric_limits<double>::max();
double maxPB= std::numeric_limits<double>::min();
double avgDB = 0.0;
double avgCS = 0.0;
double avgPB = 0.0;
int* baseArray;
int* permuteArray;
vector<int>** origClustering;
ofstream trackFile;
ofstream trackOff;
int currentGen = 0;
double origScale;
int minimumClusSize = 2;

// to compare two triples containing the item id, cluster id and distance between item and cluster
int compare(const void *p1, const void *p2){
  DistItemCluster *elem1 = (DistItemCluster *)p1;
  DistItemCluster *elem2 = (DistItemCluster *)p2;
  if(elem1->distance < elem2->distance)
    return -1;
  else if(elem1->distance > elem2->distance)
    return 1;
  else
    return 0;
}

int compareSize(const void *p1, const void *p2){
  ClusterSort *elem1 = (ClusterSort *)p1;
  ClusterSort *elem2 = (ClusterSort *)p2;
  if(elem1->size > elem2->size)
    return -1;
  else if(elem1->size < elem2->size)
    return 1;
  else
    return 0;
}


DEMain::DEMain(int dim, Item** items, int itemSize, int validityIndex, Parameters param) {
  // TODO Auto-generated constructor stub
  popObject = new Population(param.maxNumClusters, dim, param.popScaleFactor);
  strategy = stRand1Bin; //DE algo used
  numGenerations = param.gen; //number of generations
  crossoverProbability = param.CrMax - param.CrMin; //crossover probability
  origScale = param.FScale;
  scale = param.FScale; //factoring value
  popSize = dim * param.popScaleFactor; //population size
  popScaleFactor = param.popScaleFactor; //scaling factor for size of population
  maxNumClusters = param.maxNumClusters; //maximum # of clusters
  numFeatures = dim; //number of features
  numItems = itemSize; //total # of items to cluster
  minNumClusters = param.minNumClusters; // minimum # of clusters
  activationThreshold = param.threshold; //threshold value to activate a cluster's centroid
  indexForFitness = validityIndex; //which clustering validity index used
  itemsArray = items; //data structure to hold all items
  baseArray = new int[popSize];
  permuteArray = new int[popSize];
  for(int i = 0; i < popSize; i++){
    baseArray[i]=i;
  }
  //clusters = new vector<int>*[maxNumClusters]; //holds actual clusters for every pop members(calcFitness)
  distItem = new double*[numItems]; //jagged array storing dist between all items(used for calcCSIndex and calcPBIndex)
  for (int i = 0; i < numItems - 1; i++) {
    distItem[i] = new double[numItems - 1 - i];
  }
  int size = maxNumClusters * numItems;
  nearestNeighborTriples = new DistItemCluster[size]; //used for reshuffling(reshuffle & calcFitness)
  ItemUsed = new bool[numItems](); //used in reshuffle
  ClusFull = new bool[maxNumClusters](); //used in reshuffle
  avgArr = new double[maxNumClusters]; //used in calcDB
  sumArr = new double[numFeatures]; // used in calcCS
  ItemCounter = new int[numItems]; //used in reshuffle
  newClustCenters = new double*[maxNumClusters]; //used in calcCS
  isReplaceOrg = new bool[popSize];
  numClasses = param.numClasses;
  numRepeat = param.numRepeat;
  scalingArr = new double[numFeatures];
  for (int count = 0; count < maxNumClusters; count++) {
    //clusters[count] = new vector<int>;
    newClustCenters[count] = new double[dim];
  }
}

DEMain::~DEMain() {
  // TODO Auto-generated destructor stub
  //	 delete popObject;
  for (int i = 0; i < numItems; i++){
    delete [] itemsArray[i];
    delete [] distItem[i];
  }
  for(int i = 0; i < maxNumClusters; ++i) {
    //delete clusters[i];
    delete newClustCenters[i];
  }
  delete [] itemsArray;
  delete [] distItem;
  // delete [] clusters;
  delete [] nearestNeighborTriples;
  delete [] ItemUsed;
  delete [] ClusFull;
  delete [] avgArr;
  delete [] scalingArr;
  delete [] sumArr;
  delete [] ItemCounter;
  delete [] isReplaceOrg;
  delete [] newClustCenters;
  delete [] baseArray;
}

//permute base vectors
void DEMain:: permuteBaseArray(){
  int counter = 0;
  int psize = popSize;
  int rand, temp;
  bool repeat = true;
  for(int i = 0; i < popSize; i++){
    permuteArray[i]=-1;
  }
  while(counter < popSize){
    while(repeat){
      rand = uniformInRange(0, psize -1);
      if(rand != counter) repeat = false;
    }
    permuteArray[counter] = baseArray[rand];
    //swap last index value with rand's value and reduce size of array being searched
    temp = baseArray[psize - 1];
    baseArray[psize-1] = baseArray[rand];
    baseArray[rand] = temp;
    psize -= 1;
    counter++;
  }
}

/*
 * This method calculates and stores the distance between all items
 * Used while calculating CS index and standard deviation between distances for Point Biserial index
 */
void DEMain::calcDistBtwnItems(double min[], double max[]) {
  double scaleFactor = 0.0;
  double maxFeatVal = max[0];
  for(int f = 1; f <numFeatures; f++){
    if(max[f] > maxFeatVal){
      maxFeatVal = max[f];//finding the maximum feature value
    }
  }
  for(int f = 0; f <numFeatures; f++){
    scalingArr[f] = maxFeatVal/max[f] * uniformInRange(0.7,1.0); //computing the scale factor for each attribute
  }
  for (int i = 0; i < numItems - 1; i++) {
    for (int j = i + 1; j < numItems; j++) {
      distItem[i][j - i - 1] = dist(itemsArray[i]->items, itemsArray[j]->items);     
    }  //end j for
  }  // end i for

}


/* This method takes an input two arrays min and max
 * that hold the minimum and maximum value of each attribute
 * Purpose of this method is to setup the initial population
 * return type is void
 * no output expected
 */
void DEMain::setup(double min[], double max[]) {
  //initialize chromosomes for the first time
  for (int p = 0; p < popSize; p++) {
    bool isValid = false; //determines whether our individual is valid or not based on how many clusters it forms
    double fitn = 0.0;
    Individual* temp = new Individual(maxNumClusters, numFeatures); //instantiating new individual
    while (!isValid) { //continue this till a valid individual found
      int numActiveClus = 0; //# of active cluster centers for the individual
      for (int c = 0; c < maxNumClusters; c++) {
	if(p <= 0.33*popSize) temp->activationThreshold[c] = uniformInRange(0.9,1.0);
	else if(p > 0.33*popSize && p <= 0.66*popSize) temp->activationThreshold[c] = uniformInRange(0.2,1.0);
	else temp->activationThreshold[c] = uniform01();
	if (temp->activationThreshold[c] > activationThreshold) { // based on threshold, making a centroid active or inactive
	  temp->active[c] = true;
	  numActiveClus++;
	} else
	  temp->active[c] = false;
	for (int f = 0; f < numFeatures; f++) {
	  temp->clusCenter[c][f] = uniformInRange(min[f], max[f]); //randomly creating centroid to be a cluster center			     
	}//end for f
      }//end for c
      assert(numActiveClus <= maxNumClusters); //check to determine if active clusters are more than kmax
      temp->numActiveClusters = numActiveClus; //assigning this as active cluster count in the individual object
      //code to check kmin = 2
      if (temp->numActiveClusters < minNumClusters) { //if #of active clusters is less than min #of clusters possible, then activate a cluster centroid forcefully
	int num = temp->numActiveClusters;
	int activateClusNum = uniformInRange(minNumClusters, maxNumClusters); //instead of just activating 2 random centroids, activating based on a coin toss
	//while (num < minNumClusters) {
	while (num < activateClusNum) {
	  int i = uniformInRange(0, maxNumClusters - 1);
	  if (!temp->active[i]) {
	    temp->activationThreshold[i] = uniformInRange(activationThreshold, 1.0);
	    temp->active[i] = true;
	    temp->numActiveClusters++;
	    num++;
	  }
	}
      }//end if
      computeClustering(temp);
      cleanIndividual(temp, min, max);
      fitn = calcFitness(temp); //calculate fitness of individual generated for initial population
      if (fitn != -1) { //returned -1 in case of active clusters falling below kmin
	isValid = true;
      }
    } //end while
    assert(fitn != 0.0);
    temp->rawFitness = fitn; //assign fitness
    popObject->org[p] = temp;
  }//end for p
  int bestInd = 0;
  double fitness = popObject->org[0]->rawFitness;
  for (int p = 1; p < popSize; p++) { //finding the best individual in the population
    if (fitness < popObject->org[p]->rawFitness) {
      bestInd = p;
      fitness = popObject->org[p]->rawFitness;
    }
  }
  popObject->bestOrgIndex = bestInd; //updating the index of the best individual found

}

/*
 * Input parameters: pointers to arrays that hold
 * cluster center and item features
 * This method finds the distance between cluster center and an item
 * returns the distance calculated
 */
//double DEMain::dist(double* x, double* y) {
//  double Sum = 0.0;
//  for (int f = 0; f < numFeatures; f++) {
//    Sum = Sum + pow(scalingArr[f] * (x[f] - y[f]), 2.0);
//  }
//  return sqrt(Sum);
//}

double DEMain::dist(double* x, double* y) {
  double Sum = 0.0;
  for (int f = 0; f < numFeatures; f++) {
    Sum = Sum + pow((x[f] - y[f]), 2.0);
  }
  return sqrt(Sum);
}


/*
 * calculate and return DB index
 * Parameters : Individual
 */
double DEMain::calcDBIndex(Individual* org) {
  fill_n(avgArr, maxNumClusters, 0.0);
  double sum;
  //this section of code is to compute centroids of the clusters to compute DB index
  /*	for (int c = 0; c < maxNumClusters; c++) {
	for (int f = 0; f < numFeatures; f++) {
	newClustCenters[c][f] = 0.0; //clearing old centroids
	}
	}
	for (int c = 0; c < maxNumClusters; c++) {
	fill_n(sumArr, numFeatures, 0);//clearing sums stored in array
	if (org->active[c]) {
	for (vector<int>::size_type j = 0; j != clusters[c]->size(); j++) {//go through all items in cluster
	int a = clusters[c]->at(j);//find item index stored in clusters
	double* tempItem = itemsArray[a]->items;//array to hold features corresponding to item in cluster
	for (int f = 0; f < numFeatures; f++) {
	sumArr[f] += tempItem[f];//this holds the sum of all features for all items in a single cluster to compute average later
	}
	}
	for (int f = 0; f < numFeatures; f++) {
	newClustCenters[c][f] = sumArr[f] / clusters[c]->size(); //finding new centroids for the cluster
	}
	}
	}*/
  for (int c = 0; c < maxNumClusters; c++) {
    if (org->active[c]) {
      sum = 0.0;
      for (vector<int>::size_type i = 0; i != org->clusters[c]->size(); i++) {
	int a = org->clusters[c]->at(i);
	sum += dist(itemsArray[a]->items, org->clusCenter[c]);
	//sum += dist(itemsArray[a]->items, newClustCenters[c]);
      } //end for
      avgArr[c] = (sum / org->clusters[c]->size()); //finding the intra cluster distance for all active clusters
    } //end if
  } //end outer for
  sum = 0.0;
  for (int c1 = 0; c1 < maxNumClusters; c1++) {
    double maxValue = 0.0;
    for (int c2 = c1; c2 < maxNumClusters; c2++) {
      if (c1 != c2 && org->active[c1] && org->active[c2]) {
	double temp = avgArr[c1] + avgArr[c2];
	//trackOff << "values of s_i and s_j " << avgArr[c1] << " " << avgArr[c2] << endl;
	double distanceVal = dist(org->clusCenter[c1], org->clusCenter[c2]);
	//trackOff << "distance between centroid i and j " << distanceVal << endl;
	temp /=  distanceVal; //finding R =(S_i+S_j)/d_i,j
	//  temp /= dist(newClustCenters[c1], newClustCenters[c2]);
	if (temp > maxValue)
	  maxValue = temp;
      }
      //trackOff << "Rmax calc " << maxValue << endl;
    }
    sum += maxValue;		  //finding sum(Rmax)
  }
  double avg = sum / org->numActiveClusters;
  return avg;
}

/*
 * calculate and return CS index for an organism
 * Parameters : Individual
 */
double DEMain::calcCSIndex(Individual* org) {
  /*  for (int c = 0; c < maxNumClusters; c++) {
      for (int f = 0; f < numFeatures; f++) {
      newClustCenters[c][f] = 0.0; //clearing old centroids
      }
      }*/
  double finalIntraSum = 0.0;
  double finalInterSum = 0.0;
  double csVal = 0.0;
  for (int c = 0; c < maxNumClusters; c++) {
    double intraClusSum = 0.0;
    double tempIntraDist;
    //fill_n(sumArr, numFeatures, 0);
    if (org->active[c]) {
      for (vector<int>::size_type i1 = 0; i1 != org->clusters[c]->size(); i1++) {
	int a = org->clusters[c]->at(i1);
	double maxIntraDist = numeric_limits<double>::min();
	//	double* tempItem = itemsArray[a]->items;
	//	for (int f = 0; f < numFeatures; f++) {
	//	  sumArr[f] += tempItem[f]; //to compute centroids
	//	}
	for (vector<int>::size_type i2 = 0; i2 != org->clusters[c]->size();
	     i2++) { //finding max distance between items in a cluster
	  if (i2 != i1) {
	    int b = org->clusters[c]->at(i2);
	    if (b < a) {
	      tempIntraDist = distItem[b][a - (b + 1)];
	    } else {
	      tempIntraDist = distItem[a][b - (a + 1)];
	    }
	    if (tempIntraDist > maxIntraDist) {
	      maxIntraDist = tempIntraDist;
	    }

	  }
	} //end for i2
	intraClusSum += maxIntraDist; //assigning intra cluster distance for this cluster
      } // end for i1
      finalIntraSum += (intraClusSum / org->clusters[c]->size()); //updating the average intra cluster sum for all active centroids
      //    for (int f = 0; f < numFeatures; f++) {
      //	newClustCenters[c][f] = sumArr[f] / clusters[c]->size(); //finding new centroids
      //    }
    } //endif
  } // end for c
  double interClusSum = 0.0;
  double tempInterDist;
  for (int c1 = 0; c1 < maxNumClusters; c1++) { //trying to find min distance between two clusters
    double minInterDist = numeric_limits<double>::max();
    if (org->active[c1]) {
      for (int c2 = 0; c2 < maxNumClusters; c2++) {
	if (c1 != c2 && org->active[c2]) {
	  //tempInterDist = dist(newClustCenters[c1], newClustCenters[c2]);
	  tempInterDist = dist(org->clusCenter[c1], org->clusCenter[c2]);
	  if (tempInterDist < minInterDist) {
	    minInterDist = tempInterDist;
	  }
	} //end outer if
      } //end for c2
      interClusSum += minInterDist; //assigning minimum distance between clusters
    }
  } //end for c1
  finalInterSum = interClusSum;
  csVal = finalIntraSum / finalInterSum; //formula for CS index
  return csVal;
}

/*
 * computing the standard deviation of all possible distances between items
 */
double DEMain::calcSD(){
  double sum = 0.0;
  for (int i = 0; i < numItems - 1; i++) {
    for (int j = i + 1; j < numItems ; j++) {
      sum += distItem[i][j - i - 1];
    } //end j for
  } // end i for
  double num = numItems*(numItems - 1)/2;
  double mean = sum/num;
  double sumSD = 0.0;
  for(int i = 0; i < numItems - 1; i++){
    for (int j = i + 1; j < numItems ; j++) {
      sumSD += pow((distItem[i][j - i - 1] - mean), 2.0);
    }
  }
  double SD = sqrt(sumSD/num);
  return SD;
}

/*calculates and returns the Point biserial distance
 * Parameter : Individual
 */
double DEMain::calcPBIndex(Individual* org) {
  double sd = calcSD();
  //calculate average intra group distance
  double t = numItems * (numItems - 1) / 2;
  double pbIndex;
  int countIntraGroup = 0; // this is w_d
  int countInterGroup = 0; //this is b_d
  double intraClusAvgDist = 0.0;
  double interClusAvgDist = 0.0;
  for (int c1 = 0; c1 < maxNumClusters; c1++) {
    if (org->active[c1]) {
      double sumIntraCluster = 0.0;
      double sumInterCluster = 0.0;
      int n = org->clusters[c1]->size();
      countIntraGroup += (n * (n - 1)) / 2;
      countInterGroup += (n * (numItems - n)) / 2;
      for (vector<int>::size_type i1 = 0; i1 != org->clusters[c1]->size(); i1++) {
	int a = org->clusters[c1]->at(i1);
	//find distance between elements in separate clusters; go through all elements in other clusters
	for (int cl = 0; cl < c1; cl++) {
	  if (org->active[cl]) {
	    for (vector<int>::size_type ind = 0;
		 ind != org->clusters[cl]->size(); ind++) {
	      int b = org->clusters[cl]->at(ind);
	      if (b < a) {
		sumInterCluster += distItem[b][a - (b + 1)];
	      } else {
		sumInterCluster += distItem[a][b - (a + 1)];
	      }
	    }
	  }
	}
	for (int cl = c1 + 1; cl < maxNumClusters; cl++) {
	  if (org->active[cl]) {
	    for (vector<int>::size_type ind = 0;
		 ind != org->clusters[cl]->size(); ind++) {
	      int b = org->clusters[cl]->at(ind);
	      if (b < a) {
		sumInterCluster += distItem[b][a - (b + 1)];
	      } else {
		sumInterCluster += distItem[a][b - (a + 1)];
	      }
	    }
	  }
	}
	//finding distance between items in the same cluster
	for (vector<int>::size_type i2 = 0; i2 != org->clusters[c1]->size();
	     i2++) {
	  if (i2 != i1) {
	    int b = org->clusters[c1]->at(i2);
	    if (b < a) {
	      sumIntraCluster += distItem[b][a - (b + 1)];
	    } else {
	      sumIntraCluster += distItem[a][b - (a + 1)];
	    }

	  }
	} //end for k
      } //end for j
      intraClusAvgDist += sumIntraCluster / countIntraGroup; //value of d_w
      interClusAvgDist += sumInterCluster / countInterGroup; //value of d_b
    } //end if
  } //end for i
  double totalSums = (interClusAvgDist - intraClusAvgDist);
  double t_val = t * t;
  double sqrtVal = sqrt((countIntraGroup * countInterGroup) / t_val);
  pbIndex = totalSums * sqrtVal / sd;
  return pbIndex;
}

/*
 * This method cleans up the chromosome by deactivating the tiny clusters and moving their centroids to a better position
 * It also redistributes items from the tiny clusters to other larger clusters when required
 */
void DEMain:: cleanIndividual(Individual* org, double min[], double max[]) {
  double temp_dist, min_dist;
  for (int c = 0; c < maxNumClusters; c++) { //inactivating a cluster that's empty
    if (org->active[c]) {
      if (org->clusters[c]->size() < minimumClusSize) {
	for (int f = 0; f < numFeatures; f++) { //since an inactive cluster can also contribute to crossover,
	  //we reinitialize the centroid since we know that the existing centroid is not good as it's empty
	  //reinitialize only those features which are close to the boundary
	  double gap = max[f] - min[f];
	  if (org->clusCenter[c][f] >= min[f] && org->clusCenter[c][f] <= (min[f] + 0.1 * gap)) {
	    org->clusCenter[c][f] = uniformInRange(min[f], (min[f] + (0.5 * gap)));//check if 1.5*min[f] < max[f]
	  } else if (org->clusCenter[c][f] <= max[f] && org->clusCenter[c][f] >= (max[f] - 0.1 * gap)) {
	    org->clusCenter[c][f] = uniformInRange((max[f] - (0.5 * gap)), max[f]);
	  }
	}//end for
	if (!org->clusters[c]->empty()) {
	  for (vector<int>::size_type j = 0; j != org->clusters[c]->size();j++) {
	    //go through all items in cluster
	    int addItemIndex = org->clusters[c]->at(j);
	    temp_dist = 0;
	    min_dist = numeric_limits<double>::max();
	    int clusCtr = 0;
	    int minInd = -1;
	    while (clusCtr < maxNumClusters) {
	      if (clusCtr != c && org->active[clusCtr]) {
		temp_dist = dist(org->clusCenter[clusCtr], itemsArray[addItemIndex]->items);
		if (temp_dist < min_dist) {
		  temp_dist = min_dist;
		  minInd = clusCtr;
		}
	      }
	      clusCtr++;
	    }
	    assert(minInd != -1);
	    org->clusters[minInd]->push_back(addItemIndex);
	  }   //end for
	  org->clusters[c]->clear(); // after element has been moved, cluster with single item made empty
	}
	org->active[c] = false; // deactivating the empty cluster
	org->activationThreshold[c] = uniformInRange(0.0,activationThreshold);
	org->numActiveClusters--;
      } //end if size
    } //end org->active
  } //end for
}

/*
 * Input parameters : Pointer to chromosome, index of population's chromosome,
 * bool parameter to know whether function called during initial setup or during DE
 * This method computes the clustering of a chromosome and then calls the fitness function inside
 * return  type : double; returns the fitness computed
 */
void DEMain::computeClustering(Individual* org) {
  int min_index = -1;
  double temp_dist = 0;
  int itemIndex = 0;
  while (itemIndex < numItems) { //form clusters
    min_index = -1;
    double min = numeric_limits<double>::max();
    for (int clusIndex = 0; clusIndex < maxNumClusters; clusIndex++) {
      if (org->active[clusIndex]) {
	temp_dist = dist(org->clusCenter[clusIndex],itemsArray[itemIndex]->items);
	if (temp_dist < min) { //finding the closest cluster center to a particular item
	  min = temp_dist;
	  min_index = clusIndex;
	}
      }
    }
    assert(min_index != -1);
    org->clusters[min_index]->push_back(itemIndex); //adding item to its closest cluster center
    itemIndex++;
  } // end while
}

/*
 * Input parameters : Pointer to parent and offspring, index of parent, generation number, min-max value for each feature
 * This method determines whether a parent is to be replaced by an offspring or not
 */
Individual* DEMain::replacement(Individual* org, double min[], double max[]) {
  computeClustering(org);
  cleanIndividual(org, min, max);
  double defaultFit = calcFitness(org);
  org->rawFitness = defaultFit;
  ClusterSort* objClus = new ClusterSort[maxNumClusters];
  for (int c = 0; c < maxNumClusters; c++) {
    if (org->active[c]) {
      objClus[c].size = org->clusters[c]->size();
    } else {
      objClus[c].size = 0;
    }
    objClus[c].clusIndex = c;
  }
  qsort(objClus, maxNumClusters, sizeof(ClusterSort), compareSize);
  if (objClus[0].size > 0.4 * numItems && org->numActiveClusters != maxNumClusters) {
    //create copy of original individual
    Individual* orgDup = new Individual(*org);
    centroidAddition(orgDup, objClus);
    computeClustering(orgDup);
    cleanIndividual(orgDup, min, max);
    double newFitness = calcDBIndex(orgDup);
    orgDup->rawFitness = newFitness;
    //cout << "new fit " << newFitness << " old fit " << defaultFit << endl;
    if (newFitness > defaultFit) {
      //delete org;
      return orgDup;
    }
    else {
      delete orgDup;
      return org;
    }
  }
  else{
    return org;
  }
}

/*
 * This method computes and adds centroid to the chromosome passed
 */
void DEMain::centroidAddition(Individual* org, ClusterSort* objClus) {
  double maxDist;
  double avgDistCal = 0;
  double sdSum = 0;
  maxDist = dist(org->clusCenter[objClus[0].clusIndex], org->clusCenter[objClus[1].clusIndex]);
  for (vector<int>::size_type j = 0; j != objClus[0].size; j++) {
    int a = org->clusters[objClus[0].clusIndex]->at(j); //find item index stored in clusters
    avgDistCal += dist(itemsArray[a]->items, org->clusCenter[objClus[0].clusIndex]);
  }
  avgDistCal /= objClus[0].size;
  for (vector<int>::size_type j = 0; j != objClus[0].size; j++) {
    int a = org->clusters[objClus[0].clusIndex]->at(j); //find item index stored in clusters
    sdSum += pow((dist(itemsArray[a]->items, org->clusCenter[objClus[0].clusIndex]) - avgDistCal), 2.0);
  }
  sdSum /= objClus[0].size;
  sdSum = sqrt(sdSum);
  // cout << "avg dist " << avgDistCal << " sdSum " << sdSum << endl;
  double newScale;
  //need to check if the centroid is too close or far from the larger cluster
  if ((maxDist / 2) < (avgDistCal + sdSum)) {
    //centroid needs to be close
    newScale = (avgDistCal) / maxDist;
  } else {
    //to calculate new centroid
    newScale = (avgDistCal + sdSum) / maxDist;
  }
  double* newCentroid = new double[numFeatures];
  //    double* midCentroid = new double[numFeatures];
  for (int f = 0; f < numFeatures; f++) {
    newCentroid[f] = org->clusCenter[objClus[0].clusIndex][f] + newScale*(org->clusCenter[objClus[1].clusIndex][f]- org->clusCenter[objClus[0].clusIndex][f]);
    //      midCentroid[f] = (org->clusCenter[objClus[1].clusIndex][f] + org->clusCenter[objClus[0].clusIndex][f])/2;
    //      cout << "new centroid " << newCentroid[f] << " ";
  }
  //    cout << endl;
  double min = numeric_limits<double>::max();
  int minIndex = -1;
  //    cout << "active centroids " << org->numActiveClusters << endl;
  for (int c = 0; c < maxNumClusters; c++) {
    if (org->active[c] == false) {
      //	cout << "false" << endl;
      //find index of centroid that's closest to new centroid
      double tempDist = dist(org->clusCenter[c], newCentroid);
      //	cout << "temp dist after finding centroid " << tempDist << endl;
      if (tempDist < min) {
	min = tempDist;
	minIndex = c;
      }
    }
  }		//end for
  //replace centroid at minIndex with the new centroid. recompute clustering, clean up and calc. fitness
  assert(minIndex != -1);
  org->active[minIndex] = true;
  org->activationThreshold[minIndex] = uniformInRange(activationThreshold, 1.0);
  for (int f = 0; f < numFeatures; f++) {
    org->clusCenter[minIndex][f] = newCentroid[f];
  }
}

/*
 * Input parameters : Pointer to chromosome, index of population's chromosome,
 * bool parameter to know whether function called during initial setup or during DE
 * This method calculates the fitness of a chromosome and returns it
 * return  type : double
 */
double DEMain::calcFitness(Individual* org) {
  double fit = 0.0;
  if (org->numActiveClusters >= minNumClusters) { //if after deactivating cluster centers above, an individual's active cluster centers fall below 2(kmin), assign it lowest fitness
    //based on the validity index find DB, CS or PB index for the individual to compute the fitness
    if (indexForFitness == 1) {
      //minimization problem
      double dBValue = calcDBIndex(org);
      fit = 1 / dBValue;
    } //end isDB if
    else if (indexForFitness == 2) {
      //calculate CS index; minimization problem
      double csVal = calcCSIndex(org);
      fit = 1 / csVal;
    } //end else
    else {
      double pbVal = calcPBIndex(org);
      fit = pbVal;//maximization problem
    }
    return fit * 100;
  } else {
    return -1;
  }
	
}

/*
 * Input parameters : indexes for chromosomes
 * This method makes sure that we get unique indexes for performing crossover
 * return type: void
 */

void DEMain::selectSamples(int orgIndex, int &s1, int &s2, int &s3) { // finding unique indexes to perform crossover
  //  do {
  //    s1 = uniformInRange(0, popSize-1);
  //  } while (s1 == orgIndex);
  s1 = permuteArray[orgIndex];
  do {
    s2 = uniformInRange(0, popSize-1);
  } while ((s2 == orgIndex) || (s2 == s1));

  do {
    s3 = uniformInRange(0, popSize-1);
  } while ((s3 == orgIndex) || (s3 == s2) || (s3 == s1));
  return;
}

/*
 * Input parameters: index for chromosome chosen and the generation number
 * This method performs crossover to create an offspring
 * returns pointer to offspring created
 */
Individual* DEMain::crossover(int orgIndex, double genNum, double min[], double max[]) {
  double fit = 0.0;
  double f_scale,d;
  int s1, s2, s3;
  bool isExploration = true;
  if(genNum >= 0.9*numGenerations){
    isExploration = false;
  }
  double crossoverProb = crossoverProbability * ((numGenerations - genNum) / numGenerations);//based on formula in paper
  if(genNum == currentGen+100){
    scale -= 0.1;
    currentGen = genNum;
    //     cout << scale << " and gen " << genNum << endl;
  }
  if(isExploration) {
    d = 0.5;
    f_scale = scale + (d*(uniform01() - 0.5)); //based on formula in paper
  }
  else { 
    d = 0.001;
    crossoverProb = uniformInRange(0.8,0.98);
  }//need to discuss

  selectSamples(orgIndex, s1, s2, s3);//to get unique individuals to perform DE crossover on
  Individual* child = new Individual(maxNumClusters, numFeatures);
  int counterActiveClus = 0;//keeps track of active cluster centers
  for (int c = 0; c < maxNumClusters; c++) {
    bool change = uniform01() < crossoverProb? true : false;
    if(uniform01() < crossoverProb){
      child->activationThreshold[c] =  popObject->org[s1]->activationThreshold[c] +
	f_scale*(popObject->org[s2]->activationThreshold[c] -  popObject->org[s3]->activationThreshold[c]);    
    }
    else{
      child->activationThreshold[c] =  popObject->org[orgIndex]->activationThreshold[c];               
    }
    for (int f = 0; f < numFeatures; f++) {//binomial crossover computation for centroids
      assert(popObject->org[orgIndex] != NULL);
      //jitter
      if(!isExploration) { 
	f_scale = scale + (d*(uniform01() - 0.5));
      }
      if(change) {
	child->clusCenter[c][f] = popObject->org[s1]->clusCenter[c][f]
	  + f_scale*(popObject->org[s2]->clusCenter[c][f]- popObject->org[s3]->clusCenter[c][f]);
	if (child->clusCenter[c][f] < min[f]) {
	  child->clusCenter[c][f] =  popObject->org[s1]->clusCenter[c][f]
	    + abs( f_scale*(popObject->org[s2]->clusCenter[c][f]- popObject->org[s3]->clusCenter[c][f]) );
	}
	else if (child->clusCenter[c][f] > max[f]){//if a feature for cluster center computed goes out of its bounds, assign it a random value
	  child->clusCenter[c][f] =  popObject->org[s1]->clusCenter[c][f]
	    - abs( f_scale*(popObject->org[s2]->clusCenter[c][f]- popObject->org[s3]->clusCenter[c][f]));
	}
      } else {	
	child->clusCenter[c][f] = popObject->org[orgIndex]->clusCenter[c][f];
      }
      if ((child->clusCenter[c][f] < min[f]) || (child->clusCenter[c][f] > max[f])){
	child->clusCenter[c][f] = uniformInRange(min[f], max[f]);
      } 
     
    }//for i
    if (child->activationThreshold[c] > 1 || child->activationThreshold[c] < 0)//if a threshold value after crossover goes out of its bounds, recompute it
      child->activationThreshold[c] = uniform01();
    if (child->activationThreshold[c] > activationThreshold) {
      child->active[c] = true;
      counterActiveClus++;
    } else
      child->active[c] = false;
  }//for j
  assert(counterActiveClus <= maxNumClusters);
 
  child->numActiveClusters = counterActiveClus;
  if(child->numActiveClusters < minNumClusters){//if #of active cluster centers below kmin, forcefully make more clusters active till kmin satisfied
    int num = child->numActiveClusters;
    // int randNumClusters = uniformInRange(1, maxNumClusters/2);
    // while(num < randNumClusters) {
    while (num < minNumClusters) {
      int i = uniformInRange(0, maxNumClusters - 1);
      if(!child->active[i]){
	child->activationThreshold[i] = uniformInRange(activationThreshold, 1.0);
	child->active[i] = true;
	child->numActiveClusters++;
	num++;
      }
    }
  }
  return child;

}

/*
 * Reinitialize variables in between cycle
 */
void DEMain::restart(){
  scale = origScale;
  currentGen = 0;
  minimumClusSize = 2;
}

void DEMain::initializePopCycle(Individual* temp, double min[], double max[]){
  int numActiveClus = 0; //# of active cluster centers for the individual
  double fitn;
  for (int c = 0; c < maxNumClusters; c++) {
    temp->activationThreshold[c] = uniform01();
    if (temp->activationThreshold[c] > activationThreshold) { // based on threshold, making a centroid active or inactive
      temp->active[c] = true;
      numActiveClus++;
    } else
      temp->active[c] = false;
    //int randomVal = uniformInRange(0, numItems-1);
    for (int f = 0; f < numFeatures; f++) {
      // temp->clusCenter[c][f] = itemsArray[randomVal]->items[f];
      temp->clusCenter[c][f] = uniformInRange(min[f], max[f]); //randomly creating centroid to be a cluster center
    } //end for f
  } //end for c
  assert(numActiveClus <= maxNumClusters); //check to determine if active clusters are more than kmax
  temp->numActiveClusters = numActiveClus; //assigning this as active cluster count in the individual object
  //code to check kmin = 2
  if (temp->numActiveClusters < minNumClusters) { //if #of active clusters is less than min #of clusters possible, then activate a cluster centroid forcefully
    int num = temp->numActiveClusters;
    int activateClusNum = uniformInRange(minNumClusters, maxNumClusters); //instead of just activating 2 random centroids, activating based on a coin toss
    //while (num < minNumClusters) {
    while (num < activateClusNum) {
      int i = uniformInRange(0, maxNumClusters - 1);
      if (!temp->active[i]) {
	temp->activationThreshold[i] = uniformInRange(
						      activationThreshold, 1.0);
	temp->active[i] = true;
	temp->numActiveClusters++;
	num++;
      }
    }
  }		  //end if
  computeClustering(temp);
  cleanIndividual(temp, min, max);
  fitn = calcFitness(temp); //calculate fitness of individual generated for initial population
  temp->rawFitness = fitn;
}

/*
 *during initial cycles, let crossover restart
 *towards later stages, based on the avg/best change in population, add random chromosomes
 */
void DEMain::perturbPop(int cycle, double min[], double max[]){
  restart();
  double fitn = 0.0;
  if(cycle > 0 && cycle < 0.4*numRepeat){
    //currently replacing a small portion of population irrespective of avg/best
    //replace 20% of population randomly
    int repl = 0;
    while(repl < 0.2*popSize){
      int indexRepl = uniformInRange(0,popSize-1);
      if(indexRepl != popObject->bestOrgIndex){
	//create new offspring and replace
	Individual* temp = new Individual(maxNumClusters, numFeatures); //instantiating new individual
	initializePopCycle(temp, min, max);
	delete popObject->org[indexRepl];//could archive these
	popObject->org[indexRepl] = temp;
	repl++;
      }      
    }
  }
  else if(cycle >= 0.4*numRepeat && cycle < 0.9*numRepeat){
    //currently replacing a small portion of population irrespective of avg/best
    //replace 10% of population randomly
    int repl =0;
    while(repl < 0.1*popSize){
      int indexRepl = uniformInRange(0,popSize-1);
      if(indexRepl != popObject->bestOrgIndex){
	//create new offspring and replace
	Individual* temp = new Individual(maxNumClusters, numFeatures); //instantiating new individual
	initializePopCycle(temp, min, max);
	delete popObject->org[indexRepl];
	popObject->org[indexRepl] = temp;
	repl++;
      }
      
    }
  }
  else if(cycle >= 0.9*numRepeat && cycle < numRepeat){
    //replace fewer elements from population
    //no changes in the last cycle
    int repl =0;
    while(repl < 0.05*popSize){
      int indexRepl = uniformInRange(0,popSize-1);
      if(indexRepl != popObject->bestOrgIndex){
	//create new offspring and replace
	Individual* temp = new Individual(maxNumClusters, numFeatures); //instantiating new individual
	initializePopCycle(temp, min, max);
	delete popObject->org[indexRepl];
	popObject->org[indexRepl] = temp;
	repl++;
      }//end if 
    }//end while
  }//end else if

}


/*
 * This method runs the DE algorithm
 */
void DEMain::run(double min[], double max[], string filename) {
  double g = 0;
  fill_n(isReplaceOrg, popSize, false);
  double fitness;
  try {
    trackOff.open("record.txt");
    trackOff << setiosflags(ios::left);
    trackOff << setw(5) << "Gen" << setw(5) << "Avg DB" << "|" << setw(10) << "Max DB"  << "|" << setw(10) << "Min DB" << endl;
    for(int cycle = 0; cycle < numRepeat; cycle++){
      while (g < numGenerations) {    //till max generation reached
	permuteBaseArray();
	minDB = numeric_limits<double>::max();
	maxDB = numeric_limits<double>::min();
	/* minCS = numeric_limits<double>::max();
	   maxCS = numeric_limits<double>::min(); 
	   minPB = numeric_limits<double>::max();
	   maxPB = numeric_limits<double>::min();*/ 
	avgDB = 0.0, avgCS = 0.0, avgPB = 0.0;
	Population* newpop = new Population(maxNumClusters, numFeatures, popScaleFactor);
	for (int p = 0; p < popSize; p++) {
	  Individual *child, *offspring;
	  child = crossover(p, g, min, max); //generate an offspring my performing DE crossover
	  offspring = replacement(child, min, max);
	  double avg;
	  if (popObject->org[p]->rawFitness <= offspring->rawFitness) { //if offspring better than parent, replace the parent
	    isReplaceOrg[p] = true;
	    newpop->org[p] = offspring;
	    avg = 100*(1/offspring->rawFitness);
	    avgDB += avg;
	    if (minDB > avg) {//keeping track of the minimum and maximum DB index encountered for all individuals in all generations
	      minDB = avg;
	    }
	    if (maxDB < avg) {
	      maxDB = avg;
	    }	  
	  } else { //otherwise delete the offspring
	    isReplaceOrg[p] = false;
	    delete offspring;
	    newpop->org[p] = popObject->org[p];
	    avg = 100*(1/popObject->org[p]->rawFitness);
	    avgDB += avg;
	    if (minDB > avg) {//keeping track of the minimum and maximum DB index encountered for all individuals in all generations
	      minDB = avg;
	    }
	    if (maxDB < avg) {
	      maxDB = avg;
	    }
	  }
	}
	trackOff << setw(5) << g << setw(5) << avgDB << "|" << setw(10) << maxDB  << "|" << setw(10) << minDB << endl;
	for (int p = 0; p < popSize; p++) { //based on the parents that are supposed to be replaced, freeing the space occupied in memory
	  if (isReplaceOrg[p]) {
	    delete popObject->org[p];
	  }
	}
	delete[] popObject->org;
	//delete popObject;
	popObject = newpop; //new population generated becomes the current population
	int bestInd = 0;
	fitness = popObject->org[0]->rawFitness;
	for (int p = 1; p < popSize; p++) { //keeping track of the best member of the population
	  if (fitness < popObject->org[p]->rawFitness) {
	    bestInd = p;
	    fitness = popObject->org[p]->rawFitness;
	  }
	}
	popObject->bestOrgIndex = bestInd;
	g++; //increment generation number
      }//end while
      //cause perturbation in population
      g=0;
      perturbPop(cycle,min,max);
    }
    cout << endl;
    //cout << "Out of bounds " << countOOB << endl;
    cout << "Stopped at generation " << g << endl;
    //find worst chromosome
    int worstInd = 0;
    fitness = popObject->org[0]->rawFitness;
    for (int p = 1; p < popSize; p++) { //for testing purposes
      //finding the worst member in the population at the end
      if (fitness > popObject->org[p]->rawFitness) {
	worstInd = p;
	fitness = popObject->org[p]->rawFitness;
      }
    }
    trackOff.close();
    report(popObject->bestOrgIndex, worstInd, filename);
  } catch (exception& e) {
    cerr << e.what() << endl;
  }
}

/*
 * Input parameter : index for the best chromosome of the population
 * This method outputs the clusters for the best chromosome in the population to a file
 * return type : void
 */
void DEMain::report(int bestPopIndex, int worstInd, string filename) {
  ofstream outputFile;
  ofstream centroids;
  centroids.open("bestCentroids");
  outputFile.open(filename.c_str());
  outputFile << "Parameters used to train this dataset were :" << endl;
  outputFile << "------------------------------------------------" << endl;
  outputFile << "Scale Factor : " << scale << endl;
  outputFile << "Crossover probability : " << crossoverProbability << endl;
  outputFile << "Maximum generations : " << numGenerations << endl;
  outputFile << "Population scaling factor : " << popScaleFactor << endl;
  outputFile << "Activation Threshold : " << activationThreshold << endl;
  outputFile << "Number of cycles : " << numRepeat << endl;
  outputFile << "------------------------------------------------" << endl;
  outputFile << endl;
  outputFile << "The final clusters obtained are:" << endl;
  centroids << "Centroids for best cluster obtained are : " << endl;
  int clusIndex = -1;

  int* clusClass = new int[numClasses+1];
  for(int c = 0; c < numClasses+1; c++){
    clusClass[c] = 0;
  }
  Individual* org = popObject->org[bestPopIndex];

  //original clustering for the data set to compute nmi
  origClustering = new vector<int>*[numClasses+1];
  for (int count = 0; count <= numClasses; count++) {
    origClustering[count] = new vector<int>;
  }
  for(int i =0; i < numItems; i++) {
    int index = itemsArray[i]->typeClass;
    origClustering[index]->push_back(i);
  }
  double valNMI = MI(bestPopIndex, 0, true);
  outputFile << "NMI for best cluster is " << valNMI << endl;

  double* arr;
  int activeCount = 0;
  for (int c = 0; c < maxNumClusters; c++) {
    if (org->active[c]) {			//prints out the best output
      centroids << "Centroid at index " << c << " of size " << org->clusters[c]->size() << endl;
      for(int f=0;f < numFeatures;f++){
	centroids << org->clusCenter[c][f] << " " ;
      }
      centroids <<endl;
      activeCount++;
      for (vector<int>::size_type i = 0; i != org->clusters[c]->size(); i++) {
	int itemIndex = org->clusters[c]->at(i);
	clusClass[itemsArray[itemIndex]->typeClass]++;
      }
      outputFile << "Elements of cluster " << c << " :" << endl;
      outputFile << "Total # of items in cluster " << org->clusters[c]->size() << endl;
      for(int cl = 0; cl < numClasses+1; cl++){
	if(clusClass[cl] != 0){
	  outputFile << "Class " << cl << " : " << clusClass[cl] << " number of items";
	  outputFile << endl;
	}
      }//end printing for
    }//end if
    for(int cl = 0; cl < numClasses+1; cl++){
      clusClass[cl] = 0;
    }
  }

  outputFile << "Total number of clusters obtained : " << activeCount << endl;


  org = popObject->org[worstInd];
	
  outputFile << "-------------------------------------------" << endl;
  outputFile << "Worst fitness clusters are as follows" << endl;
  outputFile << endl;
  centroids << "Centroids for worst chromosome are : " << endl;
  activeCount = 0;
  for (int c = 0; c < maxNumClusters; c++) {//prints out the worst output after changes are made
    if (org->active[c]) {
      activeCount++;
      centroids << "Centroid at index " << c << " of size " << org->clusters[c]->size() << endl;
      for(int f=0;f < numFeatures;f++){
	centroids << org->clusCenter[c][f] << " " ;
      }
      centroids <<endl;
      for (vector<int>::size_type i = 0; i != org->clusters[c]->size(); i++) {
	int itemInd = org->clusters[c]->at(i);
	clusClass[itemsArray[itemInd]->typeClass]++;
      }
      outputFile << "Elements of cluster " << c << " :" << endl;
      outputFile << "Total # of items in cluster " << org->clusters[c]->size() << endl;
      for(int cl = 0; cl < numClasses+1; cl++){
	if(clusClass[cl] != 0){
	  outputFile << "Class " << cl << " : " << clusClass[cl] << " number of items";
	  outputFile << endl;
	}
      }//end printing for
    }//end if
    for(int cl = 0; cl < numClasses+1; cl++){
      clusClass[cl] = 0;
    }
  }
  outputFile << "Total number of clusters obtained : " << activeCount << endl;
  valNMI = MI(worstInd, 0, true);
  outputFile << "NMI for worst cluster is " << valNMI << endl;

  outputFile << "DB index for worst " << 1.0/popObject->org[worstInd]->rawFitness <<endl;
  if (indexForFitness == 1) {
    outputFile << "Min DB index is " << minDB << " Max DB index is "
	       << maxDB << endl;
  } else if (indexForFitness == 2) {
    outputFile << "Min CS index is " << minCS << " Max CS index is "
	       << maxCS << endl;
  } else {
    outputFile << "Min PB index is " << minPB << " Max PB index is "
	       << maxPB << endl;
  }
  outputFile.close();
  cout << "Result saved in file.";

  delete[] arr;
}

/*
 * This method computes the mutual information based scores
 * It is used to compare our result to the ground truth or to compare two clusters
 * Returns value b/w [0,1] 0 - purely independent clusters; 1 - same clusters 
 */
double DEMain::MI(int popInd1, int popInd2, bool isFinal) {
  double H1 = 0.0;
  double H2 = 0.0;
  double mi = 0.0;
  double emi = 0.0;
  double p1, p2, p_12, mi_log;
  double itemSize = numItems;
  Individual* org1 = popObject->org[popInd1];
  double size;
  bool* isSorted;// = new bool[maxNumClusters]();

  for (int c1 = 0; c1 < maxNumClusters; c1++) {   
    if (org1->active[c1]) {
      sort(org1->clusters[c1]->begin(), org1->clusters[c1]->end());
      size = org1->clusters[c1]->size();
      p1 = size / itemSize;
      //	cout << "P1 is " << p1 << " ";
      H1 += p1 * (log2(p1));//H_u = \sum(P(i)log(P(i))) to calc entropy in clustering for 1
      if(!isFinal) {   
	isSorted = new bool[maxNumClusters]();
	Individual* org2 = popObject->org[popInd2];
	for (int c2 = 0; c2 < maxNumClusters; c2++) {
	  if (org2->active[c2]) {
	    vector<int> clus3;
	    if(!isSorted[c2]) {
	      sort(org2->clusters[c2]->begin(), origClustering[c2]->end());
	      isSorted[c2] = true;
	    }
	    set_intersection(org1->clusters[c1]->begin(), org1->clusters[c1]->end(),
			     origClustering[c2]->begin(), origClustering[c2]->end(),
			     back_inserter(clus3));
	    p_12 = clus3.size() / numItems;
	    p2 = origClustering[c2]->size() / numItems;
	    mi_log = p_12 / (p1 * p2);	     
	    mi += p_12 * (log2(mi_log));
	  }//end if
	}//end for c2
      }//end if isFinal 
      else{
	int it1 = 0, it2 = 1;
	isSorted = new bool[numClasses+1]();
	for (int c2 = 1; c2 < numClasses+1; c2++) {
	  //    if(c2 != 4) {
	  vector<int> clus3;
	  if(!isSorted[c2]) {	   
	    sort(origClustering[c2]->begin(), origClustering[c2]->end());
	    isSorted[c2] = true;
	  }
	  //clusters[c1]->begin(), clusters[c1]->end(),
	  //origClustering[c2]->begin(), origClustering[c2]->end(),
	  //back_inserter(clus3));
	  double size2 = origClustering[c2]->size(); 
	  double ctr = 0;
	  while(it1 < size && it2 < size2) {
	    if(org1->clusters[c1]->at(it1) == origClustering[c2]->at(it2)){
	      ctr++;
	      it1++; it2++;
	    }
	    else if(org1->clusters[c1]->at(it1) > origClustering[c2]->at(it2)) it2++;
	    else it1++;
	  }
	  if(ctr != 0) {
	    p_12 = ctr / itemSize;
	    p2 = size2 / itemSize;
	    mi_log = p_12 / (p1 * p2);
	    mi += p_12 * (log2(mi_log));      
	  }
	  // }//end if c2!=4
	}//end for c2
      }//end else
    }//end if
      
    if (!isFinal && popObject->org[popInd2]->active[c1]) {
      p2 = origClustering[c1]->size() / numItems;
      H2 += p2 * (log2(p2));//H_v = \sum(P(j)log(P(j))) to calc entropy in clustering for 2
    }
     
  }//end for
  if(isFinal){
    for(int c1 =1; c1 < numClasses+1; c1++) {
      double size3 = origClustering[c1]->size();
      p2 = size3 / itemSize;
      if(p2!=0)
	H2 += p2 * (log2(p2));//H_v = \sum(P(j)log(P(j))) to calc entropy in clustering for 2 
    }
  }
  cout << "MI " << mi << " " << " H1 " << H1 << " h2 " << H2 << endl;
  double nmi = mi/sqrt(H1*H2);
  // double max = (H1>H2) ? H1 : H2;
  // double ami = (mi - emi)/(max - emi);
  delete [] isSorted;
  return nmi;
}

/*
 *Rand index's value ranges between [0,1]
 *Higher the value, more similar the clusters
 */

double DEMain::randIndex(int popInd1, bool isARI){
  int a=0; //number of pairs of points with the same label in C and  assigned to the same cluster in K
  int b=0; //the number of pairs with the same label, but in different clusters
  int c=0; //number of pairs in the same cluster, but with different class labels
  int d=0; //number of pairs with different label and different cluster
  double rInd = 0;
  int* trackerArray = new int[numItems];
  for(int c = 0; c < maxNumClusters; c++){
    if(popObject->org[popInd1]->active[c]){
      for (vector<int>::size_type i = 0; i != popObject->org[popInd1]->clusters[c]->size(); i++) {
	int a = popObject->org[popInd1]->clusters[c]->at(i);
	trackerArray[a] = c;
      }
    }
  }
  int clusIndex1, clusIndex2;
  for (int i = 0; i < numItems-1; i++) {
    for(int j = i+1; j < numItems; j++) {
      clusIndex1 = trackerArray[i];
      clusIndex2 = trackerArray[j];
      if (itemsArray[i]->typeClass == itemsArray[j]->typeClass) {
	if(clusIndex1 == clusIndex2) a++; //points have same label and assigned to same cluster
	else b++; //points have same label but assigned to diff cluster
      }
      else{
	if(clusIndex1 == clusIndex2) c++; //points have diff label but assigned to same class
	else d++; //points have diff label and assigned to diff class
      }
    }
  }//forming clusters
  assert((a+b+c+d) == numItems);
  if(!isARI) {
    rInd = (a+d)/numItems;
  }
  else{
    double ac = a + c;
    double ab = a + b;
    double num = a - ((ac*ab)/numItems);
    double den = ((ac+ab)/2) - ((ac*ab)/numItems);
    rInd = num/den;
  }

  return rInd;
}
