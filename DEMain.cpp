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
#include <sstream>
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
using namespace std::placeholders;

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
int* permuteArray;
vector<int>** origClustering;
ofstream trackFile;
ofstream trackOff;
int lastUpdatedGen = 0;
int minimumClusSize = 2;
double origScale;
bool isBestInd = false;
int indCentAdded;
double sd;
int counterTwoCent = 0;
int counterThreeCent = 0;
double gCProb = 0.0;
double gFScale = 0.0;
bool trackChange[20];
double probWorse = 0.2;
bool acceptWorse = false;
bool printTrackArr =false;

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

bool compareDist(const Clustering& obj1, const Clustering& obj2){
  return (obj1.distance < obj2.distance)? true: false;
}

bool compareItems(const Clustering& obj1, const Clustering& obj2){
  return (obj1.itemIndex < obj2.itemIndex) ? true: false;
}

DEMain::DEMain(int dim, Item** items, int itemSize, int validityIndex, Parameters param) {
  // TODO Auto-generated constructor stub
  popObject = new Population(param.maxNumClusters, dim, param.popScaleFactor);
  strategy = stRand1Bin; //DE algo used
  numGenerations = param.gen; //number of generations
  crossoverProbability = param.CrMax - param.CrMin; //crossover probability
  scale = param.FScale; //factoring value
  origScale = scale;
  popSize = dim * param.popScaleFactor; //population size
  popScaleFactor = param.popScaleFactor; //scaling factor for size of population
  maxNumClusters = param.maxNumClusters; //maximum # of clusters
  numFeatures = dim; //number of features
  numItems = itemSize; //total # of items to cluster
  minNumClusters = param.minNumClusters; // minimum # of clusters
  activationThreshold = param.threshold; //threshold value to activate a cluster's centroid
  indexForFitness = validityIndex; //which clustering validity index used
  itemsArray = items; //data structure to hold all items
  permuteArray = new int[popSize];
  for(int i = 0; i < popSize; i++){
    permuteArray[i]=i;
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
  sd = calcSD();//should only be called once.
}

DEMain::~DEMain() {
  // TODO Auto-generated destructor stub
  //	 delete popObject;
 // for (int i = 0; i < numItems; i++){
 //   delete [] itemsArray[i];
 //   delete [] distItem[i];
 // }
 // for(int i = 0; i < maxNumClusters; ++i) {
 //   //delete clusters[i];
 //   delete newClustCenters[i];
 // }
 // delete [] itemsArray;
 // delete [] distItem;
 // // delete [] clusters;
 // delete [] nearestNeighborTriples;
 // delete [] ItemUsed;
 // delete [] ClusFull;
 // delete [] avgArr;
 // delete [] scalingArr;
 // delete [] sumArr;
 // delete [] ItemCounter;
 // delete [] isReplaceOrg;
 // delete [] newClustCenters;
 // delete [] permuteArray;
}

//generating a permuted array for selection of base donor
void DEMain:: permuteBaseArray(){
  int rand, temp;
  for( int p = 0; p < popSize; p++){
    rand = uniformInRange(0, (popSize - p -1));
    temp = permuteArray[p];
    permuteArray[p] = permuteArray[p+rand];
    permuteArray[p+rand] = temp;
  }
}

//generating a permuted array for selection of items inside a cluster
int* DEMain:: permuteItems(Individual* org, int index){
  int rand, temp;
  int sizeArr = org->clusters[index]->size();
  int* arr = new int[sizeArr];
  for (vector<int>::size_type i = 0; i != sizeArr; i++) {
    arr[i] = (org->clusters[index]->at(i)).itemIndex;
  }
  for( int p = 0; p < sizeArr; p++){
    rand = uniformInRange(0, (sizeArr - p -1));
    temp = arr[p];
    arr[p] = arr[p+rand];
    arr[p+rand] = temp;
  }
  return arr;
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

double DEMain::findDist(int a, int b){
  double dist = 0.0;
  if (b < a) {
    dist = distItem[b][a - (b + 1)];
  } else {
    dist = distItem[a][b - (a + 1)];
  }
  //  cout << dist;
  return dist;
}

int DEMain::callBinarySearch(Individual* org, double d, int size, int index){
  int start = 0;
  int end = size -1;
  //  cout << "disst " << d << endl;
  while(start < end) {
    int mid = (start + end)/2;
    double distance = (org->clusters[index]->at(mid)).distance;
    if(d < distance) 
      end = mid - 1;
    else if(d > distance) 
      start = mid + 1;
    else 
      return mid;
  }

  if(start == end)  return end;
  else return start;
  // return -1;
}

int DEMain::callBinarySearch(DistItemCluster* obj, double dist, int size){
  int start = 0;
  int end = size -1;

  while(start < end) {
    int mid = (start + end)/2;

    if(dist < obj[mid].distance) 
      end = mid - 1;
    else if(dist > obj[mid].distance) 
      start = mid + 1;
    else 
      return mid;
  }
  return start;
  //return -1;
} 

/*
 * This method calculates and stores the distance between all items
 * Used while calculating CS index and standard deviation between distances for Point Biserial index
 */
void DEMain::calcDistBtwnItems(double min[], double max[]) {
  double maxFeatVal = max[0];
  for(int f = 1; f <numFeatures; f++){
    if(max[f] > maxFeatVal){
      maxFeatVal = max[f];//finding the maximum value of all features
    }
  }
  for(int f = 0; f <numFeatures; f++){
    scalingArr[f] = maxFeatVal/max[f] * uniformInRange(0.7,1.0); //computing the scale factor for each attribute(need to change)
    //the way it's done now can modify the max value to be smaller. Maybe keep the max value unchanged
  }
  for (int i = 0; i < numItems - 1; i++) {
    for (int j = i + 1; j < numItems; j++) {
      distItem[i][j - i - 1] = dist(itemsArray[i]->items, itemsArray[j]->items);     
    }  //end j for
  }  // end i for

}

bool DEMain::validityCheck(Individual* org){
  //if chromosome has less than minimum # of active centroids
  //if a centroid forms a tiny cluster
  if(org->numActiveClusters < minNumClusters) return false;
  else{
    for(int c = 0; c < maxNumClusters; c++){
      if(org->clusters[c]->size() < minimumClusSize){
	return false;
      }
    }
    return true;
  }
}

void DEMain::printClusters(Individual* org){
  int* clusClass = new int[numClasses+1];
  for(int c = 0; c < numClasses+1; c++){
    clusClass[c] = 0;
  }
  int activeCount = 0;
  for (int c = 0; c < maxNumClusters; c++) {
    if (org->active[c]) {			//prints out the best output
      trackFile << "Centroid [" << c << "][" << org->clusters[c]->size() << "] \t";
      for(int f=0;f < numFeatures;f++){
	trackFile << org->clusCenter[c][f] << " " ;
      }
      trackFile << endl;
      activeCount++;
      for (vector<int>::size_type i = 0; i != org->clusters[c]->size(); i++) {
	int itemIndex = (org->clusters[c]->at(i)).itemIndex;
	clusClass[itemsArray[itemIndex]->typeClass]++;
      }
      trackFile << "Cluster[" << c << "] \t";
      //trackFile << "Total # of items in cluster " << org->clusters[c]->size() << endl;
      trackFile << "[" ;
      for(int cl = 0; cl < numClasses+1; cl++){
	if(clusClass[cl] != 0){
	  trackFile <<  clusClass[cl] << ",";
	}
	else trackFile << "," ;
      }//end printing for
      trackFile << "]" << endl;
    }//end if
    for(int cl = 0; cl < numClasses+1; cl++){
      clusClass[cl] = 0;
    }
  }
  trackFile << "Total[" << activeCount << "]" << endl;
  trackFile << "Fitness[" << org->rawFitness << "]" << endl;
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
//	if(p <= 0.33*popSize) auto coinToss = bind(uniformInRange(0.9, 1.0));
//	else if(p > 0.33*popSize && p <= 0.66*popSize) auto coinToss = bind(uniformInRange(0.3, 1.0));
//	else auto coinToss = bind(uniform01);
//	temp->activationThreshold[c] = coinToss();
	if(p <= 0.33*popSize) temp->activationThreshold[c] = uniformInRange(0.9, 1.0);
	else if(p > 0.33*popSize && p <= 0.66*popSize) temp->activationThreshold[c] = uniformInRange(0.3, 1.0);
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
	int numClustersToActivate = uniformInRange(minNumClusters, maxNumClusters); //instead of just activating 2 random centroids, activating based on a coin toss
	//while (num < minNumClusters) {
	while (num < numClustersToActivate) {
	  int i = uniformInRange(0, maxNumClusters - 1);
	  if (!temp->active[i]) {
	    temp->activationThreshold[i] = uniformInRange(activationThreshold, 1.0);
	    temp->active[i] = true;
	    //	    temp->numActiveClusters++;
	    num++;
	  }
	}
	temp->numActiveClusters = num;
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

}//end setup

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
       
	sum += (org->clusters[c]->at(i)).distance;
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
	int a = (org->clusters[c]->at(i1)).itemIndex;
	double maxIntraDist = numeric_limits<double>::min();
	//	double* tempItem = itemsArray[a]->items;
	//	for (int f = 0; f < numFeatures; f++) {
	//	  sumArr[f] += tempItem[f]; //to compute centroids
	//	}
	for (vector<int>::size_type i2 = 0; i2 != org->clusters[c]->size();
	     i2++) { //finding max distance between items in a cluster
	  if (i2 != i1) {
	    
	    int b =  org->clusters[c]->at(i2).itemIndex;
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
       
	int a = (org->clusters[c1]->at(i1)).itemIndex;
	//find distance between elements in separate clusters; go through all elements in other clusters
	for (int cl = 0; cl < c1; cl++) {
	  if (org->active[cl]) {
	    for (vector<int>::size_type ind = 0;
		 ind != org->clusters[cl]->size(); ind++) {
	      
	      int b = (org->clusters[cl]->at(ind)).itemIndex;
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
	      
	      int b = (org->clusters[cl]->at(ind)).itemIndex; 
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
	    
	    int b = (org->clusters[c1]->at(i2)).itemIndex;
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
	    org->clusCenter[c][f] = uniformInRange((min[f]+0.1*gap), (min[f] + (0.5 * gap)));//check if 1.5*min[f] < max[f]
	  } else if (org->clusCenter[c][f] <= max[f] && org->clusCenter[c][f] >= (max[f] - (0.1 * gap))) {
	    org->clusCenter[c][f] = uniformInRange((max[f] - (0.5 * gap)), (max[f]-(0.1*gap)));
	  }
	}//end for
	if (!org->clusters[c]->empty()) {
	  for (vector<int>::size_type j = 0; j != org->clusters[c]->size();j++) {
	    //go through all items in cluster
	    
	    temp_dist = 0;
	    min_dist = numeric_limits<double>::max();
	    int clusCtr = 0;
	    int minInd = -1;
	    while (clusCtr < maxNumClusters) {
	      if (clusCtr != c && org->active[clusCtr]) {
		temp_dist = (org->clusters[c]->at(j)).distance;
		if (temp_dist < min_dist) {
		  temp_dist = min_dist;
		  minInd = clusCtr;
		}
	      }
	      clusCtr++;
	    }
	    assert(minInd != -1);
	    Clustering obj((org->clusters[c]->at(j)).distance, (org->clusters[c]->at(j)).itemIndex);
	    org->clusters[minInd]->push_back(obj);
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
 * Input parameters : Pointer to chromosome
 * This method computes the clustering of a chromosome and then calls the fitness function inside
 * return  type : double; returns the fitness computed
 */
void DEMain::computeClustering(Individual* org) {
  int min_index = -1;
  double temp_dist = 0;
  double min;
  int itemIndex = 0;
  for(int c = 0; c <maxNumClusters; c++){
    org->clusters[c]->clear();  
  }
  while (itemIndex < numItems) { //form clusters
    min_index = -1;
    min = numeric_limits<double>::max();
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
    Clustering obj(min,itemIndex);
    org->clusters[min_index]->push_back(obj); //adding item to its closest cluster center
    itemIndex++;
  } // end while
}

/*
 * Input parameters : Pointer to offspring, min-max value for each feature
 * This method determines whether a parent is to be replaced by an offspring or not
 */
Individual* DEMain::replacement(Individual* org, double min[], double max[]) {
  acceptWorse = false;
  counterTwoCent = 0;
  counterThreeCent = 0;
  computeClustering(org);
  cleanIndividual(org, min, max);
  double defaultFit = calcFitness(org);
  org->rawFitness = defaultFit;
  if(isBestInd)
    printClusters(org);
  ClusterSort* objClus = new ClusterSort[maxNumClusters];
  for (int c = 0; c < maxNumClusters; c++) {
    if (org->active[c]) {
      objClus[c].size = org->clusters[c]->size();
    } else {
      objClus[c].size = 0;
    }
    objClus[c].clusIndex = c;
  }
  qsort(objClus, maxNumClusters, sizeof(ClusterSort), compareSize);//to get the largest cluster in the chromosome
  if (objClus[0].size > 0.4 * numItems) {
    //create copy of original individual
    Individual* orgDup = new Individual(*org);
    if(uniform01() > 0.5) {   //with a probability decide if in between centroid addition or inside cluster centroid addition
      centroidAddition(orgDup, objClus);
      computeClustering(orgDup);
      cleanIndividual(orgDup, min, max);
      //check if there still exists a huge cluster
      int maxIndexDup = -1;
      double maxSize = numeric_limits<double>::min();
      for (int c = 0; c < maxNumClusters; c++) {
	if (orgDup->active[c] && orgDup->clusters[c]->size() > 0.4*numItems) {
	  if(maxSize < orgDup->clusters[c]->size()){
	    maxSize = orgDup->clusters[c]->size();
	    maxIndexDup = c;
	  }
	}
      }
      if(maxIndexDup != -1){ 
	//add centroids inside the large cluster if a large cluster still exists
	addOneCentroid(orgDup, maxIndexDup); 
	if(isBestInd){
	  trackFile << "Centroid added in between and inside" << endl;
	  if(counterTwoCent > 0) trackFile << "Second centroid added" << endl;
	  if(counterThreeCent > 0) trackFile << "Third centroid added" << endl;
	}
	computeClustering(orgDup);
	cleanIndividual(orgDup, min, max);
	//	cout << "outside addCent" << endl;
      }
    }
    else{
     
      addOneCentroid(orgDup, objClus[0].clusIndex);
      if(isBestInd){
	trackFile << "Centroid added inside alone" << endl;
	if(counterTwoCent > 0) trackFile << "Second centroid added" << endl;
	if(counterThreeCent > 0) trackFile << "Third centroid added" << endl;
      }
      computeClustering(orgDup);
      cleanIndividual(orgDup, min, max);
      //      cout << "only addcent " << endl;
    }
    double newFitness;   
    if(orgDup->numActiveClusters > minNumClusters){
      newFitness = calcDBIndex(orgDup);
      orgDup->rawFitness = 100 *(1/newFitness);
    }
    else{
      orgDup->rawFitness = -1;
      //newFitness = -1;
    }
    if(isBestInd){
      trackFile << "////////////////////////////////////////////////////" << endl;
      printClusters(orgDup);
    }
    if(uniform01() < probWorse) acceptWorse = true;
    if (orgDup->rawFitness > org->rawFitness) {
      delete org;
      return orgDup;
    }
    else {
      if(!acceptWorse)
	{     
	  delete orgDup;
	  return org;
	}
      else{
	delete org;
	return orgDup;
      }
    }
  }
  else{
    return org;
  }
}

/*
 * This method takes in the offspring that has the large cluster we want to split and the index of the centroid corresponding to the largest cluster
 * Output: offspring has upto 3 new centroids at the end
 */
void DEMain::addOneCentroid(Individual* org, int index){
  int size = org->clusters[index]->size();
  sort(org->clusters[index]->begin(), org->clusters[index]->end(), compareDist);
  int itemIndex;
  double maxDist = (org->clusters[index]->at(size - 1)).distance;
  int startRange = callBinarySearch(org,0.7*maxDist, size, index);
  int endRange = size - 1;
  int newCent = -1;
  DistItemCluster* obj2 = new DistItemCluster[size];
  double calcDist;
  //cout << startRange << " startRange" << endl;
  if(startRange != endRange && startRange != 0) {
    int c = uniformInRange(startRange, endRange);//find a random item in the 70-100% range
    //cout << "random c" << c << endl;
    itemIndex = (org->clusters[index]->at(c)).itemIndex;
    maxDist = numeric_limits<double>::min();
    double distMax = numeric_limits<double>::min(); 
    //cout << "c1 is " << itemIndex << endl;
    for (vector<int>::size_type i = 0; i != size; i++) {//form array of distances from the new centroid
      
      int a = (org->clusters[index]->at(i)).itemIndex;
      if(itemIndex != a)
	calcDist = findDist(itemIndex, a);
      else calcDist = 0;//distance of item from itself is 0
      //cout << "calcDist " << calcDist << endl;
      obj2[i].distance = calcDist;
      obj2[i].itemIndex = a;
      if(maxDist < calcDist) maxDist = calcDist;
    }
    //cout << "maxDist from c1 " << maxDist << endl;
    for(int i = startRange; i <= endRange; i++){
      int orgClusInd = org->clusters[index]->at(i).itemIndex;
      if(obj2[i].itemIndex == orgClusInd){
	if(obj2[i].distance >= 0.7*maxDist && obj2[i].distance <= maxDist){//check if there exists an item in the intersection
	  if(distMax < obj2[i].distance){
	    distMax = obj2[i].distance;
	    newCent = obj2[i].itemIndex;
	  }
	}
      } 
    }//end for
  }
  else{
    itemIndex = ( org->clusters[index]->at(endRange)).itemIndex;
  }
  //  cout << "centroid 1 " << itemIndex << " centroid 2 " << newCent << endl;
  //add this as centroid in the chromosome 
  double min = numeric_limits<double>::max();
  int minIndex = -1;
  //cout << "Checkpoint" << endl;
  if(org->numActiveClusters != maxNumClusters){
    for (int c = 0; c < maxNumClusters; c++) {
      if (org->active[c] == false) {
	//find index of centroid that's closest to new centroid
	double tempDist = dist(org->clusCenter[c], itemsArray[itemIndex]->items);
	if (tempDist < min) {
	  min = tempDist;
	  minIndex = c;
	}
      }
    }//end for
    //replace centroid at minIndex with the new centroid. recompute clustering, clean up and calc. fitness
    assert(minIndex != -1);
  }
  org->active[minIndex] = true;
  org->numActiveClusters++;
  org->activationThreshold[minIndex] = uniformInRange(activationThreshold, 1.0);
  for (int f = 0; f < numFeatures; f++) {
    org->clusCenter[minIndex][f] = itemsArray[itemIndex]->items[f];
  }
  //form arrays to pass to addTwoCentroid
  //cout << "checkpoint" <<endl;
  if(newCent != -1){
    // cout << "second centroid added" << endl;
    counterTwoCent++;
    min = numeric_limits<double>::max();
    minIndex = -1;
    if(org->numActiveClusters != maxNumClusters){
      for (int c = 0; c < maxNumClusters; c++) {
	if (org->active[c] == false) {
	  //find index of centroid that's closest to new centroid
	  double tempDist = dist(org->clusCenter[c], itemsArray[newCent]->items);
	  if (tempDist < min) {
	    min = tempDist;
	    minIndex = c;
	  }
	}
      }//end for
      //replace centroid at minIndex with the new centroid. recompute clustering, clean up and calc. fitness
      assert(minIndex != -1);
    }
    org->active[minIndex] = true;
    org->numActiveClusters++;
    org->activationThreshold[minIndex] = uniformInRange(activationThreshold, 1.0);
    for (int f = 0; f < numFeatures; f++) {
      org->clusCenter[minIndex][f] = itemsArray[newCent]->items[f];
    }
    //deactivating the original centroid
    org->active[index] = false;
    org->numActiveClusters--;
    org->activationThreshold[index] = uniformInRange(0.0,activationThreshold);
    //  addThreeCentroid(org, itemIndex, newCent, index, obj2, maxDist);
    //cout << " checkpoint after c3 " << endl;
  }//end if newCent
}

/*
 * Input : offspring, index of c1, index of c2, index of original centroid, sorted array of items distant from c1, maxDist of an item from c1
 * Output : if a new centroid added, old centroid is replaced with the new one or old centroid stays as it is
 */
void DEMain::addThreeCentroid(Individual* org, int ind1, int ind2, int origInd, DistItemCluster* obj2, double maxDist){
  qsort(obj2,org->clusters[origInd]->size(), sizeof(DistItemCluster), compare);
  //select 70-100% 
  int startRange = callBinarySearch(obj2,0.7*maxDist, org->clusters[origInd]->size());
  int endRange =  org->clusters[origInd]->size()- 1;
  double calcDist;
  if(startRange != endRange && startRange != 0) {
    DistItemCluster* obj3 = new DistItemCluster[org->clusters[origInd]->size()];
    double maxDist1 = numeric_limits<double>::min();
    for (vector<int>::size_type i = 0; i != org->clusters[origInd]->size(); i++) {
      int a = (org->clusters[origInd]->at(i)).itemIndex;
      if(ind2 != a)
	calcDist = findDist(ind2, a);
      else calcDist = 0.0;
      obj3[i].distance = calcDist;
      obj3[i].itemIndex = a;
      if(maxDist1 < calcDist) maxDist1 = calcDist;
    }
    int newCent = -1;
    double distMax = numeric_limits<double>::min();
    for(int i = startRange; i <= endRange; i++){
      if(obj3[i].itemIndex == obj2[i].itemIndex){
	if(obj3[i].distance >= 0.1*maxDist1 && obj3[i].distance <= maxDist1){
	  if(distMax < obj3[i].distance){
	    distMax = obj3[i].distance;
	    newCent = obj3[i].itemIndex;
	  }
	}
    }
    }//end for
    //replace old centroid with third centroid
    if(newCent != -1){
      counterThreeCent++;
      for (int f = 0; f < numFeatures; f++) {
	org->clusCenter[origInd][f] = itemsArray[newCent]->items[f];
      }
    }
    else{
      //move original centroid away with some jittering
      for (int f = 0; f < numFeatures; f++) {
	org->clusCenter[origInd][f] *= gFScale;
      }
    }
    //delete [] obj3;
  }

  //  delete [] obj2;
}

/*
 * This method computes and adds centroid to the chromosome passed
 */
void DEMain::centroidInsertion(Individual* org, int c1, int c2) {
  bool isOneCentroid = false;
  int closerCentInd = 0; // index of the new centroid(closer to larger cluster) after it replaces the one its closest to in the actual chromosome
  double maxDist;
  double avgDistCal1 = 0, avgDistCal2 = 0;
  double sdSum1 = 0, sdSum2 = 0;
  maxDist = dist(org->clusCenter[c1], org->clusCenter[c2]);
  //computing the average distance and s.d for the largest cluster center and the selected subset
  for (vector<int>::size_type j = 0; j != org->clusters[c1]->size(); j++) {

    avgDistCal1 += (org->clusters[c1]->at(j)).distance;
  }
  avgDistCal1 /= org->clusters[c1]->size();
  for (vector<int>::size_type j = 0; j != org->clusters[c1]->size(); j++) {

    sdSum1 += pow((org->clusters[c1]->at(j).distance - avgDistCal1), 2.0);
  }
  sdSum1 /= org->clusters[c1]->size();
  sdSum1 = sqrt(sdSum1);
 for (vector<int>::size_type j = 0; j != org->clusters[c2]->size(); j++) {

    avgDistCal2 += org->clusters[c2]->at(j).distance;
  }
  avgDistCal2 /= org->clusters[c2]->size();
  for (vector<int>::size_type j = 0; j != org->clusters[c2]->size(); j++) {

   sdSum2 += pow((org->clusters[c2]->at(j).distance - avgDistCal2), 2.0);
  }
  sdSum2 /= org->clusters[c2]->size();
  sdSum2 = sqrt(sdSum2);
  //cout << avgDistCal1 << " " << avgDistCal2 << " " << sdSum1 << " " << sdSum2 << endl;
  double newScale, newScale1, newScale2;
  double* newCentroid = new double[numFeatures];
  double* newCentroid1 = new double[numFeatures];
  double* newCentroid2 = new double[numFeatures];
  //need to check if the centroid is too close or far from the larger cluster
  if ((0.5 * maxDist) < (avgDistCal1 + sdSum1) && (0.5 * maxDist) < (avgDistCal2 + sdSum2)) {
    //add one centroid with some jitter when both clusters are overlapping in between
    for (int f = 0; f < numFeatures; f++) {
      newScale = 0.5 + (uniformInRange(0.0, 0.1) - 0.05);
      newCentroid[f] = org->clusCenter[c1][f] + newScale*(org->clusCenter[c2][f]- org->clusCenter[c1][f]);
    }
    isOneCentroid = true;
  } else if((0.5 * maxDist) > (avgDistCal1 + sdSum1) && (0.5 * maxDist) > (avgDistCal2 + sdSum2)){
    //introduce two new centroids if both clusters are very far away such that the introduced centroids are closer to periphery of the clusters
    newScale1 = (avgDistCal1 + 0.5*sdSum1)/maxDist;
    newScale2 = (avgDistCal2 + 0.5*sdSum2)/maxDist;
    for (int f = 0; f < numFeatures; f++) {
      newScale1 *= 1+(uniformInRange(0.0, 0.1) - 0.05);
      newScale2 *= 1+(uniformInRange(0.0, 0.1) - 0.05);
      newCentroid1[f] = org->clusCenter[c1][f] + newScale1*(org->clusCenter[c2][f]- org->clusCenter[c1][f]);
      newCentroid2[f] = org->clusCenter[c2][f] + newScale2*(org->clusCenter[c1][f]- org->clusCenter[c2][f]);
    }   
    isOneCentroid = false;
  }
  else if((0.5 * maxDist) > (avgDistCal1 + sdSum1) && (0.5 * maxDist) < (avgDistCal2 + sdSum2)){
    //inside c2 and outside c1
    //introduce two new centroids;
    newScale1 = (avgDistCal1 + 0.5*sdSum1)/maxDist;
    newScale2 = ((avgDistCal2 + 0.5*sdSum2) > (0.5*maxDist)) ? (avgDistCal2 + 0.5*sdSum2)/maxDist: 0.5;
    closerCentInd = 2;
    for (int f = 0; f < numFeatures; f++) {
      newScale1 *= 1+(uniformInRange(0.0, 0.1) - 0.05);
      newScale2 *= 1+(uniformInRange(0.0, 0.1) - 0.05);
      newCentroid1[f] = org->clusCenter[c1][f] + newScale1*(org->clusCenter[c2][f]- org->clusCenter[c1][f]);
      newCentroid2[f] = org->clusCenter[c2][f] + newScale2*(org->clusCenter[c1][f]- org->clusCenter[c2][f]);
    }   
    isOneCentroid = false;
  }
  else if((0.5 * maxDist) < (avgDistCal1 + sdSum1) && (0.5 * maxDist) > (avgDistCal2 + sdSum2)){
    //inside c1 and outside c2
    //introduce two new centroids
    newScale2 = (avgDistCal2 + 0.5*sdSum2)/maxDist;
    newScale1 = ((avgDistCal1 + 0.5*sdSum1) > (0.5*maxDist)) ? (avgDistCal1 + 0.5*sdSum1)/maxDist: (0.5);
    closerCentInd = 1;
    for (int f = 0; f < numFeatures; f++) {
      newScale1 *= 1+(uniformInRange(0.0, 0.1) - 0.05);
      newScale2 *= 1+(uniformInRange(0.0, 0.1) - 0.05);
      newCentroid1[f] = org->clusCenter[c1][f] + newScale1*(org->clusCenter[c2][f]- org->clusCenter[c1][f]);
      newCentroid2[f] = org->clusCenter[c2][f] + newScale2*(org->clusCenter[c1][f]- org->clusCenter[c2][f]);
    }
    isOneCentroid = false;
  }
  double min = numeric_limits<double>::max();
  int minIndex = -1;
  if(isOneCentroid){
    if(org->numActiveClusters != maxNumClusters){
      for (int c = 0; c < maxNumClusters; c++) {
	if (org->active[c] == false) {
	  //find index of centroid that's closest to new centroid
	  double tempDist = dist(org->clusCenter[c], newCentroid);
	  if (tempDist < min) {
	  min = tempDist;
	  minIndex = c;
	  }
	}
      }//end for
    //replace centroid at minIndex with the new centroid. recompute clustering, clean up and calc. fitness
      assert(minIndex != -1);
    }
    org->active[minIndex] = true;
    org->numActiveClusters++;
    org->activationThreshold[minIndex] = uniformInRange(activationThreshold, 1.0);
    for (int f = 0; f < numFeatures; f++) {
      org->clusCenter[minIndex][f] = newCentroid[f];
    }
  }
  else{
    //find min for both centroids
    int minIndex1 = -1, minIndex2 = -1;
    double min1 = numeric_limits<double>::max();
    if(org->numActiveClusters != maxNumClusters) {
      for (int c = 0; c < maxNumClusters; c++) {
	if (org->active[c] == false) {
	  //find index of centroid that's closest to new centroid
	  double tempDist = dist(org->clusCenter[c], newCentroid1);
	if (tempDist < min1) {
	  min1 = tempDist;
	  minIndex1 = c;
	}
	}
      }//end for
      //replace centroid at minIndex with the new centroid. recompute clustering, clean up and calc. fitness
      assert(minIndex1 != -1);
    }
    org->active[minIndex1] = true;
    org->numActiveClusters++;
    org->activationThreshold[minIndex1] = uniformInRange(activationThreshold, 1.0);
    double min2 = numeric_limits<double>::max();
    if(org->numActiveClusters != maxNumClusters) {
      //cout << "-------------" << endl;
      for (int c = 0; c < maxNumClusters; c++) {
	if (org->active[c] == false) {
	  //find index of centroid that's closest to new centroid
	  //cout << "size of c1 " << org->clusters[c1]->size() << " and size of c2 " << org->clusters[c2]->size() << endl; 
	  double tempDist = dist(org->clusCenter[c], newCentroid2);
	  //cout << "temp dist " << tempDist << "min " << min2 << endl;
	  if (tempDist < min2) {
	    min2 = tempDist;
	    minIndex2 = c;
	    //cout << minIndex2 << endl;
	  }
	}
      }//end for
      //replace centroid at minIndex with the new centroid. recompute clustering, clean up and calc. fitness
      assert(minIndex2 != -1);
    }
    org->active[minIndex2] = true;
    org->numActiveClusters++;
    org->activationThreshold[minIndex2] = uniformInRange(activationThreshold, 1.0);
    assert(minIndex1 != minIndex2);
    for (int f = 0; f < numFeatures; f++) {
      org->clusCenter[minIndex1][f] = newCentroid1[f];
    }
    for (int f = 0; f < numFeatures; f++) {
      org->clusCenter[minIndex2][f] = newCentroid2[f];
    }
  }
  //if(closerCentInd == 1) indCentAdded = minIndex1;
  //else if(closerCentInd == 2) indCentAdded = minIndex2;
  delete[] newCentroid;
  delete[] newCentroid1;
  delete[] newCentroid2;
}

/*
 * Input : offspring and object containing clusters in sorted order
 * This method finds the subset of other cluster centers that are large enough to be considered along with the largest cluster to add new centroids in between
 */
void DEMain::centroidAddition(Individual* org, ClusterSort* objClus){
  bool* clusSubSet = new bool[maxNumClusters];
  int maxSize = objClus[0].size;
  clusSubSet[0] = true; 
  for(int c = 1; c <= maxNumClusters/5; c++){
    //if(objClus[c].size >= (0.5*maxSize) && uniform01() <= 0.5)
    if(uniform01() <= 0.5 && org->active[objClus[c].clusIndex] && org->clusters[objClus[c].clusIndex]->size() > 0.25*maxSize)
      clusSubSet[c] = true;
    else
      clusSubSet[c] = false;
  }
  for(int c = 1; c <= maxNumClusters/5; c++){
    if(clusSubSet[c]){
      centroidInsertion(org, objClus[0].clusIndex, objClus[c].clusIndex);
      //cout << "Inserted centroids" << endl;
    }
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
  //change to: random displacement with wraparound
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
  printTrackArr = false;
  double featureFScale,d;
  double thresFScale = 0.5;
  int s1, s2, s3;
  bool isExploration = true;
  if(genNum >= 0.9*numGenerations){
    isExploration = false;
    probWorse = 0.1;
  }
  double crossoverProb = crossoverProbability * ((numGenerations - genNum) / numGenerations);//based on formula in paper
  int epoch = numGenerations/6;
  if(genNum == lastUpdatedGen+epoch){
    scale -= 0.1;
    lastUpdatedGen = genNum;
    //     //cout << scale << " and gen " << genNum << endl;
  }
  if(isExploration) {
    d = 0.5;
    featureFScale = scale + (d*(uniform01() - 0.5)); //based on formula in paper
  }
  else { 
    d = 0.001;
    crossoverProb = uniformInRange(0.8,0.98);
  }//need to discuss

  selectSamples(orgIndex, s1, s2, s3);//to get unique individuals to perform DE crossover on
  //cout << "selectSamples finished" << endl;
  Individual* child = new Individual(maxNumClusters, numFeatures);
  int counterActiveClus = 0;//keeps track of active cluster centers
  for (int c = 0; c < maxNumClusters; c++) {
    bool change = uniform01() < crossoverProb? true : false;
    trackChange[c] = change;
    if(uniform01() < crossoverProb){
      child->activationThreshold[c] =  popObject->org[s1]->activationThreshold[c] +
	thresFScale*(popObject->org[s2]->activationThreshold[c] -  popObject->org[s3]->activationThreshold[c]);    
    }
    else{
      child->activationThreshold[c] =  popObject->org[orgIndex]->activationThreshold[c];               
    }
    for (int f = 0; f < numFeatures; f++) {//binomial crossover computation for centroids
      assert(popObject->org[orgIndex] != NULL);
      //jitter
      if(!isExploration) { 
	featureFScale = scale + (d*(uniform01() - 0.5));
      }
      if(change ==1) {
	child->clusCenter[c][f] = popObject->org[s1]->clusCenter[c][f]
	  + featureFScale*(popObject->org[s2]->clusCenter[c][f]- popObject->org[s3]->clusCenter[c][f]);
	if (child->clusCenter[c][f] < min[f]) {
	  child->clusCenter[c][f] =  popObject->org[s1]->clusCenter[c][f]
	    + abs( featureFScale*(popObject->org[s2]->clusCenter[c][f]- popObject->org[s3]->clusCenter[c][f]) );
	}
	else if (child->clusCenter[c][f] > max[f]){//if a feature for cluster center computed goes out of its bounds, assign it a random value
	  child->clusCenter[c][f] =  popObject->org[s1]->clusCenter[c][f]
	    - abs( featureFScale*(popObject->org[s2]->clusCenter[c][f]- popObject->org[s3]->clusCenter[c][f]));
	}
      } else {	
	child->clusCenter[c][f] = popObject->org[orgIndex]->clusCenter[c][f];
      }
      if ((child->clusCenter[c][f] < min[f]) || (child->clusCenter[c][f] > max[f])){
	child->clusCenter[c][f] = uniformInRange(min[f], max[f]);
      } 
     
    }//for i

    //reporting
   

    if (child->activationThreshold[c] > 1 || child->activationThreshold[c] < 0)//if a threshold value after crossover goes out of its bounds, recompute it
      child->activationThreshold[c] = uniform01();
    if (child->activationThreshold[c] > activationThreshold) {
      child->active[c] = true;
      counterActiveClus++;
    } else
      child->active[c] = false;
  }//for j
  assert(counterActiveClus <= maxNumClusters);
  gCProb = crossoverProb;
  gFScale = featureFScale;
  child->numActiveClusters = counterActiveClus;
  if(child->numActiveClusters < minNumClusters){//if #of active cluster centers below kmin, forcefully make more clusters active till kmin satisfied
    int num = child->numActiveClusters;
    // int randNumClusters = uniformInRange(1, maxNumClusters/2);
    // while(num < randNumClusters) {
    while (num < minNumClusters) {//we don't want to activate too many centroids because we want the threshold to kick in
      int i = uniformInRange(0, maxNumClusters - 1);
      if(!child->active[i]){
	child->activationThreshold[i] = uniformInRange(activationThreshold, 1.0);
	child->active[i] = true;
	child->numActiveClusters++;
	num++;
      }
    }
  }

 if(*child == *(popObject->org[orgIndex])){
      //print s1, s2, s3
 
      //print s1, s2, s3
      cout << "Generation " << genNum << endl;
      cout << "s1 : ";
      for(int c =0; c < maxNumClusters; c++){
	for(int f = 0; f < numFeatures; f++)
	  cout << popObject->org[s1]->clusCenter[c][f] << " ";
      }
      cout << endl;
      cout << "s2 : ";
      for(int c =0; c < maxNumClusters; c++){
	for(int f = 0; f < numFeatures; f++)
	  cout << popObject->org[s2]->clusCenter[c][f] << " ";
      }
      cout << endl;
      cout << "s3 : ";
      for(int c =0; c < maxNumClusters; c++){
	for(int f = 0; f < numFeatures; f++)
	  cout << popObject->org[s3]->clusCenter[c][f] << " ";
      }
 
 }
  return child;

}

/*
 * Reinitialize variables in between cycle
 */
void DEMain::restart(){
  scale = origScale;
  isBestInd = false;
  lastUpdatedGen = 0;
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
	temp->activationThreshold[i] = uniformInRange(activationThreshold, 1.0);
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
  if(cycle > 0 && cycle < 0.4*numRepeat){
    //currently replacing a small portion of population irrespective of avg/best
    //replace 20% of population randomly
    for(int repl = 0; repl < 0.2*popSize; repl++){
      int indexRepl = permuteArray[repl];
      if(indexRepl != popObject->bestOrgIndex){
	//modify individual
	initializePopCycle(popObject->org[indexRepl], min, max);
      }      
    }
  }
  else if(cycle >= 0.4*numRepeat && cycle < 0.9*numRepeat){
    //currently replacing a small portion of population irrespective of avg/best
    //replace 10% of population randomly
    for(int repl = 0; repl < 0.1*popSize; repl++){
      int indexRepl = permuteArray[repl];
      if(indexRepl != popObject->bestOrgIndex){
	//modify individual
	initializePopCycle(popObject->org[indexRepl], min, max);
      }
      
    }
  }
  else if(cycle >= 0.9*numRepeat && cycle < numRepeat){
    //replace fewer elements from population
    //no changes in the last cycle
    for(int repl = 0; repl < 0.05*popSize; repl++){
      int indexRepl = permuteArray[repl];
      if(indexRepl != popObject->bestOrgIndex){
	//modify individual
	initializePopCycle(popObject->org[indexRepl], min, max);
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
  ofstream trackBest;
  try {
    ostringstream strs, nf;
    strs << numRepeat;
    nf << numFeatures;
    string lr = strs.str();
    string fn = nf.str();
    trackFile.open("trackCentAddition_" + lr + "_" + fn + ".txt");
    trackOff.open("record_wine3.txt");
    trackBest.open("BestChromosme_wine3.txt");
    trackBest << "Gen 0 : Index : " << popObject->bestOrgIndex << " with fitness " << popObject->org[popObject->bestOrgIndex]->rawFitness << endl; 
    trackOff << setiosflags(ios::left);
    trackOff << setw(5) << "Gen" << setw(5) << "Avg DB" << "|" << setw(10) << "Max DB"  << "|" << setw(10) << "Min DB" << endl;
    for(int cycle = 0; cycle < numRepeat; cycle++){
      while (g < numGenerations) {    //till max generation reached
	permuteBaseArray();
	//cout << "permute finished" << endl;
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
	  //cout << "crossover finished" << endl;
	  if(p == 1) {
	    trackFile << "-------------------------------------------------" << endl;
	    trackFile << "For individual at index " << p << " at generation " << g << " and cycle " << cycle << endl;
	    trackFile << "Crossover prob " << gCProb << " and F scale " << gFScale << endl;
	    trackFile << "Change array : ";
	    if(printTrackArr) {
	      for(int i = 0; i < maxNumClusters; i++){
		trackFile << trackChange[i] << " ";
	      }
	    }
	    trackFile << endl;
	    isBestInd = true;
	  }
	  else{isBestInd = false;}
	  offspring = replacement(child, min, max);
	  //cout << "replacement finished" << endl;
	  double avg;
	  if (popObject->org[p]->rawFitness <= offspring->rawFitness) { //if offspring better than parent, replace the parent
	    isReplaceOrg[p] = true;
	    if(isBestInd) {trackFile << "Offspring added and parent's fitness " << popObject->org[p]->rawFitness << endl;}
	    newpop->org[p] = offspring;
	    avg = 100*(1/offspring->rawFitness);
	    avgDB += avg;
	    if (minDB > avg) {//keeping track of the minimum and maximum DB index encountered for all individuals in all generations
	      minDB = avg;
	    }
	    if (maxDB < avg) {
	      maxDB = avg;
	    }	  
	  
	}
	else { //otherwise delete the offspring
	  if(!acceptWorse) {
	    isReplaceOrg[p] = false;
	    if(isBestInd) {trackFile << "Offspring discarded and parent's fitness " << popObject->org[p]->rawFitness << endl;}
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
	  else{
	    isReplaceOrg[p] = true;
	    if(isBestInd) {trackFile << "Offspring added and parent's fitness " << popObject->org[p]->rawFitness << endl;}
	    newpop->org[p] = offspring;
	    avg = 100*(1/offspring->rawFitness);
	    avgDB += avg;
	    if (minDB > avg) {//keeping track of the minimum and maximum DB index encountered for all individuals in all generations
	      minDB = avg;
	    }
	    if (maxDB < avg) {
	      maxDB = avg;
	    }
	  }
	}//end else
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
      trackBest << "Gen " << g+1 << " : Index : " << popObject->bestOrgIndex << " with fitness " << popObject->org[bestInd]->rawFitness << endl; 
	g++; //increment generation number
    }//end while
    //cause perturbation in population
    g=0;
    perturbPop(cycle,min,max);
  }
  //cout << endl;
  ////cout << "Out of bounds " << countOOB << endl;
  //cout << "Stopped at generation " << g << endl;
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
  trackFile.close();
  trackBest.close();
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
  centroids << "Centroids for best individual obtained are : " << endl;
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
	int itemIndex = org->clusters[c]->at(i).itemIndex;
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

	int itemInd =  org->clusters[c]->at(i).itemIndex;
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
  //cout << "Result saved in file.";

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
      sort(org1->clusters[c1]->begin(), org1->clusters[c1]->end(), compareItems);
      size = org1->clusters[c1]->size();
      p1 = size / itemSize;
      //	//cout << "P1 is " << p1 << " ";
      H1 += p1 * (log2(p1));//H_u = \sum(P(i)log(P(i))) to calc entropy in clustering for 1
      if(!isFinal) {   
	isSorted = new bool[maxNumClusters]();
	Individual* org2 = popObject->org[popInd2];
	for (int c2 = 0; c2 < maxNumClusters; c2++) {
	  if (org2->active[c2]) {
	    vector<int> clus3;
	    if(!isSorted[c2]) {
	      sort(org2->clusters[c2]->begin(), org2->clusters[c2]->end(), compareItems);
	      isSorted[c2] = true;
	    }
	    int sizeClus =org1->clusters[c1]->size();
	    vector<int>* myClus = new vector<int>[sizeClus];
	    for (vector<int>::size_type i = 0; i != sizeClus; i++) {
	      myClus->push_back((*(org1->clusters[c1]))[i].itemIndex);
	    }
	    set_intersection(myClus->begin(), myClus->end(),
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
	    int itemInd = org1->clusters[c1]->at(it1).itemIndex;
	    if(itemInd == origClustering[c2]->at(it2)){
	      ctr++;
	      it1++; it2++;
	    }
	    else if(itemInd > origClustering[c2]->at(it2)) it2++;
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
  //cout << "MI " << mi << " " << " H1 " << H1 << " h2 " << H2 << endl;
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
	int a = popObject->org[popInd1]->clusters[c]->at(i).itemIndex;
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
