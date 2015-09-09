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
vector<int>** origClustering;
ofstream trackFile;
ofstream trackOff;

int countOOB = 0;

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


DEMain::DEMain(int dim, int** tracebackArr, Item** items, int itemSize, int validityIndex, Parameters param) {
	// TODO Auto-generated constructor stub
	popObject = new Population(param.maxNumClusters, dim, param.popScaleFactor);
	strategy = stRand1Bin; //DE algo used
	numGenerations = param.gen; //number of generations
	crossoverProbability = param.CrMax - param.CrMin; //crossover probability
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

	//scratch space
	trackerArray = tracebackArr; //holds cluster indexes corresponding to every item and pop member(calcFitness and report)
	clusters = new vector<int>*[maxNumClusters]; //holds actual clusters for every pop members(calcFitness)
	distItem = new double*[numItems]; //jagged array storing dist between all items(used for calcCSIndex and calcPBIndex)
	for (int i = 0; i < numItems - 1; i++) {
		distItem[i] = new double[numItems - 1 - i];
	}
	offspringArray = new int[numItems]; //column that replaces worse parent's column from tracker(calcFitness)
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
	scalingArr = new double[numFeatures];
	for (int count = 0; count < maxNumClusters; count++) {
		clusters[count] = new vector<int>;
		newClustCenters[count] = new double[dim];
	}
}

DEMain::~DEMain() {
	// TODO Auto-generated destructor stub
	 delete popObject;
	  for (int i = 0; i < numItems; i++){
	    delete [] trackerArray[i];
	    delete [] itemsArray[i];
	    delete [] distItem[i];
	  }
	  for(int i = 0; i < maxNumClusters; ++i) {
	    delete clusters[i];
	    delete newClustCenters[i];
	  }
	  delete [] trackerArray;
	  delete [] itemsArray;
	  delete [] distItem;
	  delete [] clusters;
	  delete [] offspringArray;
	  delete [] nearestNeighborTriples;
	  delete [] ItemUsed;
	  delete [] ClusFull;
	  delete [] avgArr;
	  delete [] scalingArr;
	  delete [] sumArr;
	  delete [] ItemCounter;
	  delete [] isReplaceOrg;
	  delete [] newClustCenters;
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
  for(int f = 1; f <numFeatures; f++){
   scalingArr[f] = maxFeatVal/max[f] * uniformInRange(0.7,1.0); //computing the scale factor for each attribute
  }
  for (int i = 0; i < numItems - 1; i++) {
    for (int j = i + 1; j < numItems; j++) {
      distItem[i][j - i - 1] = dist(itemsArray[i]->items, itemsArray[j]->items);     
    }  //end j for
  }  // end i for

}

void DEMain::printClusters(int popIndex) {
  int clusIndex;
  for (int c = 0; c < maxNumClusters; ++c) {
    clusters[c]->clear();
  }
  Individual* org = popObject->org[popIndex];
  for (int i = 0; i < numItems; i++) {
    clusIndex = trackerArray[i][popIndex];
    if (org->active[clusIndex]) {
      clusters[clusIndex]->push_back(i);
    }
  }//forming clusters
  double* arr;
  int activeCount = 0;
  for (int c = 0; c < maxNumClusters; c++) {
    if (org->active[c]) {			//prints out the best output
      activeCount++;
      trackOff << "-------------------------------------------" << endl;
      trackOff << "Elements of cluster " << c << " :" << endl;
      for (vector<int>::size_type i = 0; i != clusters[c]->size(); i++) {
	int itemIndex = clusters[c]->at(i);
	arr = itemsArray[itemIndex]->items;
	for (int f = 0; f < numFeatures; f++) {
	  trackOff << arr[f] << " ";
	}
	trackOff << itemsArray[itemIndex]->typeClass;
	trackOff << endl;
      }
    }
  }
  trackOff << "Total number of clusters obtained for index " << popIndex << " is " <<  activeCount << endl;
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
				}//end for f
			}//end for c
			assert(numActiveClus <= maxNumClusters); //check to determine if active clusters are more than kmax
			temp->numActiveClusters = numActiveClus; //assigning this as active cluster count in the individual object
			//code to check kmin = 2
			if (temp->numActiveClusters < minNumClusters) { //if #of active clusters is less than min #of clusters possible, then activate a cluster centroid forcefully
				int num = temp->numActiveClusters;
				int activateClusNum = uniformInRange(1, maxNumClusters/2); //instead of just activating 2 random centroids, activating based on a coin toss
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
			fitn = calcFitness(temp, p, true, -1, min, max); //calculate fitness of individual generated for initial population
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
double DEMain::dist(double* x, double* y) {
  double Sum = 0.0;
  for (int f = 0; f < numFeatures; f++) {
      Sum = Sum + pow(scalingArr[f] * (x[f] - y[f]), 2.0);
  }
  return sqrt(Sum);
}

/*
 * Input parameters : A pointer to a chromosome, struct array that holds distance corresponding
 * to item and cluster center, size of the struct array and index of population's chromosome
 * This method reshuffles items from largest cluster into different empty active cluster centers of an individual based on a threshold value
 * return type : void
 */
void DEMain::reshuffleValid(Individual* org, int orgIndex, bool isInitial, bool isExploitation, int cIndex) {
  int fixSize = 0.1*(clusters[cIndex]->size()); //at least one-tenth of the size of largest cluster
  int maxSize = 1.1 * fixSize; //don't want gap between minsize and maxsize to be large
  int numWeakClus = 0;
  int numSmallClus = 0;
  bool* ItemUsed = new bool[numItems]();
  bool* ClusFull = new bool[maxNumClusters]();
  bool* smallClus = new bool[maxNumClusters]();
  bool smallClusFull = false;
  bool* isWithinAvg = new bool[maxNumClusters];
  bool* isSDCalc = new bool[maxNumClusters]();
  double* avgDistCal = new double[maxNumClusters];
  int ctrSmallFull = 0;
  bool isSmallClusFull = false;
  int itemIndex = 0;
  //  cout << "Min size for invalid cluster is " << fixSize << " & max size is " << maxSize << endl;
  for(int c= 0; c < maxNumClusters; c++){
    if(org->active[c]){
      if (clusters[c]->empty() || clusters[c]->size() == 1) {
	numWeakClus++; //find number of invalid clusters
	//	arrSize[c] = fixSize;
	smallClus[c]=true;
      }
    }//end active if
  }//end for

  int remSize = 0.67*(clusters[cIndex]->size()); //approx. items we can reshuffle
  int allowedClus = remSize/fixSize;

  if(numWeakClus > allowedClus) {
    int numRepeat = numWeakClus - allowedClus;
    int ctr = 0;
    while (ctr == numRepeat) {
      int i = uniformInRange(0, maxNumClusters - 1);
      if(org->active[i] && smallClus[i]){
	org->activationThreshold[i] = uniformInRange(0.0,activationThreshold);
	org->active[i] = false;
	org->numActiveClusters--;
	smallClus[i] = false;
	ctr++;
      }
    }//end while
    numWeakClus = allowedClus;
  }//if extreme case 

  for (int c = 0; c < maxNumClusters; ++c) {
    clusters[c]->clear();
    isWithinAvg[c] = true;
    avgDistCal[c] = 0.0;
  }
  double tempDist = 0.0;
  if(numWeakClus > 0) {
  //form heaps for smaller clusters and then form large clusters

    Heap** objHeap = new Heap*[maxNumClusters];
    for (int clusIndex = 0; clusIndex < maxNumClusters; clusIndex++) {
      if(org->active[clusIndex] && smallClus[clusIndex]){	
	objHeap[clusIndex] = new Heap(numItems);
	for(int i = 0; i < numItems; i++) {
	  tempDist = dist(org->clusCenter[clusIndex],itemsArray[i]->items);
	  objHeap[clusIndex]->Enqueue(tempDist, i); 
	}//end inner for
      }//end if
    }//end for
    //after heap is made, we need to form clusters
    for (int clusIndex = 0; clusIndex < maxNumClusters; clusIndex++) {
      if(smallClus[clusIndex]){	
	cout<<	objHeap[clusIndex]->getNumElements() << endl;
      }
    }
    HeapItem *temp;
    while(!isSmallClusFull) {
      int selIndex = -1;
      int itemInd;
      
      if(ctrSmallFull < numWeakClus) {
	double minDist = numeric_limits<double>::max();
	for (int c1 = 0; c1 < maxNumClusters; c1++) {
	  if(org->active[c1] && smallClus[c1]){
	    if(objHeap[c1]->elements[0].getKey() < minDist) {
	      minDist = objHeap[c1]->elements[0].getKey();
	      selIndex = c1;
	    }
	  }//if active
	}//c1 for
	assert(selIndex != -1);
	temp = objHeap[selIndex]->Dequeue();
	itemInd = temp->getData();
	if(clusters[selIndex]->size() < fixSize){
	  clusters[selIndex]->push_back(itemInd);
	  ItemUsed[itemInd] = true; 
	  if (isInitial) {
	    trackerArray[itemInd][orgIndex] = selIndex; //mainitaining tracking array
	  } else {
	    offspringArray[itemInd] = selIndex;
	  }
	}
	else{
	  if(isExploitation) {
	    if(clusters[selIndex]->size() < maxSize && isWithinAvg[selIndex]){
	      //calc avg dist for clus                                                                                                                                                                                        
	      double sdSum = 0;
	      if(!isSDCalc[selIndex]) {
		for (vector<int>::size_type j = 0; j != clusters[selIndex]->size(); j++){
		  int a = clusters[selIndex]->at(j);//find item index stored in clusters                                                                                                                                       
		  avgDistCal[selIndex] += dist(itemsArray[a]->items, org->clusCenter[selIndex]);
		}
		avgDistCal[selIndex] /= clusters[selIndex]->size();
		for (vector<int>::size_type j = 0; j != clusters[selIndex]->size(); j++){
		  int a = clusters[selIndex]->at(j);//find item index stored in clusters                                                                                                                                       
		  sdSum += pow((dist(itemsArray[a]->items, org->clusCenter[selIndex]) - avgDistCal[selIndex]),2.0);
		}
		isSDCalc[selIndex] =true;
		sdSum /= clusters[selIndex]->size();
		sdSum = sqrt(sdSum);
	      }
	      if(avgDistCal[selIndex] !=0 && objHeap[selIndex]->elements[0].getKey() <= avgDistCal[selIndex] + sdSum) {//determine thres 
		avgDistCal[selIndex] = (clusters[selIndex]->size()*avgDistCal[selIndex]) +  objHeap[selIndex]->elements[0].getKey();
		clusters[selIndex]->push_back(itemInd); //add item to cluster                                                            
		ItemUsed[itemInd]=true;
		//update avg dist                                                                                                                                                                                             
		avgDistCal[selIndex] /= clusters[selIndex]->size();
		if (isInitial) {
		  trackerArray[itemInd][orgIndex] = selIndex; //mainitaining tracking array
		} else {
		  offspringArray[itemInd] = selIndex;
		}
	      }
	      else{
		isWithinAvg[selIndex] = false;
	      }
	    }
	    else{
	      if(!ClusFull[selIndex]) {
		ClusFull[selIndex] = true;
		ctrSmallFull++;
	      }
	    }
	  }//end exploitation if
	  else{
	    if(!ClusFull[selIndex]) {
	      ClusFull[selIndex] = true;
	      ctrSmallFull++;
	    }
	  }
	}//end else
      }
      else{
	isSmallClusFull = true;
    }
    }//end while
  }
  //for testing
  cout << "Min size " << fixSize << endl;
  for(int c = 0; c < maxNumClusters; c++){
    if(org->active[c] && smallClus[c]){
      cout << clusters[c]->size() << " ";
    }
  }
  cout << endl;
  //form knn for remaining elements
  while (itemIndex < numItems) { //form clusters
    if(!ItemUsed[itemIndex]) {
      int min_index = -1;
      double min = numeric_limits<double>::max();
      for (int clusIndex = 0; clusIndex < maxNumClusters; clusIndex++) {
	if (org->active[clusIndex] && !smallClus[clusIndex]) {
	  tempDist = dist(org->clusCenter[clusIndex],itemsArray[itemIndex]->items);
	  if (tempDist < min) { //finding the closest cluster center to a particular item
	    min = tempDist;
	    min_index = clusIndex;
	  }
	}
      }
      assert(min_index != -1);
      clusters[min_index]->push_back(itemIndex); //adding item to its closest cluster center
      if (isInitial) {
	trackerArray[itemIndex][orgIndex] = min_index; //mainitaining tracking array
      } else {
	offspringArray[orgIndex] = min_index;
      }
    }//end if
    itemIndex++;
  }
  //if still any invalid clusters, deactivate them
  for(int c= 0; c < maxNumClusters; c++){
    if(org->active[c]){
      if (clusters[c]->size() < 2) {
	org->activationThreshold[c] = uniformInRange(0.0,activationThreshold);
	org->active[c] = false;
	org->numActiveClusters--;
      }
    }
  }
}

/*
 * Input parameters : A pointer to a chromosome, struct array that holds distance corresponding
 * to item and cluster center, size of the struct array and index of population's chromosome
 * This method reshuffles items from largest cluster into different empty active cluster centers of an individual based on a threshold value
 * return type : void
 *
void DEMain::reshuffleValid(Individual* org, int numTriplesArray, int orgIndex, bool isInitial, bool isExploitation, int cIndex) {
  int fixSize = 0.1*(clusters[cIndex]->size()); //at least one-tenth of the size of largest cluster
  int maxSize = 1.1 * fixSize; //don't want gap between minsize and maxsize to be large
  bool isSmallClusFull = false;
  int numWeakClus = 0;
  int numSmallClus = 0;
  bool* ItemUsed = new bool[numItems]();
  bool* smallClus = new bool[maxNumClusters]();
  //  cout << "Min size for invalid cluster is " << fixSize << " & max size is " << maxSize << endl;
  for(int c= 0; c < maxNumClusters; c++){
    if(org->active[c]){
      if (clusters[c]->empty() || clusters[c]->size() == 1) {
	numWeakClus++; //find number of invalid clusters
	//	arrSize[c] = fixSize;
	smallClus[c]=true;
      }
      else if(clusters[c]->size <= fixSize){
	numSmallClus++;//number of clusters that are already smaller than the new min size we want to have
      }
    }//end active if
  }//end for
  int remSize = 2/3*(clusters[cIndex]->size()); //approx. items we can reshuffle
  int allowedClus = remSize/fixSize;
  if(numEmptyClus >= allowedClus && numSmallClus != 0 && numSmallClus < numEmptyClus) {
    //this means that the current active small clusters might have a higher prob of getting discarded
    //in this extreme case, we probabilistically deactivate the remaining empty ones and reshuffle the item in the single sized one
    if(uniform01() < 0.5){
      for (int c = 0; c < maxNumClusters; c++) { //inactivating a cluster that's empty
	if (org->active[c]) {
	  if (clusters[c]->empty()) {
	    for (int f = 0; f < numFeatures; f++) { //since an inactive cluster can also contribute to crossover,
	      //we reinitialize the centroid since we know that the existing centroid is not good as it's empty
	      //reinitialize only those features which are close to the boundary
	      double gap = max[f]-min[f];
	      if(org->clusCenter[c][f] >= min[f] && org->clusCenter[c][f] <= min[f] + 0.1*gap) 
		org->clusCenter[c][f] = uniformInRange(min[f], (min[f] + (0.5*gap)));//check if 1.5*min[f] < max[f]
	      else if (org->clusCenter[c][f] <= max[f] && org->clusCenter[c][f] >= max[f] - 0.1*gap)
		org->clusCenter[c][f] = uniformInRange((max[f] - (0.5*gap)), max[f]);
	    }
	    org->activationThreshold = uniformInRange(0,activationThreshold);
	    org->active[c] = false;
	    org->numActiveClusters--;
	  }//end if empty
	  else if(clusters[c]->size == 1){
	  int addItemIndex = clusters[c]->at(0);
	  double temp_dist = 0;
	  double min = numeric_limits<double>::max();
	  int clusCtr = 0;
	  int minInd = -1;
	  while (clusCtr < maxNumClusters) {
	    if (clusCtr != c && org->active[clusCtr]) {
	      temp_dist = dist(org->clusCenter[clusCtr],itemsArray[addItemIndex]->items);
	      if (temp_dist < min) {
		temp_dist = min;
		minInd = clusCtr;
	      }
	    }
	    clusCtr++;
	  }
	  assert(minInd != -1);
	  clusters[minInd]->push_back(addItemIndex);
	  clusters[c]->clear(); // after element has been moved, cluster with single item made empty
	  org->active[c] = false; // deactivating the empty cluster
	  org->numActiveClusters--;
	  org->activationThreshold = uniformInRange(0,activationThreshold);
	  }//end elseif
	}
      }// end if active
    }//end if 0.5
    else{
      //keep only a few more active before forming heaps
      int keepActive = numEmptyClus - numSmallClus;
      int numLargeClus = org->numActiveClusters - numEmptyClus;
      while (org->numActiveClusters == (keepActive+numLargeClus)) {
	int i = uniformInRange(0, maxNumClusters - 1);
	if(org->active[i] && smallClus[i]){
	  org->activationThreshold[i] = uniformInRange(0,activationThreshold);
	  org->active[i] = false;
	  org->numActiveClusters--;
	}
      }//end while
      //form heaps for smaller clusters and then form large clusters
      double tempDist = 0.0;
      for (int clusIndex = 0; clusIndex < maxNumClusters; clusIndex++) {
	if(smallClus[clusIndex]){	
	  Heap *objHeap = new Heap(numItems);
	  for(int i = 0; i < numItems; i++) {
	    tempDist = dist(org->clusCenter[clusIndex],itemsArray[i]->items);
	    objHeap->Enqueue(dist, i); 
	    ItemUsed[i]=true;
	  }
	}//end if
      }//end for
    }//end else
  }//if extreme case 
 
  }*/
/*
 * Input parameters : A pointer to a chromosome, struct array that holds distance corresponding
 * to item and cluster center, size of the struct array and index of population's chromosome
 * This method reshuffles items equally into different active cluster centers of an individual
 * return type : void
 */
void DEMain::reshuffle(Individual* org, int numTriplesArray, int orgIndex, bool isInitial) {
	for (int c = 0; c < maxNumClusters; ++c) {
		clusters[c]->clear();
	}
	int fix_size = numItems / org->numActiveClusters; //maximum # of items that a cluster can hold
	if (fix_size < 2) {
		cout << numItems << " " << org->numActiveClusters;
	}
	assert(fix_size >= 2); //to assert that each cluster has at least 2 items
	int ctr = 0; //counter that goes through all elements in the triples array
	int numFullClusters = 0; //determines whether all clusters have reached their max capacity

	for (int i = 0; i < numItems; i++) //initializing counter array corresponding to how many times an item has appeared in the triples array
		ItemCounter[i] = 0;
	fill_n(ItemUsed, numItems, false); //initializing bool array such that no item has been added to cluster yet
	fill_n(ClusFull, maxNumClusters, false); // initializing bool array such that no cluster is full yet
	while (ctr < numTriplesArray) { //check for total # of items
		int itemInd = nearestNeighborTriples[ctr].itemIndex; //retrieving index of item
		int clusInd = nearestNeighborTriples[ctr].clustIndex; // retrieving index of cluster (s.t we have the closest item and cluster present
		if (!ItemUsed[itemInd]) { //only if the item has not already been placed in a cluster
			if (org->active[clusInd]) { //if retrieved cluster is valid still (could be deactivated if the cluster had 0/1 item allocated while computing fitness
				ItemCounter[itemInd]++;
				if (numFullClusters != org->numActiveClusters) { //if all clusters aren't full yet
					if (clusters[clusInd]->size() < fix_size) { //if the selected cluster isn't full yet
						clusters[clusInd]->push_back(itemInd); //add item to cluster
						ItemUsed[itemInd] = true; //mark that item as already added to a cluster
						if (isInitial) { //updating the arrays to trace back later
							trackerArray[itemInd][orgIndex] = clusInd;
						} else {
							offspringArray[itemInd] = clusInd;
						}
					} else { //if the cluster has reached its max capacity
						if (ItemCounter[itemInd] == org->numActiveClusters) { //check if a particular item wasn't added to any cluster due to the cluster being full
							clusters[clusInd]->push_back(itemInd); //added to the cluster its farthest from; otherwise the item wouldn't be added at all (special use case)
							ItemUsed[itemInd] = true;
						}
						if (!ClusFull[clusInd]) { //since the cluster has reached its max capacity, mark it as full
							ClusFull[clusInd] = true;
							numFullClusters++; //increment the counter that keeps track of #of clusters full
						}
					}
				} else { //if all clusters are full but still item remains push it in anyway
					clusters[clusInd]->push_back(itemInd);
					ItemUsed[itemInd] = true;
					if (isInitial) { //update the arrays to trace back
						trackerArray[itemInd][orgIndex] = clusInd;
					} else {
						offspringArray[itemInd] = clusInd;
					}
				}
			}
		}
		ctr++; //go through next item in triples array
	} //end while

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
	    for (vector<int>::size_type i = 0; i != clusters[c]->size(); i++) {
	      int a = clusters[c]->at(i);
	      sum += dist(itemsArray[a]->items, org->clusCenter[c]);
	      //sum += dist(itemsArray[a]->items, newClustCenters[c]);
	    } //end for
	    avgArr[c] = (sum / clusters[c]->size()); //finding the intra cluster distance for all active clusters
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
	/*	avgDB += avg;
	//trackOff << "DB index " << avg << endl;
	if (minDB > avg) {//keeping track of the minimum and maximum DB index encountered for all individuals in all generations
	  minDB = avg;
	}
	if (maxDB < avg) {
	  maxDB = avg;
	  }*/
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
      for (vector<int>::size_type i1 = 0; i1 != clusters[c]->size(); i1++) {
	int a = clusters[c]->at(i1);
	double maxIntraDist = numeric_limits<double>::min();
	//	double* tempItem = itemsArray[a]->items;
//	for (int f = 0; f < numFeatures; f++) {
//	  sumArr[f] += tempItem[f]; //to compute centroids
//	}
	for (vector<int>::size_type i2 = 0; i2 != clusters[c]->size();
	     i2++) { //finding max distance between items in a cluster
	  if (i2 != i1) {
	    int b = clusters[c]->at(i2);
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
      finalIntraSum += (intraClusSum / clusters[c]->size()); //updating the average intra cluster sum for all active centroids
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
  avgCS += csVal;
  if (minCS > csVal) { //updating the maximum and minimum CS index calculated for all individuals in every generation
    minCS = csVal;
  }
  if (maxCS < csVal) {
    maxCS = csVal;
  }
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
			int n = clusters[c1]->size();
			countIntraGroup += (n * (n - 1)) / 2;
			countInterGroup += (n * (numItems - n)) / 2;
			for (vector<int>::size_type i1 = 0; i1 != clusters[c1]->size(); i1++) {
				int a = clusters[c1]->at(i1);
				//find distance between elements in separate clusters; go through all elements in other clusters
				for (int cl = 0; cl < c1; cl++) {
					if (org->active[cl]) {
						for (vector<int>::size_type ind = 0;
								ind != clusters[cl]->size(); ind++) {
							int b = clusters[cl]->at(ind);
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
								ind != clusters[cl]->size(); ind++) {
							int b = clusters[cl]->at(ind);
							if (b < a) {
								sumInterCluster += distItem[b][a - (b + 1)];
							} else {
								sumInterCluster += distItem[a][b - (a + 1)];
							}
						}
					}
				}
				//finding distance between items in the same cluster
				for (vector<int>::size_type i2 = 0; i2 != clusters[c1]->size();
						i2++) {
					if (i2 != i1) {
						int b = clusters[c1]->at(i2);
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
	avgPB += pbIndex;
	//Point biserial : (d_b*d_w)(sqrt(w_d*b_d/t^2))/sd
	if (minPB > pbIndex) { //finding the min and max PB index found over all generations for all individuals
		minPB = pbIndex;
	}
	if (maxPB < pbIndex) {
		maxPB = pbIndex;
	}

	return pbIndex;
}

/*
 * Input parameters : Pointer to chromosome, index of population's chromosome,
 * bool parameter to know whether function called during initial setup or during DE
 * This method computes the clustering of a chromosome and then calls the fitness function inside
 * return  type : double; returns the fitness computed
 */
void DEMain::computeClustering(Individual* org, int popIndex, bool isInitial,
		int genNum, double min[], double max[]) {
  int min_index = -1;
  double temp_dist = 0;
  int itemIndex = 0;
  int ctrKnn = 0;
  bool isAdvExploration = false;
  //clearing the scratch space
  for (int c = 0; c < maxNumClusters; ++c) {
    clusters[c]->clear();
  }
  if(genNum >= 0.9*numGenerations){
    isAdvExploration = true;
  }
  
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
    /*    if(min_index == -1) {
      cout << genNum << " " << popIndex << " " << itemIndex << endl;
      cout << org->numActiveClusters << endl;
      for (int clusIndex = 0; clusIndex < maxNumClusters; clusIndex++) {
	if (org->active[clusIndex]) {
	  temp_dist = dist(org->clusCenter[clusIndex],itemsArray[itemIndex]->items);
	  cout << temp_dist << endl;
	}
      }
      }*/
    assert(min_index != -1);
    clusters[min_index]->push_back(itemIndex); //adding item to its closest cluster center
    if (isInitial) {
      trackerArray[itemIndex][popIndex] = min_index; //mainitaining tracking array
    } else {
      offspringArray[itemIndex] = min_index;
    }
    itemIndex++;
  } // end while
  
  //find out index of largest centroid
  int size = -1;
  int indexMax = -1;
  for(int c= 0; c <maxNumClusters; c++){
    if(org->active[c]){
      int clusSize = clusters[c]->size();
      if(size < clusSize){
	size = clusSize;
	indexMax = c;
      }
    }
  }//end for
  assert(indexMax != -1);
  //inactivating empty clusters
  for (int c = 0; c < maxNumClusters; c++) { //inactivating a cluster that's empty
    if (org->active[c]) {
      if(clusters[c]->size() < 2) {
	if (clusters[c]->empty()) {
	  for (int f = 0; f < numFeatures; f++) { //since an inactive cluster can also contribute to crossover,
	    //we reinitialize the centroid since we know that the existing centroid is not good as it's empty
	    //reinitialize only those features which are close to the boundary
	    double gap = max[f]-min[f];
	    if(org->clusCenter[c][f] >= min[f] && org->clusCenter[c][f] <= (min[f] + 0.1*gap)) {
	      org->clusCenter[c][f] = uniformInRange(min[f], (min[f] + (0.5*gap)));//check if 1.5*min[f] < max[f]
	    }
	    else if (org->clusCenter[c][f] <= max[f] && org->clusCenter[c][f] >= (max[f] - 0.1*gap)){
	      org->clusCenter[c][f] = uniformInRange((max[f] - (0.5*gap)), max[f]);
	    }
	  }
	  if(uniform01() < 0.5) {
	    org->activationThreshold[c] = uniformInRange(0.0, activationThreshold);
	    org->active[c] = false;
	    org->numActiveClusters--;
	  }
	}
	reshuffleValid(org, popIndex, isInitial, isAdvExploration, indexMax);
      }//end if size
    }//end org->active
  }//end for
  
}


/*
 * Input parameters : Pointer to chromosome, index of population's chromosome,
 * bool parameter to know whether function called during initial setup or during DE
 * This method calculates the fitness of a chromosome and returns it
 * return  type : double
 */
double DEMain::calcFitness(Individual* org, int popIndex, bool isInitial, int genNum, double min[], double max[]) {
	double fit = 0.0;
	computeClustering(org, popIndex, isInitial, genNum, min, max);
	//find out index of largest centroid
	int size = -1;
	int indexMax = -1;
	for(int c= 0; c <maxNumClusters; c++){
	  if(org->active[c]){
	    int clusSize = clusters[c]->size();
	    if(size < clusSize){
	      size = clusSize;
	      indexMax = c;
	    }
	  }
	}
	assert(indexMax != -1);
	//recomputing the centroid
	if(genNum >= 0.9*numGenerations){
	  // cout << "index " << indexMax << endl;
//	  for (int c = 0; c < maxNumClusters; c++) {
//	    for (int f = 0; f < numFeatures; f++) {
//	      newClustCenters[c][f] = 0.0; //clearing old centroids
//	    }
//	  }
	  for (int c = 0; c < maxNumClusters; c++) {
	    fill_n(sumArr, numFeatures, 0); //clearing sums stored in array
	    int clusterSize = clusters[c]->size();
	    int maxClusSize = clusters[indexMax]->size();
	    //	    cout << "clus size " << clusterSize << " maxClus " << maxClusSize << " and " << 0.2 * maxClusSize << endl;
	    if (org->active[c] && maxClusSize != 0 && clusterSize >= 0.2*(maxClusSize)) {
	      for (vector<int>::size_type j = 0; j != clusters[c]->size(); j++) { //go through all items in cluster
		int a = clusters[c]->at(j); //find item index stored in clusters
		double* tempItem = itemsArray[a]->items; //array to hold features corresponding to item in cluster
		for (int f = 0; f < numFeatures; f++) {
		  sumArr[f] += tempItem[f]; //this holds the sum of all features for all items in a single cluster to compute average later
		}
	      }
	      for (int f = 0; f < numFeatures; f++) {
		//newClustCenters[c][f] = sumArr[f] / clusters[c]->size(); //finding new centroids for the cluster
		//org->clusCenter[c][f] = newClustCenters[c][f]; //updating the cluster center
		//	cout << sumArr[f] << " size of cluster " << clusters[c]->size() << " ";
		org->clusCenter[c][f] = (sumArr[f]/clusters[c]->size());
		//		cout << "chromosome " << org->clusCenter[c][f] << " ";
	      } 
	      //cout<<endl;
	    }//end if
	  }
	  // cout << genNum << endl;
	  computeClustering(org, popIndex, isInitial, genNum, min, max); //recomputing clustering with new centroids
	}
//	for(int c = 0; c < maxNumClusters; c++) {
//	  if(org->active[c] && clusters[c]->empty()){
//	    org->active[c] = false;
//	    org->numActiveClusters--;
//	  }
//	}
	if (org->numActiveClusters >= minNumClusters) { //if after deactivating cluster centers above, an individual's active cluster centers fall below 2(kmin), assign it lowest fitness
		//based on the validity index find DB, CS or PB index for the individual to compute the fitness
	 
	  if (indexForFitness == 1) {
	    //minimization problem
	    double dBValue = calcDBIndex(org);
	    fit = 1 / dBValue;
	    //	    double csVal = calcCSIndex(org);

	    //	    double pbVal = calcPBIndex(org);
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
  do {
      s1 = uniformInRange(0, popSize-1);
  } while (s1 == orgIndex);

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
  //cout << "crossover method called" << endl;
  double fit = 0.0;
  int s1, s2, s3;
  double crossoverProb = crossoverProbability * ((numGenerations - genNum) / numGenerations);//based on formula in paper
  double f_scale = scale * (1 + uniform01()); //based on formula in paper
  //double f_scale = uniformInRange(scale, 0.8);//temp change for testing
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
    int valChange;
    for (int f = 0; f < numFeatures; f++) {//binomial crossover computation for centroids
      assert(popObject->org[orgIndex] != NULL);
      if(change) {
 //  if (uniform01() < cr_prob) {
	child->clusCenter[c][f] = popObject->org[s1]->clusCenter[c][f]
	  + f_scale*(popObject->org[s2]->clusCenter[c][f]- popObject->org[s3]->clusCenter[c][f]);
	if (child->clusCenter[c][f] < min[f]) {
	  //valChange =min[f] - child->clusCenter[c][f];
	  //child->clusCenter[c][f] = min[f] + valChange;

	  child->clusCenter[c][f] =  popObject->org[s1]->clusCenter[c][f]
	  + abs( f_scale*(popObject->org[s2]->clusCenter[c][f]- popObject->org[s3]->clusCenter[c][f]) );
	}
	else if (child->clusCenter[c][f] > max[f]){//if a feature for cluster center computed goes out of its bounds, assign it a random value
	  //valChange = child->clusCenter[c][f] - max[f];
	  //child->clusCenter[c][f] = max[f] - valChange;	 
	  child->clusCenter[c][f] =  popObject->org[s1]->clusCenter[c][f]
	   - abs( f_scale*(popObject->org[s2]->clusCenter[c][f]- popObject->org[s3]->clusCenter[c][f]));
	}
      } else {	
	child->clusCenter[c][f] = popObject->org[orgIndex]->clusCenter[c][f];
      }
      if ((child->clusCenter[c][f] < min[f]) || (child->clusCenter[c][f] > max[f])){
	countOOB++;
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
    int randNumClusters = uniformInRange(1, maxNumClusters/2);
    while(num < randNumClusters) {
    //while (num < minNumClusters) {
      int i = uniformInRange(0, maxNumClusters - 1);
      if(!child->active[i]){
	child->activationThreshold[i] = uniformInRange(activationThreshold, 1.0);
	child->active[i] = true;
	child->numActiveClusters++;
	num++;
      }
    }
  }
  fit = calcFitness(child, orgIndex, false, genNum, min, max);
  child->rawFitness = fit;
  /* if(orgIndex == popObject->bestOrgIndex){
    trackOff << "[[[[[[[[[[[Generation " << genNum << " ]]]]]]]]]]]" << endl;
    trackOff << "s1 s2 s3 " << s1<< " " << s2 << " " << s3 << endl;
    trackOff << "[[[[[[[[[[s1]]]]]]]]]]" <<endl;
    trackOff << *(popObject->org[s1]) << endl;
    printClusters(s1);
    trackOff << "[[[[[[[[[[s2]]]]]]]]]]" <<endl;
    trackOff << *(popObject->org[s2]) << endl;
    printClusters(s2);
    trackOff << "[[[[[[[[[[s3]]]]]]]]]]" <<endl;
    trackOff << *(popObject->org[s3]) << endl;
    printClusters(s3);
    trackOff << "[[[[[[[[[[parent]]]]]]]]]]" <<endl;
    trackOff << *(popObject->org[orgIndex]) << endl;
    printClusters(orgIndex);
    trackOff<< "[[[[[[[[[offspring]]]]]]]]]" << endl;
    trackOff << *child << endl;
    }*/
  return child;

}

/*
 * This method runs the DE algorithm
 */
void DEMain::run(double min[], double max[], string filename) {
	double g = 0;
	fill_n(isReplaceOrg, popSize, false);
	double fitness;
	try {
	  //	trackFile.open("clusters.txt");
		trackOff.open("record.txt");
		trackOff << setiosflags(ios::left);
		//		trackOff << setw(5) << "Gen" << setw(5) << "Avg DB" << "|" << setw(5) << "Avg CS" << "|" << setw(5) << "Avg PB" << "|" << setw(10) << "Max DB"
		//	 << "|" << setw(5) << "Max CS" << "|" << setw(5) << "Max PB" << "|" << setw(10) << "Min DB" <<"|" << setw(5) << "Min CS" << "|" << setw(5) << "Min PB" << endl;
		trackOff << setw(5) << "Gen" << setw(5) << "Avg DB" << "|" << setw(10) << "Max DB"  << "|" << setw(10) << "Min DB" << endl;
		while (g < numGenerations) {    //till max generation reached
		  minDB = numeric_limits<double>::max();
		  maxDB = numeric_limits<double>::min();
		  /* minCS = numeric_limits<double>::max();
		  maxCS = numeric_limits<double>::min(); 
		  minPB = numeric_limits<double>::max();
		  maxPB = numeric_limits<double>::min();*/ 
		  avgDB = 0.0, avgCS = 0.0, avgPB = 0.0;
		  Population* newpop = new Population(maxNumClusters, numFeatures, popScaleFactor);
		  for (int p = 0; p < popSize; p++) {
				Individual *offspring;
				//				trackOff << "[[[[[[[[ Generation " << g << " ]]]]]]]]" << endl;
				offspring = crossover(p, g, min, max); //generate an offspring my performing DE crossover

		//	if(offspring->numActiveClusters == numClasses){
		//	  int* clusClass = new int[numClasses+1];
		//	 
		//	  cout << "----------------------------------" <<endl;
		//	  for (int c = 0; c < maxNumClusters; c++) {
		//	    if (offspring->active[c]) {			//prints out the clusters
		//	      for(int c = 0; c < numClasses+1; c++){
		//		clusClass[c] = 0;
		//	      }
		//	      for (vector<int>::size_type i = 0; i != clusters[c]->size(); i++) {
		//		int itemIndex = clusters[c]->at(i);
		//		clusClass[itemsArray[itemIndex]->typeClass]++;
		//	      }
		//	      cout << "Elements of cluster " << c << " :" << endl;
		//	      cout << "Total # of items in cluster " << clusters[c]->size() << endl;
		//	      for(int cl = 0; cl < numClasses+1; cl++){
		//		if(clusClass[cl] != 0){
		//		  cout << "Class " << cl << " : " << clusClass[cl] << " number of items";
		//		  cout << endl;
		//		}
		//	      }
		//	    
		//	    }//end active if
		//	  }//end for printing clusters
		//	  cout << "DB index is " << 100*(1/offspring->rawFitness) << endl;
		//	}//end if 
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
				  for (int i = 0; i < numItems; i++) {
				    trackerArray[i][p] = offspringArray[i]; //updating the parent chromosome replaced with new cluster centers of offspring
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

		  //		  trackOff << setw(5) << g << setw(5) << avgDB << "|" << setw(5) << avgCS << "|" << setw(5) << avgPB << "|" << setw(10) << maxDB
		  //	   << "|" << setw(5) << maxCS << "|" << setw(5) << maxPB <<"|" << setw(10) << minDB <<"|" << setw(5) << minCS << "|" << setw(5) << minPB << endl;
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
		  //	for (int p = 0; p < popSize; p++) { //for testing purposes (keeping track of the fitness change as the algorithm progresses)
		  //		trackFile << setiosflags(ios::left);
		  //		trackFile << setw(5) << g << setw(12)
		  //				<< popObject->org[p]->rawFitness << setw(5) << p << endl;
		  //	}
		  g++; //increment generation number
		}
		cout << endl;
		//		cout << "Out of bounds " << countOOB << endl;
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
		//	trackFile.close();
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
	outputFile.open(filename.c_str());
	outputFile << "The final clusters obtained are:" << endl;
	int clusIndex = -1;
	for (int c = 0; c < maxNumClusters; ++c) {
		clusters[c]->clear();
	}
	int* clusClass = new int[numClasses+1];
	for(int c = 0; c < numClasses+1; c++){
	  clusClass[c] = 0;
	}
	Individual* org = popObject->org[bestPopIndex];
	for (int i = 0; i < numItems; i++) {
		clusIndex = trackerArray[i][bestPopIndex];
		if (org->active[clusIndex]) {
			clusters[clusIndex]->push_back(i);
		}
	}
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
	    activeCount++;
	    for (vector<int>::size_type i = 0; i != clusters[c]->size(); i++) {
	      int itemIndex = clusters[c]->at(i);
	      clusClass[itemsArray[itemIndex]->typeClass]++;
	     
	     // arr = itemsArray[itemIndex]->items;
	     // for (int f = 0; f < numFeatures; f++) {
	     //	outputFile << arr[f] << " ";
	     // }
	     // outputFile << itemsArray[itemIndex]->typeClass;
	     // outputFile << endl;
	    }
	    outputFile << "Elements of cluster " << c << " :" << endl;
	    outputFile << "Total # of items in cluster " << clusters[c]->size() << endl;
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

	for (int c = 0; c < maxNumClusters; ++c) {
	  clusters[c]->clear();
	}
	org = popObject->org[worstInd];
	for (int i = 0; i < numItems; i++) {
	  clusIndex = trackerArray[i][worstInd];
	  if (org->active[clusIndex]) {
	    clusters[clusIndex]->push_back(i);
	  }
	}
	

	outputFile << "-------------------------------------------" << endl;
	outputFile << "Worst fitness clusters are as follows" << endl;
	outputFile << endl;
	activeCount = 0;
	for (int c = 0; c < maxNumClusters; c++) {//prints out the worst output after changes are made
	  if (org->active[c]) {
	    activeCount++;
	    for (vector<int>::size_type i = 0; i != clusters[c]->size(); i++) {
	      int itemInd = clusters[c]->at(i);
	      clusClass[itemsArray[itemInd]->typeClass]++;
//	      arr = itemsArray[itemInd]->items;
//	      for (int f = 0; f < numFeatures; f++) {
//		outputFile << arr[f] << " ";
//	      }
//	      outputFile << itemsArray[itemInd]->typeClass;
//	      outputFile << endl;
	    }
	    outputFile << "Elements of cluster " << c << " :" << endl;
	    outputFile << "Total # of items in cluster " << clusters[c]->size() << endl;
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
	for (int c = 0; c < maxNumClusters; ++c) {
		clusters[c]->clear();
	}
	for (int i = 0; i < numItems; i++) {
	  int min_index = -1;
	  double min = numeric_limits<double>::max();
	  for (int c = 0; c < maxNumClusters; c++) {
	    if (org->active[c]) {
	      double temp_dist = dist(org->clusCenter[c], itemsArray[i]->items);
	      if (temp_dist < min) {
		min = temp_dist;
		min_index = c;
	      }
	    }
	  }
	  assert(min_index != -1);
	  clusters[min_index]->push_back(i);
	}

	/*outputFile << "-------------------------------------------" << endl;
	outputFile << "Worst fitness cluster without tweaking" << endl;
	for (int c = 0; c < maxNumClusters; c++) {//prints out the worst output w/o making any change like reshuffling or deactivating clusters.
		if (org->active[c] && !clusters[c]->empty()) {
			outputFile
					<< "------------------------------------------------------------"
					<< endl;
			outputFile << "Elements of cluster : " << c << endl;
			for (vector<int>::size_type j = 0; j != clusters[c]->size(); j++) {
				int itemInd = clusters[c]->at(j);
				arr = itemsArray[itemInd]->items;
				for (int f = 0; f < numFeatures; f++) {
					outputFile << arr[f] << " ";
				}
				outputFile << itemsArray[itemInd]->typeClass;
				outputFile << endl;
			}
		}
	}
	valNMI = MI(clusters, origClustering, worstInd, 0, true);
	outputFile << "NMI for worst cluster w/o reshuffling is " << valNMI << endl;*/
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
	sort(clusters[c1]->begin(), clusters[c1]->end());
	size = clusters[c1]->size();
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
		sort(clusters[c2]->begin(), origClustering[c2]->end());
		isSorted[c2] = true;
	      }
	      set_intersection(clusters[c1]->begin(), clusters[c1]->end(),
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
		if(clusters[c1]->at(it1) == origClustering[c2]->at(it2)){
		  ctr++;
		  it1++; it2++;
		}
		else if(clusters[c1]->at(it1) > origClustering[c2]->at(it2)) it2++;
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
  int clusIndex1, clusIndex2;
  for (int i = 0; i < numItems-1; i++) {
    for(int j = i+1; j < numItems; j++) {
      clusIndex1 = trackerArray[i][popInd1];
      clusIndex2 = trackerArray[j][popInd1];
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
