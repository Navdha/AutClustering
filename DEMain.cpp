/*
 * DE.cpp
 *
 *  Created on: Jul 8, 2015
 *      Author: Navs
 */

#include "DEMain.h"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <cassert>
#include <random>     // needed for the mersenne twister random number generator
#include <functional> // needed for bind
//#include <utility>
#include <algorithm>
#include <stdexcept>

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

/*int uniformInRange(int from, int to) { // generate a random uniformly in the range [from, to]

  return  (rand() % (to - from+1) + from);
  }
  double uniformInRange(double from, double to) { // generate a random uniformly in the range [from, to]
  return  uniformInRange((int)(from*10000), (int)(to*10000))/10000;
  }

  double uniform01(){
  return uniformInRange(0.0, 100.0)/100.0;
  }*/
struct compare_first_only {
  template<typename T1, typename T2>
  bool operator()(const std::pair<T1, T2>& p1, const std::pair<T1, T2>& p2) {
    return p1.first < p2.first;
  }
};

ofstream trackFile;
ofstream traceF;
//int changeCounter =0;
double minDB = std::numeric_limits<double>::max();
double maxDB= std::numeric_limits<double>::min();
int compare(const void *p1, const void *p2){
  Dist_IC *elem1 = (Dist_IC *)p1;
  Dist_IC *elem2 = (Dist_IC *)p2;
  if(elem1->distance < elem2->distance)
    return -1;
  else if(elem1->distance > elem2->distance)
    return 1;
  else
    return 0;
}


DEMain::DEMain(int kmax, int dim, int gen, int** placeholder, Item** items,
	       int itemSize, bool calcDB) {
  // TODO Auto-generated constructor stub
  //cout << "Constructor called from DEMain class" << endl;
  p = new Population(kmax, dim);
  strategy = stRand1Bin;
  generations = gen;
  probability = 0.5;
  scale = 0.5;
  pSize = dim * 10;
  this->kmax = kmax;
  this->dim = dim;
  numItems = itemSize;
  tracker = placeholder;
  attr = items;
  clusters = new vector<int>*[kmax];
  isDB = calcDB;
  distItem = new double*[numItems];
  offspring_arr = new int[numItems];
  for (int count = 0; count < kmax; count++)
    {
      clusters[count] = new vector<int>;
    }
  int size = kmax*numItems;
  knn = new Dist_IC[size];
}

DEMain::~DEMain() {
  // TODO Auto-generated destructor stub
  //	cout << "DEMAIN destructor" << endl;
  delete p;
  for (int i = 0; i < numItems; i++){
    delete [] tracker[i];
    delete [] attr[i];
    delete [] distItem[i];
  }
  for(int i = 0; i < kmax; ++i) {
    delete clusters[i];
  }
  delete[] tracker;
  delete[] attr;
  delete[] distItem;
  delete [] clusters;
  delete [] offspring_arr;
  delete [] knn;
}

/*
 * This method calculates and stores the distance between all items
 * Used to calculate CS index
 */
void DEMain::calcDistBtwnItems(){
     for (int i = 0; i < numItems; i++) {
	distItem[i] = new double[numItems - i];
	for (int j = i + 1; j < numItems; j++) {
	  distItem[i][j - i - 1] = dist(attr[i]->items, attr[j]->items);
	}//end j for
   }// end i for
}

/* This method takes an input two arrays min and max
 * that hold the minimum and maximum value of each attribute
 * Purpose of this method is to setup the initial population
 * return type is void
 * no output expected
 */
void DEMain::setup(double min[], double max[]) {
  //initialize chromosomes for the first time
 // cout<< "Setup method called" <<endl;

  for (int i = 0; i < pSize; i++) {
    bool isValid = false;
    double fitn = 0.0;
    Individual* temp = new Individual(kmax, dim);
    while(!isValid) {
      int ctr_act = 0;
      for (int j = 0; j < kmax; j++) {
	temp->threshold[j] = uniform01();
	//	cout << temp->threshold[j];
	if (temp->threshold[j] > 0.5) {
	  temp->active[j] = true;
	  ctr_act++;
	} else
	  temp->active[j] = false;
	for (int k = 0; k < dim; k++) {
	  temp->clusCenter[j][k] = uniformInRange(min[k], max[k]);
	  //cout << temp->clusCenter[j][k] << " ";
	}
      }
      assert(ctr_act <= kmax);
      temp->active_ctr = ctr_act;
      //code to check kmin = 2
      if(temp->active_ctr < 2){
	int num = temp->active_ctr;
	while (num < 2) {
	  int i = uniformInRange(0, kmax - 1);
	  if(!temp->active[i]){
	    temp->threshold[i] = uniformInRange(0.5, 1.0);
	    temp->active[i] = true;
	    temp->active_ctr++;
	    num++;
	  }
	}
      }
      fitn = calcFitness(temp, i, true, -1);
      if( fitn != -1){
	isValid = true;
      }
    } //end while
    assert(fitn != 0.0);
    temp->setFitness(fitn);
    p->chromosome[i] = temp;
  }
  int bestInd = 0;
  double fitness = p->chromosome[0]->rawFitness;
  for (int k = 1; k < pSize; k++) {
    if (fitness < p->chromosome[k]->rawFitness) {
      bestInd = k;
      fitness = p->chromosome[k]->rawFitness;
    }
  }
  p->bestChromosomeIndex = bestInd;

}


/*
 * Input parameters: pointers to arrays that hold
 * cluster center and item features
 * This method finds the distance between cluster center and an item
 * returns the distance calculated
 */
double DEMain::dist(double* x, double* y) {
  double Sum = 0.0;
  double distance = 0.0;
  for (int i = 0; i < dim; i++) {
    Sum = Sum + pow((x[i] - y[i]), 2.0);
    //		distance = sqrt(Sum);
  }
  distance = sqrt(Sum);
  return distance;
}

/*double* DEMain::avgDist(Individual* org) {
  cout << "avgDist function called" << endl;
  double* temp = new double[kmax];
  //double *tempArr = &temp;
  double* d2;
  vector<int>* c = new vector<int>;
  //int ind = -1;
  for (int i = 0; i < kmax; i++) {
  if (org->active[i]) {
  double sum = 0.0;
  delete [] d2;
  d2 = org->clusCenter[i];
  //			delete [] c;
  c = org->clusters[i];
  //for(vector<int>::iterator it = c[i]->begin(); it != c[i]->end(); ++it) {
  //for (vector<int>::size_type j = 0; j != size; j++) {
  for (std::vector<int>::const_iterator j = c->begin(); j != c->end();
  ++j) {
  sum += dist(attr[*j]->items, d2);
  }
  temp[i] = sum / (c->size());
  }
  }

  return temp;
  }*/

/*int DEMain :: compare(const void *p1, const void *p2) {
  Dist_IC *elem1 = (Dist_IC *)p1;
  Dist_IC *elem2 = (Dist_IC *)p2;
  if(elem1->distance < elem2->distance)
  return -1;
  else if(elem1->distance > elem2->distance)
  return 1;
  else
  return 0;
  }*/

/*
 * Input parameters : A pointer to a chromosome, struct array that holds distance corresponding
 * to item and cluster center, size of the struct array and index of population's chromosome
 * This method reshuffles items equally into different active cluster centers of an individual
 * return type : void
 */
void DEMain :: reshuffle(Individual* org, int size, int indpop, bool isInitial, int initialActive){//need to think
//  cout << "reshuffle method called" <<endl;
  for(int i = 0; i < kmax; ++i) {
    clusters[i]->clear();
  }
  int fix_size = numItems/org->active_ctr;
  if(fix_size < 2){
    cout << numItems << " " << org->active_ctr;
  }
  assert(fix_size >= 2);
  int ctr = 0;
  int numFullClusters = 0;
  bool* ItemUsed = new bool[numItems]();
  bool* ClusFull = new bool[kmax]();
  int* ItemCounter = new int[numItems];
  for(int i = 0; i < numItems; i++) ItemCounter[i] = 0;
  while (ctr < size) { //check for total # of items
    int itemInd = knn[ctr].itemIndex;
    int clusInd = knn[ctr].clustIndex;
    if (!ItemUsed[itemInd]) {
	ItemCounter[itemInd]++;
      if (org->active[clusInd]) {
	if (numFullClusters != org->active_ctr) {
	  if (clusters[clusInd]->size() < fix_size) {
	    clusters[clusInd]->push_back(itemInd);
	    ItemUsed[itemInd] = true;
	    if (isInitial) {
	      tracker[itemInd][indpop] = clusInd;
	    } else {
	      offspring_arr[itemInd] = clusInd;
	    }
	  }
	  else {
	    if(!ClusFull[clusInd]) {
	      ClusFull[clusInd] = true;
	      numFullClusters++;
		if(ItemCounter[itemInd] == initialActive){
		    clusters[clusInd]->push_back(itemInd);
		    ItemUsed[itemInd] = true;
		    cout << "Item index affected " << itemInd << endl;
		}
	    }
	  }
	}
	else {
	  clusters[clusInd]->push_back(itemInd);
	  ItemUsed[itemInd] = true;
	  if (isInitial) {
	    tracker[itemInd][indpop] = clusInd;
	  } else {
	    offspring_arr[itemInd] = clusInd;
	  }
	}
      }
    }
    ctr++;
  }//end while
  int totalCount = 0;
  for(int i = 0; i < kmax; i++){
  if(org->active[i]){
   totalCount += clusters[i]->size();
  }
  }
cout << "knn size " << size << endl;
  cout << totalCount << endl;
int notUsed;
for(int i = 0 ; i < numItems; i++){
if(!ItemUsed[i]){
notUsed = i; 
cout << "Item index not used yet " <<i << endl;
} 
}
for(int i = 0; i < size; i++) {
if(notUsed == knn[i].itemIndex) { cout << "Index of item not used in knn array " << i << endl;
}
}
cout << "Number of active cluster centers initially " << initialActive << endl;
cout << "Number of times an item is encountered in knn " << ItemCounter[notUsed] << endl; 
  assert(totalCount == numItems);
  delete [] ItemUsed;
  delete [] ClusFull;
}

/*
 * Input parameters : Pointer to chromosome, index of population's chromosome,
 * bool parameter to know whether function called during initial setup or during DE
 * This method calculates the fitness of a chromosome and returns it
 * return  type : double
 */
double DEMain::calcFitness(Individual* org, int index, bool isInitial, int genNum) {// not using index right now
 // cout << "calcFitness method called" << endl;
  double fit = 0.0;
  double maxValue = 0.0;
  double sum = 0.0;
  //	double eps = 0.5;
  int min_index = -1;	
//  int str_size = numItems * org->active_ctr;
  //cout << "value of str_size " << str_size << endl;
  double temp_dist;
  int ctr = 0;
  int vals = 0;
  for(int i = 0; i < kmax; ++i) {
    clusters[i]->clear();
  }
 int initialActive = org->active_ctr;
  while (ctr < numItems) { //form clusters
    min_index = -1;
    double min = numeric_limits<double>::max();
    for (int i = 0; i < kmax; i++) {
      if (org->active[i]) {
	temp_dist = dist(org->clusCenter[i], attr[ctr]->items);
	knn[vals].distance = temp_dist;
	knn[vals].clustIndex = i;
	knn[vals].itemIndex = ctr;
	if (temp_dist < min) {
	  //cout << "and current min = " << min << "current temp dist = " << temp_dist << endl;				
	  min = temp_dist;
	  min_index = i;
	}
	vals++;
      }

    }
   // cout << "# of active clusters before reshuffle after knn  " << org->active_ctr << endl;
  //  cout << "value of vals " << vals << endl;
//    assert(str_size == vals); 
    assert(min_index != -1);
    clusters[min_index]->push_back(ctr);
    if(isInitial) {
      tracker[ctr][index] = min_index;
    }
    else {
      offspring_arr[ctr] = min_index;
    }
    ctr++;

  }// end while
 cout << "value of vals " << vals << endl;
  //	trackFile.open("clusters.txt", ofstream::app);
  //	trackFile << "Number of clusters and size " << endl;
  //	trackFile << org->active_ctr << " " ;
  int tempArr[org->active_ctr];
  int tempC = 0;
  for (int i = 0; i < kmax; i++) {
    if (org->active[i]) {
      tempArr[tempC]= clusters[i]->size();
      tempC++;	
      //trackFile << clusters[i]->size() << " ";
      //	cout << "size of cluster " << i << " is " << clusters[i]->size() << endl;
    }
  }	
 
  for (int i = 0; i < kmax; i++) {
    if (org->active[i]) {
      if (clusters[i]->empty()) {
	org->active[i] = false;
	org->active_ctr--;
      }
    }

  }
  for (int i = 0; i < kmax; i++) {
    if (org->active[i]) {
      if (clusters[i]->size() == 1) {
	if (uniform01() < 0.5) {
	  int a = clusters[i]->at(0);
	  temp_dist = 0;
	  double min = numeric_limits<double>::max();
	  ctr = 0;
	  int minInd = -1;
	  while (ctr < kmax) {
	    if (ctr != i && org->active[ctr]) {
	      temp_dist = dist(org->clusCenter[ctr], attr[a]->items);
	      if (temp_dist < min) {
		temp_dist = min;
		minInd = ctr;
	      }
	    }
	    ctr++;
	  }
	  assert(minInd != -1);
	  clusters[minInd]->push_back(a);
	  clusters[i]->clear();
	  org->active[i] = false;
	  org->active_ctr--;
	} else {
	  //reshuffle items in clusters
	  qsort(knn, vals, sizeof(Dist_IC), compare);
	  reshuffle(org, vals, index, isInitial, initialActive);
	  break;
	}//random if-else end
      }//size if
    }//active if
  }//end for
  double avg = 0.0;
  if(org->active_ctr > 1) {
  if(isDB) {   
    double avgArr[kmax];
    for (int i = 0; i < kmax; i++) {
      if (org->active[i]) {
	sum =0.0;
	for (vector<int>::size_type j = 0; j != clusters[i]->size(); j++){
	  int a = clusters[i]->at(j);
	  sum += dist(attr[a]->items, org->clusCenter[i]);			
	}//end for
	avgArr[i] = sqrt(sum / clusters[i]->size());
      }//end if
		
    }//end outer for
    sum = 0.0;	
    for (int i = 0; i < kmax; i++) {
      maxValue = 0.0;
      for (int j = 0; j < kmax; j++) {
	if (i != j && org->active[i] && org->active[j]) {
	  double temp = avgArr[i] + avgArr[j];
	  temp /= dist(org->clusCenter[i], org->clusCenter[j]);
	  if (temp > maxValue)
	    maxValue = temp;
	}
      }
      sum += maxValue;
    }
    double avg = sum / org->active_ctr;
    //cout << "DB Index is " << avg;
    if(minDB > avg){
      minDB = avg;
    }
    if(maxDB < avg){
      maxDB = avg;
    }
    fit = (1 / avg);// + eps);
  }//end isDB if
  else {
  calcDistBtwnItems();
  double** newClustCenters = new double*[kmax];
  double finalIntraSum = 0.0;
  double finalInterSum = 0.0;
  double csVal = 0.0;
  for(int i = 0; i < kmax; i++){
	double sumArr[dim];
	double maxIntraDist = numeric_limits<double>::min();
	double intraClusSum = 0.0;
	double tempIntraDist;
	fill_n(sumArr, dim, 0);
	if(org->active[i]){
	   for (vector<int>::size_type j = 0; j != clusters[i]->size(); j++) {
		int a = clusters[i]->at(j);
		double* tempItem = attr[a]->items;
		for(int k = 0; k < dim; k++){
		    sumArr[k] += tempItem[k];
		}
		for(vector<int>::size_type k = 0; k != clusters[i]->size(); k++) {
		    if (k != j) {
			int b = clusters[i]->at(k);
			if (b < a) {
			    tempIntraDist = distItem[b][a - (b + 1)];
			} else {
			    tempIntraDist = distItem[a][b - (a + 1)];
			}
			if (tempIntraDist > maxIntraDist) {
			    maxIntraDist = tempIntraDist;
			}
		    }
		}//end for k
		intraClusSum += maxIntraDist;
	}// end for j
	finalIntraSum += (intraClusSum / clusters[i]->size());
	newClustCenters[i] = new double[dim];
	for(int m = 0; m < dim; m++){
	    newClustCenters[i][m] = sumArr[m]/clusters[i]->size(); //finding new centroids
	}
  }//endif
 }// end for i
  double interClusSum;
  for (int i = 0; i < kmax; i++) {
	double minInterDist = numeric_limits<double>::max();
	interClusSum = 0.0;
	double tempInterDist;
	for (int j = 0; j < kmax; j++) {
	     if (i != j && org->active[i] && org->active[j]) {
		tempInterDist = dist(newClustCenters[i], newClustCenters[j]);
		if (tempInterDist < minInterDist){
		    tempInterDist = minInterDist;
		}
	     }//end outer if
	}//end for j
	interClusSum += tempInterDist;
   }//end for i
   finalInterSum += interClusSum;
   csVal = finalIntraSum/finalInterSum;
   fit = 1/csVal;
  }//end else
    trackFile.open("clusters.txt", ofstream::app);
    trackFile << setiosflags(ios::left) ;
    trackFile << setw(5) << genNum << setw(12) << fit*100 << setw(5) << org->active_ctr << setw(5) <<tempC ;
    for(int i = 0; i < kmax; i++) {
      if(org->active[i]){ trackFile << setw(5) << clusters[i]->size() ;}
    }       
    trackFile << setw(5) << "|" ;
    for(int i = 0; i <tempC; i++){trackFile << setw(5) <<  tempArr[i] ;}
    trackFile << setw(5) << index << endl;
    trackFile.close();

    return fit*100;
  }
  else {
    cout << "Invalid chromosome" << endl;
    return -1;
  }

}

/*
 * Input parameters : indexes for chromosomes
 * This method makes sure that we get unique indexes for performing crossover
 * return type: void
 */

void DEMain::selectSamples(int org, int *s1, int *s2, int *s3) {
  if (s1) {
    do {
      *s1 = uniformInRange(0, pSize-1);
    } while (*s1 == org);
  }

  if (s2) {
    do {
      *s2 = uniformInRange(0, pSize-1);
    } while ((*s2 == org) || (*s2 == *s1));
  }

  if (s3) {
    do {
      *s3 = uniformInRange(0, pSize-1);
    } while ((*s3 == org) || (*s3 == *s2) || (*s3 == *s1));
  }
 // cout << "selectSamples called" << endl;
  return;
}

/*
 * Input parameters: index for chromosome chosen and the generation number
 * This method performs crossover to create an offspring
 * returns pointer to offspring created
 */
Individual* DEMain::crossover(int org, int gen, double min[], double max[]) {
  //cout << "crossover method called" << endl;
  int s1, s2, s3;
  double cr_prob = probability * ((generations - gen) / generations);
//  scale = 0.1;  //temporary change
  double f_scale = scale * (1+uniform01());
  selectSamples(org, &s1, &s2, &s3);
  Individual* child = new Individual(kmax, dim);
  int counter = 0;
  for (int j = 0; j < kmax; j++) {
    if(uniform01() < cr_prob){
      child->threshold[j] =  p->chromosome[s1]->threshold[j] +
	f_scale*(p->chromosome[s2]->threshold[j] -  p->chromosome[s3]->threshold[j]);
    }
    else{
      child->threshold[j] =  p->chromosome[org]->threshold[j];
    }
    bool change = uniform01() < cr_prob? true : false;
    for (int i = 0; i < dim; i++) {
      assert(p->chromosome[org] != NULL);
      if (change) {
	child->clusCenter[j][i] = p->chromosome[s1]->clusCenter[j][i]
	  + f_scale*(p->chromosome[s2]->clusCenter[j][i]- p->chromosome[s3]->clusCenter[j][i]);

	//cout << p->chromosome[s1]->clusCenter[j][i] << " " << p->chromosome[s2]->clusCenter[j][i] <<  " " << p->chromosome[s3]->clusCenter[j][i] << " " << child->clusCenter[j][i] << endl;
      } else {
	child->clusCenter[j][i] = p->chromosome[org]->clusCenter[j][i];

      }

      if ((child->clusCenter[j][i] < min[i]) || (child->clusCenter[j][i] > max[i])){
	child->clusCenter[j][i] = uniformInRange(min[i], max[i]);
      }
    }//for i
    if (child->threshold[j] > 1 || child->threshold[j] < 0)
      child->threshold[j] = uniform01();
    if (child->threshold[j] > 0.5) {
      child->active[j] = true;
      counter++;
    } else
      child->active[j] = false;
  }//for j
  assert(counter <= kmax);
  child->active_ctr = counter;
  if(child->active_ctr < 2){
    int num = child->active_ctr;
    while (num < 2) {
      int i = uniformInRange(0, kmax - 1);
      if(!child->active[i]){
	child->threshold[i] = uniformInRange(0.5, 1.0);
	child->active[i] = true;
	child->active_ctr++;
	num++;
      }
    }
  }
 
  traceF << "Parent index : " << org << " s1, s2, s3 index " << s1 << " " << s2 << " " << s3 <<endl;
  return child;

}

/*
 * This method runs the DE algorithm
 */
void DEMain::run(double min[], double max[]) {
 // cout << "run method called" << endl;
  int i = 0;
  bool * new_pop = new bool[pSize];
  double fitness;
  try {
    traceF.open("popEval.txt", ofstream::app);

    while (i < generations) {
      Population* newpop = new Population(kmax, dim);
      traceF << "For generation " << i << endl;
      for (int c = 0; c < pSize; c++) {
	//cout << c << " iteration" << endl;
	Individual *offspring;
	offspring = crossover(c, i, min, max);
	fitness = calcFitness(offspring, c, false, i);
	offspring->setFitness(fitness);
	if (p->chromosome[c]->rawFitness <= offspring->rawFitness) {
	  traceF << endl;
	  traceF << "Parent replaced" << endl;
	  new_pop[c] = true;
//	  cout << "offspring added" << endl;
	  newpop->chromosome[c] = offspring;
	  for (int d = 0; d < numItems; d++) {
	    tracker[d][c] = offspring_arr[d]; //updating the parent chromosome replaced with new cluster centers of offspring
	  }
	} else {
	  traceF << endl;
	  traceF << "Parent not replaced" << endl;
	  new_pop[c] = false;
//	  cout << "offspring discarded" << endl;
	  delete offspring;
	  newpop->chromosome[c] = p->chromosome[c];
	  //delete [] offspring_arr;
	}
      }
     // cout << "Generation " << i << " completed" << endl;
      //assert(newpop != NULL);
      for (int c = 0; c < pSize; c++) {
	if (new_pop[c]) {
	  delete p->chromosome[c];
	}
      }
      delete [] p->chromosome;
      p = newpop;
      int bestInd = 0;
      fitness = p->chromosome[0]->rawFitness;
      for (int k = 1; k < pSize; k++) {
	if (fitness < p->chromosome[k]->rawFitness) {
	  bestInd = k;
	  fitness = p->chromosome[k]->rawFitness;
	}
      }
      p->bestChromosomeIndex = bestInd;
      //	trackFile.open("clusters.txt", ofstream::app);
      //	trackFile << "Generation " << i << " finished" << endl;
      //	trackFile << "Best chromosome index is " << bestInd << " and fitness is " << fitness << endl;
      //	trackFile << "--------------------------------------" << endl;
      traceF << "-----------------------------------------------" << endl;
      //	trackFile.close();
      i++;
    }
    traceF.close();
    cout << endl;
    cout << "Stopped at generation " << i << endl;
    //find best chromosome
    /*	int bestInd = 0;
	double fitness = p->chromosome[0]->rawFitness;
	for (int k = 1; k < pSize; k++) {
	if (fitness < p->chromosome[k]->rawFitness) {
	bestInd = k;
	fitness = p->chromosome[k]->rawFitness;
	}
	}*/
    //assert(bestInd != -1);
    report(p->bestChromosomeIndex);
  }
  catch (exception& e) {
    cerr << e.what() << endl;
  }
}

/*
 * Input parameter : index for the best chromosome of the population
 * This method outputs the clusters for the best chromosome in the population to a file
 * return type : void
 */
void DEMain::report(int index) {
  ofstream outputFile;
  outputFile.open("data_1000.txt");
  outputFile << "The final clusters obtained are:" << endl;
  int clus_index = -1;
  for(int i = 0; i < kmax; ++i) {
    clusters[i]->clear();
  }
  Individual* org = p->chromosome[index];
  for(int i = 0; i < numItems; i++){
    clus_index = tracker[i][index];
    if(org->active[clus_index]){
      clusters[clus_index]->push_back(i);
    }
  }
	
  double* arr;
  int activeCount = 0;
  for (int k = 0; k < kmax; k++) {
    if (org->active[k]) {
      activeCount++;
      outputFile << "Elements of cluster " << k << " :" << endl;
      for (vector<int>::size_type j = 0; j != clusters[k]->size(); j++){
	int a = clusters[k]->at(j);
	arr = attr[a]->items;				
	for (int m = 0; m < dim; m++) {
	  outputFile << arr[m] << " ";
	}
	outputFile << attr[a]->typeClass;
	outputFile << endl;
      }
    }
  }
  /*	for (int i = 0; i < kmax; ++i) {
	clusters[i]->clear();
	}
	for(int i = 0; i < numItems; i++) {
	int min_index = -1;
	double min = numeric_limits<double>::max();
	for (int j = 0; j < kmax; j++) {
	if (org->active[j]) {
	double temp_dist = dist(org->clusCenter[j], attr[i]->items);
	if (temp_dist < min) {
	min = temp_dist;
	min_index = j;
	}
	}
	}
	assert(min_index != -1);
	clusters[min_index]->push_back(i);
	}
	for (int k = 0; k < kmax; k++) {
	if (org->active[k]) {
	outputFile << "------------------------------------------------------------" << endl;
	outputFile << "Elements of cluster : " << k << endl;
	//			for (std::vector<int>::const_iterator j = clusters[k]->begin();j != clusters[k]->end(); ++j) {
	for (vector<int>::size_type j = 0; j != clusters[k]->size(); j++){
	int a = clusters[k]->at(j);
	arr = attr[a]->items;				
	for (int m = 0; m < dim; m++) {
	outputFile << arr[m] << " ";
	}
	outputFile << attr[a]->typeClass;
	outputFile << endl;
	}
	}
	}*/
  outputFile << "Total number of clusters obtained : " << activeCount
	     << endl;
  outputFile << "Min DB is " << minDB << " Max DB is " << maxDB << endl;
  outputFile.close();
  cout << "Result saved in file.";
  delete [] arr;
}
