/*
 * DE.cpp
 *
 *  Created on: Jul 8, 2015
 *      Author: Navs
 */

#include "DEMain.h"
#include "Parameters.h"
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
#include <iterator>
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

struct compare_first_only {
  template<typename T1, typename T2>
  bool operator()(const std::pair<T1, T2>& p1, const std::pair<T1, T2>& p2) {
    return p1.first < p2.first;
  }
};

ofstream trackFile; //used to record the clusters size etc
//int changeCounter =0;
double minDB = std::numeric_limits<double>::max();
double maxDB= std::numeric_limits<double>::min();
double minCS = std::numeric_limits<double>::max();
double maxCS= std::numeric_limits<double>::min();
double minPB = std::numeric_limits<double>::max();
double maxPB= std::numeric_limits<double>::min();

int compare(const void *p1, const void *p2){ // to compare two triples containing the item id, cluster id and distance between item and cluster
  Dist_IC *elem1 = (Dist_IC *)p1;
  Dist_IC *elem2 = (Dist_IC *)p2;
  if(elem1->distance < elem2->distance)
    return -1;
  else if(elem1->distance > elem2->distance)
    return 1;
  else
    return 0;
}


DEMain::DEMain(int dim, int** placeholder, Item** items, int itemSize, int validityIndex, Parameters param) {
  // TODO Auto-generated constructor stub
  //cout << "Constructor called from DEMain class" << endl;
  p = new Population(param.kmax, dim);
  strategy = stRand1Bin; //DE algo used
  generations = param.gen; //number of generations
  probability = param.CrMax - param.CrMin; //crossover probability
  scale = param.FScale; //factoring value
  pSize = dim * 10; //population size
  kmax = param.kmax; //maximum # of clusters
  this->dim = dim; //number of features
  numItems = itemSize; //total # of items to cluster
  kmin = param.kmin; // minimum # of clusters
  thresholdVal = param.threshold; //threshold value to activate a cluster's centroid
  indexForFit = validityIndex; //which clustering validity index used
  attr = items; //data structure to hold all items

  //scratch space
  tracker = placeholder; //holds cluster indexes corresponding to every item and pop member(calcFitness and report)
  clusters = new vector<int>*[kmax]; //holds actual clusters for every pop members(calcFitness)
  distItem = new double*[numItems]; //jagged array storing dist between all items(used for calcCSIndex and calcPBIndex)
  for (int i = 0; i < numItems -1; i++) {
    distItem[i] = new double[numItems- 1 - i];
  } 
  offspring_arr = new int[numItems];//column that replaces worse parent's column from tracker(calcFitness)
  int size = kmax*numItems;
  knn = new Dist_IC[size];//used for reshuffling(reshuffle & calcFitness)
  ItemUsed = new bool[numItems](); //used in reshuffle
  ClusFull = new bool[kmax](); //used in reshuffle
  avgArr = new double[kmax];//used in calcDB
  sumArr = new double[dim];// used in calcCS
  ItemCounter = new int[numItems]; //used in reshuffle
  newClustCenters = new double*[kmax];//used in calcCS
  new_pop = new bool[pSize];
  for (int count = 0; count < kmax; count++)
      {
        clusters[count] = new vector<int>;
        newClustCenters[count] = new double[dim];
      }
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
    delete newClustCenters[i];
  }
  delete [] tracker;
  delete [] attr;
  delete [] distItem;
  delete [] clusters;
  delete [] offspring_arr;
  delete [] knn;
  delete [] ItemUsed;
  delete [] ClusFull;
  delete [] avgArr;
  delete [] sumArr;
  delete [] ItemCounter;
  delete [] new_pop;
  delete [] newClustCenters;
}

/*
 * This method calculates and stores the distance between all items
 * Used while calculating CS index and standard deviation between distances for Point Biserial index
 */
void DEMain::calcDistBtwnItems() {
	for (int i = 0; i < numItems - 1; i++) {
		for (int j = i + 1; j < numItems; j++) {
			distItem[i][j - i - 1] = dist(attr[i]->items, attr[j]->items);
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
	// cout<< "Setup method called" <<endl;
	for (int i = 0; i < pSize; i++) {
		bool isValid = false; //determines whether our individual is valid or not based on how many clusters it forms
		double fitn = 0.0;
		Individual* temp = new Individual(kmax, dim); //instantiating new individual
		while (!isValid) {//continue this till a valid individual found
			int ctr_act = 0; //# of active cluster centers for the individual
			for (int j = 0; j < kmax; j++) {
				temp->threshold[j] = uniform01();
				if (temp->threshold[j] > thresholdVal) { // based on threshold, making a centroid active or inactive
					temp->active[j] = true;
					ctr_act++;
				} else
					temp->active[j] = false;
				//				 int randomVal = uniformInRange(0, numItems-1);
				for (int k = 0; k < dim; k++) {
				  // temp->clusCenter[j][k] = attr[randomVal]->items[k];
				   temp->clusCenter[j][k] = uniformInRange(min[k], max[k]); //randomly creating centroid to be a cluster center
				}
			}
			assert(ctr_act <= kmax); //check to determine if active clusters are more than kmax
			temp->active_ctr = ctr_act; //assigning this as active cluster count in the individual object
			//code to check kmin = 2
			if (temp->active_ctr < kmin) {//if #of active clusters is less than min #of clusters possible, then activate a cluster centroid forcefully
				int num = temp->active_ctr;
				while (num < kmin) {
					int i = uniformInRange(0, kmax - 1);
					if (!temp->active[i]) {
						temp->threshold[i] = uniformInRange(thresholdVal, 1.0);
						temp->active[i] = true;
						temp->active_ctr++;
						num++;
					}
				}
			}
			fitn = calcFitness(temp, i, true, -1); //calculate fitness of individual generated for initial population
			if (fitn != -1) { //returned -1 in case of active clusters falling below kmin
				isValid = true;
			}
		} //end while
		assert(fitn != 0.0);
		temp->setFitness(fitn); //assign fitness
		p->chromosome[i] = temp;
	}
	int bestInd = 0;
	double fitness = p->chromosome[0]->rawFitness;
	for (int k = 1; k < pSize; k++) { //finding the best individual in the population
		if (fitness < p->chromosome[k]->rawFitness) {
			bestInd = k;
			fitness = p->chromosome[k]->rawFitness;
		}
	}
	p->bestChromosomeIndex = bestInd;//updating the index of the best individual found

}


/*
 * Input parameters: pointers to arrays that hold
 * cluster center and item features
 * This method finds the distance between cluster center and an item
 * returns the distance calculated
 */
double DEMain::dist(double* x, double* y) {
  double Sum = 0.0;
  for (int i = 0; i < dim; i++) {
    Sum = Sum + pow((x[i] - y[i]), 2.0);
  }
  return sqrt(Sum);
}

/*
 * Input parameters : A pointer to a chromosome, struct array that holds distance corresponding
 * to item and cluster center, size of the struct array and index of population's chromosome
 * This method reshuffles items equally into different active cluster centers of an individual
 * return type : void
 */
void DEMain::reshuffle(Individual* org, int size, int indpop, bool isInitial) { //need to think
//  cout << "reshuffle method called" <<endl;
	for (int i = 0; i < kmax; ++i) {
		clusters[i]->clear();
	}
	int fix_size = numItems / org->active_ctr; //maximum # of items that a cluster can hold
	if (fix_size < 2) {
		cout << numItems << " " << org->active_ctr;
	}
	assert(fix_size >= 2);//to assert that each cluster has at least 2 items
	int ctr = 0; //counter that goes through all elements in the triples array
	int numFullClusters = 0;//determines whether all clusters have reached their max capacity

	for (int i = 0; i < numItems; i++) //initializing counter array corresponding to how many times an item has appeared in the triples array
		ItemCounter[i] = 0;
	fill_n(ItemUsed, numItems, false); //initializing bool array such that no item has been added to cluster yet
	fill_n(ClusFull, kmax, false); // initializing bool array such that no cluster is full yet
	while (ctr < size) { //check for total # of items
		int itemInd = knn[ctr].itemIndex; //retrieving index of item
		int clusInd = knn[ctr].clustIndex; // retrieving index of cluster (s.t we have the closest item and cluster present
		if (!ItemUsed[itemInd]) {//only if the item has not already been placed in a cluster
			if (org->active[clusInd]) {//if retrieved cluster is valid still (could be deactivated if the cluster had 0/1 item allocated while computing fitness
				ItemCounter[itemInd]++;
				if (numFullClusters != org->active_ctr) {//if all clusters aren't full yet
					if (clusters[clusInd]->size() < fix_size) {//if the selected cluster isn't full yet
						clusters[clusInd]->push_back(itemInd);//add item to cluster
						ItemUsed[itemInd] = true;//mark that item as already added to a cluster
						if (isInitial) {//updating the arrays to trace back later
							tracker[itemInd][indpop] = clusInd;
						} else {
							offspring_arr[itemInd] = clusInd;
						}
					} else {//if the cluster has reached its max capacity
						if (ItemCounter[itemInd] == org->active_ctr) {//check if a particular item wasn't added to any cluster due to the cluster being full
							clusters[clusInd]->push_back(itemInd);//added to the cluster its farthest from; otherwise the item wouldn't be added at all (special use case)
							ItemUsed[itemInd] = true;
						}
						if (!ClusFull[clusInd]) {//since the cluster has reached its max capacity, mark it as full
							ClusFull[clusInd] = true;
							numFullClusters++;//increment the counter that keeps track of #of clusters full
						}
					}
				} else {//if all clusters are full but still item remains push it in anyway
					clusters[clusInd]->push_back(itemInd);
					ItemUsed[itemInd] = true;
					if (isInitial) {//update the arrays to trace back
						tracker[itemInd][indpop] = clusInd;
					} else {
						offspring_arr[itemInd] = clusInd;
					}
				}
			}
		}
		ctr++;//go through next item in triples array
	} //end while

}

/*
 * calculate and return DB index
 * Parameters : Individual
 */
double DEMain::calcDBIndex(Individual* org) {
	fill_n(avgArr, kmax, 0.0);
	double sum;
	for (int i = 0; i < kmax; i++) {
		for (int m = 0; m < dim; m++) {
			newClustCenters[i][m] = 0.0; //clearing old centroids
		}
	}
	//this section of code is to compute centroids of the clusters to compute DB index
	for (int i = 0; i < kmax; i++) {
		fill_n(sumArr, dim, 0);//clearing sums stored in array
		if (org->active[i]) {
			for (vector<int>::size_type j = 0; j != clusters[i]->size(); j++) {//go through all items in cluster
				int a = clusters[i]->at(j);//find item index stored in clusters
				double* tempItem = attr[a]->items;//array to hold features corresponding to item in cluster
				for (int k = 0; k < dim; k++) {
					sumArr[k] += tempItem[k];//this holds the sum of all features for all items in a single cluster to compute average later
				}
			}
			for (int m = 0; m < dim; m++) {
				newClustCenters[i][m] = sumArr[m] / clusters[i]->size(); //finding new centroids for the cluster
			}
		}
	}
	for (int i = 0; i < kmax; i++) {
		if (org->active[i]) {
			sum = 0.0;
			for (vector<int>::size_type j = 0; j != clusters[i]->size(); j++) {
				int a = clusters[i]->at(j);
				//sum += dist(attr[a]->items, org->clusCenter[i]);
				 sum += dist(attr[a]->items, newClustCenters[i]);
			} //end for
			avgArr[i] = sqrt(sum / clusters[i]->size()); //finding the intra cluster distance for all active clusters
		} //end if

	} //end outer for
	sum = 0.0;
	for (int i = 0; i < kmax; i++) {
		double maxValue = 0.0;
		for (int j = 0; j < kmax; j++) {
			if (i != j && org->active[i] && org->active[j]) {
				double temp = avgArr[i] + avgArr[j];
				//temp /= dist(org->clusCenter[i], org->clusCenter[j]); //finding R =(S_i+S_j)/d_i,j
				temp /= dist(newClustCenters[i], newClustCenters[j]);
				if (temp > maxValue)
					maxValue = temp;
			}
		}
		sum += maxValue;		  //finding sum(Rmax)
	}
	double avg = sum / org->active_ctr;
	if (minDB > avg) {//keeping track of the minimum and maximum DB index encountered for all individuals in all generations
		minDB = avg;
	}
	if (maxDB < avg) {
		maxDB = avg;
	}
	return avg;
}

/*
 * calculate and return CS index for an organism
 * Parameters : Individual
 */
double DEMain::calcCSIndex(Individual* org){
	for (int i = 0; i < kmax; i++) {
		for (int m = 0; m < dim; m++) {
			newClustCenters[i][m] = 0.0; //clearing old centroids
		}
	}
	  double finalIntraSum = 0.0;
	  double finalInterSum = 0.0;
	  double csVal = 0.0;
	  for(int i = 0; i < kmax; i++){
		double intraClusSum = 0.0;
		double tempIntraDist;
		fill_n(sumArr, dim, 0);
		if(org->active[i]){
		   for (vector<int>::size_type j = 0; j != clusters[i]->size(); j++) {
			int a = clusters[i]->at(j);
			double maxIntraDist = numeric_limits<double>::min();
			double* tempItem = attr[a]->items;
			for(int k = 0; k < dim; k++){
			    sumArr[k] += tempItem[k]; //to compute centroids
			}
			for(vector<int>::size_type k = 0; k != clusters[i]->size(); k++) {//finding max distance between items in a cluster
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
			intraClusSum += maxIntraDist; //assigning intra cluster distance for this cluster
		}// end for j
		finalIntraSum += (intraClusSum / clusters[i]->size()); //updating the average intra cluster sum for all active centroids
		for(int m = 0; m < dim; m++){
		    newClustCenters[i][m] = sumArr[m]/clusters[i]->size(); //finding new centroids
		}
	  }//endif
	 }// end for i
	  double interClusSum = 0.0;
	  double tempInterDist;
	  for (int i = 0; i < kmax; i++) {//trying to find min distance between two clusters
		double minInterDist = numeric_limits<double>::max();
		if(org->active[i]){
		for (int j = 0; j < kmax; j++) {
		     if (i != j && org->active[j]) {
			tempInterDist = dist(newClustCenters[i], newClustCenters[j]);
			if (tempInterDist < minInterDist){
			    minInterDist = tempInterDist;
			}
		     }//end outer if
		}//end for j
		interClusSum += minInterDist; //assigning minimum distance between clusters
		}
	 }//end for i
	   finalInterSum = interClusSum;
	   csVal = finalIntraSum/finalInterSum; //formula for CS index
	   if(minCS > csVal){//updating the maximum and minimum CS index calculated for all individuals in every generation
	      minCS = csVal;
	    }
	    if(maxCS < csVal){
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
	for (int i = 0; i < kmax; i++) {
		if(org->active[i]){
			double sumIntraCluster = 0.0;
			double sumInterCluster = 0.0;
			int n = clusters[i]->size();
			countIntraGroup += (n * (n - 1)) / 2;
			countInterGroup += (n * (numItems - n)) / 2;
			for (vector<int>::size_type j = 0; j != clusters[i]->size(); j++) {
				int a = clusters[i]->at(j);
				for (int l = 0; l < i; l++) {
					if (org->active[l]) {
						for (vector<int>::size_type ind = 0; ind != clusters[l]->size(); ind++) {
							int b = clusters[l]->at(ind);
							if (b < a) {
								sumInterCluster += distItem[b][a - (b + 1)];
							} else {
								sumInterCluster += distItem[a][b - (a + 1)];
							}
						}
					}
				}
				for (int l = i + 1; l < kmax; l++) {
					if (org->active[l]) {
						for (vector<int>::size_type ind = 0; ind != clusters[l]->size(); ind++) {
							int b = clusters[l]->at(ind);
							if (b < a) {
								sumInterCluster += distItem[b][a - (b + 1)];
							} else {
								sumInterCluster += distItem[a][b - (a + 1)];
							}
						}
					}
				}
				for (vector<int>::size_type k = 0; k != clusters[i]->size();
						k++) {
					if (k != j) {
						int b = clusters[i]->at(k);
						if (b < a) {
							sumIntraCluster += distItem[b][a - (b + 1)];
						} else {
							sumIntraCluster += distItem[a][b - (a + 1)];
						}

					}
				} //end for k
			}//end for j
			intraClusAvgDist += sumIntraCluster/countIntraGroup; //value of d_w
			interClusAvgDist += sumInterCluster/countInterGroup; //value of d_b
		}//end if
	}//end for i
	double totalSums = (interClusAvgDist - intraClusAvgDist);
	double t_val = t * t;
	double sqrtVal = sqrt((countIntraGroup * countInterGroup) / t_val);
	pbIndex = totalSums * sqrtVal / sd;
	//Point biserial : (d_b*d_w)(sqrt(w_d*b_d/t^2))/sd
	if (minPB > pbIndex) {//finding the min and max PB index found over all generations for all individuals
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
 * This method calculates the fitness of a chromosome and returns it
 * return  type : double
 */
double DEMain::calcFitness(Individual* org, int index, bool isInitial, int genNum) {// not using index right now
 // cout << "calcFitness method called" << endl;
	double fit = 0.0;
	int min_index = -1;
	double temp_dist;
	int itemctr = 0;
	int vals = 0;
	//clearing the scratch space
	for (int i = 0; i < kmax; ++i) {
		clusters[i]->clear();
	}

	//	int initialActive = org->active_ctr;
	while (itemctr < numItems) { //form clusters
		min_index = -1;
		double min = numeric_limits<double>::max();
		for (int i = 0; i < kmax; i++) {
			if (org->active[i]) {
				temp_dist = dist(org->clusCenter[i], attr[itemctr]->items);
				knn[vals].distance = temp_dist; //storing the triples set to knn array
				knn[vals].clustIndex = i;
				knn[vals].itemIndex = itemctr;
				if (temp_dist < min) { //finding the closest cluster center to a particular item
					min = temp_dist;
					min_index = i;
				}
				vals++;
			}

		}

		assert(min_index != -1);
		clusters[min_index]->push_back(itemctr); //adding item to its closest cluster center
		if (isInitial) {
			tracker[itemctr][index] = min_index; //mainitaining tracking array
		} else {
			offspring_arr[itemctr] = min_index;
		}
		itemctr++;

	} // end while

	//to print the initial cluster size before reshuffling/inactivating/moving clusters
	/* int tempArr[org->active_ctr];
	 int tempC = 0;
	 for (int i = 0; i < kmax; i++) {
	 if (org->active[i]) {
	 tempArr[tempC]= clusters[i]->size();
	 tempC++;
	 }
	 }*/	
	//inactivating empty clusters
	for (int i = 0; i < kmax; i++) {//inactivating a cluster that's empty
		if (org->active[i]) {
			if (clusters[i]->empty()) {
				org->active[i] = false;
				org->active_ctr--;
			}
		}

	}
	bool doReshuffle = uniform01() < thresholdVal ? false : true;
	for (int i = 0; i < kmax; i++) {
		if (org->active[i]) {
			if (clusters[i]->size() == 1) {
				//moving a single item in a cluster to its closest cluster center with a probability
				if (doReshuffle) {
					int a = clusters[i]->at(0);
					temp_dist = 0;
					double min = numeric_limits<double>::max();
					int ctr = 0;
					int minInd = -1;
					while (ctr < kmax) {
						if (ctr != i && org->active[ctr]) {
							temp_dist = dist(org->clusCenter[ctr],
									attr[a]->items);
							if (temp_dist < min) {
								temp_dist = min;
								minInd = ctr;
							}
						}
						ctr++;
					}
					assert(minInd != -1);
					clusters[minInd]->push_back(a);
					clusters[i]->clear(); // after element has been moved, cluster with single item made empty
					org->active[i] = false; // deactivating the empty cluster
					org->active_ctr--;
				} else {
					//reshuffle items in clusters with probability
					qsort(knn, vals, sizeof(Dist_IC), compare);
					reshuffle(org, vals, index, isInitial);
					break;
				} //random if-else end
			} //size if
		} //active if
	} //end for
	if (org->active_ctr >= kmin) {//if after deactivating cluster centers above, an individual's active cluster centers fall below 2(kmin), assign it lowest fitness
		//based on the validity index find DB, CS or PB index for the individual to compute the fitness
		if (indexForFit == 1) {
			double dBValue = calcDBIndex(org);
			fit = 1 / dBValue;
		} //end isDB if
		else if (indexForFit == 2) {
			//calculate CS index
			double csVal = calcCSIndex(org);
			fit = 1 / csVal;
		} //end else
		else {
			double pbVal = calcPBIndex(org);
			fit = pbVal;
		}
		/* trackFile.open("clusters.txt", ofstream::app);
		 trackFile << setiosflags(ios::left) ;
		 trackFile << setw(5) << genNum << setw(12) << fit*100 << setw(5) << org->active_ctr << setw(5) <<tempC ;
		 for(int i = 0; i < kmax; i++) {
		 if(org->active[i]){ trackFile << setw(5) << clusters[i]->size() ;}
		 }
		 trackFile << setw(5) << "|" ;
		 for(int i = 0; i <tempC; i++){trackFile << setw(5) <<  tempArr[i] ;}
		 trackFile << setw(5) << "|" << index << endl;
		 trackFile.close();*/
		return fit * 100;
	} else {
		cout << "Invalid chromosome" << endl;
		return -1;
	}

}

/*
 * Input parameters : indexes for chromosomes
 * This method makes sure that we get unique indexes for performing crossover
 * return type: void
 */

void DEMain::selectSamples(int org, int &s1, int &s2, int &s3) { // finding unique indexes to perform crossover
  if (s1) {
    do {
      s1 = uniformInRange(0, pSize-1);
    } while (s1 == org);
  }

  if (s2) {
    do {
      s2 = uniformInRange(0, pSize-1);
    } while ((s2 == org) || (s2 == s1));
  }

  if (s3) {
    do {
      s3 = uniformInRange(0, pSize-1);
    } while ((s3 == org) || (s3 == s2) || (s3 == s1));
  }
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
  double cr_prob = probability * ((generations - gen) / generations);//based on formula in paper
  //  double f_scale = scale * (1+uniform01()); //based on formula in paper
  double f_scale = uniformInRange(scale, 0.8);//temp change for testing
  selectSamples(org, s1, s2, s3);//to get unique individuals to perform DE crossover on
  Individual* child = new Individual(kmax, dim);
  int counter = 0;//keeps track of active cluster centers
  for (int j = 0; j < kmax; j++) {
    bool change = uniform01() < cr_prob? true : false;
    if(change) {//binomial crossover
      // if(uniform01() < cr_prob){
      child->threshold[j] =  p->chromosome[s1]->threshold[j] +
	f_scale*(p->chromosome[s2]->threshold[j] -  p->chromosome[s3]->threshold[j]);
    }
    else{
      child->threshold[j] =  p->chromosome[org]->threshold[j];
    }
    
    for (int i = 0; i < dim; i++) {//binomial crossover computation for every feature
      assert(p->chromosome[org] != NULL);
      if (change) {
	child->clusCenter[j][i] = p->chromosome[s1]->clusCenter[j][i]
	  + f_scale*(p->chromosome[s2]->clusCenter[j][i]- p->chromosome[s3]->clusCenter[j][i]);
      } else {
	child->clusCenter[j][i] = p->chromosome[org]->clusCenter[j][i];

      }

      if ((child->clusCenter[j][i] < min[i]) || (child->clusCenter[j][i] > max[i])){//if a feature for cluster center computed goes out of its bounds, assign it a random value
	child->clusCenter[j][i] = uniformInRange(min[i], max[i]);
      }
    }//for i
    if (child->threshold[j] > 1 || child->threshold[j] < 0)//if a threshold value after crossover goes out of its bounds, recompute it
      child->threshold[j] = uniform01();
    if (child->threshold[j] > thresholdVal) {
      child->active[j] = true;
      counter++;
    } else
      child->active[j] = false;
  }//for j
  assert(counter <= kmax);
  child->active_ctr = counter;
  if(child->active_ctr < kmin){//if #of active cluster centers below kmin, forcefully make more clusters active till kmin satisfied
    int num = child->active_ctr;
    while (num < kmin) {
      int i = uniformInRange(0, kmax - 1);
      if(!child->active[i]){
	child->threshold[i] = uniformInRange(thresholdVal, 1.0);
	child->active[i] = true;
	child->active_ctr++;
	num++;
      }
    }
  }

  return child;

}

/*
 * This method runs the DE algorithm
 */
void DEMain::run(double min[], double max[], string filename) {
 // cout << "run method called" << endl;
	int i = 0;
	fill_n(new_pop, pSize, false);
	double fitness;
	try {
	  trackFile.open("clusters.txt");
	  while (i < generations) {//till max generation reached
			Population* newpop = new Population(kmax, dim);
			for (int c = 0; c < pSize; c++) {
				Individual *offspring;
				offspring = crossover(c, i, min, max);//generate an offspring my performing DE crossover
				fitness = calcFitness(offspring, c, false, i);//compute its fitness
				offspring->setFitness(fitness);
				if (p->chromosome[c]->rawFitness <= offspring->rawFitness) {//if offspring better than parent, replace the parent
				        new_pop[c] = true;
					newpop->chromosome[c] = offspring;
					for (int d = 0; d < numItems; d++) {
						tracker[d][c] = offspring_arr[d]; //updating the parent chromosome replaced with new cluster centers of offspring
					}
				} else {//otherwise delete the offspring
					new_pop[c] = false;
					delete offspring;
					newpop->chromosome[c] = p->chromosome[c];
				}
			}
			for (int c = 0; c < pSize; c++) {//based on the parents that are supposed to be replaced, freeing the space occupied in memory
				if (new_pop[c]) {
					delete p->chromosome[c];
				}
			}
			delete[] p->chromosome;
			p = newpop;//new population generated becomes the current population
			int bestInd = 0;
			fitness = p->chromosome[0]->rawFitness;
			for (int k = 1; k < pSize; k++) { //keeping track of the best member of the population
			  if (fitness < p->chromosome[k]->rawFitness) {
					bestInd = k;
					fitness = p->chromosome[k]->rawFitness;
				}
			}
			p->bestChromosomeIndex = bestInd;
			for (int k = 0; k < pSize; k++) {//for testing purposes (keeping track of the fitness change as the algorithm progresses)
			  trackFile << setiosflags(ios::left) ;
			  trackFile << setw(5) << i << setw(12) << p->chromosome[k]->rawFitness << setw(5) <<  k << endl;
			}
			i++;
		}
		cout << endl;
		cout << "Stopped at generation " << i << endl;
		//find worst chromosome
		int worstInd = 0;
		double fitness = p->chromosome[0]->rawFitness;
		
		for (int k = 1; k < pSize; k++) {//for testing purposes
//finding the worst member in the population at the end
			if (fitness > p->chromosome[k]->rawFitness) {
				worstInd = k;
				fitness = p->chromosome[k]->rawFitness;
			}
		}
		trackFile.close();
		report(p->bestChromosomeIndex, worstInd, filename);
	} catch (exception& e) {
		cerr << e.what() << endl;
	}
}

/*
 * Input parameter : index for the best chromosome of the population
 * This method outputs the clusters for the best chromosome in the population to a file
 * return type : void
 */
void DEMain::report(int index, int worstInd, string filename) {
  ofstream outputFile;
  outputFile.open(filename.c_str());
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
    if (org->active[k]) {//prints out the best output
      activeCount++;
      outputFile << "-------------------------------------------" << endl;
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

  for(int i = 0; i < kmax; ++i) {
    clusters[i]->clear();
  }

  org = p->chromosome[worstInd];
  for(int i = 0; i < numItems; i++){
    clus_index = tracker[i][worstInd];
    if(org->active[clus_index]){
      clusters[clus_index]->push_back(i);
    }
  }
  outputFile << "-------------------------------------------" << endl;
  outputFile << "Worst fitness clusters are as follows" << endl;
  outputFile << endl;
  activeCount = 0;
  for (int k = 0; k < kmax; k++) {//prints out the worst output after changes are made
    if (org->active[k]) {
      activeCount++;
      outputFile << "-------------------------------------------" << endl;
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
      outputFile << endl;
    }
  }

  	for (int i = 0; i < kmax; ++i) {
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
	outputFile << "-------------------------------------------" << endl;
	outputFile << "Worst fitness cluster without tweaking" << endl;
	for (int k = 0; k < kmax; k++) {//prints out the worst output w/o making any change like reshuffling or deactivating clusters.
	  if (org->active[k] && !clusters[k]->empty()) {
	outputFile << "------------------------------------------------------------" << endl;
	outputFile << "Elements of cluster : " << k << endl;
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
  outputFile << "Total number of clusters obtained : " << activeCount
	     << endl;
  if(indexForFit == 1){
   outputFile << "Min DB index is " << minDB << " Max DB index is " << maxDB << endl;
  }
  else if(indexForFit == 2){
	  outputFile << "Min CS index is " << minCS << " Max CS index is " << maxCS << endl;
  }
  else{
	  outputFile << "Min PB index is " << minPB << " Max PB index is " << maxPB << endl;
  }
  outputFile.close();
  cout << "Result saved in file.";
  delete [] arr;
}
