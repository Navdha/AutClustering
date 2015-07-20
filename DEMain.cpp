/*
 * DE.cpp
 *
 *  Created on: Jul 8, 2015
 *      Author: Navs
 */

#include "DEMain.h"
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

struct compare_first_only {
    template<typename T1, typename T2>
    bool operator()(const std::pair<T1, T2>& p1, const std::pair<T1, T2>& p2) {
        return p1.first < p2.first;
    }
};

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

int* offspring_arr;

DEMain::DEMain(int kmax, int dim, int gen, int** placeholder, Item** items,
		int itemSize) {
	// TODO Auto-generated constructor stub
	cout << "Constructor called from DEMain class" << endl;
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
}

DEMain::~DEMain() {
	// TODO Auto-generated destructor stub
	delete p;
	for (int i = 0; i < numItems; i++)
		delete tracker[i];

	delete[] tracker;
	for (int i = 0; i < numItems; i++)
		delete attr[i];

	delete[] attr;
}

/* This method takes an input two arrays min and max
 * that hold the minimum and maximum value of each attribute
 * Purpose of this method is to setup the initial population
 * return type is void
 * no output expected
 */
void DEMain::setup(double min[], double max[]) {
	//initialize chromosomes for the first time
	cout<< "Setup method called" <<endl;
	for (int i = 0; i < pSize; i++) {
		Individual* temp = new Individual(kmax, dim);
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
			//	cout << min[k] << " " << max[k] << endl;
				temp->clusCenter[j][k] = uniformInRange(min[k], max[k]);
				//cout << temp->clusCenter[j][k] << " ";
			}
		}
		temp->active_ctr = ctr_act;
		cout << endl;
		double fitn = calcFitness(temp, i, true);
		temp->setFitness(fitn);
		//delete p->chromosome[i];
		p->chromosome[i] = temp;
		//delete temp;

	}


}


/*
 * Input parameters: pointers to arrays that hold
 * cluster center and item features
 * This method finds the distance between cluster center and an item
 * returns the distance calculated
 */
double DEMain::dist(double* x, double* y) {
	double Sum;
	double distance;
	for (int i = 0; i < dim; i++) {
		Sum = Sum + pow((x[i] - y[i]), 2.0);
		distance = sqrt(Sum);
	}
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
void DEMain :: reshuffle(Individual* org, Dist_IC *nearestDist, int size, int ind){//need to think
	cout << "reshuffle method called" <<endl;

	for(int i = 0; i < kmax; ++i) {
		  delete org->clusters[i];
		}
	delete [] org->clusters;
	int fix_size = numItems/kmax;
	org->clusters = new vector<int>*[kmax];
	for (int count = 0; count < kmax; count++)
		{
			org->clusters[count] = new vector<int>;
		}
	int ctr = 0;

	while (ctr < size) { //form clusters

		int i = 0;
		//sort array
		//for(int i = 0; i < kmax; i++){
			//if (org->active[i]) {
		//std::stable_sort(nearestDist[ctr].begin(), nearestDist[ctr].end(), compare_first_only());

		while(i < kmax) {
			int index = nearestDist[ctr].clustIndex;
			if (org->active[index] && org->clusters[index]->size() < fix_size) {
				org->clusters[index]->push_back(nearestDist[ctr].itemIndex);
				tracker[nearestDist[ctr].itemIndex][ind] = index;
				break;
			}
			else {
				i++;
			}
		}
			//}
		//}


		//tracker[ctr][index] = min_index;
		ctr++;
	}
	//delete org->clusters;
}

/*
 * Input parameters : Pointer to chromosome, index of population's chromosome,
 * bool parameter to know whether function called during initial setup or during DE
 * This method calculates the fitness of a chromosome and returns it
 * return  type : double
 */
double DEMain::calcFitness(Individual* org, int index, bool isInitial) {// not using index right now
	cout << "calcFitness method called" << endl;
	double fit = 0.0;
	double maxValue = 0.0;
	double sum = 0.0;
	double eps = 0.5;
	//double* avgArr;
	int str_size = numItems * (org->active_ctr);
	cout << str_size << endl;
	Dist_IC* knn = new Dist_IC [str_size];
	//std::vector<std::pair<double, int> > nearest_dist;
	//std::vector<std::vector< std::pair<double, int>>> all_nearest_dist;
	/* for each item read from csv, find distance of each item from the active cluster center
	 * using Euclidean distance. The least distant cluster center would be selected.
	 * As such form k clusters. check if each cluster is valid. if not, reshuffle.
	 * Finally, find fitness of each cluster using DB index.
	 */
	int ctr = 0;
	int vals = 0;


	//nearest_dist = new std::vector<std::pair<double, int> >*[numItems];
	while (ctr < numItems && vals < str_size) { //form clusters
		int min_index = -1;
		double min = numeric_limits<double>::max();
		//all_nearest_dist.push_back(nearest_dist);
		for (int i = 0; i < kmax; i++) {
			//nearest_dist = new std::vector<std::pair<double, int> >;
			if (org->active[i]) {
				double temp_dist = dist(org->clusCenter[i], attr[ctr]->items);
				knn[vals].distance = temp_dist;
				knn[vals].clustIndex = i;
				knn[vals].itemIndex = ctr;
				//cout << temp_dist << " " ;
				//all_nearest_dist[ctr].push_back(std::make_pair(temp_dist, i));
				if (temp_dist < min) {
					min = temp_dist;
					min_index = i;
				}
				vals++;
			}
			else{
				//all_nearest_dist[i].push_back(std::make_pair(-1.0,-1));
				//knn[vals].distance = -1.0;
				//knn[vals].clustIndex = -1;
				//knn[vals].itemIndex = -1;
			}

		}
		assert(min_index != -1);
		org->clusters[min_index]->push_back(ctr);
		if(isInitial) {
			tracker[ctr][index] = min_index;
		}
		else {
			offspring_arr = new int[numItems];
			offspring_arr[ctr] = min_index;
		}
		ctr++;

	}
	//check if each cluster is valid
	for (int i = 0; i < kmax; i++) {
		if (org->active[i]) {
			if (org->clusters[i]->size() < 2) {
				//org->setValid(false);
				//reshuffle items in clusters
				qsort(knn, str_size, sizeof(Dist_IC), compare);
			//	for (int n=0; n<str_size; n++)
				    // printf ("Sorted dist = %f item = %d  cluster_ind = %d \n",knn[n].distance, knn[n].itemIndex, knn[n].clustIndex);
				reshuffle(org, knn, str_size, index);
			}
		}
	}
	double avgArr[kmax];
	for (int i = 0; i < kmax; i++) {
		if (org->active[i]) {
			for (std::vector<int>::const_iterator j = org->clusters[i]->begin(); j != org->clusters[i]->end(); ++j) {
				sum += dist(attr[*j]->items, org->clusCenter[i]);
			}
			avgArr[i] = sum / (org->clusters[i]->size());
		}
	}
	//avgArr = avgDist(org);
	cout << "BAck to calcFitness" << endl;
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
	if (kmax == 1) {
		cout << "cluster number is 1, the value will be 0" << endl;
	}
	double avg = sum / org->active_ctr;
	fit = 1 / (avg + eps);
	//delete org->clusters;
	delete [] knn;
	return fit;
}

/*
 * Input parameters : indexes for chromosomes
 * This method makes sure that we get unique indexes for performing crossover
 * return type: void
 */

void DEMain::selectSamples(int org, int *s1, int *s2, int *s3) {
	if (s1) {
		do {
			*s1 = uniformInRange(0, pSize);
		} while (*s1 == org);
	}

	if (s2) {
		do {
			*s2 = uniformInRange(0, pSize);
		} while ((*s2 == org) || (*s2 == *s1));
	}

	if (s3) {
		do {
			*s3 = uniformInRange(0, pSize);
		} while ((*s3 == org) || (*s3 == *s2) || (*s3 == *s1));
	}
	cout << "selectSamples called" << endl;
	return;
}

/*
 * Input parameters: index for chromosome chosen and the generation number
 * This method performs crossover to create an offspring
 * returns pointer to offspring created
 */
Individual* DEMain::crossover(int org, int gen) {
	cout << "crossover method called" << endl;
	int s1, s2, s3;
	double cr_prob = probability * ((generations - gen) / generations);
	double f_scale = scale * uniform01();
	//cout << cr_prob << " " << f_scale << " prob & scale factor" << endl;
	selectSamples(org, &s1, &s2, &s3);
	//cout << "s1, s2, s3 " << s1 << " " << s2 << " " << s3 << endl;
	Individual* child = new Individual(kmax, dim);
	int counter = 0;
	for (int j = 0; j < kmax; j++) {
		child->threshold[j] = uniform01();
		if (child->threshold[j] > 0.5){
			child->active[j] = true;
			counter++;
		}
		else
			child->active[j] = false;

		for (int i = 0; i < dim; i++) {
//			if(child->active[j])
			if (uniform01() < cr_prob) {
				child->clusCenter[j][i] =
						p->chromosome[s1]->clusCenter[j][i]	+ f_scale*(abs(p->chromosome[s2]->clusCenter[j][i]- p->chromosome[s3]->clusCenter[j][i]));
				//cout << p->chromosome[s1]->clusCenter[j][i] << " " << p->chromosome[s2]->clusCenter[j][i] <<  " " << p->chromosome[s3]->clusCenter[j][i] << " " << child->clusCenter[j][i] << endl;
			} else
				child->clusCenter[j][i] = p->chromosome[org]->clusCenter[j][i];
		}
	}
	child->active_ctr = counter;
	return child;

}

/*
 * This method runs the DE algorithm
 */
void DEMain::run() {
	cout << "run method called" << endl;
	int i = 0;
	Population* newpop = new Population(kmax, dim);
	bool * new_pop = new bool[pSize];
	try {
	while (i < generations) {
		for (int c = 0; c < pSize; c++) {
			cout << c << " iteration" << endl;
			Individual *offspring;
			offspring = crossover(c, i);
			double fitness = calcFitness(offspring, c, false);
			offspring->setFitness(fitness);
			if (p->chromosome[c]->rawFitness <= offspring->rawFitness) {
				//cout << "Good offspring" << endl;
				new_pop[c] = true;
				cout << "offspring added" << endl;
				newpop->chromosome[c] = offspring;
				//delete p->chromosome[c];
				for (int d = 0; d < numItems; d++) {
					tracker[d][c] = offspring_arr[d];
				}
			} else {
				new_pop[c] = false;
				//cout << "Bad offspring" << endl;
				cout << "offspring discarded" << endl;
				delete offspring;
				newpop->chromosome[c] = p->chromosome[c];
				delete [] offspring_arr;
			}
		}
		cout << "Generation " << i << " completed" << endl;
		//assert(newpop != NULL);
		for (int c = 0; c < pSize; c++) {
			if (new_pop[c]) {
				delete p->chromosome[c];
			} else {

			}
		}
		delete p;
		p = newpop;
		i++;
	}
	cout << endl;
	cout << "Stopped at generation " << i << endl;
	//find best chromosome
	int bestInd = -1;
	double fitness = p->chromosome[0]->rawFitness;
	for (int k = 1; k < pSize; k++) {
		if (fitness < p->chromosome[k]->rawFitness) {
			bestInd = k + 1;
			fitness = p->chromosome[k]->rawFitness;
		}
	}
	assert(bestInd != -1);
	report(bestInd);
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
	outputFile.open("data.txt");
	outputFile << "The final clusters obtained are:" << endl;
	Individual* org = p->chromosome[index];
	vector<int>* tempClust;
	double* arr;
	outputFile << "Total number of clusters obtained : " << org->active_ctr
			<< endl;
	for (int k = 0; k < kmax; k++) {
		if (org->active[k]) {
			delete tempClust;
			tempClust = org->clusters[k];
			outputFile << "Elements of cluster : " << endl;
			for (std::vector<int>::const_iterator j = tempClust->begin();
					j != tempClust->end(); ++j) {
				arr = attr[*j]->items;
				for (int m = 0; m < dim; m++) {
					outputFile << arr[m] << " ";
				}
				outputFile << endl;
			}
		}
	}
	outputFile.close();
	cout << "Result saved in file.";
}
