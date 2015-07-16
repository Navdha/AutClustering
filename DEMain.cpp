/*
 * DE.cpp
 *
 *  Created on: Jul 8, 2015
 *      Author: Navs
 */

#include "DEMain.h"
#include <iostream>
#include <cmath>
#include <limits>
#include <cassert>
#include <random>     // needed for the mersenne twister random number generator
#include <functional>  // needed for bind

using namespace std;

random_device seed;
mt19937 mersenne(seed());  //Mersenne Twister random number generator
int uniformInRange(int from, int to)
{// generate a random uniformly in the range [from, to]
  uniform_int_distribution<int> uniformFT(from, to);
  return uniformFT(mersenne);
}
double uniformInRange(double from, double to)
{// generate a random uniformly in the range [from, to]
  uniform_real_distribution<double> uniformFT(from, to);
  return uniformFT(mersenne);
}
//uniform01() generates a random double uniformly in the range [0,1)
auto uniform01 = bind(uniform_real_distribution<double>(0,1), mt19937(seed()));

DEMain::DEMain(int kmax, int dim, int gen, int** placeholder, const Item** items, int itemSize) {
	// TODO Auto-generated constructor stub
	p = new Population(kmax, dim);
	strategy = stRand1Bin;
	generations = gen;
	probability = 0.5;
	scale = 0.5;
	pSize = dim*10;
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

	delete [] tracker;
	for (int i = 0; i < numItems; i++)
		 delete attr[i];

	delete [] attr;
}

void DEMain :: setup(double *min,double *max,int deStrategy,double diffScale,double crossoverProb){
	//initialize chromosomes for the first time
	for(int i = 0; i < pSize; i++){
		Individual* temp;
		int ctr_act = 0;
		temp = new Individual(kmax, dim);
		for(int j =0; j < kmax; j++){
			temp->threshold[j] = uniform01();
			if(temp->threshold[j] > 0.5) {
				temp->active[j] = true;
				ctr_act++;
			}
			else
				temp->active[j] = false;
			temp->active_ctr = ctr_act;
			for(int k = 0; k < dim; k++) {
				temp->clusCenter[j][k] =uniformInRange(min[k],max[k]);
			}
		}
		double fitn = calcFitness(temp, i);
		temp->setFitness(fitn);
		//delete p->chromosome[i];
		p->chromosome[i] = temp;
	}
}

double DEMain:: dist(double* x, double* y)
{
	double Sum;
	double distance;
	for(int i=0;i<dim;i++)
	{
		Sum = Sum + pow((x[i]-y[i]),2.0);
		distance = sqrt(Sum);
	}
	return distance;
}

double* DEMain::avgDist(Individual* org){
	double* temp = new double[kmax];
	//double *tempArr = &temp;
	double* d2;
	//int ind = -1;
	for(int i=0;i<kmax;i++){
		double sum = 0.0;
		delete d2;
		d2 = org->clusCenter[i];
		vector<int>* c = org->clusters[i];
		//for(vector<int>::iterator it = c[i]->begin(); it != c[i]->end(); ++it) {
		for(vector<int>::size_type j = 0; j != c->size(); j++) {
			sum += dist(attr[j]->items, d2);
		}
		temp[i] = sum/(c->size());
	}

	return temp;
}

double DEMain::calcFitness(Individual* org, int index){
	double fit = 0.0;
	double maxValue=0.0;
	double sum=0.0;
	double eps = 0.5;
	double* avgArr;
	/* for each item read from csv, find distance of each item from the active cluster center
	 * using Euclidean distance. The least distant cluster center would be selected.
	 * As such form k clusters. check if each cluster is valid. if not, reshuffle.
	 * Finally, find fitness of each cluster using DB index.
	 */
	int ctr = 0;
	double min = numeric_limits<double>::max();
	int min_index = -1;
	while(ctr < numItems){ //form clusters
	for(int i = 0; i < kmax; i++){
		if(org->active[i]){
			double temp_dist = dist(org->clusCenter[i], attr[ctr]->items);
			if(temp_dist < min){
				min = temp_dist;
				min_index = i;
			}
		}
	}
	assert(min_index != -1);
	org->clusters[min_index]->push_back(ctr);
	tracker[ctr][index] = min_index;
	ctr++;
	}
	//check if each cluster is valid
	for(int i = 0; i < kmax; i++){
		if(org->active[i]){
			if(org->clusters[i]->size() < 2){
				//org->setValid(false);
				//reshuffle items in clusters
			}
		}
	}
	avgArr = avgDist(org);
	for (int i = 0; i < kmax; i++) {
		for (int j = i + 1; j < kmax; j++) {
			double temp = avgArr[i] + avgArr[j];
			temp /= dist(org->clusCenter[i], org->clusCenter[j]);
			if (temp > maxValue)
				maxValue = temp;
		}
		sum += maxValue;
	}
	if(kmax == 1) {cout << "cluster number is 1, the value will be 0" << endl;}
	double avg = sum/kmax;
	fit = 1/(avg + eps);
	return fit;
}

void DEMain::selectSamples(int org, int *s1, int *s2, int *s3){
	if (s1)
		{
			do
			{
				*s1 = uniformInRange(0, pSize);
			}
			while (*s1 == org);
		}

		if (s2)
		{
			do
			{
				*s2 = uniformInRange(0, pSize);
			}
			while ((*s2 == org) || (*s2 == *s1));
		}

		if (s3)
		{
			do
			{
				*s3 = uniformInRange(0, pSize);
			}
			while ((*s3 == org) || (*s3 == *s2) || (*s3 == *s1));
		}
		return;
}

Individual* DEMain::crossover(int org, int gen){
	int s1, s2, s3;
	double cr_prob = probability*((generations - gen)/generations);
	double f_scale = scale*uniform01();
	selectSamples(org, &s1, &s2, &s3);
	Individual* child = new Individual(kmax, dim);
	for(int j=0;j<kmax;j++){
		child->threshold[j] = uniform01();
		if(child->threshold[j] > 0.5)
			child->active[j] = true;
		else
			child->active[j] = false;
		for(int i=0;i<dim;i++){
			if(uniform01() < cr_prob){
				child->clusCenter[j][i] = p->chromosome[s1]->clusCenter[j][i] + f_scale*(abs(p->chromosome[s2]->clusCenter[j][i] - p->chromosome[s3]->clusCenter[j][i]));
			}
			else
				child->clusCenter[j][i] = p->chromosome[org]->clusCenter[j][i];
		}
	}
	return child;

}

void DEMain::run(){
	long i = 0;
	Individual* offspring;
	Population* newpop;
	bool new_pop[pSize]; //is this really saving time? since we have two for loops
	while (i < generations){
		for(int c= 0; c < pSize; c++){
			offspring = crossover(c, i);
			double fitness = calcFitness(offspring, c);
			offspring->setFitness(fitness);
			if(p->chromosome[c]->rawFitness <= offspring->rawFitness){
				new_pop[c] = true;
			}
			else
				new_pop[c] = false;
		}
		newpop = new Population(kmax, dim);
		assert(newpop != NULL);
		for(int c = 0; c < pSize; c++){
			if(new_pop[c]){
				newpop->chromosome[c] = offspring;
			}
			else{
				newpop->chromosome[c] = p->chromosome[c];
			}
		}
		delete p;
		p = newpop;

	}
}
