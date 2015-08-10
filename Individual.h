/*
 * Individual.h
 *
 *  Created on: Jul 7, 2015
 *      Author: Navs
 */

#ifndef INDIVIDUAL_H_
#define INDIVIDUAL_H_
#include <vector>
#include <iostream>

using namespace std;
class Individual {
public:
	Individual(int kmax, int dim);
	~Individual();
	Individual(const Individual& obj);
	double rawFitness;
	bool valid;
	double** clusCenter;
	double* threshold;
	bool* active;
	int active_ctr;
	int k;
	int dimension;
	//vector<int>** clusters;

	bool isValid();
	bool getValid();
	void setValid(bool valid);
	double getFitness();
	void setFitness(double rawFitness);
	//bool operator<=(const Individual& right);
	friend ostream& operator<<(ostream& o, const Individual& org);

};

#endif /* INDIVIDUAL_H_ */
