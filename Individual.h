/*
 * Individual.h
 *
 *  Created on: Jul 7, 2015
 *      Author: Navs
 */

#ifndef INDIVIDUAL_H_
#define INDIVIDUAL_H_

#define MAX_ROW 20
#define DIM 20
/*class Item{
public:
	Item(int ftr){
		this->feature = ftr;
	}
	int feature;
	double items[feature];
};*/
class Individual {
public:
	Individual();
	virtual ~Individual();
	bool isValid();
	bool getValid();
	void setValid(bool valid);
	double getFitness();
	void setFitness(double rawFitness);
private:
	double rawFitness;
	bool valid;
	double clusCenter[MAX_ROW][DIM];
	double threshold[MAX_ROW];
};

#endif /* INDIVIDUAL_H_ */
