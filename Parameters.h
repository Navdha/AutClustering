/*
 * Parameters.h
 *
 *  Created on: Aug 6, 2015
 *      Author: Navdha Sah
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

class Parameters {
public:
	double CrMax;
	double CrMin;
	double FScale;
	double threshold;
	int kmax;
	int kmin;
	long gen;
	Parameters(double CrMaximum, double CrMinimum, double FScaleProb, double threshVal, int kmaxVal, int kminVal, long numGen);
    ~Parameters();
};


#endif /* PARAMETERS_H_ */
