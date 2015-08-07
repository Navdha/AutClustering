/*
 * Parameters.cpp
 *
 *  Created on: Aug 6, 2015
 *      Author: Navdha Sah
 */

#include "Parameters.h"

Parameters::Parameters(double CrMaximum, double CrMinimum, double FScaleProb,
		double threshVal, int kmaxVal, int kminVal, long numGen) {
	// TODO Auto-generated constructor stub
	kmax = kmaxVal;
	kmin = kminVal;
	threshold = threshVal;
	CrMax = CrMaximum;
	CrMin = CrMinimum;
	FScale = FScaleProb;
	gen = numGen;
}

Parameters::~Parameters() {
	// TODO Auto-generated destructor stub
}
