/*
 * Parameters.cpp
 *
 *  Created on: Aug 6, 2015
 *      Author: Navdha Sah
 */

#include "Parameters.h"

Parameters::Parameters(double CrMaximum, double CrMinimum, double FScaleProb,
		       double threshVal,  int maxClusters, int minClusters, int popScale, double numGen, int numC) {
  // TODO Auto-generated constructor stub
  maxNumClusters = maxClusters;
  minNumClusters = minClusters;
  threshold = threshVal;
  CrMax = CrMaximum;
  CrMin = CrMinimum;
  FScale = FScaleProb;
  popScaleFactor = popScale;
  gen = numGen;
  numClasses = numC;
}

Parameters::~Parameters() {
  // TODO Auto-generated destructor stub
}
