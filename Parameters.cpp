/*
 * Parameters.cpp
 *
 *  Created on: Aug 6, 2015
 *      Author: Navdha Sah
 */

#include "Parameters.h"

Parameters::Parameters(){
  //default constructor
}

Parameters::Parameters(double CrMaximum, double CrMinimum, double FScaleProb,
		       double threshVal,  int maxClusters, int minClusters, int popScale, double numGen, int numC, int numR) {
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
  numRepeat = numR;
}

Parameters& Parameters::operator=(Parameters& obj){
  maxNumClusters = obj.maxNumClusters;
  minNumClusters = obj.minNumClusters;
  threshold = obj.threshold;
  CrMax = obj.CrMax;
  CrMin = obj.CrMin;
  FScale = obj.FScale;
  popScaleFactor = obj.popScaleFactor;
  gen = obj.gen;
  numClasses = obj.numClasses;
  numRepeat = obj.numRepeat;
  return *this;
}

Parameters::~Parameters() {
  // TODO Auto-generated destructor stub
}
