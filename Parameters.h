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
  int maxNumClusters;
  int minNumClusters;
  int popScaleFactor;
  int numClasses;
  double gen;
  Parameters(double CrMaximum, double CrMinimum, double FScaleProb, double threshVal, int maxClusters, int minClusters, int popScale, double numGen, int numC);
  ~Parameters();
};


#endif /* PARAMETERS_H_ */
