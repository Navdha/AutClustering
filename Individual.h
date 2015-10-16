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

class Clustering {
 public:
  double distance;
  int itemIndex;
  Clustering(double dist, int item);
  Clustering(const Clustering& obj);
  //  bool operator<(const Clustering& obj);
  bool operator==(const Clustering& obj);
};

class Individual {
 public:
  Individual(int kmax, int dim);
  Individual(const Individual& org);
  bool operator==(const Individual& org);
  ~Individual();

  double rawFitness;
  double** clusCenter;  // array holding the centroids of active clusters
  double* activationThreshold;
  bool* active;
  int numActiveClusters; // number of active clusters
  int maxNumClusters;          // maximum number of clusters
  int numFeatures;  // number of features
  vector<Clustering>** clusters;


  friend ostream& operator<<(ostream& o, const Individual& org);

};//end Individual

#endif /* INDIVIDUAL_H_ */
