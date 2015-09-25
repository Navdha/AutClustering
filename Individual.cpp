/*
 * Individual.cpp
 *
 *  Created on: Jul 7, 2015
 *      Author: Navs
 */

#include "Individual.h"
#include <stdexcept>

Individual::Individual(int kmax, int dim) {
  // TODO Auto-generated constructor stub
  //	cout << "Individual class constructor called." << endl;
  rawFitness = 0.0;
  maxNumClusters = kmax;
  numFeatures = dim;
  clusCenter = new double* [maxNumClusters];
  clusters = new vector<int>*[maxNumClusters];
  for (int count = 0; count < maxNumClusters; count++)
    {
      clusCenter[count] = new double[numFeatures];
      clusters[count] = new vector<int>;
    }
  activationThreshold =  new double[maxNumClusters];
  active = new bool[maxNumClusters];
  numActiveClusters = 0;
}

Individual::Individual(const Individual& org){
  //copy constructor
  rawFitness = org.rawFitness;
  maxNumClusters = org.maxNumClusters;
  numFeatures = org.numFeatures;
  clusCenter = new double*[maxNumClusters];
  activationThreshold = new double[maxNumClusters];
  clusters = new vector<int>*[maxNumClusters];
  active = new bool[maxNumClusters];
  numActiveClusters = org.numActiveClusters;
  for(int count = 0; count < maxNumClusters; count++){
    clusCenter[count] = new double[numFeatures];
    activationThreshold[count] = org.activationThreshold[count];
    active[count] = org.active[count];
    clusters[count] = org.clusters[count];
    for(int i = 0; i < numFeatures; i++){
      clusCenter[count][i] = org.clusCenter[count][i];
    }
  }
}

ostream& operator<<(ostream& o, const Individual& org)
{//overloaded insertion operator
  o << "#active clusters: " << org.numActiveClusters << endl;
  o << "fitness: " << org.rawFitness << endl;

  for (int i = 0; i < org.maxNumClusters; i++){
    if (org.active[i]){
      o << "cluster " << i << " centroid: (" << org.clusCenter[i][0];
      for (int j = 1; j < org.numFeatures; j++){
	o << ", " << org.clusCenter[i][j];
      }
      o << ")" << endl;
    }//end if
  }//end for i
  return o;
}//operator<<

Individual::~Individual() {
  // TODO Auto-generated destructor stub
  //cout << "Destructor called" << endl;
  for(int i = 0; i < maxNumClusters; ++i) {
    delete [] clusCenter[i];
    delete clusters[i];
  }
  delete [] clusCenter;
  delete [] active;
  delete [] activationThreshold;
  delete [] clusters;

}

