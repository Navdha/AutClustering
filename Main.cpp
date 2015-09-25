/*
 * Main.cpp
 *
 *  Created on: Jul 14, 2015
 *      Author: Navdha Sah
 */
#include "Population.h"
#include "Item.h"
#include "Parameters.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include <cstdlib>
#include <stdio.h>
#include <limits>
#include <ctime>
#include "DEMain.h"

using namespace std;

string exec(const char* cmd) {
  FILE* pipe = popen(cmd, "r");
  if (!pipe) return "ERROR";
  char buffer[128];
  std::string result = "";
  while(!feof(pipe)) {
    if(fgets(buffer, 128, pipe) != NULL)
      result += buffer;
  }
  pclose(pipe);
  return result;
}

string get_date()
{
  time_t now;
  char the_date[80];
  the_date[0] = '\0';

  now = time(0);

  if (now != -1)
    {
      strftime(the_date, 80, "%d_%m_%Y.%X", localtime(&now));
    }

  return string(the_date);
}

int main(int argc, char** argv){
  string ip_file, ip, lines, numCommas, buffer, item;
  int val, dim, counter = 0;
  Item** objects;
  int** track;
  srand(time(NULL));
  if (argc != 13){
    cerr << "Usage: " << argv[0]
	 << " filename CrMax CrMin scaleFactor thresholdVal maxNumClusters minNumClusters popScalingFactor numGenerations validityIndex numOfClasses cycleRepetition"
	 << endl;
    exit(1);
  }
  ip_file = argv[1];
  double CrMaximum = atof(argv[2]);
  double CrMinimum = atof(argv[3]);
  double FScaleProb = atof(argv[4]);
  double threshVal = atof(argv[5]);
  int maxNumClusters = atoi(argv[6]);
  int minNumClusters = atoi(argv[7]);
  int popScaleFactor = atoi(argv[8]);
  double gen = atoi(argv[9]);
  int validityIndex = atoi(argv[10]);
  int numClasses = atoi(argv[11]);
  int numRepeat = atoi(argv[12]);
  string str = get_date();
  char *cstr = new char[str.length() + 1];
  strcpy(cstr, str.c_str());
  
  char* filename = argv[1];
  // for (int i = 2; i < argc-1; i++) {
    strcat(filename, "_");
    strcat(cstr, "_");
    //   strcat(filename, argv[4]);
    //}
  if (validityIndex == 1) {
    strcat(filename, "_DB");
  } else if (validityIndex == 2) {
    strcat(filename, "_CS");
  } else {
    strcat(filename, "_PB");
  }
  strcat(filename, ".txt");
  string resultFileName(filename);
  delete [] cstr;
  ip = "wc -l " + ip_file; // find the number of lines in csv file that determines the number of items to cluster.
  lines = exec(ip.c_str());
  if(!lines.empty()){
    val = atoi(lines.c_str());
  }
  ip = "head -1 " + ip_file + " | grep -o ',' | wc -l"; // find the number of commas to determine the number of features.
  numCommas = exec(ip.c_str()); //rename lines to numCommas
  if(!numCommas.empty()){
    dim = atoi(numCommas.c_str());
  }
  cout << val << endl;
  int numFeatures = dim;
	if (val > 0) {
		int popSize = popScaleFactor * (numFeatures);
		double min[numFeatures];
		double max[numFeatures];
		for (int i = 0; i < numFeatures; i++) {
			min[i] = std::numeric_limits<double>::max();
			max[i] = std::numeric_limits<double>::min();
		}
		objects = new Item*[val];
		ifstream input_stream(ip_file.c_str());
		while (getline(input_stream, buffer, '\n')) {
			istringstream in(buffer);
			objects[counter] = new Item(numFeatures);
			for (int i = 0; i < numFeatures; i++) {
				getline(in, item, ',');
				double value = atof(item.c_str());
				objects[counter]->items[i] = value;
				if (min[i] > value) {
					min[i] = value;
				}
				if (max[i] < value) {
					max[i] = value;
				}
			}
			getline(in, item, ',');
			objects[counter]->typeClass = atoi(item.c_str());
			counter++;
		}
			cout << "The arrays are " <<endl;
			for(int i = 0; i < numFeatures; i++){
						cout << min[i] << " ";
	
					}
			cout <<endl;
			for(int i = 0; i < numFeatures; i++){
			cout << max[i] << " ";
			}
	
		Parameters param(CrMaximum, CrMinimum, FScaleProb, threshVal,
				 maxNumClusters, minNumClusters, popScaleFactor, gen, numClasses, numRepeat);
		DEMain obj(numFeatures, objects, val, validityIndex, param);
		obj.calcDistBtwnItems(min, max);
		obj.setup(min, max);
		obj.run(min, max, resultFileName);
	}
}



