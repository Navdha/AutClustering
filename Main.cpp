/*
 * Main.cpp
 *
 *  Created on: Jul 14, 2015
 *      Author: Navdha Sah
 */
#include "Population.h"
#include "Item.h"
#include "DEMain.h"
#include "Parameters.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include <cstdlib>
#include <stdio.h>
#include <limits>
#include <ctime>

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

int main(int argc, char** argv){
	string ip_file, ip, lines, buffer, item;
	int val, dim, counter = 0;
	Item** objects;
	int** track;
	srand(time(NULL));
	//cin >> ip_file;
	 if (argc != 10){
	    cerr << "Usage: " << argv[0]
	         << " filename CrMax CrMin scaleFactor thresholdVal kmax kmin numGenerations validityIndex"
	         << endl;
	    exit(1);
	  }
	ip_file = argv[1];
	double CrMaximum = atof(argv[2]);
	double CrMinimum = atof(argv[3]);
	double FScaleProb = atof(argv[4]);
	double threshVal = atof(argv[5]);
	int kmax = atoi(argv[6]);
	int kmin = atoi(argv[7]);
	int gen = atoi(argv[8]);
	int validityIndex = atoi(argv[9]);

	char* filename = argv[1];
	for (int i = 2; i < argc-1; i++) {
		strcat(filename, "_");
		strcat(filename, argv[i]);
	}
	if (validityIndex == 1) {
		strcat(filename, "_DB");
	} else if (validityIndex == 2) {
		strcat(filename, "_CS");
	} else {
		strcat(filename, "_PB");
	}
	strcat(filename, ".txt");
	string resultFileName(filename);

	ip = "wc -l " + ip_file; // find the number of lines in csv file that determines the number of items to cluster.
	lines = exec(ip.c_str());
	if(!lines.empty()){
		val = atoi(lines.c_str());
	}
	ip = "head -1 " + ip_file + " | grep -o ',' | wc -l"; // find the number of commas to determine the number of features.
	lines = exec(ip.c_str()); //rename lines to numCommas
	if(!lines.empty()){
		dim = atoi(lines.c_str());
	}
	cout << val << endl;
	//	int numFeatures = dim -1;
	int numFeatures = dim;
	if(val > 0){
		int col = 10*(numFeatures);//rename to popSize
		double min [numFeatures];
		double max [numFeatures];
		for(int i = 0; i < numFeatures; i++){
			min[i] = std::numeric_limits<double>::max();
			max[i] = std::numeric_limits<double>::min();
		}
		objects = new Item*[val];
		track=new int* [val];
		for(int i = 0; i < val; i++){
			track[i] = new int[col];
		}
		ifstream input_stream(ip_file.c_str());
		while (getline(input_stream, buffer, '\n'))
			{
				istringstream in(buffer);
				objects[counter] = new Item(numFeatures);
				//	getline(in, item, ',');
				//	objects[counter]->typeClass = atoi(item.c_str());
				for(int i = 0; i < numFeatures; i++)
				     {
						getline(in, item, ',');
						double value = atof(item.c_str());
						objects[counter]->items[i] = value;
						if(min[i] > value){
							min[i] = value;
						}
						if(max[i] < value){
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

		Parameters param(CrMaximum, CrMinimum, FScaleProb, threshVal, kmax, kmin, gen);

		DEMain obj(numFeatures, track, objects, val, validityIndex, param);
		obj.calcDistBtwnItems();
		obj.setup(min, max);
		obj.run(min, max, resultFileName);
	}
}



