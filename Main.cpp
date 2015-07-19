/*
 * Main.cpp
 *
 *  Created on: Jul 14, 2015
 *      Author: Navdha Sah
 */
#include "Population.h"
#include "Item.h"
#include "DEMain.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <stdio.h>
#include <limits>

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

int main(){
	const int kmax = 20;
	const int gen = 10;
	string ip_file, ip, lines, buffer, item;
	int val, dim, counter = 0;
	Item** objects;
	int** track;
	cout << "Enter your file name";
	//cin >> ip_file;
	 ip_file = "glass.csv";
	ip = "wc -l " + ip_file; // find the number of lines in csv file that determines the number of items to cluster.
	lines = exec(ip.c_str());
	if(!lines.empty()){
		val = atoi(lines.c_str());
	}
	ip = "head -1 " + ip_file + " | grep -o ',' | wc -l"; // find the number of commas to determine the number of features.
	lines = exec(ip.c_str());
	if(!lines.empty()){
		dim = atoi(lines.c_str());
	}
	if(val > 0){
		int col = 10*(dim-2);
		double min [dim-2];
		double max [dim-2];
		for(int i = 0; i < dim-2; i++){
			min[i] = std::numeric_limits<double>::max();
			max[i] = std::numeric_limits<double>::min();
		}
		objects = new Item*[val];
		track=new int* [val];
		for(int i = 0; i < val; i++){
			track[i] = new int[col];
			track[i][0] = i;
		}
		ifstream input_stream(ip_file.c_str());
		while (getline(input_stream, buffer, '\n'))
			{
				istringstream in(buffer);
				objects[counter] = new Item(dim-2);
				getline(in, item, ',');
				for(int i = 0; i < dim-2; i++)
				     {
						getline(in, item, ',');
						cout << item << " ";
						double value = atoi(item.c_str());
						objects[counter]->items[i] = value;
						if(min[i] > value){
							min[i] = value;
						}
						if(max[i] < value){
							max[i] = value;
						}
						getline(in, item, ',');
						objects[counter]->typeClass = atoi(item.c_str());
				     }
				cout << endl;
				counter++;
			}
		DEMain obj(kmax, dim-2, gen, track, objects, val);
		obj.setup(min, max);
		obj.run();
	}
}



