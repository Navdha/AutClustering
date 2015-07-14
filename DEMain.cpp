/*
 * DE.cpp
 *
 *  Created on: Jul 8, 2015
 *      Author: Navs
 */

#include "DEMain.h"
#include "Population.h"
#include "Item.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <stdio.h>

using namespace std;
DEMain::DEMain() {
	// TODO Auto-generated constructor stub

}

DEMain::~DEMain() {
	// TODO Auto-generated destructor stub
}

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
	string ip_file, ip, lines, buffer, item;
	int val, dim, counter = 0;
	Item** objects;
	cout << "Enter your file name";
	cin >> ip_file;
	ip = "wc -l " + ip_file;
	//ip = "wc -l glass.csv";
	lines = exec(ip.c_str());
	if(!lines.empty()){
		val = atoi(lines.c_str());
	}
	ip = "head -1 " + ip_file + " | grep -o ',' | wc -l";
	lines = exec(ip.c_str());
	if(!lines.empty()){
		dim = atoi(lines.c_str());
	}
	if(val > 0){
		objects = new Item*[val];

		ifstream input_stream(ip_file.c_str());
		while (getline(input_stream, buffer, '\n'))
			{
				istringstream in(buffer);
				objects[counter] = new Item(dim-1);
				getline(in, item, ',');
				for(int i = 0; i < dim-1; i++)
				     {
						getline(in, item, ',');
						cout << item << " ";
						objects[counter]->items[i] = atoi(item.c_str());
						getline(in, item, ',');
						objects[counter]->typeClass = atoi(item.c_str());
				     }
				cout << endl;
				counter++;
			}


	}
}
