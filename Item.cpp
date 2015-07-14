/*
 * Item.cpp
 *
 *  Created on: Jul 10, 2015
 *      Author: Navs
 */

#include "Item.h"
#include <iostream>

using namespace std;

Item::Item(int dim) {
	// TODO Auto-generated constructor stub
	//cout << "Constructor called" << endl;
 this->dim = dim;
 items = new double[dim];
 typeClass = 0;
}

Item::~Item() {
	// TODO Auto-generated destructor stub
}

