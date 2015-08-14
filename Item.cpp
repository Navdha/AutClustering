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
  numFeatures = dim;
  items = new double[numFeatures];
  typeClass = 0;
}

Item::~Item() {
  // TODO Auto-generated destructor stub
}

