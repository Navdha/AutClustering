/*
 * Item.h
 *
 *  Created on: Jul 10, 2015
 *      Author: Navs
 */

#ifndef ITEM_H_
#define ITEM_H_

class Item {
 public:
  Item(int dim);
  ~Item();
  int numFeatures;
  double* items;
  int typeClass;
};

#endif /* ITEM_H_ */
