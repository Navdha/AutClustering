
#include <iostream>
#include <cassert>
#include "Heap.h"

using namespace std;

HeapItem::HeapItem()
{
     distKey = 0.0;
     itemInd = 0.0;
}
HeapItem::HeapItem(double dist, int itemIndex)
{
     distKey = dist;
     itemInd = itemIndex;
}
HeapItem::~HeapItem()
{
}
// Return item priority
double HeapItem::getKey()
{
     return distKey;
}

// Set the priority key value
void HeapItem::setKey(double key)
{
     distKey = key;
}

// Return data item
int HeapItem::getData()
{
     return itemInd;
}

// Set the data item value
void HeapItem::setData(int i)
{
     itemInd = i;
}

////////////////////////////////////////////////////////
//   Class Heap implementation
////////////////////////////////////////////////////////

Heap::Heap(int size)
{
     // Create heap of given size
     elements = new HeapItem[size];
     numElements = 0;
     heapLength = size;
}

Heap::~Heap()
{
     delete[] elements;
}


void Heap::ReheapDown(int root, int bottom)
{
     int minChild;
     int rightChild;
     int leftChild;
     HeapItem temp;

     leftChild = root * 2 + 1;          // Get index of root's left child
     rightChild = root * 2 + 2;          // Get index of root's right child

     // Check base case in recursive calls.  If leftChild's index is less
     // than or equal to the bottom index we have not finished recursively 
     // reheaping.
     if(leftChild <= bottom)               
     {
       if(leftChild == bottom)          // If this root has no right child then       {
	 minChild = leftChild;     //     leftChild must hold min key
     
       else
	 {     // Get the one lowest in the tree (highest index in the array)
	   if(elements[leftChild].getKey() >= elements[rightChild].getKey())
	     minChild = rightChild;
	   else
	     minChild = leftChild;
	 }
       if(elements[root].getKey() > elements[minChild].getKey())
	 {
	   // Swap these two elements
	   temp = elements[root];
	   elements[root] = elements[minChild];
	   elements[minChild] = temp;
	   // Make recursive call till reheaping completed
	   ReheapDown(minChild, bottom);
       }
     }
}



void Heap::ReheapUp(int root, int bottom)
{
     int parent;
     HeapItem temp;

     // Check base case in recursive calls.  If bottom's index is greater
     // than the root index we have not finished recursively reheaping.
     if(bottom > root)
     {
          parent = (bottom -1) / 2;
          if(elements[parent].getKey() > elements[bottom].getKey())
          {
               // Swap these two elements
               temp = elements[parent];
               elements[parent] = elements[bottom];
               elements[bottom] = temp;
               // Make recursive call till reheaping completed
               ReheapUp(root, parent);
          }
     }
}

// Add an item to the heap

bool Heap::Enqueue(HeapItem *item)
{
     if(numElements < heapLength)
     {
          elements[numElements] = *item; // Copy item into array
          ReheapUp(0, numElements);
          numElements++;
          return true;
     }
     return false;
}

// Add an item to the heap

void Heap::Enqueue(double dist, int ind)
{
     bool retVal;
     HeapItem *temp = new HeapItem(dist, ind);
     retVal = Enqueue(temp);
     delete temp;  // Delete this dynamically created one
     assert(retVal==true);
}

// Get item at the root
HeapItem *Heap::Dequeue()
{
  HeapItem *temp = new HeapItem(elements[0].getKey(), elements[0].getData());
    numElements--;
    // Copy last item into root
    elements[0] = elements[numElements];
    // Reheap the tree
    ReheapDown(0, numElements - 1);
    if(numElements == 0) return NULL;
    else
      return temp;
}

// Return number of elements in the heap
int Heap::getNumElements()
{
     return numElements;
}
