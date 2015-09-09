
#ifndef HEAP_H
#define HEAP_H

class HeapItem
{
     private:
          double distKey;                         // Distance between centroid and item
          int itemInd; 
          
     public:
          HeapItem();                              // Default constructor
HeapItem(double dist, int itemIndex);     // Constructor
          ~HeapItem();                         // Destructor
          double getKey();                         // Return item priority
          void setKey(double key);               // Set the priority key value
          int getData();                    // Return data item
	  void setData(int itemInd);          // Set the data item value
};


class Heap
{
     private:
          int          numElements;              // Number of elements in the heap
          int          heapLength;               // Size of the array

     public:
          HeapItem     *elements;                 // Pointer to dynamically allocated array
          Heap(int size);                     // Parameterized constructor
          ~Heap();                                 // Destructor
          void ReheapDown(int root, int bottom);   // Reheap after removing item
          void ReheapUp(int root, int bottom);     // Reheap after inserting item
          bool Enqueue(HeapItem *item);            // Add an item to the heap
	  void Enqueue(double dist, int itemInd);      // Add an item to the heap
          HeapItem *Dequeue();                     // Get item at the root
          int getNumElements();                    // Return number of elements in the heap
};

#endif
