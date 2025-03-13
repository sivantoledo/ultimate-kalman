#ifndef CONCURRENT_SET_H
#define CONCURRENT_SET_H

#include <stdint.h>
#include <stdlib.h>
 
#include "parallel.h"
#include "concurrent_set.h"

typedef struct concurrent_set_st {
  int            size;
  void**         pointers;          
  spin_mutex_t** locks; // pointers to locks, to allow locks of any type
  void   (*foreach)(void*);
} concurrent_set_t;

concurrent_set_t* concurrent_set_create (int capacity, void (*foreach)(void*));
void              concurrent_set_free   (concurrent_set_t* set);
void              concurrent_set_insert (concurrent_set_t* set, void* element);
void              concurrent_set_foreach(concurrent_set_t* set);

#endif
