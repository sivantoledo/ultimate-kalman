#ifndef CONCURRENT_SET_H
#define CONCURRENT_SET_H

#include <stdint.h>
#include <stdlib.h>
 
#include "concurrent_set.h"

struct concurrent_set_st;
typedef struct concurrent_set_st concurrent_set_t;

concurrent_set_t* concurrent_set_create (int capacity, void (*foreach)(void*));
void              concurrent_set_free   (concurrent_set_t* set);
void              concurrent_set_insert (concurrent_set_t* set, void* element);
void              concurrent_set_foreach(concurrent_set_t* set);

#endif
