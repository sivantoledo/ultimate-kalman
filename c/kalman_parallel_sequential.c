#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

void foreach_in_range(void** array, int length, size_t n, void (*func)(void**, int, size_t, size_t)) {
	(*func)( array, length, 0, n );
}

void foreach_in_range_two(void** array1, void** array2, int length, size_t n, void (*func)(void**, void**, int, size_t, size_t)) {
	(*func)( array1, array2, length, 0, n );
}
