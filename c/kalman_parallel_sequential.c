#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

void foreach_in_range(void (*func)(void**, int, size_t, size_t), void** array, int length, size_t n) {
	(*func)( array, length, 0, n );
}

void foreach_in_range_two(void (*func)(void**, void**, int, size_t, size_t), void** array1, void** array2, int length, size_t n) {
	(*func)( array1, array2, length, 0, n );
}
