#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

void foreach_step_in_range(void** step_pointers, int length, size_t n, void (*func)(void**, int, size_t, size_t)) {
	(*func)( step_pointers, length, 0, n );
}
