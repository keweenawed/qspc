#include <stdint.h>

/* Computes all positive integer divisors. Returns the quantity found.
 *   value: The integer to find the divisors of.
 *   list: The divisors are written here. The caller must provide enough
 *     memory to prevent overflow. */
int64_t QSPC_integer_divisors(int64_t value, int64_t *list)
{
	int64_t index1 = 0;

	/* TODO: Use more efficient algorithm and memoization. */
	for (int64_t index2 = 1; index2 <= value / 2; ++index2) {
		if (value % index2 == 0) {
			list[index1++] = index2;
		}
	}

	list[index1++] = value;

	return index1;
}

