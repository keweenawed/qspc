#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include "qspc.h"

/* Returns the positive greatest common divisor of value1 and value2. */
static int64_t gcd_pair(int64_t value1, int64_t value2)
{
	if (value2 < 0) value2 = -value2;

	if (value1 == 0) return value2;

	if (value1 < 0) value1 = -value1;

	return gcd_pair(value2 % value1, value2);
}

/* Returns the positive greatest common divisor of all the given values.
 *   pattern: An array of length modulus.
 *   modulus: Must be a positive integer. */
int64_t QSPC_pattern_gcd(int64_t *pattern, int64_t modulus)
{
	int64_t gcd = gcd_pair(modulus, pattern[0]);

	for (int64_t index = 1; index < modulus; ++index) {
		gcd = gcd_pair(pattern[index], gcd);

		if (gcd == 1) break;
	}

	return gcd;
}

int64_t *QSPC_divisor_cache[QSPC_COEFFICIENT_BOUND];
int64_t QSPC_divisor_values[QSPC_COEFFICIENT_BOUND];

/* Sets *divisors to point to the array of divisors of the provided value.
 * Returns the number of divisors in this array. */
int64_t QSPC_divisors(int64_t value, int64_t **divisors)
{
	*divisors = QSPC_divisor_cache[value];

	return QSPC_divisor_values[value];
}

/* Computes and stores the divisors of every integer between 0 and
 * QSPC_COEFFICIENT_BOUND. Called at program initialization. */
void QSPC_generate_divisors(void)
{
	int64_t buffer[QSPC_COEFFICIENT_BOUND];
	int64_t num_divisors;
	size_t length;

	for (int64_t index1 = 1; index1 < QSPC_COEFFICIENT_BOUND; ++index1) {
		num_divisors = 0;

		/* Brute force compute the divisors. This is only done once,
		 * so there is little need for a more efficient method. */
		for (int64_t index2 = 1; index2 <= index1 / 2; ++index2) {
			if (index1 % index2 == 0) {
				buffer[num_divisors++] = index2;
			}
		}

		buffer[num_divisors++] = index1;
		length = (size_t)num_divisors  * sizeof(int64_t);
		QSPC_divisor_cache[index1] = malloc(length);

		for (int64_t index2 = 0; index2 < num_divisors; ++index2) {
			QSPC_divisor_cache[index1][index2] = buffer[index2];
		}

		QSPC_divisor_values[index1] = num_divisors;
	}
}

/* Frees up the list of divisors. Not important now, but may be useful if
 * the scope of this project becomes large enough. */
void QSPC_delete_divisors(void)
{
	for (int64_t index = 1; index < QSPC_COEFFICIENT_BOUND; ++index)
		free(QSPC_divisor_cache[index]);
}

/* Helper function for QSPC_find_pattern. Checks a particular pattern length.
 * Returns true if the pattern exists, and false otherwise.
 *   powers: The list of powers of the factored series.
 *   period: The pattern length to check. */
static bool check_period(int64_t *powers, int64_t period)
{
	for (int64_t index1 = 0;; ++index1) {
		for (int64_t index2 = 1; index2 <= period; ++index2) {

			if (period * index1 + index2
			    == QSPC_COEFFICIENT_BOUND) return true;

			if (powers[index2] != powers[period * index1
			    + index2]) return false;
		}
	}
}

/* Looks for a repeating pattern in the powers of a factored series. Returns
 * the length of the pattern if it exists, or 0 otherwise.
 *   powers: The list of powers of the factored series.
 *   pattern: If a pattern is found, the sequence is written here. */
int64_t QSPC_find_pattern(int64_t *powers, int64_t *pattern)
{
	for (int64_t index1 = 1; index1 <= QSPC_PATTERN_BOUND; ++index1) {
		if (!check_period(powers, index1)) continue;

		for (int64_t index2 = 0; index2 < index1; ++index2) {
			pattern[index2] = powers[index2 + 1];
		}

		return index1;
	}

	return 0;
}

