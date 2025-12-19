#include <stdbool.h>
#include <stdint.h>
#include "qspc.h"

/* Returns the greatest common divisor of value1 and value2. */
static int64_t gcd_pair(int64_t value1, int64_t value2)
{
	if (value1 == 0) return value2;

	return gcd_pair(value2 % value1, value2);
}

/* Returns the greatest common divisor of all the given values.
 *   pattern: An array of length modulus.
 *   modulus: Must be a positive integer. */
int64_t QSPC_pattern_gcd(int64_t *pattern, int64_t modulus)
{
	int64_t gcd = gcd_pair(modulus, pattern[0]);

	for (int64_t index = 1; index < modulus; ++index) {
		gcd = gcd_pair(pattern[index], gcd);

		if (gcd < 0) gcd = -gcd;

		if (gcd == 1) break;
	}

	return gcd;
}

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

