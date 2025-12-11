#include <stdint.h>

/* Computes the Cauchy product of two truncated series.
 *   series1: Coefficients of the first series.
 *   series2: Coefficients of the second series.
 *   result: Where the coefficients of the product is written.
 *   bound: The length of each of these arrays. */
static void truncated_product(int64_t *series1, int64_t *series2,
			      int64_t *result, uint64_t bound)
{
	for (uint64_t index1 = 0; index1 < bound; ++index1) {
		result[index1] = 0;

		for (uint64_t index2 = 0; index2 <= index1; ++index2) {
			result[index1] += series1[index1 - index2]
					* series2[index2];
		}
	}
}

/* Computes the coefficients of a q-Pochhammer symbol in the numerator
 * of the form $(\pm q^a; q^b)_n$.
 *   dilation1: The value a, which can be 0.
 *   dilation2: The value b, which must be at least 1. 
 *   factors: The value n, which can be 0. 
 *   sign: Set to 1 or -1 to give the sign in front of $q^a$.
 *   result: The array of coefficients the product is written to.
 *   bound: The length of this array. The product is truncated to fit, and
 *     zeroes are padded if the result is smaller. */
static void expand_q_pochhammer_num(uint64_t dilation1, uint64_t dilation2,
				    uint64_t factors, int64_t sign,
				    int64_t *result, uint64_t bound)
{
	int64_t buffer[bound];

	result[0] = 1;

	for (uint64_t index = 1; index < bound; ++index) result[index] = 0;

	for (uint64_t index1 = 0; index1 < factors; ++index1) {
		uint64_t offset = index1 * dilation2 + dilation1;

		if (offset >= bound) break;

		/* Multiplication here is simple enough to avoid needing
		 * to use truncated_product. */
		for (uint64_t index = 0; index < bound - offset; ++index)
			buffer[index] = result[index];

		for (uint64_t index = 0; index < bound - offset; ++index)
			result[index + offset] -= sign * buffer[index];
	}
}

/* Computes the coefficients of a q-Pochhammer symbol in the denominator.
 * of the form $(\pm q^a; q^b)_n^{-1}$.
 *   dilation1: The value a, which must be at least 1.
 *   dilation2: The value b, which must be at least 1. 
 *   factors: The value n, which can be 0. 
 *   sign: Set to 1 or -1 to give the sign in front of $q^a$.
 *   result: The array of coefficients the product is written to.
 *   bound: The length of this array. The product is truncated to fit. */
static void expand_q_pochhammer_den(uint64_t dilation1,
				    uint64_t dilation2,
				    uint64_t factors, int64_t sign,
				    int64_t *result, uint64_t bound)
{
	result[0] = 1;

	for (uint64_t index = 1; index < bound; ++index) result[index] = 0;

	for (uint64_t index1 = 0; index1 < factors; ++index1) {
		int64_t buffer1[bound];
		int64_t buffer2[bound];

		for (uint64_t index2 = 0; index2 < bound; ++index2) {
			buffer1[index2] = 0;
			buffer2[index2] = result[index2];
		}

		if (sign == 1) {
			for (uint64_t index2 = 0; index2 < bound; index2
			     += dilation2 * index1 + dilation1)
				buffer1[index2] = 1;
		} else {
			int64_t flip = 1;

			for (uint64_t index2 = 0; index2 < bound;
			     index2 += dilation2 * index1 + dilation1) {
				flip = (flip == 1) ? -1 : 1;
				buffer1[index2] = flip;
			}
		}

		truncated_product(buffer1, buffer2, result, bound);
	}
}

/* Returns the degree of the q-Multinomial coefficient.
 *   top: The parameter on the top.
 *   bottom: An array of the parameters on the bottom.
 *   length: The number of parameters on the bottom. */
static inline uint64_t q_multinomial_degree(uint64_t top, uint64_t *bottom,
					    uint64_t length)
{
	uint64_t degree = top * (top + 1) / 2;

	// check if parameters make sense?

	for (uint64_t index = 0; index < length; ++index) {
		degree -= bottom[index] * (bottom[index] + 1) / 2;
	}

	return degree;
}

/* Computes the q-Multinomial coefficient.
 *   top: The parameter on the top in the q-Multinomial.
 *   bottom: An array of the parameters on the bottom in the q-Multinomial.
 *   length: The number of parameters on the bottom.
 *   result: Where the coefficients of the result are written.
 *   bound: The length of this array. If the result is smaller, zeroes are
 *     padded at the end, and otherwise the result is truncated to fit. */
static void expand_q_multinomial(uint64_t top, uint64_t *bottom,
				 uint64_t length, int64_t *result,
				 uint64_t bound)
{
	uint64_t degree = q_multinomial_degree(top, bottom, length);
	uint64_t adj_bound = (degree <= bound) ? degree + 1 : bound;
	int64_t buffer1[adj_bound];
	int64_t buffer2[adj_bound];
	uint64_t bottom_sum = 0;

	/* If the bottom parameters do not exactly sum to top, the result
	 * is taken by convention to be 0. */
	for (uint64_t index = 0; index < length; ++index)
		bottom_sum += bottom[index];

	if (bottom_sum != top) {
		for (uint64_t index = 0; index < bound; ++index)
			result[index] = 0;

		return;
	}

	result[0] = 1;

	for (uint64_t index = 1; index < bound; ++index) result[index] = 0;

	for (uint64_t index1 = 0; index1 < length; ++index1) {
		for (uint64_t index2 = 0; index2 < adj_bound; ++index2)
			buffer2[index2] = result[index2];

		expand_q_pochhammer_den(1, 1, bottom[index1], 1,
					 buffer1, adj_bound);
		truncated_product(buffer1, buffer2, result, adj_bound);
	}

	for (uint64_t index = 0; index < adj_bound; ++index)
		buffer2[index] = result[index];

	expand_q_pochhammer_num(1, 1, top, 1, buffer1, adj_bound);
	truncated_product(buffer1, buffer2, result, adj_bound);
}

/* Computes the q-Binomial coefficient.
 *   top: The parameter on the top of the q-Binomial.
 *   bottom: The parameter on the bottom of the q-Binomial.
 *   result: The array of coefficients the result is written to.
 *   bound: The length of this array. */
static inline void expand_q_binomial(uint64_t top, uint64_t bottom,
				     int64_t *result, uint64_t bound)
{
	if (bottom > top) {
		for (uint64_t index = 0; index < bound; ++index)
			result[index] = 0;
	} else {
		uint64_t parameters[2] = {bottom, top - bottom};

		expand_q_multinomial(top, parameters, 2, result, bound);
	}
}

extern uint64_t QSPC_integer_divisors(int64_t, int64_t *);

/* Uniquely factors a truncated series with constant term 1 into a product of
 * geometric series so that when expanded, the coefficients match up to the
 * bound. This product takes the form $\prod_{k=1}^n \frac{1}{(1-q^k)^{a_k}}$.
 *   series: The series to be factored.
 *   powers: The list of geometric series powers $a_k$. The first term $a_0$
 *     is taken to be 0 for convenience.
 *   bound: The length of the series array. The highest power coefficient that
 *     is guaranteed to match is that of $q^n$, where $n$ is one less than
 *     the value of bound. */
void QSPC_find_product_form(int64_t *series, int64_t *powers, uint64_t bound)
{
	/* Assume bound is at least 2. */
	powers[0] = 0;
	powers[1] = series[1];

	/* This algorithm is derived from observations in the book The Theory
	 * of Partitions by George Andrews. */
	for (int64_t index1 = 1; index1 < (int64_t)bound; ++index1) {
		int64_t power = 0;
		int64_t length;

		/* This is a very lazy estimate for the number of divisors. */
		int64_t divisors[bound];

		for (int64_t index2 = 1; index2 < index1; ++index2) {
			length = QSPC_integer_divisors(index2, divisors);

			for (int64_t index3 = 0; index3 < length; ++index3) {
				power -= series[index1 - index2]
				      * divisors[index3]
				      * powers[divisors[index3]];
			}
		}

		length = QSPC_integer_divisors(index1, divisors);

		for (int64_t index2 = 0; index2 < length - 1; ++index2) {
			power -= divisors[index2] * powers[divisors[index2]];
		}

		power /= index1;
		power += series[index1];
		powers[index1] = power;
	}
}

/* Stores data to describe a particular q-series. */
struct QSPC_series
{
	/* TODO: Support for q-multinomials and multiple summation indices. */ 

	/* Set to 1 to attach a alternating sign, and otherwise set to 0. */
	int64_t sign_flip;

	/* Number of q-Pochhammer symbols on the numerator. */
	uint64_t num_qps;

	/* Number of q-Pochhammer symbols on the denominator. */
	uint64_t den_qps;

	/* Each q-Pochhammer symbol on the numerator will take the form
	 * $(\pm q^a; q^b_{cn+d}$. These arrays encode these values
	 * respectively. */
	int64_t *num_signs;
	int64_t *num_dil_1;
	int64_t *num_dil_2;
	int64_t *num_fac_deg_1;
	int64_t *num_fac_deg_0;

	/* Same as above, but for the denominator. */
	int64_t *den_signs;
	int64_t *den_dil_1;
	int64_t *den_dil_2;
	int64_t *den_fac_deg_1;
	int64_t *den_fac_deg_0;

	/* These values encode the power $q^{(an^2 + bn)/c}$. */
	int64_t pow_num_deg_2;
	int64_t pow_num_deg_1;
	int64_t pow_den;
};

/* Helper function for build_series. Finds the contributions to the series
 * given by a particular summation index.
 *   parameters: The parameters that encode the series. 
 *   result: The array the coefficients of the terms are written to.
 *   bound: The length of this array.
 *   summation_index: The index of the term being computed. */
static void build_series_term(struct QSPC_series *parameters, int64_t *result,
			      uint64_t bound, int64_t summation_index)
{
	int64_t buffer1[bound];
	int64_t buffer2[bound];

	result[0] = 1;

	for (uint64_t index = 1; index < bound; ++index) result[index] = 0;

	for (uint64_t index1 = 0; index1 < parameters->num_qps; ++index1) {
		for (uint64_t index2 = 0; index2 < bound; ++index2)
			buffer1[index2] = result[index2];

		expand_q_pochhammer_num(parameters->num_dil_1[index1],
					parameters->num_dil_2[index1],
					parameters->num_fac_deg_1[index1]
					* summation_index
					+ parameters->num_fac_deg_0[index1],
					-1, buffer2, bound);
		truncated_product(buffer1, buffer2, result, bound);
	}

	for (uint64_t index1 = 0; index1 < parameters->den_qps; ++index1) {
		for (uint64_t index2 = 0; index2 < bound; ++index2)
			buffer1[index2] = result[index2];

		expand_q_pochhammer_den(parameters->den_dil_1[index1],
					parameters->den_dil_2[index1],
					parameters->den_fac_deg_1[index1]
					* summation_index
					+ parameters->den_fac_deg_0[index1],
					1, buffer2, bound);
		truncated_product(buffer1, buffer2, result, bound);
	}
}

/* Computes the truncated coefficients of a q-series series.
 *   parameters: The parameters that encode the series. 
 *   result: The array the coefficients of the terms are written to.
 *   bound: The length of this array, and the number of coefficients found. */
static void build_series(struct QSPC_series *parameters, int64_t *result,
			 uint64_t bound)
{

	for (uint64_t index = 0; index < bound; ++index) result[index] = 0;

	for (uint64_t index1 = 0;; ++index1) {
		uint64_t offset = (parameters->pow_num_deg_2 * index1 * index1
				 + parameters->pow_num_deg_1 * index1)
				 / parameters->pow_den;
		int64_t flip;

		/* This assumes that the power at least weakly grows with the
		 * summation index. If this is not the case, this can get
		 * stuck in an infinite loop. */
		if (offset >= bound) return;

		int64_t buffer[bound - offset];

		build_series_term(parameters, buffer, bound - offset, index1);

		if (parameters->sign_flip && (index1 % 2) == 1) {
			flip = -1;
		} else {
			flip = 1;
		}

		for (uint64_t index2 = 0; index2 < bound - offset; ++index2)
			result[index2 + offset] += flip * buffer[index2];
	}
}

