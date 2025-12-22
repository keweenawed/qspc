
/* Largest number of cached parameters allowed at any given time. */
#define QSPC_JOB_QUEUE_MAX 10

/* Number of parameters to cache per thread. */
#define QSPC_JOB_CACHE_SIZE 10

/* The number of threads to use. */
#define QSPC_NUM_THREADS 4

/* Maximum values that the coefficients of the powers on q-series can take. */
#define QSPC_MAX_POWER_DEG_1 4
#define QSPC_MAX_POWER_DEG_2 4

/* Maximum values the coefficients on q-Pochhammer subscripts can take. */
#define QSPC_MAX_FAC_DEG_0 3
#define QSPC_MAX_FAC_DEG_1 3

/* Maximum values the diliations on q-Pochhammer symbols can take. */
#define QSPC_MAX_DIL_1 3
#define QSPC_MAX_DIL_2 3

/* The maximum number of q-Pochhammer symbols to allow on the numerator or
 * denominator of a q-series. */
#define QSPC_MAX_NUM_QPS 1

/* The parameters for a particular q-series are encoded in an array of
 * integers with this length. The first 4 * QSPC_MAX_NUM_QPS entries
 * give the numerator q-Pochhammer symbols $(q^a; q^b)_{cn+d}$:
 *   index + 0   c
 *   index + 1   d
 *   index + 2   a
 *   index + 3   b
 * Here 0 <= index < QSPC_MAX_NUM_QPS. Given a particular index, if c = 0
 * the whole expression is taken to equal 1, and the other parameters must
 * be set to 0, as well as with all larger indices. The denominator is encoded
 * the same way directly following this, and the last 4 entries are:
 *   QSPC_PARAMETER_LENGTH - 4   The degree 2 term for the leading power.
 *   QSPC_PARAMETER_LENGTH - 3   The degree 1 term.
 *   QSPC_PARAMETER_LENGTH - 2   The denominator under both of these terms.
 *   QSPC_PARAMETER_LENGTH - 1   Set to -1 to give an alternating sign (-1)^n
 *                               and otherwise set to 1. */
#define QSPC_PARAMETER_LENGTH (8 * QSPC_MAX_NUM_QPS + 4)

/* The number of terms to compute for each q-series. Larger values are likely
 * to result in integer overflow without using a big integer library. */
#define QSPC_COEFFICIENT_BOUND 100

/* The largest pattern length to check for in a factored q-series.*/
#define QSPC_PATTERN_BOUND 20

/* Returns the number of q-Pochhammer symbols in the numerator of a q-series.
 *  parameters: The parameters that encode the series. */
static inline int64_t QSPC_num_qps(int64_t *parameters)
{
	int64_t length = 0;

	do {
		if (parameters[4 * length + 0] == 0) break;

		++length;
	} while (length < QSPC_MAX_NUM_QPS);

	return length;
}

/* Returns the number of q-Pochhammer symbols in the denominator of a
 * q-series.
 *  parameters: The parameters that encode the series. */
static inline int64_t QSPC_den_qps(int64_t *parameters)
{
	int64_t length = 0;

	do {
		if (parameters[4 * QSPC_MAX_NUM_QPS + 4 * length + 0] == 0)
			break;

		++length;
	} while (length < QSPC_MAX_NUM_QPS);

	return length;
}

