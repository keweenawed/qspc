#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include "qspc.h"

extern void QSPC_report_identity(int64_t *, int64_t *, int64_t);
extern int64_t QSPC_find_pattern(int64_t *, int64_t *);
extern void QSPC_find_product_form(int64_t *, int64_t *, int64_t);
extern void QSPC_build_series(int64_t *, int64_t *, int64_t);
extern int64_t QSPC_pattern_gcd(int64_t *, int64_t);

extern pthread_mutex_t QSPC_print_lock;

/* Helper function for work_recursive_step. Given a combination of parameters,
 * this function attempts to factor the q-series, and if successful, has
 * the result printed out.
 *   parameters: The combination of parameters to try. */
static void try_combination(int64_t *parameters)
{
	int64_t buffer1[QSPC_COEFFICIENT_BOUND];
	int64_t buffer2[QSPC_COEFFICIENT_BOUND];
	int64_t buffer3[QSPC_PATTERN_BOUND];
	int64_t period;

	QSPC_build_series(parameters, buffer1, QSPC_COEFFICIENT_BOUND);
	QSPC_find_product_form(buffer1, buffer2, QSPC_COEFFICIENT_BOUND);
	period = QSPC_find_pattern(buffer2, buffer3);

	if (period == 0) return;

	/* Throw out any dilated results since these are redundant. */
	if (QSPC_pattern_gcd(buffer3, period) != 1) return;

	QSPC_report_identity(parameters, buffer3, period);
}

/* Recursively generates all combinations of allowed series parameters, and
 * checks if they are candidates for identities.
 *   parameters: The series parameters, which are treated by the function as
 *     a list of loop indices.
 *   depth: Tracks the recursion depth.
 *   identity: Value that identifies which thread is calling. Used to
 *     divide the work. */
static void work_recursive_step(int64_t *parameters, int64_t depth,
				int64_t identity)
{
	int64_t offset;
	bool first_step;

	switch (depth) {
	case QSPC_PARAMETER_LENGTH - 2:
		parameters[QSPC_PARAMETER_LENGTH - 2] = 1;
		parameters[QSPC_PARAMETER_LENGTH - 1] = 1;
		try_combination(parameters);
	
		/* Try the same but with an alternating sign. */
		parameters[QSPC_PARAMETER_LENGTH - 1] = -1;
		try_combination(parameters);

		/* If both power coefficients are odd, we can try to find an
		 * identity with both of them divided by 2. */
		if (parameters[QSPC_PARAMETER_LENGTH - 4] % 2 == 1 &&
		    parameters[QSPC_PARAMETER_LENGTH - 3] % 2 == 1) {
			parameters[QSPC_PARAMETER_LENGTH - 2] = 2;
			parameters[QSPC_PARAMETER_LENGTH - 1] = 1;
			try_combination(parameters);
			parameters[QSPC_PARAMETER_LENGTH - 1] = -1;
			try_combination(parameters);
		}

		return;

	case QSPC_PARAMETER_LENGTH - 3:

		/* This is an arbitrary choice, but the threads need to
		 * delegate work somewhere. TODO work out a queue system? */
		parameters[depth] = identity;
		work_recursive_step(parameters, depth + 1, identity);

		//for (parameters[depth] = 0; parameters[depth]
		//     < QSPC_MAX_POWER_DEG_1; ++parameters[depth]) {
		//	work_recursive_step(parameters, depth + 1, identity);
		//}

		return;
	case QSPC_PARAMETER_LENGTH - 4:
		for (parameters[depth] = 1; parameters[depth]
		     < QSPC_MAX_POWER_DEG_2; ++parameters[depth]) {
			work_recursive_step(parameters, depth + 1, identity);
		}

		return;
	}

	first_step = false;

	if (depth < 4 * QSPC_MAX_NUM_QPS) {
		offset = 0;

		if (depth < 3) first_step = true;
	} else {
		offset = 4 * QSPC_MAX_NUM_QPS;

		if (depth < 4 * QSPC_MAX_NUM_QPS + 3) first_step = true;
	}

	switch (depth % 4) {
	case 0:
		parameters[depth] = 0;
		work_recursive_step(parameters, offset + 4 * QSPC_MAX_NUM_QPS,
				    identity);

		for (parameters[depth] = 1; parameters[depth]
		     <= QSPC_MAX_FAC_DEG_1; ++parameters[depth]) {
			if (!first_step && parameters[depth - 4] <
			    parameters[depth]) return;

			work_recursive_step(parameters, depth + 1, identity);
		}

		break;
	case 1:
		for (parameters[depth] = 0; parameters[depth]
		     <= QSPC_MAX_FAC_DEG_0; ++parameters[depth]) {
			work_recursive_step(parameters, depth + 1, identity);
		}

		break;
	case 2:
		for (parameters[depth] = 1; parameters[depth]
		     <= QSPC_MAX_DIL_1; ++parameters[depth]) {
			work_recursive_step(parameters, depth + 1, identity);
		}

		break;

	case 3:
		for (parameters[depth] = 1; parameters[depth]
		     <= QSPC_MAX_DIL_2; ++parameters[depth]) {
			work_recursive_step(parameters, depth + 1, identity);
		}
	}
}

/* Entry point for each thread.
 *   identity: A value between 0 and QSPC_NUM_THREADS - 1 that uniquely
 *     identifies the thread. */
static void *worker_thread(void *identity)
{
	int64_t parameters[QSPC_PARAMETER_LENGTH];

	work_recursive_step(parameters, 0, (int64_t)identity);
	pthread_exit(NULL);
}

int main(int argc, char **argv)
{
	pthread_t threads[QSPC_NUM_THREADS];

	(void)argc;
	(void)argv;

	printf("\\documentclass[10pt]{article}\n");
	printf("\\usepackage{amsmath}\n");
	printf("\\usepackage[margin=0.1in]{geometry}\n\\begin{document}\n");

	pthread_mutex_init(&QSPC_print_lock, NULL);

	for (int64_t index = 0; index < QSPC_NUM_THREADS; ++index) {
		pthread_create(&threads[index], NULL, worker_thread,
			       (void *)(index));
	}

	for (int64_t index = 0; index < QSPC_NUM_THREADS; ++index) {
		pthread_join(threads[index], NULL);
	}

	pthread_mutex_destroy(&QSPC_print_lock);
	printf("\\end{document}\n");

	return 0;
}

