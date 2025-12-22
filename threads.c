#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "qspc.h"

extern void QSPC_report_identity(int64_t *, int64_t *, int64_t);
extern int64_t QSPC_find_pattern(int64_t *, int64_t *);
extern void QSPC_find_product_form(int64_t *, int64_t *, int64_t);
extern void QSPC_build_series(int64_t *, int64_t *, int64_t);
extern int64_t QSPC_pattern_gcd(int64_t *, int64_t);
extern void QSPC_generate_divisors(void);
extern void QSPC_delete_divisors(void);

extern pthread_mutex_t QSPC_print_lock;

/* The main thread generates a (FILO) queue of parameters for the worker
 * threads to try. This takes the form of a linked list. */
struct QSPC_job_entry
{
	/* Number of parameter combinations. */
	int64_t entries;

	/* Points to the next linked list entry. */
	struct QSPC_job_entry *next;

	/* An array of parameter combinations. */
	int64_t parameters[QSPC_JOB_CACHE_SIZE][QSPC_PARAMETER_LENGTH];
};

/* Points to the front of the queue. */
static struct QSPC_job_entry *QSPC_job_queue;

/* Number of elements in the queue. */
static volatile int64_t QSPC_job_queue_length;

/* Set to true on initialization, and then to false when the main thread is
 * finished generating parameter combinations, and is only waiting for the
 * worker threads to finish. */
static volatile bool QSPC_keep_working;

/* Lock to modify the queue.  */
static pthread_mutex_t QSPC_job_lock;

/* Lock and conditional variable for the main thread. */
static pthread_mutex_t QSPC_generator_lock;
static pthread_cond_t QSPC_generator_cond;

/* Helper function for work_recursive_step. Adds a parameter combination to
 * the queue, handling the multithreading logic. */
static void submit_parameters(int64_t *parameters)
{
	/* By design this pointer cannot be NULL. If this condition is met,
	 * the element on the front of the queue is full. */
	if (QSPC_job_queue->entries == QSPC_JOB_CACHE_SIZE) {
		struct QSPC_job_entry *job;

		/* If the queue is completely full, wait until it is emptied.
		 * One of the worker threads will then signal to continue. */
		if (QSPC_job_queue_length == QSPC_JOB_QUEUE_MAX) {
			pthread_mutex_unlock(&QSPC_job_lock);
			pthread_cond_wait(&QSPC_generator_cond,
					  &QSPC_generator_lock);
		}

		/* Prepend a new queue element. */
		job = malloc(sizeof(struct QSPC_job_entry));
		job->entries = 0;
		job->next = QSPC_job_queue;
		QSPC_job_queue = job;
		++QSPC_job_queue_length;
	}

	/* Add the parameter combination to the element at the front of the
	 * queue, which is guaranteed to exist and not be full here. */
	for (int64_t index = 0; index < QSPC_PARAMETER_LENGTH; ++index) {
		QSPC_job_queue->parameters[QSPC_job_queue->entries][index]
			= parameters[index];
	}

	++QSPC_job_queue->entries;
}

/* Recursively generates all combinations of allowed series parameters, and
 * checks if they are candidates for identities.
 *   parameters: The series parameters, which are treated by the function as
 *     a list of loop indices.
 *   depth: Tracks the recursion depth. */
static void work_recursive_step(int64_t *parameters, int64_t depth)
{
	int64_t offset;
	bool first_step;

	switch (depth) {

	/* The furthest depth, where the combinations are finished. */
	case QSPC_PARAMETER_LENGTH - 2:
		parameters[QSPC_PARAMETER_LENGTH - 2] = 1;
		parameters[QSPC_PARAMETER_LENGTH - 1] = 1;
		submit_parameters(parameters);
	
		/* Try the same but with an alternating sign. */
		parameters[QSPC_PARAMETER_LENGTH - 1] = -1;
		submit_parameters(parameters);

		/* If both power coefficients are odd, we can try to find an
		 * identity with both of them divided by 2. */
		if (parameters[QSPC_PARAMETER_LENGTH - 4] % 2 == 1 &&
		    parameters[QSPC_PARAMETER_LENGTH - 3] % 2 == 1) {
			parameters[QSPC_PARAMETER_LENGTH - 2] = 2;
			parameters[QSPC_PARAMETER_LENGTH - 1] = 1;
			submit_parameters(parameters);
			parameters[QSPC_PARAMETER_LENGTH - 1] = -1;
			submit_parameters(parameters);
		}

		return;

	case QSPC_PARAMETER_LENGTH - 3:
		for (parameters[depth] = 0; parameters[depth]
		     < QSPC_MAX_POWER_DEG_1; ++parameters[depth]) {
			work_recursive_step(parameters, depth + 1);
		}

		return;
	case QSPC_PARAMETER_LENGTH - 4:
		for (parameters[depth] = 1; parameters[depth]
		     < QSPC_MAX_POWER_DEG_2; ++parameters[depth]) {
			work_recursive_step(parameters, depth + 1);
		}

		return;
	}

	switch (depth % 4) {
	case 0:
		/* Here we take a number of arbitrary extra steps to reduce
		 * the number of duplicate results. */
		first_step = false;

		if (depth < 4 * QSPC_MAX_NUM_QPS) {
			offset = 0;

			if (depth < 3) first_step = true;
		} else {
			offset = 4 * QSPC_MAX_NUM_QPS;

			if (depth < 4 * QSPC_MAX_NUM_QPS + 3)
				first_step = true;
		}

		/* If the degree 1 parameter on the subscript of the
		 * q-Pochhammer symbol is 0, skip the rest of the
		 * parameters since the whole symbol is taken to equal 1. */
		parameters[depth] = 0;
		work_recursive_step(parameters, offset + 4
				    * QSPC_MAX_NUM_QPS);

		for (parameters[depth] = 1; parameters[depth]
		     <= QSPC_MAX_FAC_DEG_1; ++parameters[depth]) {

			/* Weakly order the q-Pochhammer symbols. */
			if (!first_step && parameters[depth - 4] <
			    parameters[depth]) return;

			work_recursive_step(parameters, depth + 1);
		}

		break;
	case 1:
		for (parameters[depth] = 0; parameters[depth]
		     <= QSPC_MAX_FAC_DEG_0; ++parameters[depth]) {
			work_recursive_step(parameters, depth + 1);
		}

		break;
	case 2:
		for (parameters[depth] = 1; parameters[depth]
		     <= QSPC_MAX_DIL_1; ++parameters[depth]) {
			work_recursive_step(parameters, depth + 1);
		}

		break;

	case 3:
		for (parameters[depth] = 1; parameters[depth]
		     <= QSPC_MAX_DIL_2; ++parameters[depth]) {
			work_recursive_step(parameters, depth + 1);
		}
	}
}

/* Given a combination of parameters, this function generates the q-series
 * and then attempts to factor it. If successful, the identity is printed. */
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

/* Entry point for each worker thread. */
static void *worker_thread(void *argument)
{
	struct QSPC_job_entry *job;

	(void)argument;

	for (;;) {
		pthread_mutex_lock(&QSPC_job_lock);

		if (QSPC_job_queue_length == 0) {
			if (!QSPC_keep_working) {
				pthread_mutex_unlock(&QSPC_job_lock);
				pthread_exit(NULL);
			}

			/* If the queue is empty and QSPC_keep_working is true
			 * then the main thread is waiting to be woken up. */
			pthread_cond_signal(&QSPC_generator_cond);
			continue;
		}

		/* Take a job off the queue. */
		job = QSPC_job_queue;
		QSPC_job_queue = QSPC_job_queue->next;
		--QSPC_job_queue_length;
		pthread_mutex_unlock(&QSPC_job_lock);

		/* Here the actual work is done. Check every combination. */
		for (int64_t index = 0; index < job->entries; ++index) {
			try_combination(job->parameters[index]);
		}

		free(job);
	}
}

int main(int argc, char **argv)
{
	pthread_t threads[QSPC_NUM_THREADS];
	int64_t parameters[QSPC_PARAMETER_LENGTH];

	/* Prevent compiler warnings for unused variables. */
	(void)argc;
	(void)argv;

	QSPC_keep_working = true;

	/* Create a permanent list of divisors for QSPC_find_product_form. */
	QSPC_generate_divisors();

	pthread_mutex_init(&QSPC_print_lock, NULL);
	pthread_mutex_init(&QSPC_job_lock, NULL);
	pthread_mutex_init(&QSPC_generator_lock, NULL);
	pthread_cond_init(&QSPC_generator_cond, NULL);

	/* The multithreading logic requires these to start locked. */
	pthread_mutex_lock(&QSPC_job_lock);
	pthread_mutex_lock(&QSPC_generator_lock);

	/* Header for the LaTeX file. */
	printf("\\documentclass[10pt]{article}\n");
	printf("\\usepackage{amsmath}\n");
	printf("\\usepackage[margin=0.1in]{geometry}\n\\begin{document}\n");

	for (int64_t index = 0; index < QSPC_NUM_THREADS; ++index) {
		pthread_create(&threads[index], NULL, worker_thread, NULL);
	}

	/* Start generating the job queue. */
	QSPC_job_queue = malloc(sizeof(struct QSPC_job_entry));
	QSPC_job_queue->entries = 0;
	QSPC_job_queue_length = 1;
	work_recursive_step(parameters, 0);

	/* At this point, the work is nearly done. Wait for each thread to
	 * finish the queue, and then clean up. */
	QSPC_keep_working = false;
	pthread_mutex_unlock(&QSPC_job_lock);
	pthread_mutex_unlock(&QSPC_generator_lock);

	for (int64_t index = 0; index < QSPC_NUM_THREADS; ++index) {
		pthread_join(threads[index], NULL);
	}

	QSPC_delete_divisors();
	pthread_mutex_destroy(&QSPC_print_lock);
	pthread_mutex_destroy(&QSPC_job_lock);
	pthread_mutex_destroy(&QSPC_generator_lock);
	pthread_cond_destroy(&QSPC_generator_cond);

	printf("\\end{document}\n");

	return 0;
}

