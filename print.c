#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include "qspc.h"

pthread_mutex_t QSPC_print_lock;

/* Helper function for QSPC_report_identity that nicely prints powers of q. */
static inline void print_power(int64_t power)
{
	if (power == 0) {
		printf(" 1 ");
	} else if (power == 1) {
		printf(" q ");
	} else {
		printf(" q^{%lld} ", power);
	}
}

/* Prints out a sum-product identity formatted in LaTeX.
 *   parameters: The series parameters.
 *   signature: The pattern of powers for the product.
 *   modulus: The length of signature. */
void QSPC_report_identity(int64_t *parameters, int64_t *signature,
			  int64_t modulus)
{
	bool product_frac = false;
	bool numerator_empty = true;
	int64_t num_qps = QSPC_num_qps(parameters);
	int64_t den_qps = QSPC_den_qps(parameters);

	/* This function needs to be thread safe. */
	pthread_mutex_lock(&QSPC_print_lock);

	printf("\\begin{equation}\n");

	/* Arbitrary choice to help equations fit on the page. */
	if (modulus >= 10) {
		printf("\\begin{aligned}\n&");
	}

	for (int64_t index = 0; index < modulus; ++index) {
		if (signature[index] > 0) {
			product_frac = true;
			break;
		}
	}

	if (product_frac) printf("\\frac{");

	for (int64_t index = 0; index < modulus; ++index) {
		if (signature[index] >= 0) continue;

		printf("(");
		print_power(index + 1);
		printf("; q^{%lld})_\\infty ", modulus);

		if (signature[index] != -1) printf("^{%lld}",
		    -signature[index]);

		numerator_empty = false;
	}

	if (product_frac) {
		if (numerator_empty) printf("1");

		printf("}{");

		for (int64_t index = 0; index < modulus; ++index) {
			if (signature[index] <= 0) continue;

			printf("(");
			print_power(index + 1);
			printf("; q^{%lld})_\\infty ", modulus);

			if (signature[index] != 1)
				printf("^{%lld}", signature[index]);
		}

		printf("}");
	}

	if (modulus == 1 && signature[0] == 0) printf("1");

	if (modulus >= 10) printf("\\\\&");

	printf(" = \\sum_{n=0}^\\infty ");

	if (den_qps != 0) printf("\\frac{");

	if (parameters[QSPC_PARAMETER_LENGTH - 1] == -1) printf("(-1)^n");

	printf("q^{");

	if (parameters[QSPC_PARAMETER_LENGTH - 2] != 1 &&
	    parameters[QSPC_PARAMETER_LENGTH - 3] != 0) printf("(");

	if (parameters[QSPC_PARAMETER_LENGTH - 4] == 1) {
		printf("n^2");
	} else if (parameters[QSPC_PARAMETER_LENGTH - 4] > 1) {
		printf("%lld n^2", parameters[QSPC_PARAMETER_LENGTH - 4]);
	}

	if (parameters[QSPC_PARAMETER_LENGTH - 4] != 0 &&
	    parameters[QSPC_PARAMETER_LENGTH - 3] != 0)
		printf("+");

	if (parameters[QSPC_PARAMETER_LENGTH - 3] == 1) {
		printf("n");
	} else if (parameters[QSPC_PARAMETER_LENGTH - 3] > 1) {
		printf("%lld n", parameters[QSPC_PARAMETER_LENGTH - 3]);
	}

	if (parameters[QSPC_PARAMETER_LENGTH - 2] != 1) {
		if (parameters[QSPC_PARAMETER_LENGTH - 3] != 0) printf(")");

		printf("/%lld", parameters[QSPC_PARAMETER_LENGTH - 2]);
	}

	printf("}");

	for (int64_t index = 0; index < num_qps; ++index) {
		printf("(-");

		print_power(parameters[4 * index + 2]);
		printf(";");
		print_power(parameters[4 * index + 3]);
		printf(")");

		if (parameters[4 * index + 0] == 1) {
			printf("_{n");
		} else {
			printf("_{%lld n", parameters[4 * index + 0]);
		}

		if (parameters[4 * index + 1] != 0)
			printf(" + %lld", parameters[4 * index + 1]);

		printf("}");
	}

	if (den_qps == 0) {
		if (modulus >= 10) {
			printf("\n\\end{aligned}");
		}

		printf("\n\\end{equation}\n\n");
		pthread_mutex_unlock(&QSPC_print_lock);
		return;
	}

	printf("}{");

	for (int64_t index = 0; index < den_qps; ++index) {
		printf("(");

		print_power(parameters[4 * QSPC_MAX_NUM_QPS + 4 * index + 2]);
		printf(";");
		print_power(parameters[4 * QSPC_MAX_NUM_QPS + 4 * index + 3]);
		printf(")");

		if (parameters[4 * QSPC_MAX_NUM_QPS + 4 * index + 0] == 1) {
			printf("_{n");
		} else {
			printf("_{%lld n", parameters[4 * QSPC_MAX_NUM_QPS
						      + 4 * index + 0]);
		}

		if (parameters[4 * QSPC_MAX_NUM_QPS + 4 * index + 1] != 0)
			printf(" + %lld", parameters[4 * QSPC_MAX_NUM_QPS
						     + 4 * index + 1]);

		printf("}");
	}

	if (modulus >= 10) {
		printf("}\n\\end{aligned}");
		printf("\n\\end{equation}\n\n");
	} else {
		printf("}\n\\end{equation}\n\n");
	}

	pthread_mutex_unlock(&QSPC_print_lock);
}

