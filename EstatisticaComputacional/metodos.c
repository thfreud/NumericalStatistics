#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_multimin.h>
#include <time.h>
#include "metodos.h"
#define M 4294967295
#define a 1664525
#define d 1013904223

/* Implementing a multiplicative congruential generator for a sequencies of pseudo random number */
/* Constants for the generator */

/*function to allocate memory and inicialize it*/
double* allocation(int size) {
	double* memory = (double*)malloc(size * sizeof(double));
	for (int i = 0; i < size; i++) {
		memory[i] = 0;
	}
	return memory;
}

/*Pseudorandom number generator with parameters M, a, d*/
double* generator(size_t seed, int size) {
	double* memory = allocation(size);
	for (int i = 0; i < size; i++) {
		memory[i] = (double)(seed = (seed * a + d) % M) / M;
	}

	return memory;

}

/* This function implements a weibull random variable generator using the cumulative probability function inverse transform method */

double* rWeibull(size_t seed, double size, double b, double c) {
	//srand(1990);
	double* sample = (double*)malloc((size + 1) * sizeof(double));
	double* u = generator(seed, size + 1);
	for (size_t i = 0; i < size + 1; i++) {
		if (i == 0) {
			sample[i] = size;
		}
		else
		{
			//double u = (double)rand() / RAND_MAX;

			sample[i] = b * pow((-log(1 - u[i])), (double)1 / c);
		}


	}
	free(u);
	return sample;
}

/*This is the log-likelihood function of a weibull distribution with two parameters,
and here the first and second functions are used.
*/
double my_f(const gsl_vector* v, void* params)
{

	double c, b;
	double* p = (double*)params;


	c = gsl_vector_get(v, 0);
	b = gsl_vector_get(v, 1);

	double total = 0;
	for (size_t i = 0; i < p[0]; i++) {
		total += log(c) + (c - 1) * log(p[i + 1]) - c * log(b) - pow(p[i + 1] / b, c);
	}
	return -1 * total;

}

/*Here we declare the gradient.*/
void my_df(const gsl_vector* v, void* params, gsl_vector* df)
{
	double c, b;
	double* p = (double*)params;



	c = gsl_vector_get(v, 0);
	b = gsl_vector_get(v, 1);

	double sum_t1 = 0;
	double sum_t2 = 0;

	for (size_t i = 0; i < p[0]; i++) {
		sum_t2 += (-c / b) + (c / pow(b, c + 1) * pow(p[i + 1], c));
		sum_t1 += (1 / c) + log(p[i + 1]) - log(b) - log(p[i + 1] / b) * pow(p[i + 1] / b, c);

	}
	gsl_vector_set(df, 0, (-1) * sum_t1);
	gsl_vector_set(df, 1, (-1) * sum_t2);

}

/* Compute both f and df together. */
void my_fdf(const gsl_vector* x, void* params, double* f, gsl_vector* df)
{
	*f = my_f(x, params);
	my_df(x, params, df);
}