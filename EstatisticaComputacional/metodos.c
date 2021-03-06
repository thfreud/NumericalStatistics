#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_multimin.h>
#include <time.h>
#include <gsl/gsl_statistics.h>
#include "metodos.h"
#define M 4294967295
#define a 1664525
#define d 1013904223

/*function to allocate memory and inicialize it*/
double* allocation(int size) {
	double* memory = (double*)malloc(size * sizeof(double));
	if (memory != NULL) {
		for (int i = 0; i < size; i++) {
			memory[i] = 0;
		}
	}
	return memory;
}

/*Pseudorandom number generator with parameters M, a, d*/
double* generator(size_t *seed, int size) {
	double* memory = allocation(size);
	for (int i = 0; i < size; i++) {
		memory[i] = (double)(*seed = (*seed * a + d) % M) / M;
	}
	
	return memory;
}

/* This function implements a weibull random variable generator using the cumulative probability function inverse transform method */

double* rWeibull(double size, double b, double c, double * uniform) {
	double* sample = (double*)malloc((size + 1) * sizeof(double));
	//double* u = generator(seed, size + 1);
	for (size_t i = 0; i < size + 1; i++) {
		if (i == 0) {
			//doing this, we avoid giving another vector in the likelihood function declaration
			sample[i] = size;
		}
		else
		{
			sample[i] = b * pow((-log(uniform[i])), (double)1 / c);
		}
	}
	
	return sample;
}

/*This is the log-likelihood function of a weibull distribution with two parameters.
*/
double my_f(const gsl_vector* v, void* params)
{
	double c, b;
	double* p = (double*)params;
	
	c = gsl_vector_get(v, 0);
	b = gsl_vector_get(v, 1);
	
	double total = 0;
	for (size_t i = 0; i < p[0]; i++) {
		total += (c - 1) * log(p[i + 1])  - pow(p[i + 1] / b, c);
	}
	//we keep this part outside the loop for improvements in performance.
	total += p[0]*(log(c) - c * log(b));
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
		sum_t2 += (c / pow(b, c + 1) * pow(p[i + 1], c));
		sum_t1 += log(p[i + 1]) - log(p[i + 1] / b) * pow(p[i + 1] / b, c);

	}
	//we keep this part outside the loop for improvements in performance.
	//p[0] corresponds to the size of the sample vector.
	sum_t2 += p[0]*(-c / b);
	sum_t1 += p[0]*((1 / c) - log(b));
	gsl_vector_set(df, 0, (-1) * sum_t1);
	gsl_vector_set(df, 1, (-1) * sum_t2);
}

/* Compute both f and df together. */
void my_fdf(const gsl_vector* x, void* params, double* f, gsl_vector* df)
{
	*f = my_f(x, params);
	my_df(x, params, df);
}

/*Functions for computations of basic statistics*/

void printData(double * reg_shape, double * reg_scale, int nRep, double shape, double scale, int sample_size) {
	double mShape = gsl_stats_mean(reg_shape, 1, nRep);
	double mScale = gsl_stats_mean(reg_scale, 1, nRep);
	double varShape = gsl_stats_variance(reg_shape, 1, nRep);
	double varScale = gsl_stats_variance(reg_scale, 1, nRep);
	double largestShape = gsl_stats_max(reg_shape, 1, nRep);
	double largestScale = gsl_stats_max(reg_scale, 1, nRep);
	double skewShape = gsl_stats_skew(reg_shape, 1, nRep);
	double kurtosisShape = gsl_stats_kurtosis(reg_shape, 1, nRep);
	double skewScale = gsl_stats_skew(reg_scale, 1, nRep);
	double kurtosisScale = gsl_stats_kurtosis(reg_scale, 1, nRep);

	printf("\n--- SIMULATION RESULTS FOR %d REPETITIONS WITH SAMPLE SIZE OF %d ---\n", nRep, sample_size);
	printf("\nShape: %.5f \nScale: %.5f\n", mShape, mScale);
	printf("Variance_Shape: %.5f \nVariance_Scale: %.5f\n", varShape, varScale);
	printf("Largest_Shape: %.5f \nLargest_Scale: %.5f\n", largestShape, largestScale);
	printf("Skewness_Shape: %.5f \nSkewness_Scale: %.5f\n", skewShape, skewScale);

	/*As gsl implements the excess kurtosis approach, one has to sum 3 in order to get kurtosis.*/
	printf("Kurtosis_Shape: %.5f \nKurtosis_Scale: %.5f\n", kurtosisShape+3, kurtosisScale+3);
	printf("The bias of each estimator is: \n");
	printf("Bias_Shape: %.4f \nBias_Scale: %.4f\n", mShape - shape, mScale - scale);
	printf("The Relative Bias :\nRBias_Shape: %.4f%c \nRBias_Scale: %.4f%c",
		100 * (mShape - shape) / shape, '%', 100 * (mScale - scale) / scale, '%');
}