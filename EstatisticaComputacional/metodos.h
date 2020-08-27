#pragma once
double* allocation(int size);
double* generator(size_t *seed, int size);

/*This function implements a weibull random variable generator using the cumulative probability function inverse transform method */
double* rWeibull(double size, double b, double c);

/*This is the log-likelihood function of a weibull distribution with two parameters,
and here the first and second functions are used.
*/
double my_f(const gsl_vector* v, void* params);

/*Here we declare the gradient.*/
void my_df(const gsl_vector* v, void* params, gsl_vector* df);

/* Compute both f and df together. */
void my_fdf(const gsl_vector* x, void* params, double* f, gsl_vector* df);