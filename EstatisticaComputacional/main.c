#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_multimin.h>
#include <time.h>
#include "metodos.h"

/* main method */

int main(void)
{
    //To keep tracking the time of execution
    clock_t begin = clock();
    /*Setting Monte Carlo configuration*/
    /*Setting Weibull parameters for the samples to be generated*/
    double shape = 3; //b
    double scale = 3; //c
    int taxa_sucesso = 0;
    int sample_size = 10;
    int nRep = 100;
    double* reg_shape = allocation(nRep);
    double* reg_scale = allocation(nRep);
    size_t semente = 1990;
    /* Generating a nRep size vector of uniform random numbers necessary for generate Weibull samples */
    double* uniform = 0;
    uniform = generator(&semente, nRep * sample_size);
    int j = 0;
    for (int i = 0; i < nRep; i++) {
        size_t iter = 0;
        int status;

        const gsl_multimin_fdfminimizer_type* T;
        gsl_multimin_fdfminimizer* s;

        /* Put here the size of the sample and generate the Weibull sample according
           to our function previously declared */
        int offsetUniform = offsetUniform = (sample_size - 1)* j;
        ++j;
        double* par = 0;
        par = rWeibull(sample_size, scale, shape, (uniform + offsetUniform));
       
        gsl_vector* x;
        gsl_multimin_function_fdf my_func;

        my_func.n = 2;
        my_func.f = my_f;
        my_func.df = my_df;
        my_func.fdf = my_fdf;
        my_func.params = par;

        /* Starting point, x = (5,7) */
        //definindo as estimativas iniciais

        x = gsl_vector_alloc(2);
        gsl_vector_set(x, 0, 5);
        gsl_vector_set(x, 1, 1);


        T = gsl_multimin_fdfminimizer_vector_bfgs2;
        // T = gsl_multimin_fdfminimizer_conjugate_fr;
        s = gsl_multimin_fdfminimizer_alloc(T, 2);
        gsl_vector* grad = gsl_multimin_fdfminimizer_gradient(s);
        gsl_multimin_fdfminimizer_set(s, &my_func, x, 1e-3, 1e-5);

        do
        {
            iter++;
            status = gsl_multimin_fdfminimizer_iterate(s);


            if (status)
                break;

            status = gsl_multimin_test_gradient(s->gradient, 1e-3);

            if (status == GSL_SUCCESS) {
                ++taxa_sucesso;
                printf("Minimum found at:\n");


                printf("%5d %.5f %.5f %10.5f %.5f %.5f\n", iter,
                    gsl_vector_get(s->x, 0),
                    gsl_vector_get(s->x, 1),
                    s->f,
                    gsl_vector_get(grad, 0), gsl_vector_get(grad, 1));

                /* saving estimators*/
                /* configuring access in a two-dimensional array */
                *(reg_shape + i) = gsl_vector_get(s->x, 0);
                *(reg_scale + i) = gsl_vector_get(s->x, 1);
            }
            
        } while (status == GSL_CONTINUE && iter < 100);
        if (status != GSL_SUCCESS) {
            
            free(uniform);
            --i;
            uniform = generator(&semente, (nRep-i) * sample_size);
            j = 0;
            
            continue;
        }
        gsl_multimin_fdfminimizer_free(s);
        gsl_vector_free(x);
        free(par);

        /* end loop */
    }

    
    /* computing the mean value of the estimators */
    //double mShape = cum_shape / nRep;
    //double mScale = cum_scale / nRep;
    double mShape = avg(reg_shape, nRep);
    double mScale = avg(reg_scale, nRep);
    double varShape = variance(reg_shape, nRep, mShape);
    double varScale = variance(reg_scale, nRep, mScale);
    printf("\nSIMULATION RESULTS FOR %d REPETITIONS WITH SAMPLE SIZE OF %d\n",nRep,sample_size);
    printf("\nShape: %.5f Scale: %.5f\n", mShape, mScale);
    printf("Variance_Shape: %.5f Variance_Scale: %.5f\n", varShape, varScale);
    printf("The bias of each estimator is: \n");
    printf("Bias_Shape: %.4f \nBias_Scale: %.4f\n", mShape - shape, mScale - scale);
    printf("The Relative Bias :\nRBias_Shape: %.4f%c \nRBias_Scale: %.4f%c",
        100 * (mShape - shape) / shape,'%', 100 * (mScale - scale) / scale,'%');

    printf("\nShape \t Scale \n");
    for (int i = 0; i < nRep; i++) {
        printf("%.2f \t%.2f\n",reg_shape[i], reg_scale[i]);
    }
    free(reg_shape);
    free(reg_scale);
    free(uniform);
    printf("\nTaxa de Sucesso: %d", taxa_sucesso);
    clock_t end = clock();
    
    double TIME_E = ((double)(end - begin)) / CLOCKS_PER_SEC;
    printf("\nThe time of execution was of %.3f seconds.",  TIME_E);
    
    return 0;
}