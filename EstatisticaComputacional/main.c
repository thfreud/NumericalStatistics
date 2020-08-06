/*Incluir um programa para estimativa de viés usando bootstrap*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_multimin.h>
#include <time.h>
#include "metodos.h"

/* main method */

int main(void)
{
    /* Setting Monte Carlo configuration*/
    int nRep = 10;
    double* param = allocation(nRep * 2);
    double size = 0;
    size_t semente = 1990;
    
    for (int i = 0; i < nRep; i++) {

        size_t iter = 0;
        int status;

        const gsl_multimin_fdfminimizer_type* T;
        gsl_multimin_fdfminimizer* s;

        /* Put here the size of the sample and generate the Weibull sample according
        to our function previously declared*/
        size += 800;
        double* par = rWeibull(semente, size, 3, 3);


        //printf("%.2f\n", par[0]);

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
        gsl_multimin_fdfminimizer_set(s, &my_func, x, 1e-2, 1e-6);

        do
        {
            iter++;
            status = gsl_multimin_fdfminimizer_iterate(s);


            if (status)
                break;

            status = gsl_multimin_test_gradient(s->gradient, 1e-1);

            if (status == GSL_SUCCESS)
                printf("Minimum found at:\n");
            /*
            printf("%5d %.5f %.5f %10.5f %.5f %.5f\n", iter,
                gsl_vector_get(s->x, 0),
                gsl_vector_get(s->x, 1),
                s->f,
                gsl_vector_get(grad, 0), gsl_vector_get(grad, 1)); */
                /* saving estimators*/
               // param[i][0] = gsl_vector_get(s->x, 0);
               //param[i][1] = gsl_vector_get(s->x, 1);
                /* configuring access in a two-dimensional array */
            * (param + i * 2) = gsl_vector_get(s->x, 0);
            *(param + i * 2 + 1) = gsl_vector_get(s->x, 1);

        } while (status == GSL_CONTINUE && iter < 100);

        gsl_multimin_fdfminimizer_free(s);
        gsl_vector_free(x);
        free(par);

        /* end loop */
    }

    double shape = 0;
    double scale = 0;
    for (int i = 0; i < nRep; i++) {
        //shape += param[i][0];
        //scale += param[i][1];
        shape += *(param + i * 2);
        scale += *(param + i * 2 + 1);


    }
    /* computing the mean value of the estimators */
    double mShape = shape / nRep;
    double mScale = scale / nRep;
    printf("\n%.5f %.5f\n", mShape, mScale);
    printf("The bias of each estimator is: \n");
    printf("B_Shape: %.4f \nB_Scale: %.4f\n", mShape - 3, mScale - 3);
    printf("The Relative Bias :\nR_B_Shape: %.4f \nR_B_Scale: %.4f", 100 * (mShape - 3) / 3, 100 * (mScale - 3) / 3);


    free(param);
    return 0;
}