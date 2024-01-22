//
// Created by Fernando Cores Prado on 4/12/23.
//

#include <stdlib.h>
#include <time.h>
#include <memory.h>
#include <pthread.h>
#include "Standard_MultMat.h"
#include "Errors.h"

double elapsed_std;
static float ** result;
pthread_mutex_t mutex;


/*
* Standard Matrix multiplication with O(n^3) time complexity.
*/
float ** standardMultiplication(float ** matrixA,float ** matrixB,int n,int t)
{
    //return standardMultiplication_ijk(matrixA,matrixB,n);
    return concurrent_standardMultiplication_ikj(matrixA,matrixB,n,t);
}

/*
* Threads distribution for standard ikj Matrix.
*/
float ** concurrent_standardMultiplication_ikj(float ** matrixA,float ** matrixB,int n,int t)
{
    int i,k,spare;
    pthread_t threads[t];
    ThreadArgs args[t];
    pthread_mutex_init(&mutex, NULL);
    result = (float**)malloc(n*sizeof(float *));
    for(i;i<n;i++){
        result[i]=(float*)malloc(n*sizeof(float));
        memset(result[i],0,n*sizeof(float));
    }

    k = (n*n / t);
    spare = n*n - k*t;
    printf("n: %d, t: %d, result: %d, spare: %d\n",n,t,k,spare);

    for(i = 0; i < t; i++) {
        args[i].matrixA = matrixA;
        args[i].matrixB = matrixB;
        args[i].n = n;

        if (pthread_create(&threads[i], NULL, standardMultiplication_ikj, &(args[i])) != 0) {
            Error("[Standard]: pthread creation error\n\n");
        }
    }

    return standardMultiplication_ikj(&(args[0]));
}

/*
* Standard ijk Matrix multiplication with O(n^3) time complexity.
*/
float ** standardMultiplication_ijk(float ** matrixA,float ** matrixB,int n)
{
    struct timespec start, finish;
    float ** result;
    int i,j,k;

    clock_gettime(CLOCK_MONOTONIC, &start);

    result = (float**)malloc(n*sizeof(float *));
    for(i=0;i<n;i++){
        result[i]=(float*)malloc(n*sizeof(float));
        memset(result[i],0,n*sizeof(float));
        for(j=0;j<n;j++){
            for(k=0;k<n;k++) {
                result[i][j]=result[i][j]+(matrixA[i][k]*matrixB[k][j]);
            }
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed_std = (finish.tv_sec - start.tv_sec);
    elapsed_std += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

    return result;
}

/*
* Standard ikj Matrix multiplication with O(n^3) time complexity.
*/
float ** standardMultiplication_ikj(PtrArgs args)
{
    float ** matrixA = args -> matrixA;
    float ** matrixB = args -> matrixB;
    int n = args -> n, i = args -> i, j = args -> j, cells_n = args -> cells_n;


    struct timespec start, finish;
    //int i,j,k;
    int k,cells_calc;

    clock_gettime(CLOCK_MONOTONIC, &start);

    //result = (float**)malloc(n*sizeof(float *));
    /*for(i;i<n;i++){
        result[i]=(float*)malloc(n*sizeof(float));
        memset(result[i],0,n*sizeof(float));
        for(k=0;k<n;k++) {
            for(j=0;j<n;j++){
                result[i][j]=result[i][j]+(matrixA[i][k]*matrixB[k][j]);
            }
        }
    }*/

    while (cells_calc < cells_n)
    {
        for(k=0;k<n;k++) {
            pthread_mutex_lock(&mutex);
            result[i][j]=result[i][j]+(matrixA[i][k]*matrixB[k][j]);
            pthread_mutex_unlock(&mutex);
        }
        cells_calc++;
        if (j == n - 1)
        {
            i++;
            j = 0;
        } else {
            j++;
        }
    }
    
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed_std = (finish.tv_sec - start.tv_sec);
    elapsed_std += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

    return result;
}
