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
    concurrent_standardMultiplication_ikj(matrixA,matrixB,n,t);

    return result;
}

/*
* Threads distribution for standard ikj Matrix.
*/
void * concurrent_standardMultiplication_ikj(float ** matrixA,float ** matrixB,int n,int t)
{
    int i,k,j,spare,last_i = 0,last_j = 0;
    struct timespec start, finish;
    pthread_t threads[t];
    ThreadArgs args[t];
    pthread_mutex_init(&mutex, NULL);
    result = (float**)malloc(n*sizeof(float *));
    for(i;i<n;i++){
        result[i]=(float*)malloc(n*sizeof(float));
        memset(result[i],0,n*sizeof(float));
    }

    clock_gettime(CLOCK_MONOTONIC, &start);

    k = (n*n / t);
    spare = n*n - k*t;

    float ** _return;

    for(i = 0; i < t; i++) {
        args[i].matrixA = matrixA;
        args[i].matrixB = matrixB;
        args[i].n = n;
        args[i].i = last_i;
        args[i].j = last_j;
        if (spare > 0) {
            spare--;
            args[i].cells_n = k + 1;
        } else {
            args[i].cells_n = k;
        }
        
        // Update last position
        last_j += args[i].cells_n;
        while(last_j >= n) {
            last_j -= n;
            last_i++;
        }

        if (pthread_create(&threads[i], NULL, (void *(*) (void*)) standardMultiplication_ikj, &(args[i])) != 0) {
            Error("[Standard]: pthread creation error\n\n");
        }
    }
    
    for (i = 0; i < t; i++)
    {
        pthread_join(threads[i], (void **) NULL);
    }
    pthread_mutex_destroy(&mutex);

    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed_std = (finish.tv_sec - start.tv_sec);
    elapsed_std += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
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
void * standardMultiplication_ikj(PtrArgs args)
{
    float ** matrixA = args -> matrixA;
    float ** matrixB = args -> matrixB;
    int n = args -> n, i = args -> i, j = args -> j, cells_n = args -> cells_n;

    struct timespec start, finish;
    int k,cells_calc = 0;

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
}
