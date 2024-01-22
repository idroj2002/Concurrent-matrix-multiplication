//
// Created by Fernando Cores Prado on 4/12/23.
//

#include "Strassens_MultMat.h"
#include <time.h>
#include <stdlib.h>
#include "Matrix.h"

double elapsed_str;
int Dim2StopRecursivity = 10;
static float ** result;
static int availableThreads;
pthread_mutex_t mutex;

/*
* Wrapper function over strassensMultRec.
*/
float ** strassensMultiplication(float ** matrixA, float** matrixB,int n,int t)
{
    struct timespec start, finish;
    clock_gettime(CLOCK_MONOTONIC, &start);
    availableThreads = t - 1;
    pthread_mutex_init(&mutex, NULL);

    if (n>32)
        Dim2StopRecursivity = n/16;

    strassenMultRec(matrixA,matrixB,n,result);

    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed_str = (finish.tv_sec - start.tv_sec);
    elapsed_str += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

    return result;
}

/*
* Strassen's Multiplication algorithm using Divide and Conquer technique.
*/
void * strassensMultRec(float ** matrixA, float** matrixB,int n,float** finalResult){
    float ** result = createZeroMatrix(n);
    if(n>Dim2StopRecursivity) {
        //Consult available threads
        int k, spare;
        int currentAvailableThreads = 0;
        pthread_mutex_lock(&mutex);
        if (availableThreads > 6) {
            currentAvailableThreads = 6;
            availableThreads -= 6;
        } else if (availableThreads > 0)
        {
            currentAvailableThreads = availableThreads;
            availableThreads = 0;
        }
        pthread_mutex_unlock(&mutex);

        //Calculate operations per thread
        k = (7 / currentAvailableThreads);
        spare = 7 - k*currentAvailableThreads;

        //Create threads
        pthread_t threads[currentAvailableThreads];
        ThreadArgs args[currentAvailableThreads];

        //Divide the matrix
        float ** a11 = divide(matrixA, n, 0, 0);
        float ** a12 = divide(matrixA, n, 0, (n/2));
        float ** a21 = divide(matrixA, n, (n/2), 0);
        float ** a22 = divide(matrixA, n, (n/2), (n/2));
        float ** b11 = divide(matrixB, n, 0, 0);
        float ** b12 = divide(matrixB, n, 0, n/2);
        float ** b21 = divide(matrixB, n, n/2, 0);
        float ** b22 = divide(matrixB, n, n/2, n/2);

        //Recursive call for Divide and Conquer
        float** p1 = strassensMultRec(a11,subMatrix(b12,b22,n/2),n/2);
        float** p2 = strassensMultRec(addMatrix(a11,a12,n/2),b22,n/2);
        float** p3 = strassensMultRec(addMatrix(a21,a22,n/2),b11,n/2);
        float** p4 = strassensMultRec(a22,subMatrix(b21,b11,n/2),n/2);
        float** p5 = strassensMultRec(addMatrix(a11,a22,n/2),addMatrix(b11,b22,n/2),n/2);
        float** p6 = strassensMultRec(subMatrix(a12,a22,n/2),addMatrix(b21,b22,n/2),n/2);
        float** p7 = strassensMultRec(subMatrix(a21,a11,n/2),addMatrix(b11,b12,n/2),n/2);
        free(a11); free(a12); free(a21); free(a22);
        free(b11); free(b12); free(b21); free(b22);

        float** c11 = addMatrix(subMatrix(addMatrix(p5,p4,n/2),p2,n/2),p6,n/2);
        float** c12 = addMatrix(p1,p2,n/2);
        float** c21 = addMatrix(p3,p4,n/2);
        float** c22 = addMatrix(subMatrix(addMatrix(p5,p1,n/2),p3,n/2),p7,n/2);
        free(p5); free(p3); free(p1); free(p4);
        free(p2); free(p7); free(p6);

        //Compose the matrix
        compose(c11,result,0,0,n/2);
        compose(c12,result,0,n/2,n/2);
        compose(c21,result,n/2,0,n/2);
        compose(c22,result,n/2,n/2,n/2);

        free(c11); free(c12); free(c21); free(c22);
    }
    else {
        //This is the terminating condition for recurssion.
        //result[0][0]=matrixA[0][0]*matrixB[0][0];
        result = standardMultiplication(matrixA,matrixB, n, 1);
    }

    finalResult = result;
}

/*
* This method combines the matrix in the result matrix
*/
void compose(float** matrix,float** result,int row,int col,int n){
    int i,j,r=row,c=col;
    for(i=0;i<n;i++){
        c=col;
        for(j=0;j<n;j++){
            result[r][c]=matrix[i][j];
            c++;
        }
        r++;
    }
}

/*
* Sub-divide the matrix according to row and col specified
*/
float** divide(float ** matrix,int n, int row,int col) {
    int n_new=n/2;

    float ** array = createZeroMatrix(n_new);
    int i,j,r=row,c=col;
    for(i = 0;i < n_new; i++) {
        c=col;
        for(j = 0; j < n_new; j++) {
            array[i][j] = matrix[r][c];
            c++;
        }
        r++;
    }
    return array;
}

/*
* Add the two input matrix
*/
float** addMatrix(float** matrixA,float** matrixB,int n){
    float ** res = createZeroMatrix(n);
    int i,j;
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            res[i][j]=matrixA[i][j]+matrixB[i][j];

    return res;
}

/*
* Substract the two matrix
*/
float** subMatrix(float** matrixA,float** matrixB,int n){
    float ** res = createZeroMatrix(n);
    int i,j;
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            res[i][j]=matrixA[i][j]-matrixB[i][j];

    return res;
}

