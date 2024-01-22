//
// Created by Fernando Cores Prado on 4/12/23.
//

#ifndef MULTMAT_SEC3_STANDARD_MULTMAT_H
#define MULTMAT_SEC3_STANDARD_MULTMAT_H

extern double elapsed_std;

extern struct ThreadArgs {
    float ** matrixA;
    float ** matrixB;
    int n;
};
typedef struct ThreadArgs ThreadArgs,*PtrArgs;

float ** standardMultiplication(float ** matrixA,float ** matrixB,int n,int t);
float ** concurrent_standardMultiplication_ikj(float ** matrixA,float ** matrixB,int n,int t);
float ** standardMultiplication_ijk(float ** matrixA,float ** matrixB,int n);
float ** standardMultiplication_ikj(PtrArgs args);

#endif //MULTMAT_SEC3_STANDARD_MULTMAT_H
