//
// Created by Fernando Cores Prado on 4/12/23.
//

#ifndef MULTMAT_SEC3_STANDARD_MULTMAT_H
#define MULTMAT_SEC3_STANDARD_MULTMAT_H

extern double elapsed_std;

extern struct StandardArgs {
    float ** matrixA;
    float ** matrixB;
    int n;
    int i;
    int j;
    int cells_n;
};
typedef struct StandardArgs StandardArgs,*PtrStandardArgs;

float ** standardMultiplication(float ** matrixA,float ** matrixB,int n,int t);
void * concurrent_standardMultiplication_ikj(float ** matrixA,float ** matrixB,int n,int t);
float ** standardMultiplication_ijk(float ** matrixA,float ** matrixB,int n);
void * standardMultiplication_ikj(PtrStandardArgs args);

#endif //MULTMAT_SEC3_STANDARD_MULTMAT_H
