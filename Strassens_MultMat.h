//
// Created by Fernando Cores Prado on 4/12/23.
//

#ifndef MULTMAT_SEC3_STRASSENS_MULTMAT_H
#define MULTMAT_SEC3_STRASSENS_MULTMAT_H


extern double elapsed_str;

extern struct ThreadArgs {
    float ** matrixA;
    float ** matrixB;
    int n;
    float ** result;
};
typedef struct ThreadArgs ThreadArgs,*PtrArgs;

// Functions Prototypes
float** strassensMultiplication(float **, float **,int,int);
float** standardMultiplication(float **,float **,int,int);
void * strassensMultRec(float **, float**,int n,float**);
float** divide(float ** matrixA,int n, int row,int col);
float** addMatrix(float**,float**,int);
float** subMatrix(float**,float**,int);
void compose(float**,float**,int,int,int);


#endif //MULTMAT_SEC3_STRASSENS_MULTMAT_H
