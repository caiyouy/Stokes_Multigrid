/**
 * @file sparseMatrix.h
 * @author Caiyou Yuan
 * @brief 
 * @version 0.1
 * @date 2019-01-12
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#ifndef _MY_SPARSE_H
#define _MY_SPARSE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/// CSR sparse matrix
typedef struct _sparse_matrix
{
    /// number of rows
    int _m;
    /// number of columns
    int _n; 
    /// row start index
    int *rs;    
    /// column index
    int *ci;    
    /// element value
    double *ev; 
} sparse_matrix;

/**
 * @brief initialize pattern
 * 
 * @param A sparse_matrix pointer
 * @param m number of rows
 * @param n number of columns
 * @param K estimated upper-bound for the nonzeros
 */
void init_pattern(sparse_matrix *A,
                  int m,
                  int n,
                  int K);

/**
 * @brief compress redundant space
 * 
 * @param A sparse_matrix pointer 
 */
void compress(sparse_matrix *A);

/**
 * @brief add nonzero entry
 * 
 * @param A sparse_matrix pointer
 * @param i inserted row
 * @param j inserted column
 */
void add_nz_entry(sparse_matrix *A, int i, int j);

/**
 * @brief add nonzero value
 * 
 * @param A sparse_matrix pointer
 * @param i inserted row
 * @param j inserted column
 * @param v inserted value
 */
void add_nz_value(sparse_matrix *A, int i, int j, double v);

/**
 * @brief sparse_matrix multiply a vector 
 * 
 * @param A sparse_matrix pointer
 * @param x inputed vector pointer
 * @param b result vector pointer
 */
void mv(sparse_matrix *A,
        double *x,
        double *b);

/**
 * @brief sparse_matrix's transportation multiply a vector
 * 
 * @param A sparse_matrix pointer
 * @param x inputed vector pointer
 * @param b result vector pointer
 */
void mTv(sparse_matrix *A,
         double *x,
         double *b);

/**
 * @brief GS mixed GS Iteration For Solving Ax=b
 * 
 * @param A sparse_matrix pointer
 * @param x solution vector pointer
 * @param b right hand vector pointer
 * @param tol solution precious parameter
 * @param max_step  maximum step
 * @return int exit flag, 0 for converge, 1 for not converge 
 */
int gs(sparse_matrix *A,
        double *x,
        double *b,
        double tol,
        int max_step);

/**
 * @brief inner-product between vectors
 * 
 * @param x inputed vector pointer
 * @param y inputed vector pointer
 * @param N vecter length
 * @return double result
 */
double vv(double *x,
          double *y,
          int N);

/**
 * @brief CG Iteration For Solving Ax=b
 * 
 * @param A sparse_matrix pointer
 * @param x solution vector pointer
 * @param b right hand vector pointer
 * @param tol solution precious parameter
 * @param max_step  maximum step
 * @return int exit flag, 0 for converge, 1 for not converge 
 */
int cg(sparse_matrix *A,
        double *x,
        double *b,
        double tol,
        int max_step);

/**
 * @brief CG mixed GS Iteration For Solving Ax=b
 * 
 * @param A sparse_matrix pointer
 * @param x solution vector pointer
 * @param b right hand vector pointer
 * @param tol solution precious parameter
 * @param max_step  maximum step
 * @return int exit flag, 0 for converge, 1 for not converge 
 */
int gs_cg(sparse_matrix *A,
           double *x,
           double *b,
           double tol,
           int max_step);

/**
 * @brief Sparse Matrix Multiply
 * 
 * @param A sparse matrix
 * @param B sparse matrix
 * @param C result sparse matrix
 */
void mm(sparse_matrix *A,
        sparse_matrix *B,
        sparse_matrix *C);

/**
 * @brief Output sparse_matrix
 * 
 * @param A sparse matrix
 * @param str file name
 */
void output(sparse_matrix *A, char* str);

#endif 
