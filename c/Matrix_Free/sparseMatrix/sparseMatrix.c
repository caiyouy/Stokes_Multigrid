/**
 * @file sparseMatrix.c
 * @author Caiyou Yuan
 * @brief CSR format of sparse matrix
 * @version 0.1
 * @date 2019-01-12
 * 
 * @copyright Copyright (c) 2019
 * 
 */


#include "sparseMatrix.h"

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
                  int K)
{
    A->_n = n;
    A->_m = m;
    A->rs = (int *)malloc((m + 1) * sizeof(int));
    for (int i = 0; i <= m; ++i)
    {
        A->rs[i] = i * K;
    }
    A->ci = (int *)malloc((K * m) * sizeof(int));
    for (int i = 0; i < K * m; ++i)
    {
        A->ci[i] = n; // Marks for compressible
    }
}

/**
 * @brief compress redundant space
 * 
 * @param A sparse_matrix pointer 
 */
void compress(sparse_matrix *A)
{
    int nnz = 0;
    int m = A->_m;
    int n = A->_n;
    int K = A->rs[1];
    // Set rs
    for (int i = 0; i < m; ++i)
    {
        for (int j = A->rs[i]; j < A->rs[i + 1]; ++j)
        {
            if (A->ci[j] == n)
            {
                int nnz_old = nnz;
                nnz += j - A->rs[i];
                A->rs[i] = nnz_old; // A->rs[i] revised after nnz computation
                break;
            }
            if (j == A->rs[i + 1] - 1) // not find A->ci[j]==n at last
            {
                int nnz_old = nnz;
                nnz += j + 1 - A->rs[i];
                A->rs[i] = nnz_old;
            }
        }
    }
    A->rs[m] = nnz;
    // Set ci
    int *ci = A->ci; //copy ci address
    A->ci = (int *)malloc(nnz * sizeof(int));
    for (int i = 0; i < m; ++i)
    {
        for (int j = A->rs[i]; j < A->rs[i + 1]; ++j)
        {
            A->ci[j] = ci[i * K + j - A->rs[i]];
        }
    }
    free(ci);
    A->ev = (double *)calloc(nnz, sizeof(double));
}

/**
 * @brief add nonzero entry
 * 
 * @param A sparse_matrix pointer
 * @param i inserted row
 * @param j inserted column
 */
void add_nz_entry(sparse_matrix *A, int i, int j)
{
    int k;
    for (k = A->rs[i]; k < A->rs[i + 1]; ++k)
    {
        if (A->ci[k] == A->_n || A->ci[k] == j) // Empty or already added
        {
            A->ci[k] = j;
            break;
        }
    }
    if (k == A->rs[i + 1]) //not find place
    {
        printf("Add_nz_entry (%d,%d) fail! \n ", i, j);
    }
}

/**
 * @brief add nonzero value
 * 
 * @param A sparse_matrix pointer
 * @param i inserted row
 * @param j inserted column
 * @param v inserted value
 */
void add_nz_value(sparse_matrix *A, int i, int j, double v)
{
    int k;
    for (k = A->rs[i]; k < A->rs[i + 1]; ++k)
    {
        if (A->ci[k] == j)
        {
            A->ev[k] = A->ev[k] + v;
            break;
        }
    }
    if (k == A->rs[i + 1]) // not find
    {
        printf("Add_nz_value (%d,%d) fail! \n", i, j);
    }
}

/**
 * @brief sparse_matrix multiply a vector 
 * 
 * @param A sparse_matrix pointer
 * @param x inputed vector pointer
 * @param b result vector pointer
 */
void mv(sparse_matrix *A,
        double *x,
        double *b)
{
    for (int i = 0; i < A->_m; ++i)
    {
        b[i] = 0; // set component 0 
        for (int j = A->rs[i]; j < A->rs[i + 1]; ++j)
        {
            b[i] += A->ev[j] * x[A->ci[j]]; //accumulate
        }
    }
}

/**
 * @brief sparse_matrix's transportation multiply a vector
 * 
 * @param A sparse_matrix pointer
 * @param x inputed vector pointer
 * @param b result vector pointer
 */
void mTv(sparse_matrix *A,
         double *x,
         double *b)
{
    for (int i = 0; i < A->_n; i++)
    {
        b[i] = 0;
    }
    for (int i = 0; i < A->_m; ++i)
    {
        for (int j = A->rs[i]; j < A->rs[i + 1]; ++j)
        {
            b[A->ci[j]] += A->ev[j] * x[i];
        }
    }
}

/**
 * @brief Gauss-Seidel Iteration For Solving Ax=b
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
       int max_step)
{

    double res;
    int step = 0;
    int N = 100;
    while (1)
    {
        for (int i = 0; i < A->_m; ++i) // For every row
        {
            double ri = b[i];
            for (int j = A->rs[i] + 1; j < A->rs[i + 1]; ++j) // For every nonzero off-diagonal in the row
            {
                ri -= A->ev[j] * x[A->ci[j]];
            }
            x[i] = ri / A->ev[A->rs[i]];
        }
        step = step + 1;
        if (step % 100 == 0)
        {
            res = 0;
            for (int i = 0; i < A->_m; i++)
            {
                double tmp = b[i];
                for (int j = A->rs[i]; j < A->rs[i + 1]; ++j)
                {
                    tmp -= A->ev[j] * x[A->ci[j]];
                }
                res += tmp * tmp;
            }
            res = sqrt(res);
            if (res <= tol || step >= max_step)
                break;
        }
    }
    if (res <= tol){
        return 0;
    }
    else
        return 1;
}

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
          int N)
{
    double res = 0;
    for (int i = 0; i < N; i++)
    {
        res += x[i] * y[i];
    }
    return res;
}

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
       int max_step)
{
    int m = A->_m;
    double rho, mu, eta;
    double *r = (double *)malloc(m * sizeof(double));
    double *z = (double *)malloc(m * sizeof(double));
    double *p = (double *)malloc(m * sizeof(double));
    mv(A, x, z);
    for (int i = 0; i < m; i++) //initialization
    {
        r[i] = b[i] - z[i];
        p[i] = r[i];
        rho += r[i] * r[i];
    }
    int step = 0;
    while (1)
    {
        step++;
        if (sqrt(rho) <= tol || step >= max_step)
                break;
        mv(A, p, z);
        eta = vv(z, p, m);
        eta = rho / eta;
        for (int i = 0; i < m; i++) //update solution x and residual r
        {
            x[i] = x[i] + eta * p[i];
            r[i] = r[i] - eta * z[i];
        }
        double tmp_rho = vv(r, r, m);
        mu = tmp_rho / rho;
        rho = tmp_rho;
        for (int i = 0; i < m; i++) //update direction p
        {
            p[i] = r[i] + mu * p[i];
        }
    }
    if (sqrt(rho) <= tol){
        return 0;
    }
    else{
        printf("%1.2e\n",sqrt(rho));
        return 1;
    }
}

/**
 * @brief Sparse Matrix Multiply
 * 
 * @param A sparse matrix
 * @param B sparse matrix
 * @param C result sparse matrix
 */
void mm(sparse_matrix *A,
        sparse_matrix *B,
        sparse_matrix *C)
{
    // Attempt to multiply A and B (Determine max_nnz in the row)
    int K = 0;
    int *index = (int *)calloc(B->_n, sizeof(int));
    int *flag = (int *)calloc(B->_n, sizeof(int));
    for (int i1 = 0; i1 < A->_m; ++i1)
    {
        int nz = 0; // nnz of i-th row
        for (int j1 = A->rs[i1]; j1 < A->rs[i1 + 1]; ++j1)
        {
            int i2 = A->ci[j1];
            for (int j2 = B->rs[i2]; j2 < B->rs[i2 + 1]; ++j2)
            {
                int k = B->ci[j2];
                if (flag[k] == 0) // have not counted into the nz
                {
                    index[nz++] = k; // keep index for restore
                    flag[k] = 1;
                }
            }
        }
        if (K < nz) // K record the maximum of nnz of all row
        {
            K = nz;
        }
        // Restore index and flag
        int tmp;
        for (int i = 0; i < nz; i++)
        {
            flag[index[i]] = 0;
            index[i] = 0;
        }
    }
    init_pattern(C, A->_m, B->_n, K);
    // Add nz entry
    for (int i1 = 0; i1 < A->_m; i1++)
    {
        for (int j1 = A->rs[i1]; j1 < A->rs[i1 + 1]; j1++)
        {
            int i2 = A->ci[j1];
            for (int j2 = B->rs[i2]; j2 < B->rs[i2 + 1]; j2++)
            {
                int k = B->ci[j2];
                add_nz_entry(C, i1, k);
            }
        }
    }
    //Compress
    compress(C);
    // Add nz value
    for (int i1 = 0; i1 < A->_m; i1++)
    {
        for (int j1 = A->rs[i1]; j1 < A->rs[i1 + 1]; j1++)
        {
            int i2 = A->ci[j1];
            double val1 = A->ev[j1];
            for (int j2 = B->rs[i2]; j2 < B->rs[i2 + 1]; j2++)
            {
                int k = B->ci[j2];
                double val2 = B->ev[j2];
                add_nz_value(C, i1, k, val1 * val2);
            }
        }
    }
}

/**
 * @brief Output sparse_matrix
 * 
 * @param A sparse matrix
 * @param str file name
 */
void output(sparse_matrix *A, char *str)
{
    FILE *fp;
    char name[20];
    sprintf(name, "%s.csv", str);
    fp = fopen(name, "w");
    for (int i = 0; i < A->_m; i++)
    {
        for (int j = A->rs[i]; j < A->rs[i + 1]; j++)
        {
            fprintf(fp, "%d %d %1.12e\n", i, A->ci[j], A->ev[j]);
        }
    }
    fclose(fp);
}
