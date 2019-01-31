#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "./sparseMatrix/sparseMatrix.h"

#define PI atan(1.0)*4.0

double Func_fx(double x, double y)
{
    return -4*PI*PI*(2*cos(2*PI*x)-1)*sin(2*PI*y)+x*x;
}

double Func_fy(double x, double y)
{
    return 4*PI*PI*(2*cos(2*PI*y)-1)*sin(2*PI*x);
}

double Func_true_ux(double x, double y)
{
    return (1-cos(2*PI*x))*sin(2*PI*y);
}

double Func_true_uy(double x, double y)
{
    return -(1-cos(2*PI*y))*sin(2*PI*x);
}

double Func_true_p(double x, double y)
{
    return x*x*x/3.0-1.0/12;
}

double Func_true_dux(double x,double y)
{
    return 2*PI*(1-cos(2*PI*x))*cos(2*PI*y);
}

double Func_true_duy(double x, double y)
{
    return -2*PI*(1-cos(2*PI*y))*cos(2*PI*x);
}

void get_rhs(int N, double *rhs)
{
    double h=1.0/N;
    int disp=0;
    for(int j=0;j<N;j++)
    {
        for(int i=0;i<N-1;i++)
        {
            double tmp=0;
            //Boundary Condition
            if(j==0){
                tmp+=-1.0/h*Func_true_dux(i*h+h,0);
            }
            if(j==N-1){
                tmp+=1.0/h*Func_true_dux(i*h+h,1);
            }
            tmp+=Func_fx((i+1)*h,(j+0.5)*h);
            rhs[disp+j*(N-1)+i]=tmp;
        }
    }
    disp+=N*N-N;
    for(int j=0; j<N-1;j++)
    {
        for(int i=0;i<N;i++)
        {
            double tmp=0;
            if(i==0){
                tmp+=-1.0/h*Func_true_duy(0,(j+1)*h);
            }
            if(i==N-1){
                tmp+=1.0/h*Func_true_duy(1,(j+1)*h);
            }
            tmp+=Func_fy((i+0.5)*h,(j+1)*h);
            rhs[disp+j*N+i]=tmp;
        }
    }
    disp+=N*N-N;
    for(int j=0; j<N;j++)
    {
        for(int i=0;i<N;i++)
        {
            rhs[disp+j*N+i]=0;
        }
    }
}

void get_exact_solution(int N, double *exact)
{
    double h=1.0/N;
    int disp=0;
    for(int j=0;j<N;j++)
    {
        for(int i=0;i<N-1;i++)
        {
            exact[disp+j*(N-1)+i]=Func_true_ux((i+1)*h,(j+0.5)*h);
        }
    }
    disp+=N*N-N;
    for(int j=0;j<N-1;j++)
    {
        for(int i=0;i<N;i++)
        {
            exact[disp+j*N+i]=Func_true_uy((i+0.5)*h,(j+1)*h);
        }
    }
    disp+=N*N-N;
    for(int j=0;j<N;j++)
    {
        for(int i=0;i<N;i++)
        {
            exact[disp+j*N+i]=Func_true_p((i+0.5)*h,(j+0.5)*h);
        }
    }
}

double get_norm2(int N, double *vec)
{
    int dof=2*N*(N-1)+N*N;
    double res=0;
    for(int i=0;i<dof;i++)
    {
        res+=vec[i]*vec[i];
    }
    res=1.0/N*sqrt(res);
    return res;
}

double get_normMax(int N, double *vec)
{
    int dof=2*N*(N-1)+N*N;
    double res=fabs(vec[0]);
    for(int i=1;i<dof;i++)
    {
        double tmp=fabs(vec[i]);
        if(res<tmp){
            res=tmp;
        }
    }
    return res;
}

void output_vec(int N, double* vec)
{
    int dof=2*N*(N-1)+N*N;
    for(int i=0;i<dof;i++)
    {
        printf("%d : %lf\n",i,vec[i]);
    }
}

void get_error(int N, double* solution,double* err)
{
    int dof=2*N*(N-1)+N*N;
    double *tmp=(double *)malloc(dof*sizeof(double));
    get_exact_solution(N,tmp);
    for(int i=0;i<dof;i++)
    {
        err[i]=fabs(tmp[i]-solution[i]);
    }
}

void get_residual(int N,
                  double *solution,
                  double *rhs,
                  double *residual)
{
    double h=1.0/N;
    int tmpN=N*N-N;
    int disp=0;
    //residual ux
    for(int j=0;j<N;j++)
    {
        for(int i=0;i<N-1;i++)
        {
            double tmp=0;
            //-Laplace ux
            int index=disp+j*(N-1)+i;
            if(j==0||j==N-1){
                tmp+=solution[index]*3.0;
            }
            else{
                tmp+=solution[index]*4.0;
            }
            if(i>0){
                tmp-=solution[index-1];
            }
            if(i<N-2){
                tmp-=solution[index+1];
            }
            if(j>0){
                tmp-=solution[index-(N-1)];
            }
            if(j<N-1){
                tmp-=solution[index+(N-1)];
            }
            tmp=1.0/h/h*tmp;
            //dp/dx
            int index_p=j*N+i+2*tmpN;
            tmp+=(solution[index_p+1]-solution[index_p])/h;
            residual[index]=rhs[index]-tmp;
        }
    }
    //residual uy
    disp+=tmpN;
    for(int j=0;j<N-1;j++)
    {
        for(int i=0;i<N;i++)
        {
            double tmp=0;
            //-Laplace v
            int index=disp+j*N+i;
            if(i==0||i==N-1){
                tmp+=solution[index]*3.0;
            }
            else{
                tmp+=solution[index]*4.0;
            }
            if(i>0){
                tmp-=solution[index-1];
            }
            if(i<N-1){
                tmp-=solution[index+1];
            }
            if(j>0){
                tmp-=solution[index-N];
            }
            if(j<N-2){
                tmp-=solution[index+N];
            }
            tmp=1.0/h/h*tmp;
            //dp/dy
            int index_p=j*N+i+2*tmpN;
            tmp+=(solution[index_p+N]-solution[index_p])/h;
            residual[index]=rhs[index]-tmp;
        }
    }
    //residual p
    disp+=tmpN;
    for(int j=0;j<N;j++)
    {
        for(int i=0;i<N;i++)
        {
            double tmp=0;
            int index=disp+j*N+i;
            if(i>0){
                tmp=tmp-solution[j*(N-1)+i-1];
            }
            if(i<N-1){
                tmp=tmp+solution[j*(N-1)+i];
            }
            if(j>0){
                tmp=tmp-solution[tmpN+(j-1)*N+i];
            }
            if(j<N-1){
                tmp=tmp+solution[tmpN+j*N+i];
            }
            residual[index]=rhs[index]-tmp/h;
        }
    }
}

void test_exact()
{
    int N=4;
    int dof=N*N+2*N*(N-1);
    double *solution=(double *)malloc(dof*sizeof(double));
    double *rhs=(double *)malloc(dof*sizeof(double));
    double *residual=(double *)malloc(dof*sizeof(double));
    get_exact_solution(N,solution);
    output_vec(N,solution);
    get_rhs(N,rhs);
    get_residual(N,solution,rhs,residual);
    printf("%lf\n",get_normMax(N,residual));
}

void DGS(int N,
         double *solution,
         double *rhs)
{
    int tmpN=N*N-N;
    double h=1.0/N;
    double* ux=(double *)malloc(tmpN*sizeof(double));
    double* uy=(double *)malloc(tmpN*sizeof(double));
    double* p=(double *)malloc(N*N*sizeof(double));
    double* rhs_ux=(double *)malloc(tmpN*sizeof(double));
    double* rhs_uy=(double *)malloc(tmpN*sizeof(double));
    double* rhs_p=(double *)malloc(N*N*sizeof(double));
    memcpy(ux,solution,tmpN*sizeof(double));
    memcpy(uy,solution+tmpN,tmpN*sizeof(double));
    memcpy(p,solution+2*tmpN,N*N*sizeof(double));
    memcpy(rhs_ux,rhs,tmpN*sizeof(double));
    memcpy(rhs_uy,rhs+tmpN,tmpN*sizeof(double));
    memcpy(rhs_p,rhs+2*tmpN,N*N*sizeof(double));
    //GS for Ax*ux=fx-Bx*p
    for(int j=0;j<N;j++)
    {
        for(int i=0;i<N-1;i++)
        {
            int index=j*(N-1)+i;
            double tmp=rhs_ux[index];
            int index_p=j*N+i;
            tmp-=(p[index_p+1]-p[index_p])/h;
            tmp=tmp*h*h;
            if(i>0){
                tmp+=ux[index-1];
            }
            if(i<N-2){
                tmp+=ux[index+1];
            }
            if(j>0){
                tmp+=ux[index-(N-1)];
            }
            if(j<N-1){
                tmp+=ux[index+(N-1)];
            }
            if(j==0||j==N-1){
                ux[index]=tmp/3.0;
            }
            else{
                ux[index]=tmp/4.0;
            }
        }
    }
    //GS for Ay*uy=fy-By*p
    for(int j=0;j<N-1;j++)
    {
        for(int i=0;i<N;i++)
        {
            int index=j*N+i;
            double tmp=rhs_uy[index];
            int index_p=j*N+i;
            tmp-=(p[index_p+N]-p[index_p])/h;
            tmp=tmp*h*h;
            if(i>0){
                tmp+=uy[index-1];
            }
            if(i<N-1){
                tmp+=uy[index+1];
            }
            if(j>0){
                tmp+=uy[index-N];
            }
            if(j<N-2){
                tmp+=uy[index+N];
            }
            if(i==0||i==N-1){
                uy[index]=tmp/3.0;
            }
            else{
                uy[index]=tmp/4.0;
            }
        }
    }
    // div u=0
    double div;
    int k;
    for(int j=0;j<N;j++)
    {
        for(int i=0;i<N;i++)
        {
            div=0;
            k=0;
            if(i>0){
                div-=ux[j*(N-1)+i-1];
                k++;
            }
            if(i<N-1){
                div+=ux[j*(N-1)+i];
                k++;
            }
            if(j>0){
                div-=uy[(j-1)*N+i];
                k++;
            }
            if(j<N-1){
                div+=uy[j*N+i];
                k++;
            }
            double r=rhs_p[j*N+i]-div/h;
            double delta=r/k*h;
            // Revise ux,uy,p
            if(i>0){
                ux[j*(N-1)+i-1]-=delta;
                p[j*N+i-1]+=-delta/h;
            }
            if(i<N-1){
                ux[j*(N-1)+i]+=delta;
                p[j*N+i+1]+=-delta/h;
            }
            if(j>0){
                uy[(j-1)*N+i]-=delta;
                p[j*N+i-N]+=-delta/h;
            }
            if(j<N-1){
                uy[j*N+i]+=delta;
                p[j*N+i+N]+=-delta/h;
            }
            p[j*N+i]+=k*delta/h;
        }
    }
    //Output
    memcpy(solution,ux,tmpN*sizeof(double));
    memcpy(solution+tmpN,uy,tmpN*sizeof(double));
    memcpy(solution+2*tmpN,p,N*N*sizeof(double));
    free(ux);free(uy);free(p);
    free(rhs_ux);free(rhs_uy);free(rhs_p);
}

void test_DGS()
{
    int N=4;
    int dof=N*N+2*N*(N-1);
    double *solution=(double *)calloc(dof,sizeof(double));
    double *rhs=(double *)malloc(dof*sizeof(double));
    double *residual=(double *)malloc(dof*sizeof(double));
    get_rhs(N,rhs);
    for(int i=0;i<1;i++)
    {
        DGS(N,solution,rhs);
    }
    output_vec(N,solution);

}

void get_matrix(int N,
                sparse_matrix* A,
                sparse_matrix* B)
{
    int dof_U=2*N*(N-1);
    int dof_P=N*N;
    init_pattern(A,dof_U,dof_U,5);
    init_pattern(B,dof_U,dof_P,2);
    for(int j=0;j<N;j++)
    {
        for(int i=0;i<N-1;i++)
        {
            int index=j*(N-1)+i;
            add_nz_entry(A,index,index);
            if(i>0){
                add_nz_entry(A,index,index-1);
            }
            if(i<N-2){
                add_nz_entry(A,index,index+1);
            }
            if(j>0){
                add_nz_entry(A,index,index-(N-1));
            }
            if(j<N-1){
                add_nz_entry(A,index,index+(N-1));
            }
            add_nz_entry(B,index,j*N+i);
            add_nz_entry(B,index,j*N+i+1);
        }
    }
    int disp=N*(N-1);
    for(int j=0;j<N-1;j++)
    {
        for(int i=0;i<N;i++)
        {
            int index=j*N+i+disp;
            add_nz_entry(A,index,index);
            if(i>0){
                add_nz_entry(A,index,index-1);
            }
            if(i<N-1){
                add_nz_entry(A,index,index+1);
            }
            if(j>0){
                add_nz_entry(A,index,index-N);
            }
            if(j<N-2){
                add_nz_entry(A,index,index+N);
            }
            int index_p=j*N+i;
            add_nz_entry(B,index,index_p);
            add_nz_entry(B,index,index_p+N);
        }
    }
    compress(A);
    compress(B);
    double h=1.0/N;
    double tmp=1.0/h/h;
    for(int j=0;j<N;j++)
    {
        for(int i=0;i<N-1;i++)
        {
            int index=j*(N-1)+i;
            if(j==0||j==N-1){
                add_nz_value(A,index,index,3*tmp);
            }
            else{
                add_nz_value(A,index,index,4*tmp);
            }
            if(i>0){
                add_nz_value(A,index,index-1,-tmp);
            }
            if(i<N-2){
                add_nz_value(A,index,index+1,-tmp);
            }
            if(j>0){
                add_nz_value(A,index,index-(N-1),-tmp);
            }
            if(j<N-1){
                add_nz_value(A,index,index+(N-1),-tmp);
            }
            add_nz_value(B,index,j*N+i,-1/h);
            add_nz_value(B,index,j*N+i+1,1/h);
        }
    }
    for(int j=0;j<N-1;j++)
    {
        for(int i=0;i<N;i++)
        {
            int index=j*N+i+disp;
            if(i==0||i==N-1){
                add_nz_value(A,index,index,3*tmp);
            }
            else{
                add_nz_value(A,index,index,4*tmp);
            }
            if(i>0){
                add_nz_value(A,index,index-1,-tmp);
            }
            if(i<N-1){
                add_nz_value(A,index,index+1,-tmp);
            }
            if(j>0){
                add_nz_value(A,index,index-N,-tmp);
            }
            if(j<N-2){
                add_nz_value(A,index,index+N,-tmp);
            }
            int index_p=j*N+i;
            add_nz_value(B,index,index_p,-1/h);
            add_nz_value(B,index,index_p+N,1/h);
        }
    }
    // output(A,"A");
}

void Uzawa(int N,
           double *solution,
           double *rhs)
{
    int dof_U=2*N*(N-1);
    int dof_P=N*N;
    double* U=(double *)malloc(dof_U*sizeof(double));
    double* P=(double *)malloc(dof_P*sizeof(double));
    double* rhs_U=(double *)malloc(dof_U*sizeof(double));
    double* rhs_P=(double *)malloc(dof_P*sizeof(double));
    memcpy(U,solution,dof_U*sizeof(double));
    memcpy(P,solution+dof_U,dof_P*sizeof(double));
    memcpy(rhs_U,rhs,dof_U*sizeof(double));
    memcpy(rhs_P,rhs+dof_U,dof_P*sizeof(double));
    //Form A,B
    sparse_matrix A,B;
    get_matrix(N,&A,&B);
    //AU=F-BP
    double* tmp=(double *)malloc(dof_U*sizeof(double));
    mv(&B,P,tmp);
    for(int i=0;i<dof_U;i++){
        tmp[i]=rhs_U[i]-tmp[i];
    }
    int flag=cg(&A,U,tmp,1e-7,N*N);
    if(flag==-1){
        printf("CG not converge\n");
        return;
    }
    //p1,p2
    double *p1=(double *)malloc(dof_U*sizeof(double));
    double *p2=(double *)malloc(dof_P*sizeof(double));
    mTv(&B,U,p2);
    for(int i=0;i<dof_P;i++){
        p2[i]=p2[i]-rhs_P[i];
    }
    mv(&B,p2,tmp);
    // flag=cg(&A,p1,tmp,1e-1,N*N);
    // if(flag==1){
    //     printf("CG not converge\n");
    //     return;
    // }
    //alpha
    // double alpha=vv(p2,p2,dof_P);
    // free(tmp);
    // tmp=(double *)malloc(dof_P*sizeof(double));
    // mTv(&B,p1,tmp);
    // if(fabs(alpha)>1e-4){
    //     alpha=alpha/vv(p2,tmp,dof_P);
    // }
    double alpha=1;
    // for(int i=0;i<dof_U;i++){
    //     U[i]=U[i]-alpha*p1[i];
    // }
    for(int i=0;i<dof_P;i++){
        P[i]=P[i]+alpha*p2[i];
    }
    //Output
    memcpy(solution,U,dof_U*sizeof(double));
    memcpy(solution+dof_U,P,dof_P*sizeof(double));
}

void test_Uzawa()
{
    int N=256;
    int dof=N*N+2*N*(N-1);
    double *solution=(double *)calloc(dof,sizeof(double));
    double *rhs=(double *)malloc(dof*sizeof(double));
    double *residual=(double *)malloc(dof*sizeof(double));
    get_rhs(N,rhs);
    for(int i=0;i<5;i++)
    {
        Uzawa(N,solution,rhs);
        get_residual(N,solution,rhs,residual);
    printf("%d : %1.5e\n",i, get_normMax(N,residual));
    }
    
}

void restriction(int N,
                 double *in,
                 double *out)
{
    int n=N/2;
    int tmpN=N*N-N;
    double* ux_new=(double *)malloc(n*(n-1)*sizeof(double));
    double* uy_new=(double *)malloc(n*(n-1)*sizeof(double));
    double* p_new=(double *)malloc(n*n*sizeof(double));
    double* ux_old=(double *)malloc(tmpN*sizeof(double));
    double* uy_old=(double *)malloc(tmpN*sizeof(double));
    double* p_old=(double *)malloc(N*N*sizeof(double));
    memcpy(ux_old,in,tmpN*sizeof(double));
    memcpy(uy_old,in+tmpN,tmpN*sizeof(double));
    memcpy(p_old,in+2*tmpN,N*N*sizeof(double));
    //ux_new
    for(int j=0;j<n;j++)
    {
        for(int i=0;i<n-1;i++)
        {
            int index=j*2*(N-1)+2*i+1;
            double tmp=ux_old[index-1]+2*ux_old[index]+ux_old[index+1];
            index=index+N-1;
            tmp+=ux_old[index-1]+2*ux_old[index]+ux_old[index+1];
            ux_new[j*(n-1)+i]=tmp/8.0;
        }
    }
    //uy_new
    for(int j=0;j<n-1;j++)
    {
        for(int i=0;i<n;i++)
        {
            int index=(2*j+1)*N+2*i;
            double tmp=uy_old[index]*2+uy_old[index+1]*2;
            index=index-N;
            tmp+=uy_old[index]+uy_old[index+1];
            index=index+2*N;
            tmp+=uy_old[index]+uy_old[index+1];
            uy_new[j*n+i]=tmp/8.0;
        }
    }
    //p_new
    for(int j=0;j<n;j++)
    {
        for(int i=0;i<n;i++)
        {
            int index=j*2*N+2*i;
            double tmp=p_old[index]+p_old[index+1]+
                       p_old[index+N]+p_old[index+N+1];
            p_new[j*n+i]=tmp/4.0;
        }
    }
    //Assemble
    int tmp=n*n-n;
    memcpy(out,ux_new,tmp*sizeof(double));
    memcpy(out+tmp,uy_new,tmp*sizeof(double));
    memcpy(out+2*tmp,p_new,n*n*sizeof(double));
    free(ux_new);free(ux_old);
    free(uy_new);free(uy_old);
    free(p_new); free(p_old);
}

void test_restriction()
{
    int N=4;
    int dof=N*N+2*N*(N-1);
    double *in=(double *)malloc(dof*sizeof(double));
    double *out=(double *)malloc(dof*sizeof(double));    
    for(int i=0;i<dof;i++)
    {
        in[i]=i;
    }
    restriction(N,in,out);
    output_vec(N/2,out);
    free(in);free(out);
}

void prolongation(int N,
                  double *in,
                  double *out)
{
    int n=N*2;
    int tmpN=N*N-N;
    double* ux_new=(double *)malloc(n*(n-1)*sizeof(double));
    double* uy_new=(double *)malloc(n*(n-1)*sizeof(double));
    double* p_new=(double *)malloc(n*n*sizeof(double));
    double* ux_old=(double *)malloc(tmpN*sizeof(double));
    double* uy_old=(double *)malloc(tmpN*sizeof(double));
    double* p_old=(double *)malloc(N*N*sizeof(double));
    memcpy(ux_old,in,tmpN*sizeof(double));
    memcpy(uy_old,in+tmpN,tmpN*sizeof(double));
    memcpy(p_old,in+2*tmpN,N*N*sizeof(double));
    //ux_new
    // odd column
    for(int j=0;j<n;j=j+2)
    {
        int index_j=j/2;
        for(int i=1;i<n-1;i=i+2)
        {
            int index_i=(i-1)/2;
            int index=index_j*(N-1)+index_i;
            if(j==0){
                ux_new[i]=ux_old[index_i];
            }
            else{
                ux_new[j*(n-1)+i]=0.75*ux_old[index]+0.25*ux_old[index-(N-1)];
            }
        }
    }
    for(int j=1;j<n;j=j+2)
    {
        int index_j=(j-1)/2;
        for(int i=1;i<n-1;i=i+2)
        {
            int index_i=(i-1)/2;
            int index=index_j*(N-1)+index_i;
            if(j==n-1){
                ux_new[j*(n-1)+i]=ux_old[index];
            }
            else{
                ux_new[j*(n-1)+i]=0.75*ux_old[index]+0.25*ux_old[index+N-1];
            }
        }
    }
    // even column
    for(int j=0;j<n;j++)
    {
        for(int i=0;i<n-1;i=i+2)
        {
            int index=j*(n-1)+i;
            double tmp=0;
            if(i>0){
                tmp+=0.5*ux_new[index-1];
            }
            if(i<n-2){
                tmp+=0.5*ux_new[index+1];
            }
            ux_new[index]=tmp;
        }
    }
    //uy_new
    // odd row
    for(int j=1;j<n-1;j=j+2)
    {
        int index_j=(j-1)/2;
        for(int i=0;i<n;i=i+2)
        {
            int index_i=i/2;
            int index=index_j*N+index_i;
            if(i==0){
                uy_new[j*n+i]=uy_old[index];
            }
            else{
                uy_new[j*n+i]=0.75*uy_old[index]+0.25*uy_old[index-1];
            }
        }
        for(int i=1;i<n;i=i+2)
        {
            int index_i=(i-1)/2;
            int index=index_j*N+index_i;
            if(i==n-1){
                uy_new[j*n+i]=uy_old[index];
            }
            else{
                uy_new[j*n+i]=0.75*uy_old[index]+0.25*uy_old[index+1];
            }
        }
    }
    // even row
    for(int j=0;j<n-1;j=j+2)
    {
        for(int i=0;i<n;i++)
        {
            int index=j*n+i;
            double tmp=0;
            if(j>0){
                tmp+=0.5*uy_new[index-n];
            }
            if(j<n-2){
                tmp+=0.5*uy_new[index+n];
            }
            uy_new[index]=tmp;
        }
    }
    // p_new
    for(int j=0;j<N;j++)
    {
        for(int i=0;i<N;i++)
        {
            int index=j*2*n+2*i;
            double tmp=p_old[j*N+i];
            p_new[index]=tmp;
            p_new[index+1]=tmp;
            p_new[index+n]=tmp;
            p_new[index+n+1]=tmp;
        }
    }
    //Assemble
    int tmp=n*n-n;
    memcpy(out,ux_new,tmp*sizeof(double));
    memcpy(out+tmp,uy_new,tmp*sizeof(double));
    memcpy(out+2*tmp,p_new,n*n*sizeof(double));
    free(ux_new);free(ux_old);
    free(uy_new);free(uy_old);
    free(p_new); free(p_old);
}

void test_prolongation()
{
    int N=2;
    int dof=N*N+2*N*(N-1);
    double *in=(double *)malloc(dof*sizeof(double));
    double *out=(double *)malloc(dof*sizeof(double));    
    for(int i=0;i<dof;i++)
    {
        in[i]=i;
    }
    prolongation(N,in,out);
    output_vec(2*N,out);
    free(in);free(out);
}

void Vcycle(int N,
              double* solution,
              double* rhs, 
              int mu1, int mu2, int flag)
{
    if(N==2){
        //Solving by DGS
        for(int i=0;i<10;i++)
        {
            DGS(N,solution,rhs);
        }
        return;
    }
    void (*smoother)(int,double *,double *);
    if(flag==0){
        smoother=&DGS;
    }
    else{
        smoother=&Uzawa;
    }
    //pre-smooth
    for(int i=0;i<mu1;i++)
    {
        (*smoother)(N,solution,rhs);
    }
    int n=N/2;
    int dof_N=2*N*(N-1)+N*N;
    int dof_n=2*n*(n-1)+n*n;
    double *residual_N=(double *)malloc(dof_N*sizeof(double));
    double *solution_n=(double *)calloc(dof_n,sizeof(double));
    double *residual_n=(double *)malloc(dof_n*sizeof(double));
    double *corrector_N=(double *)malloc(dof_N*sizeof(double));
    get_residual(N,solution,rhs,residual_N);
    restriction(N,residual_N,residual_n);
    Vcycle(n,solution_n,residual_n,mu1,mu2,flag);
    prolongation(n,solution_n,corrector_N);
    for(int i=0;i<dof_N;i++)
    {
        solution[i]+=corrector_N[i];
    }
    //after-smooth
    for(int i=0;i<mu1;i++)
    {
        (*smoother)(N,solution,rhs);
    }
    free(residual_N);free(residual_n);
    free(solution_n);free(corrector_N);
}

void test_multigrid(int N)
{
    int dof=2*N*(N-1)+N*N;
    double* solution=(double *)calloc(dof,sizeof(double));
    double* rhs=(double *)malloc(dof*sizeof(double)); 
    double* residual=(double *)calloc(dof,sizeof(double));
    double* error=(double *)calloc(dof,sizeof(double));
    get_rhs(N,rhs);   
    int mu1=2;
    int mu2=2;
    int flag=0;// Flag for DGS or Uzawa
    int nIter=1;
    printf("N : %d\n",N);
    printf("Num of V-Cycles\t: Residual Norm\n");
    time_t t0=clock();
    for(int i=0;i<10;i++)
    {
        Vcycle(N,solution,rhs,mu1,mu2,flag);
        get_residual(N,solution,rhs,residual);
        printf("NO.%d\t\t: %1.5e\n",nIter,get_norm2(N,residual));
        nIter++;
    }
    get_error(N,solution,error);
    printf("Error: %1.5e\n",get_norm2(N,error));
    printf("Time : %1.5e s\n",(double)(clock()-t0)/CLOCKS_PER_SEC); 
    free(solution);free(rhs);free(residual);free(error);  
}

int main()
{
    // test_exact();
    // test_DGS();
    // test_Uzawa();
    // test_restriction();
    // test_prolongation();
    test_multigrid(512);
    return 0;
}