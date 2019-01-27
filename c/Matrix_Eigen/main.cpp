#include <iostream>
#include <cmath>
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"
#include "Eigen/IterativeLinearSolvers"
#define PI M_PI

using namespace Eigen;
using namespace std;
typedef SparseMatrix<double,RowMajor> spMat;

double Func_fx(double x, double y)
{
    return -4*PI*PI*(2*cos(2*PI*x)-1)*sin(2*PI*y)
           +x*x;
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

void get_Matrix(int N, spMat* Ax,spMat* Bx,spMat* Ay,spMat* By)
{
    double h=1.0/N;
    // Allocate size
    int tmpN=N*N-N;
    Ax->resize(tmpN,tmpN);
    Ax->reserve(5*tmpN);
    Ay->resize(tmpN,tmpN);
    Ay->reserve(5*tmpN);
    Bx->resize(tmpN,N*N);
    Bx->reserve(2*tmpN);
    By->resize(tmpN,N*N);
    By->reserve(2*tmpN);
    //(-Laplace)ux + dp/dx=fx
    for(int j=0;j<N;j++)
    {
        for(int i=0;i<N-1;i++)
        {
            int index=j*(N-1)+i;
            // Diagonal element
            if(j==0||j==N-1){
                // Ax.insert(index,index)=-4+21.0/23.0;
                Ax->insert(index,index)=-3;
            }
            else{
                Ax->insert(index,index)=-4;
            }
            // Not the down-most row
            if(j>0){
                Ax->insert(index,index-(N-1))=1;
            }
            // Not the up-most row
            if(j<N-1){
                Ax->insert(index,index+N-1)=1;
            }
            // Not the left-most col
            if(i>0){
                Ax->insert(index,index-1)=1;
            }
            // Not the right-most row
            if(i<N-2){
                Ax->insert(index,index+1)=1;
            }
            Bx->insert(index,j*N+i)=-1;
            Bx->insert(index,j*N+i+1)=1;
        }
    }
    *Ax=(*Ax)*(-1.0/h/h);
    *Bx=(*Bx)*(1.0/h);
    // (-Laplace)uy + dp/dy = fy
    for(int j=0; j<N-1;j++)
    {
        for(int i=0;i<N;i++)
        {
            int index=j*N+i;
            if(i==0||i==N-1){
                Ay->insert(index,index)=-3;
            }
            else{
                Ay->insert(index,index)=-4;
            }
            if(j>0){
                Ay->insert(index,index-N)=1;
            }
            if(j<N-2){
                Ay->insert(index,index+N)=1;
            }
            if(i>0){
                Ay->insert(index,index-1)=1;
            }
            if(i<N-1){
                Ay->insert(index,index+1)=1;
            }
            By->insert(index,j*N+i)=-1;
            By->insert(index,j*N+i+N)=1;
        }
    }
    *Ay=(*Ay)*(-1.0/h/h);
    *By=(*By)*(1.0/h);
}

void get_rhs(int N, VectorXd *res)
{
    double h=1.0/N;
    VectorXd fx(N*N-N),fy(N*N-N);
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
            fx(j*(N-1)+i)=tmp;
        }
    }
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
            fy(j*N+i)=tmp;
        }
    }
    //Output
    res->resize(2*(N-1)*N+N*N);
    *res << fx,fy,VectorXd::Zero(N*N);
}

void get_exact_solution(int N, VectorXd *exact)
{
    double h=1.0/N;
    Matrix<double,Dynamic,1> ux_exact(N*N-N),
                             uy_exact(N*N-N),
                             p_exact(N*N);
    for(int j=0;j<N;j++)
    {
        for(int i=0;i<N-1;i++)
        {
            ux_exact(j*(N-1)+i)=Func_true_ux((i+1)*h,(j+0.5)*h);
        }
    }
    for(int j=0;j<N-1;j++)
    {
        for(int i=0;i<N;i++)
        {
            uy_exact(j*N+i)=Func_true_uy((i+0.5)*h,(j+1)*h);
        }
    }
    for(int j=0;j<N;j++)
    {
        for(int i=0;i<N;i++)
        {
            p_exact(j*N+i)=Func_true_p((i+0.5)*h,(j+0.5)*h);
        }
    }
    //Output
    exact->resize(2*(N-1)*N+N*N);
    (*exact)<< ux_exact,uy_exact,p_exact;
}

void get_residual(int N,
                  spMat* Ax,spMat* Bx,spMat* Ay,spMat* By,
                  VectorXd *solution,
                  VectorXd *rhs,
                  VectorXd *res)
{
    int tmpN=N*N-N;
    VectorXd ux(tmpN),uy(tmpN),p(N*N);
    VectorXd fx(tmpN),fy(tmpN),fp(N*N);
    ux=solution->segment(0,tmpN);
    uy=solution->segment(tmpN,tmpN);
    p=solution->segment(2*tmpN,N*N);
    fx=rhs->segment(0,tmpN);
    fy=rhs->segment(tmpN,tmpN);
    fp=rhs->segment(2*tmpN,N*N);
    res->resize(2*N*(N-1)+N*N);
    *res<<fx-(*Ax)*ux-(*Bx)*p,
          fy-(*Ay)*uy-(*By)*p,
          fp-Bx->transpose()*ux-By->transpose()*uy;
    // cout<<*solution<<endl;
    // cout<<*solution<<endl;
    // cout<<res->lpNorm<Infinity>()<<endl;
}

void error(int N, VectorXd *solution)
{
    VectorXd tmp;
    get_exact_solution(N,&tmp);
    tmp=tmp-(*solution);
    // cout << "err : " << tmp <<endl;
    cout <<1.0/N*tmp.norm()<<endl;
}

void test_exact()
{
    int N=4;
    spMat Ax,Bx,Ay,By;
    VectorXd exact,rhs,res;
    get_exact_solution(N,&exact);
    get_rhs(N,&rhs);
    get_Matrix(N,&Ax,&Bx,&Ay,&By);
    get_residual(N,&Ax,&Bx,&Ay,&By,&exact,&rhs,&res);
    cout << Ax <<endl;
    cout << By <<endl;
}

void DGS(int N, spMat* Ax,spMat* Bx,spMat* Ay,spMat* By,VectorXd *rightItem,
         VectorXd* solution)
{
    int tmpN=N*N-N;
    double h=1.0/N;
    VectorXd ux(tmpN),uy(tmpN),p(N*N);
    VectorXd fx(tmpN),fy(tmpN),fp(N*N);
    ux=solution->segment(0,tmpN);
    uy=solution->segment(tmpN,tmpN);
    p=solution->segment(2*tmpN,N*N);
    fx=rightItem->segment(0,tmpN);
    fy=rightItem->segment(tmpN,tmpN);
    fp=rightItem->segment(2*tmpN,N*N);
    //GS for Ax*ux=fx-Bx*p
    VectorXd rhs(tmpN);
    rhs=fx-(*Bx)*p;
    for(int j=0;j<Ax->outerSize();++j)
    {
        double tmp=rhs[j];
        double diag=0;
        for(SparseMatrix<double,RowMajor>::InnerIterator it(*Ax,j);it;++it)
        {
            if(it.col()!=j){
                tmp-=it.value()*ux(it.col());
            }
            else{
                diag=it.value();
            }
        }
        ux(j)=tmp/diag;
    }
    //GS for Ay*uy=fy-By*p
    rhs=fy-(*By)*p;
    for(int j=0;j<Ay->outerSize();++j)
    {
        double tmp=rhs[j];
        double diag=0;
        for(SparseMatrix<double,RowMajor>::InnerIterator it(*Ay,j);it;++it)
        {
            if(it.col()!=j){
                tmp-=it.value()*uy(it.col());
            }
            else{
                diag=it.value();
            }
        }
        uy(j)=tmp/diag;
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
                div-=ux(j*(N-1)+i-1);
                k++;
            }
            if(i<N-1){
                div+=ux(j*(N-1)+i);
                k++;
            }
            if(j>0){
                div-=uy((j-1)*N+i);
                k++;
            }
            if(j<N-1){
                div+=uy(j*N+i);
                k++;
            }
            double delta=-div/k;
            // Revise ux,uy,p
            if(i>0){
                ux(j*(N-1)+i-1)-=delta;
                p(j*N+i-1)+=-delta/h;
            }
            if(i<N-1){
                ux(j*(N-1)+i)+=delta;
                p(j*N+i+1)+=-delta/h;
            }
            if(j>0){
                uy((j-1)*N+i)-=delta;
                p(j*N+i-N)+=-delta/h;
            }
            if(j<N-1){
                uy(j*N+i)+=delta;
                p(j*N+i+N)+=-delta/h;
            }
            p(j*N+i)+=k*delta/h;
        }
    }
    //Output
    *solution<<ux,uy,p;
}

void test_DGS()
{
    int N=32;
    spMat Ax,Bx,Ay,By;
    VectorXd rhs,solution,res;
    solution.resize(2*N*(N-1)+N*N);
    get_Matrix(N,&Ax,&Bx,&Ay,&By);
    get_rhs(N,&rhs);
    for(int i=0;i<3000;i++)
    {
        DGS(N,&Ax,&Bx,&Ay,&By,&rhs,&solution);
        // cout<<i<<endl;
    }
    get_residual(N,&Ax,&Bx,&Ay,&By,&solution,&rhs,&res);
    cout<< res.norm()<<endl;
    error(N,&solution);
}

void Uzawa(int N, spMat* Ax,spMat* Bx,spMat* Ay,spMat* By,VectorXd *rightItem,
           VectorXd* solution)
{
    int tmpN=N*N-N;
    double h=1.0/N;
    VectorXd ux(tmpN),uy(tmpN),p(N*N);
    VectorXd fx(tmpN),fy(tmpN),fp(N*N);
    // ux=solution->segment(0,tmpN);
    // uy=solution->segment(tmpN,tmpN);
    p=solution->segment(2*tmpN,N*N);
    fx=rightItem->segment(0,tmpN);
    fy=rightItem->segment(tmpN,tmpN);
    // fp=rightItem->segment(2*tmpN,N*N);
    //solve Ax*ux=fx-Bx*p
    VectorXd rhs(N*N-N);
    rhs=fx-(*Bx)*p;
    // ConjugateGradient<SparseMatrix<double,RowMajor> > solver;
    SparseLU<SparseMatrix<double,RowMajor> > solver;
    solver.compute(*Ax);
    ux = solver.solve(rhs);
    //solve Ax*ux=fx-Bx*p
    rhs=fy-(*By)*p;
    solver.compute(*Ay);
    uy = solver.solve(rhs);
    p=p+(Bx->transpose()*ux+By->transpose()*uy);
    //Output
    *solution<<ux,uy,p;
}

void test_Uzawa()
{
    int N=4;
    spMat Ax,Bx,Ay,By;
    VectorXd rhs,solution;
    solution.resize(2*N*(N-1)+N*N);
    get_Matrix(N,&Ax,&Bx,&Ay,&By);
    get_rhs(N,&rhs);
    for(int i=0;i<2;i++)
    {
        Uzawa(N,&Ax,&Bx,&Ay,&By,&rhs,&solution);
        VectorXd res;
        get_residual(N,&Ax,&Bx,&Ay,&By,&solution,&rhs,&res);
        cout << res.norm()<<endl;
    }        
}

VectorXd restriction(int N,VectorXd in)
{
    int n=N/2;
    VectorXd ux_new(n*(n-1));
    VectorXd uy_new(n*(n-1));
    VectorXd p_new(n*n);
    VectorXd ux_old(N*(N-1)); ux_old=in.segment(0,N*(N-1));
    VectorXd uy_old(N*(N-1)); uy_old=in.segment(N*(N-1),N*(N-1));
    VectorXd p_old(N*N); p_old=in.segment(2*N*(N-1),N*N);
    //ux_new
    for(int j=0;j<n;j++)
    {
        for(int i=0;i<n-1;i++)
        {
            int index=j*2*(N-1)+2*i+1;
            double tmp=ux_old(index-1)+2*ux_old(index)+ux_old(index+1);
            index=index+N-1;
            tmp+=ux_old(index-1)+2*ux_old(index)+ux_old(index+1);
            ux_new(j*(n-1)+i)=tmp/8.0;
        }
    }
    //uy_new
    for(int j=0;j<n-1;j++)
    {
        for(int i=0;i<n;i++)
        {
            int index=(2*j+1)*N+2*i;
            double tmp=uy_old(index)*2+uy_old(index+1)*2;
            index=index-N;
            tmp+=uy_old(index)+uy_old(index+1);
            index=index+2*N;
            tmp+=uy_old(index)+uy_old(index+1);
            uy_new(j*n+i)=tmp/8.0;
        }
    }
    //p_new
    for(int j=0;j<n;j++)
    {
        for(int i=0;i<n;i++)
        {
            int index=j*2*N+2*i;
            double tmp=p_old(index)+p_old(index+1)+
                       p_old(index+N)+p_old(index+N+1);
            p_new(j*n+i)=tmp/4.0;
        }
    }
    //Assemble
    VectorXd res(2*(n-1)*n+n*n);
    res<<ux_new,uy_new,p_new;
    return res;
}

void test_restriction()
{
    int N=4;
    VectorXd ux(N*(N-1));
    VectorXd uy(N*(N-1));
    VectorXd p(N*N);
    ux=VectorXd::LinSpaced(N*(N-1),0,N*(N-1)-1);
    uy=VectorXd::LinSpaced(N*(N-1),0,N*(N-1)-1);
    p=VectorXd::LinSpaced(N*N,0,N*N-1);
    VectorXd in(2*N*(N-1)+N*N);
    in << ux,uy,p;
    int n=2;
    std::cout << in <<std::endl;
    VectorXd res(2*(n-1)*n+n*n);
    res=restriction(N,in);
    std::cout << "res:" <<std::endl;
    std::cout << res<<std::endl;
}

VectorXd prolongation(int N,VectorXd in)
{
    int n=N*2;
    VectorXd ux_new(n*(n-1));
    VectorXd uy_new(n*(n-1));
    VectorXd p_new(n*n);
    VectorXd ux_old(N*(N-1)); ux_old=in.segment(0,N*(N-1));
    VectorXd uy_old(N*(N-1)); uy_old=in.segment(N*(N-1),N*(N-1));
    VectorXd p_old(N*N); p_old=in.segment(2*N*(N-1),N*N);
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
                ux_new(i)=0.75*ux_old(index_i);
                // ux_new(i)=ux_old(index_i)-h/4*Func_true_dux(index_i*h+h,0);
            }
            else{
                ux_new(j*(n-1)+i)=0.75*ux_old(index)+0.25*ux_old(index-(N-1));
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
                ux_new(j*(n-1)+i)=0.75*ux_old(index);
                // ux_new(j*(n-1)+i)=ux_old(index)+h/4*Func_true_dux(index_i*h+h,1);
            }
            else{
                ux_new(j*(n-1)+i)=0.75*ux_old(index)+0.25*ux_old(index+N-1);
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
                tmp+=0.5*ux_new(index-1);
            }
            if(i<n-2){
                tmp+=0.5*ux_new(index+1);
            }
            ux_new(index)=tmp;
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
                uy_new(j*n+i)=0.75*uy_old(index);
                // uy_new(j*n+i)=uy_old(index)-h/4*Func_true_duy(0,index_j*h+h);
            }
            else{
                uy_new(j*n+i)=0.75*uy_old(index)+0.25*uy_old(index-1);
            }
        }
        for(int i=1;i<n;i=i+2)
        {
            int index_i=(i-1)/2;
            int index=index_j*N+index_i;
            if(i==n-1){
                uy_new(j*n+i)=0.75*uy_old(index);
                // uy_new(j*n+i)=uy_old(index)+h/4*Func_true_duy(1,index_j*h+h);
                // cout<< "Prolongation" << h/4*Func_true_duy(1,index_j*h+h)<<endl;
            }
            else{
                uy_new(j*n+i)=0.75*uy_old(index)+0.25*uy_old(index+1);
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
                tmp+=0.5*uy_new(index-n);
            }
            if(j<n-2){
                tmp+=0.5*uy_new(index+n);
            }
            uy_new(index)=tmp;
        }
    }
    //p_new
    for(int j=0;j<N;j++)
    {
        for(int i=0;i<N;i++)
        {
            int index=j*2*n+2*i;
            double tmp=p_old(j*N+i);
            p_new(index)=tmp;
            p_new(index+1)=tmp;
            p_new(index+n)=tmp;
            p_new(index+n+1)=tmp;
        }
    }
    //Assemble
    VectorXd res(2*(n-1)*n+n*n);
    res<<ux_new,uy_new,p_new;
    return res;
}

void test_prolongation()
{
    int N=4;
    VectorXd ux(N*(N-1));
    VectorXd uy(N*(N-1));
    VectorXd p(N*N);
    ux=VectorXd::LinSpaced(N*(N-1),0,N*(N-1)-1);
    uy=VectorXd::LinSpaced(N*(N-1),0,N*(N-1)-1);
    p=VectorXd::LinSpaced(N*N,0,N*N-1);
    VectorXd in(2*N*(N-1)+N*N);
    in << ux,uy,p;
    int n=8;
    std::cout << in <<std::endl;
    VectorXd res(2*(n-1)*n+n*n);
    res=prolongation(N,in);
    std::cout << "res:" <<std::endl;
    std::cout << res<<std::endl;

}

double Vcycle(int N,
            VectorXd* rhs, VectorXd* solution,
            int mu1, int mu2, int flag)
{
    spMat Ax,Bx,Ay,By;
    VectorXd residual;
    get_Matrix(N,&Ax,&Bx,&Ay,&By);
    if(N==2){
        // direct solver
        Matrix<double,8,8,RowMajor> A;
        Matrix<double,4,4,RowMajor> tmpA;
        Matrix<double,4,4,RowMajor> tmpB;
        tmpA<< 12,-4, 0, 0,
               -4,12, 0, 0,
                0, 0,12,-4,
                0, 0,-4,12;
        tmpB<< -2, 2, 0, 0,
                0, 0,-2, 2,
               -2, 0, 2, 0,
                0,-2, 0, 2;
        A<<tmpA,tmpB,tmpB.transpose(),Matrix<double,4,4,RowMajor>::Zero(4,4);
        *solution=A.partialPivLu().solve(*rhs);
        get_residual(N,&Ax,&Bx,&Ay,&By,solution,rhs,&residual);
        cout << N<< " : " << residual.lpNorm<Infinity>()<<endl;
        return residual.lpNorm<Infinity>();
    }
    void (*smoother)(int, spMat*,spMat*,spMat*,spMat*,VectorXd*,VectorXd* );
    if(flag==0){
        smoother=&DGS;
    }
    else{
        smoother=&Uzawa;
    }
    //pre-smooth
    get_residual(N,&Ax,&Bx,&Ay,&By,solution,rhs,&residual);
    cout << " pre "<< N<< " : " << residual.lpNorm<Infinity>()<<endl;
    VectorXd rhs_new=restriction(N,residual);
    int n=N/2;
    int len=2*n*(n-1)+n*n;
    VectorXd tmp(len);
    tmp=VectorXd::Zero(len);
    Vcycle(n,&rhs_new,&tmp,mu1,mu2,flag);
    VectorXd corrector=prolongation(n,tmp);
    *solution=*solution+corrector;
    for(int i=0;i<mu2;i++)
    {
        (*smoother)(N,&Ax,&Bx,&Ay,&By,rhs,solution);
    }
    get_residual(N,&Ax,&Bx,&Ay,&By,solution,rhs,&residual);
    cout << " after "<< N<< " : " << residual.lpNorm<Infinity>()<<endl;
    return residual.lpNorm<Infinity>();
}

void test_multigrid(int N)
{
    VectorXd rhs;
    get_rhs(N,&rhs);
    VectorXd solution(2*N*(N-1)+N*N);
    int mu1=2;
    int mu2=2;
    int flag=0;
    int nIter=0;
    // while(1)
    // {
    //     double res=Vcycle(N,&rhs,&solution,mu1,mu2,flag);
    //     cout<<nIter<<": residual  "<< res <<endl;
    //     // error(N,&solution);
    //     if(res<1e-6){
    //         break;
    //     }
    //     nIter++;
    // }
    for(int i=0;i<10;i++)
    {
        double res=Vcycle(N,&rhs,&solution,mu1,mu2,flag);
        cout<<"nIter "<< nIter<<": residual  "<< res <<endl;
        nIter++;
    }
    error(N,&solution);   
}

int main(int argc, char const *argv[])
{
    // test_exact();
    // error(N,&solution);
    test_multigrid(128);
    // test_prolongation();
    // test_Uzawa();    
    // test_DGS();
    return 0;
}