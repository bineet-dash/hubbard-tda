#ifndef _TDA_HPP_INCLUDED_
#define _TDA_HPP_INCLUDED_

#include <iostream>
#include <cmath>
#include <fstream>
#include <lapacke.h>
#include <vector>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

typedef std::complex <double> cd;

extern double t;
extern double U_prime;
extern int L;
extern MatrixXd sigma;
extern MatrixXcd U;

inline double del(int a1, int a2){return (a1==a2)?1:0;}
inline cd jn(cd z){return conj(z);}

bool diagonalize(MatrixXcd& A, vector<double>& lambda, char eigenvec_choice='N')
{
  int N = A.cols();
  int LDA = A.outerStride();
  int INFO = 0;
  double* w = new  double [N];
  char Nchar = 'N';
  char Vchar = 'V';
  char Uchar = 'U';
  int LWORK = int(A.size())*4;
  __complex__ double* WORK= new __complex__ double [LWORK];
  double* RWORK = new double [3*LDA];

  zheev_( &eigenvec_choice, &Uchar, &N, reinterpret_cast <__complex__ double*> (A.data()), &LDA, w, WORK, &LWORK, RWORK, &INFO );

  lambda.clear();
  for(int i=0; i<N; i++) lambda.push_back(w[i]);

  delete[] w; delete[] RWORK; delete[] WORK;
  return INFO==0;
}

VectorXd Eigenvalues(MatrixXcd A)
{
  std::vector<double> lambda; MatrixXcd eigenvec;
  if(diagonalize(A,lambda,'N'))
 	{
		Map<ArrayXd> b(lambda.data(),lambda.size());
  	return b;
	}
}

MatrixXcd Eigenvectors(MatrixXcd A)
{
	std::vector<double> lambda;
  if(diagonalize(A,lambda,'V'))
		return A; 
}

MatrixXcd construct_h0(void)
{
    MatrixXcd Mc = MatrixXcd::Zero(2*L,2*L);
    for(int row=0;row<2*L-1; row++) Mc(row,row+1)=Mc(row+1,row)=-t;
    Mc(L-1,0)=Mc(0,L-1)=-t; //PBC
    Mc(2*L-1,L)=Mc(L,2*L-1)= -t; //PBC
    Mc(L,L-1)= Mc(L-1,L)=0;
    return Mc;
}

MatrixXcd matrixelement_sigmax(MatrixXd randsigma)
{
  MatrixXcd Mcx = MatrixXcd::Zero(2*L,2*L);
  for(int row=0; row<L; row++)
  {
    Mcx(row,row+L) = cd(randsigma(row,0),0);
    Mcx(row+L,row) = cd(randsigma(row,0),0);
  }
  return Mcx;
}

MatrixXcd matrixelement_sigmay(MatrixXd randsigma)
{
	MatrixXcd Mcy = MatrixXcd::Zero(2*L,2*L);
	for(int row=0; row<L; row++)
	{
		Mcy(row,row+L) = cd(0,-randsigma(row,1));
		Mcy(row+L,row) = cd(0,randsigma(row,1));
	}
	return Mcy;
}

MatrixXcd matrixelement_sigmaz(MatrixXd randsigma)
{
	MatrixXcd Mcz = MatrixXcd::Zero(2*L,2*L);
	for(int row=0; row<L; row++)
		Mcz(row,row)= cd(randsigma(row,2),0);
	for(int row=L; row<2*L; row++)
		Mcz(row, row)=cd(-randsigma(row-L,2),0);
	return Mcz;
}


cd matrix_elem(int i, int m, int j, int n)
{
	cd res = 0;
	// // res -= del(i,j)*del(m,n)*(h.sum());
	for(int site=0; site<L; site++)
	{
		res -= 0.25*U_prime*sigma(site,2)*(-del(m,n)*conj(U(site,j))*U(site,i) + del(i,j)*conj(U(site,m))*U(site,n)+ del(i,j)*del(m,n)*U.row(site).squaredNorm());
	}
	for(int site=0; site<L; site++)
	{
		res += 0.25*U_prime*sigma(site,2)*(-del(m,n)*conj(U(site+L,j))*U(site+L,i) + del(i,j)*conj(U(site+L,m))*U(site+L,n)+ del(i,j)*del(m,n)*U.row(site).squaredNorm());
	}

	cd interaction = 0; 
	for(int site=0; site<L; site++)
	{
		cd temp1 = conj(U(site,m))*U.row(site).dot(U.row(site+L))*U(site+L,n) - conj(U(site,m))*U(site,n)*U.row(site+L).squaredNorm() 
						+ conj(U(site+L,m))*U(site+L,n)*U.row(site).squaredNorm() - conj(U(site+L,m))*U(site,n)*U.row(site+L).dot(U.row(site))
						+ del(m,n)*(U.row(site).squaredNorm()*U.row(site+L).squaredNorm() - U.row(site+L).dot(U.row(site))*U.row(site).dot(U.row(site+L)));
		
		cd temp2 = -del(m,n)*U.row(site).squaredNorm()*conj(U(site+L,j))*U(site+L,i)+ conj(U(site,m))*U(site,n)*conj(U(site+L,j))*U(site+L,i)
						+ conj(U(site,j))*U(site,n)*conj(U(site+L,m))*U(site+L,i) - del(m,n)*conj(U(site,j))*U(site+L,i)*U.row(site).dot(U.row(site+L));

		cd temp3 = conj(U(site,m))*U(site,i)*conj(U(site+L,j))*U(site+L,n) + del(m,n)*U.row(site+L).dot(U.row(site))*U(site,i)*conj(U(site+L,j))
						-conj(U(site,j))*U(site,i)*conj(U(site+L,m))*U(site+L,n) - del(m,n)*conj(U(site,j))*U(site,i)*U.row(site+L).squaredNorm();
		
		interaction += U_prime*del(i,j)*temp1+ temp2 +temp3;
	}

	res += interaction;
	return res;
}

double get_mu(double temperature, std::vector<double> v)
{
  sort (v.begin(), v.end());
  double bisection_up_lim = v.back();
  double bisection_low_lim = v.front();

  double mu, no_of_electrons; int count=0;
  double epsilon = 0.000001;

  for(; ;)
  {
    no_of_electrons=0;
    mu = 0.5*(bisection_low_lim+bisection_up_lim);

    for(auto it = v.begin(); it!= v.end(); it++)
    {
      double fermi_func = 1/(exp((*it-mu)/temperature)+1);
      no_of_electrons += fermi_func;
    }
    if(abs(no_of_electrons-L) < epsilon)
    {
      return mu; break;
    }
    else if(no_of_electrons > L+epsilon)
    {
       if(abs(bisection_up_lim-v.front())<0.001){return mu; break;}
       else {bisection_up_lim=mu;}
    }
    else if(no_of_electrons < L-epsilon)
    {bisection_low_lim=mu;}
  }
}

#endif