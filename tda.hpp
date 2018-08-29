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
	// res -= del(i,j)*del(m,n)*(h.sum());
	for(int site=0; site<L; site++)
	{
		res -= 0.25*U_prime*sigma(site,2)*(-del(m,n)*conj(U(site,j-1))*U(site,i-1) + del(i,j)*conj(U(site,m-1))*U(site,n-1)+ del(i,j)*del(m,n)*U.row(site).squaredNorm());
	}
	for(int site=0; site<L; site++)
	{
		res += 0.25*U_prime*sigma(site,2)*(-del(m,n)*conj(U(site+L,j-1))*U(site+L,i-1) + del(i,j)*conj(U(site+L,m-1))*U(site+L,n-1)+ del(i,j)*del(m,n)*U.row(site).squaredNorm());
	}

	cd interaction = 0; 
	for(int site=0; site<L; site++)
	{
		cd temp1 = conj(U(site,m-1))*U.row(site).dot(U.row(site+L))*U(site+L,n-1) - conj(U(site,m-1))*U(site,n-1)*U.row(site+L).squaredNorm() 
						+ conj(U(site+L,m-1))*U(site+L,n-1)*U.row(site).squaredNorm() - conj(U(site+L,m-1))*U(site,n-1)*U.row(site+L).dot(U.row(site))
						+ del(m,n)*(U.row(site).squaredNorm()*U.row(site+L).squaredNorm() - U.row(site+L).dot(U.row(site))*U.row(site).dot(U.row(site+L)));
		
		cd temp2 = -del(m,n)*U.row(site).squaredNorm()*conj(U(site+L,j-1))*U(site+L,i-1)+ conj(U(site,m-1))*U(site,n-1)*conj(U(site+L,j-1))*U(site+L,i-1)
						+ conj(U(site,j-1))*U(site,n-1)*conj(U(site+L,m-1))*U(site+L,i-1) - del(m,n)*conj(U(site,j-1))*U(site+L,i-1)*U.row(site).dot(U.row(site+L));

		cd temp3 = conj(U(site,m-1))*U(site,i-1)*conj(U(site+L,j-1))*U(site+L,i-1) + del(m,n)*U.row(site+L).dot(U.row(site))*U(site,i-1)*conj(U(site+L,j-1))
						-conj(U(site,j-1))*U(site,i-1)*conj(U(site+L,m-1))*U(site+L,n-1) - del(m,n)*conj(U(site,j-1))*U(site,i-1)*U.row(site+L).squaredNorm();
		
		interaction += U_prime*del(i,j)*temp1+temp2+temp3;
	}

	res += interaction;
	return res;
}

#endif

// int main()
// {
//     U = MatrixXcd::Zero(4,4);
//     h = VectorXd::Zero(4);
//     U << -1,0,0,1,
//         0,1,-1,0,
//         1,0,0,1,
//         0,1,1,0;
//     h << 3, -1, 1, 1;
    
//     MatrixXcd H = MatrixXcd::Zero(4,4);

//     H(0,0) = matrix_elem(1,3,1,3);
//     H(0,1) = matrix_elem(1,3,1,4);
//     H(0,2) = matrix_elem(1,3,2,3);
//     H(0,3) = matrix_elem(1,3,2,4);

//     H(1,0) = matrix_elem(1,4,1,3);
//     H(1,1) = matrix_elem(1,4,1,4);
//     H(1,2) = matrix_elem(1,4,2,3);
//     H(1,3) = matrix_elem(1,4,2,4);
    
//     H(2,0) = matrix_elem(2,3,1,3);
//     H(2,1) = matrix_elem(2,3,1,4);
//     H(2,2) = matrix_elem(2,3,2,3);
//     H(2,3) = matrix_elem(2,3,2,4);
    
//     H(3,0) = matrix_elem(2,4,1,3);
//     H(3,1) = matrix_elem(2,4,1,4);
//     H(3,2) = matrix_elem(2,4,2,3);
//     H(3,3) = matrix_elem(2,4,2,4);

//     cout << H.real() << endl;
//     exit(3);

//     ComplexEigenSolver <MatrixXcd> ces;
//     ces.compute(H);
//     cout << ces.eigenvalues()<< endl;

// }