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
inline pair<int,int> mi(int index){return make_pair(int(index/L), index%L+L);}
inline double Sqr(double x){return x*x;}
inline double gs_energy(VectorXd hf_eivals) {return hf_eivals.block(0,0,hf_eivals.size()/2,1).sum();}
inline double gs_energy(vector<double> v) {return accumulate(v.begin(), v.begin()+v.size()/2, 0.00);}
inline cd filter(cd x){return (abs(x)<1e-3)?0.0:x;}
inline double filter_d(double x) {return (abs(x)<1e-3)?0.0:x;}
inline double fermi_fn(double e_minus_mu, double T) {return (isinf(exp(e_minus_mu/T)))? 0: 1/(exp(e_minus_mu/T)+1);}


#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876
double ran0(long *idum)
{
   long  k;
   double ans;

   *idum ^= MASK;
   k = (*idum)/IQ;
   *idum = IA*(*idum - k*IQ) - IR*k;
   if(*idum < 0) *idum += IM;
   ans=AM*(*idum);
   *idum ^= MASK;
   return ans;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef MASK

bool zheev_cpp(MatrixXcd& A, vector<double>& lambda, char eigenvec_choice='N')
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

bool zgeev_cpp(MatrixXcd& A, vector<double>& lambda, char eigenvec_choice='N')
{  
  int N = A.cols();
  int LDA = A.outerStride();
  int INFO = 0;
  __complex__ double* w = new __complex__ double [N];
  __complex__ double* vl;
  __complex__ double* vr;
  char Nchar = 'N';
  char Vchar = 'V';
  char Uchar = 'U';
  int LWORK = pow(2, N);
  __complex__ double* WORK= new __complex__ double [LWORK];
  double* RWORK = new double [LWORK];
  
  zgeev_(&Nchar, &eigenvec_choice, &N, reinterpret_cast <__complex__ double*> (A.data()), &LDA, w, vl, &LDA, vr, &LDA, WORK, &LWORK, RWORK, &INFO );

  lambda.clear();
  for(int i=0; i<N; i++) lambda.push_back(__real__ w[i]);

  delete[] w; delete[] RWORK; delete[] WORK;
  return INFO==0;
}

vector <double> stdEigenvalues(MatrixXcd A, bool (*diagonalization_routine)(MatrixXcd&, vector <double>&, char)=&zheev_cpp)
{
  std::vector<double> lambda; 
  if(diagonalization_routine(A,lambda,'N')) return lambda;
}

VectorXd Eigenvalues(MatrixXcd A, bool (*diagonalization_routine)(MatrixXcd&, vector <double>&, char)=&zheev_cpp)
{
  std::vector<double> lambda;
  if(diagonalization_routine(A,lambda,'N'))
 	{
		Map<ArrayXd> b(lambda.data(),lambda.size());
  	return b;
	}
}

MatrixXcd Eigenvectors(MatrixXcd A, bool (*diagonalization_routine)(MatrixXcd&, vector <double>&, char)=&zheev_cpp)
{
	std::vector<double> lambda;
  if(diagonalization_routine(A,lambda,'V')) return A; 
}

pair<MatrixXcd, vector<double>> stdEigenspectrum(MatrixXcd A, bool (*diagonalization_routine)(MatrixXcd&, vector <double>&, char)=&zheev_cpp)
{
  std::vector<double> lambda;
  if(diagonalization_routine(A,lambda,'V')) return make_pair(A,lambda);
}

pair<MatrixXcd, VectorXd> Eigenspectrum(MatrixXcd A, bool (*diagonalization_routine)(MatrixXcd&, vector <double>&, char)=&zheev_cpp)
{
  std::vector<double> lambda;
  if(diagonalization_routine(A,lambda,'V'))
 	{
    Map<ArrayXd> b(lambda.data(),lambda.size());
    return make_pair(A,b);
	}
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

cd tda_matrix_elem(int i, int j, int k, int l, double e_hf, double T, vector <double> hf)
{
  cd u_iljk = 0.0;
  for(int s=0; s<L; s++) u_iljk += conj(U(s,l)*U(s+L,i)-U(s,i)*U(s+L,l))*(U(s,k)*U(s+L,j)-U(s,j)*U(s+L,k));

  cd elem = del(i,k)*del(j,l)*(e_hf+hf[i]-hf[j]) + u_iljk*(fermi_fn(hf[k],T)-fermi_fn(hf[l],T));
  return elem;
}

MatrixXcd construct_truncated_tda(vector <pair<int,int>> excitations, double e_hf, double T, vector <double> hf)
{
  int N = excitations.size();

  MatrixXcd H_tda = MatrixXcd::Zero(N, N);
  for(int it1=0; it1<H_tda.rows(); it1++)
  {
    for(int it2=0; it2<H_tda.cols(); it2++)
    { 
      cd hij= tda_matrix_elem(excitations[it1].first, excitations[it1].second, excitations[it2].first, excitations[it2].second, T, e_hf, hf);
      H_tda(it1,it2) = hij; 
    }
  }
  return H_tda;
}


// MatrixXcd construct_tda(double e_hf)
// {
//   MatrixXcd H_tda = MatrixXcd::Zero(L*L,L*L);
//   for(int it1=0; it1<H_tda.rows(); it1++)
//   {
//     for(int it2=0; it2<H_tda.cols(); it2++)
//     cd hij= tda_matrix_elem(excitations[it1].first, excitations[it1].second, excitations[it2].first, excitations[it2].second, T, e_hf, hf);
//   }
//   return H_tda;
// }

double tda_free_energy(VectorXd tda_eivals, double e_hf, double temperature)
{
	double Z_rest = 1.00;
	for(int it=0; it< tda_eivals.size(); it++) Z_rest += exp(-(tda_eivals(it))/temperature);
	return e_hf-temperature*log(Z_rest);
}

double tda_internal_energy(VectorXd tda_eivals, double e_hf, double temperature)
{
  double num = e_hf, denom=1.0 ;
  for(int it=0; it< tda_eivals.size(); it++)
  {
    num += exp(-(tda_eivals(it))/temperature)*((tda_eivals(it))+e_hf);
    denom += exp(-(tda_eivals(it))/temperature);
  }
  return num/denom;
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

double get_mu(double temperature, VectorXd v)
{
  vector<double> stdv (v.data(),v.data()+v.size());
  return get_mu(temperature, stdv);
}


double spa_free_energy(MatrixXcd Mc, double temperature)
{
  std::vector<double> eigenvalues;
  zheev_cpp(Mc, eigenvalues, 'N');
  sort(eigenvalues.begin(),eigenvalues.end());

  double free_energy = 0; double ekt =0;
  double mu = get_mu(temperature, eigenvalues);

  for(auto it=eigenvalues.begin(); it!= eigenvalues.end(); it++)
  {
    ekt = (*it-mu)/temperature;
    if(!isinf(exp(-ekt))) free_energy += -temperature*log(1+exp(-ekt));
    else  free_energy += (*it-mu);
  }
  return free_energy+L*mu;
}

double spa_internal_energy(MatrixXcd Mc, double temperature)
{
  std::vector<double> eigenvalues;
  zheev_cpp(Mc, eigenvalues, 'N');
  sort(eigenvalues.begin(),eigenvalues.end());
  double mu = get_mu(temperature, eigenvalues);

  double internal_energy=0.0 ; double e_min = eigenvalues.front();
  for(auto it=eigenvalues.begin(); it!= eigenvalues.end(); it++)
  {
    internal_energy += (*it)/(exp((*it-mu)/temperature)+1);
  }
  return internal_energy;
}

string current_time_str(void)
{
  time_t rawtime;
  struct tm * timeinfo;
  char buffer[80];
  time (&rawtime);
  timeinfo = localtime(&rawtime);
  strftime(buffer,sizeof(buffer),"%S-%M-%I-%Y-%m-%d",timeinfo);
  string str(buffer);
  return str;
}

void spinarrangement_Mathematica_output(MatrixXd M, ofstream& fout)
{
  double lattice_separation = 1.0;
  // cout << "Enter lattice separation (a): ";
  // cin >> lattice_separation;

  fout << "Show[Graphics3D[{" << endl;
  for(int i=0; i< M.rows(); i++)
  {
    fout << "Arrow[{{" << lattice_separation*i << ", 0, 0}, {" << M(i,0)+lattice_separation*i << ","  << M(i,1) << "," << M(i,2) << "}}]";
    if(i!=M.rows()-1) fout << ",\n";
  }
  fout <<"}] ]" << endl << endl;
}

inline double pf(double e, double temperature) {return 1/(exp(-e/temperature)+1);}

cd fdot(int site, int site_conj, double temperature, VectorXd hf)
{
  cd result = 0;
  for(int i=0; i<hf.size(); i++) result += U(site,i)*conj(U(site_conj,i))*pf(hf(i),temperature);
  return result;
}

#endif



/* double tda_free_energy(VectorXd tda_eivals, double e_hf, double temperature)
{
	// Z = exp(-beta*e_min)*(1+\sum_e exp(-beta*(e-e_min)))
	double Z_rest = 1.00;
	for(int it=0; it< tda_eivals.size(); it++) Z_rest += exp(-(tda_eivals(it)-e_hf)/temperature);
	return e_hf-temperature*log(Z_rest);
} */

/* cd matrix_elem(int i, int m, int j, int n, double e_hf)
{
	cd res = del(i,j)*del(m,n)*(e_hf-L*U_prime*L/4); //change this to 0.25*U*\sum m_i^2
	for(int site=0; site<L; site++)
	{
		res -= 0.25*U_prime*sigma(site,2)*(-del(m,n)*conj(U(site,j))*U(site,i) + del(i,j)*conj(U(site,m))*U(site,n)+ del(i,j)*del(m,n)*U.row(site).squaredNorm());
	}
	for(int site=0; site<L; site++)
	{
		res += 0.25*U_prime*sigma(site,2)*(-del(m,n)*conj(U(site+L,j))*U(site+L,i) + del(i,j)*conj(U(site+L,m))*U(site+L,n)+ del(i,j)*del(m,n)*U.row(site+L).squaredNorm());
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
		
		interaction += U_prime*(del(i,j)*temp1+ temp2 +temp3);
	}

	res += interaction;
	return res;
}
 */