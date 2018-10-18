#ifndef _RPA_HPP_INCLUDED
#define _RPA_HPP_INCLUDED

#include "tda.hpp"

cd rpa_matrix_elem(int i, int j, int k, int l, double T, VectorXd hf)
{
  cd u_iljk = 0.0;
  for(int s=0; s<L; s++) u_iljk += conj(U(s,l)*U(s+L,i)-U(s,i)*U(s+L,l))*(U(s,k)*U(s+L,j)-U(s,j)*U(s+L,k));

  cd elem = del(i,k)*del(j,l)*(hf(i)-hf(j)) + u_iljk*(fermi_fn(hf(k),T)-fermi_fn(hf(l),T));
  return elem;
}

MatrixXcd construct_RPA(VectorXd spa_eivals, double T)
{
  double mu = get_mu(T, spa_eivals);
  VectorXd spa_eivals_minus_mu = spa_eivals.array()-mu;
  int N = spa_eivals.size();
  
  MatrixXcd H_rpa = MatrixXcd::Zero(N*N/2,N*N/2);

  for(int j=0; j<N/2; j++)
  {
    for(int i=N/2; i<N; i++)
    {
      int it1 = j*L+(i-L);
      for(int l=0; l<N/2; l++)
      {
        for(int k=N/2; k<N; k++)
        {
          int it2 = l*L+(k-L); // A_matrix 
          H_rpa(it1,it2) = rpa_matrix_elem(i, j, k, l, T, spa_eivals_minus_mu);
          H_rpa(it1+2*L, it2+2*L) = -conj(H_rpa(it1,it2));
        }
      }
      for(int l=N/2; l<N; l++)
      {
        for(int k=0; k<N/2; k++)
        {
          int it2 = k+(l-L)*L+2*L;
          H_rpa(it1,it2) = rpa_matrix_elem(i, j, k, l, T, spa_eivals_minus_mu);
          H_rpa(it1+2*L, it2-2*L) = -conj(H_rpa(it1,it2));
        }
      }
    }
  }
  return H_rpa;
}

#endif