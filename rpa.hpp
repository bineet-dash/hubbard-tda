#ifndef _RPA_HPP_INCLUDED
#define _RPA_HPP_INCLUDED

#include "tda.hpp"

VectorXd hf;
double temperature; 

double ne(int i, int j) //normal ordered expectation
{
  return (i==j)? 1/(exp(hf(i)/temperature)+1) : 0;
}

double ane(int i, int j) //antinormal ordered expectation
{
  return (i==j)? 1/(exp(-hf(i)/temperature)+1) : 0;
}

ofstream debugout;

cd Viljk(int i, int l, int j, int k, double T, VectorXd spa_eivals, ostream& out=debugout)
{
  int N = hf.size();
  cd result = 0.0;
  hf = spa_eivals;
  temperature = T;

  for(int s=0; s<L; s++)
  {
    for(int x1=0; x1<N; x1++)
    {
      for(int x2=0; x2<N; x2++)
      {
        for(int x3=0; x3<N; x3++)
        {
          for(int x4=0; x4<N; x4++)
          {
            // cout << x1 << x2 << x3 << x4 << endl;
            double wicks_part = ne(j,l)*(ane(i,x1)*ane(x2,x3)*ane(x4,k) + ane(i,x1)*ane(x2,k)*ne(x3,x4)+ ane(i,x3)*ne(x1,x2)*ane(x4,k)
                                        -ane(i,x3)*ne(x1,x4)*ane(x2,k) + ane(i,k)*ne(x1,x2)*ne(x3,x4) + ane(i,k)*ne(x1,x4)*ane(x2,x3))
            +ne(j,x4)*(-ane(i,k)*ne(x1,x2)*ne(x3,l)-ane(i,x1)*ane(x2,k)*ne(x3,l)+ane(i,x3)*ane(x2,k)*ne(x1,l)-ane(i,k)*ane(x2,x3)*ne(x1,l))
            +ne(j,x2)*(ane(i,x1)*ne(x3,l)*ane(x4,k)+ane(i,k)*ne(x1,x4)*ne(x3,l)-ane(i,x3)*ne(x1,l)*ane(x4,k)-ane(i,k)*ne(x1,l)*ne(x3,x4));

            out << x1 << " " << x2 << " " << x3 << " " << x4 << "  " << filter(wicks_part).real() << " " << filter(conj(U(s,x1)*U(s+L,x3))*U(s,x2)*U(s+L,x4)).real() << endl;

            result += conj(U(s,x1)*U(s+L,x3))*U(s,x2)*U(s+L,x4)*wicks_part;
          }
        }
      }
    }
  }
  return result;
}

double fermi_fn(double e_minus_mu, double T) {return (isinf(e_minus_mu/T))? 0:1/(exp(e_minus_mu/T)+1); }

cd rpa_matrix_elem(int i, int j, int k, int l, double T, VectorXd spa_eivals_minus_mu)
{
  return Viljk(i,l,j,k,T,spa_eivals_minus_mu)*(fermi_fn(k,T)-fermi_fn(l,T));//+del(m,n)*del(i,j)*(spa_eivals_minus_mu(m)-spa_eivals_minus_mu(i));
}

inline pair<int,int> rpa_pair(int index){return make_pair(int(index/(2*L)), index%(2*L));}

MatrixXcd construct_RPA(VectorXd spa_eivals, double T)
{
  double mu = get_mu(T, spa_eivals);
  VectorXd spa_eivals_minus_mu = spa_eivals.array()-mu;
  int N = spa_eivals.size();
  
  MatrixXcd H_rpa = MatrixXcd::Zero(N*N/2,N*N/2);
  // for(int it1=0; it1<H_rpa.rows(); it1++)
  // {
  //   for(int it2=0; it2<H_rpa.cols(); it2++)
  //   {
  //     // cout << it1 << " " << it2 << " started.\n";
  //     H_rpa(it1,it2) = rpa_matrix_elem(rpa_pair(it1).first, rpa_pair(it1).second, rpa_pair(it2).first, rpa_pair(it2).second, T, spa_eivals_minus_mu);
  //     // cout << it1 << " " << it2 << " done.\n";
  //   }
  // }

  for(int j=0; j<N/2; j++)
  {
    for(int i=N/2; i<N; i++)
    {
      for(int l=0; l<N/2; l++)
      {
        for(int k=N/2; k<N; k++)
        {
          int it1 = j*L+(i-L);
          int it2 = l*L+(k-L);
          H_rpa(it1,it2) = rpa_matrix_elem(i, j, k, l, T, spa_eivals_minus_mu);
          // if(filter(H_rpa(it1,it2))!= 0.0) cout << i << " " << j << " " << k << " " << l << " " << H_rpa(it1,it2) << endl;
        }
      }
    }
  }

  for(int j=N/2; j<N; j++)
  {
    for(int i=0; i<N/2; i++)
    {
      for(int l=N/2; l<N; l++)
      {
        for(int k=0; k<N/2; k++)
        {
          int it1 = i*L+(j-L)+2*L;
          int it2 = k*L+(j-L)+2*L;
          H_rpa(it1,it2) = rpa_matrix_elem(i, j, k, l, T, spa_eivals_minus_mu);
          // if(filter(H_rpa(it1,it2))!= 0.0) cout << i << " " << j << " " << k << " " << l << " " << H_rpa(it1,it2) << endl;          
        }
      }
    }
  }

  return H_rpa;
}

#endif