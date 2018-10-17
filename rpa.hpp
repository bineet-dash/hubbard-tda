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

// cd Viljk(int i, int l, int j, int k, double T, VectorXd spa_eivals, ostream& out=debugout)
// {
//   int N = hf.size();
//   cd result = 0.0;
//   hf = spa_eivals;
//   temperature = T;

//   for(int s=0; s<L; s++)
//   {
//     for(int x1=0; x1<N; x1++)
//     {
//       for(int x2=0; x2<N; x2++)
//       {
//         for(int x3=0; x3<N; x3++)
//         {
//           for(int x4=0; x4<N; x4++)
//           {
            
//             double wicks_part = ne(j,l)*(ane(i,x1)*ane(x2,x3)*ane(x4,k) + ane(i,x1)*ane(x2,k)*ne(x3,x4)+ ane(i,x3)*ne(x1,x2)*ane(x4,k)
//                                         -ane(i,x3)*ne(x1,x4)*ane(x2,k) + ane(i,k)*ne(x1,x2)*ne(x3,x4) + ane(i,k)*ne(x1,x4)*ane(x2,x3))
//             +ne(j,x4)*(-ane(i,k)*ne(x1,x2)*ne(x3,l)-ane(i,x1)*ane(x2,k)*ne(x3,l)+ane(i,x3)*ane(x2,k)*ne(x1,l)-ane(i,k)*ane(x2,x3)*ne(x1,l))
//             +ne(j,x2)*(ane(i,x1)*ne(x3,l)*ane(x4,k)+ane(i,k)*ne(x1,x4)*ne(x3,l)-ane(i,x3)*ne(x1,l)*ane(x4,k)-ane(i,k)*ne(x1,l)*ne(x3,x4));
            
//             // if(x1==x2==x3==x4==1)
//             out << x1 << x2 << x3 << x4 << endl;
//             out <<  " " <<ne(j,l)*(ane(i,x1)*ane(x2,x3)*ane(x4,k)) <<  " " <<ne(j,l)*ane(i,x1)*ane(x2,k)*ne(x3,x4) <<  " " <<ne(j,l)*ane(i,x3)*ne(x1,x2)*ane(x4,k)
//                 <<  " " <<- ne(j,l)*ane(i,x3)*ne(x1,x4)*ane(x2,k) <<  " " << ne(j,l)*ane(i,k)*ne(x1,x2)*ne(x3,x4) <<  " " <<ne(j,l)*ane(i,k)*ne(x1,x4)*ane(x2,x3)
//                 <<  " " <<ne(j,x4)*(-ane(i,k)*ne(x1,x2)*ne(x3,l)) <<  " " <<-ne(j,x4)*ane(i,x1)*ane(x2,k)*ne(x3,l) <<  " " <<ne(j,x4)*ane(i,x3)*ane(x2,k)*ne(x1,l)
//                 <<  " " <<-ne(j,x4)*ane(i,k)*ane(x2,x3)*ne(x1,l) <<  " " << ne(j,x2)*ane(i,x1)*ne(x3,l)*ane(x4,k) <<  " " <<ne(j,x2)*ane(i,k)*ne(x1,x4)*ne(x3,l) 
//                 <<  " " <<-ne(j,x2)*ane(i,x3)*ne(x1,l)*ane(x4,k) <<  " " <<-ne(j,x2)*ane(i,k)*ne(x1,l)*ne(x3,x4) <<  " " <<endl <<  " " <<endl;

//             if(filter(wicks_part)!= 0.0|| filter(conj(U(s,x1)*U(s+L,x3))*U(s,x2)*U(s+L,x4))!= 0.0)
//             // out << x1 << " " << x2 << " " << x3 << " " << x4 << "  " << filter(wicks_part).real() << " " << filter(conj(U(s,x1)*U(s+L,x3))*U(s,x2)*U(s+L,x4)).real() << endl;
//             result += conj(U(s,x1)*U(s+L,x3))*U(s,x2)*U(s+L,x4)*wicks_part;
//           }
//         }
//       }
//     }
//   }
//   return result;
// }

// cd Viljk(int i, int l, int j, int k, ostream& out=debugout)
// {
//   cd res = 0.0;
//   return res;
// }

// cd Viljk(int i, int l, int j, int k, double T, VectorXd spa_eivals, ostream& out=debugout)
// {
//   cd res = 0.0;
//   for(int s=0; s<L; s++)
//   {
//     res += conj(U(s,i)*U(s+L,l))*U(s,k)*U(s+L,j);
//   }
//   res *= -U_prime;
//   // cout << j << i << l << k << " " << res << " ";
//   return res;
// }

double fermi_fn(double e_minus_mu, double T) {return (isinf(exp(e_minus_mu/T)))? 0: 1/(exp(e_minus_mu/T)+1);}

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

  // for(int j=N/2; j<N; j++)
  // {
  //   for(int i=0; i<N/2; i++)
  //   {
  //     int it1 = i+(j-L)*L+2*L;

  //     for(int l=0; l<N/2; l++)
  //     {
  //       for(int k=N/2; k<N; k++)
  //       {
  //         int it2 = l*L+(k-L);
  //         H_rpa(it1,it2) = rpa_matrix_elem(i, j, k, l, T, spa_eivals_minus_mu);
  //         // /* if(filter(H_rpa(it1,it2))!= 0.0)  */cout << j << i << l << k << " " << it1 << it2 << " " << H_rpa(it1,it2) << endl;          
  //       }
  //     }

  //     for(int l=N/2; l<N; l++)
  //     {
  //       for(int k=0; k<N/2; k++)
  //       {
  //         int it2 = k+(l-L)*L+2*L;
  //         H_rpa(it1,it2) = rpa_matrix_elem(i, j, k, l, T, spa_eivals_minus_mu);
  //         // /* if(filter(H_rpa(it1,it2))!= 0.0)  */cout << j << i << l << k << " " << it1 << it2 << " " << H_rpa(it1,it2) << endl;          
  //       }
  //     }
  //   }
  // }

  return H_rpa;
}

#endif