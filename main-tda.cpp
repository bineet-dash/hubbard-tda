#include "tda.hpp"

double t=1; 
double U_prime=2;
int L=2;
MatrixXd sigma;
MatrixXcd U;

inline pair<int,int> mi(int index)
{
   switch (index)
   {
      case 0:
         return make_pair(0,2);
         break;
      case 1:
         return make_pair(0,3);
         break;
      case 2:
         return make_pair(1,2);
         break;
      case 3:
         return make_pair(1,3);
   
      default:
         exit(3);
         break;
   }
}

int main()
{
   sigma = MatrixXd::Zero(L,3);
   sigma(0,2) = 1; sigma(1,2) = 1;

   MatrixXcd H_spa = construct_h0() - U_prime/2*matrixelement_sigmaz(sigma);
   U = Eigenvectors(H_spa);

   MatrixXcd H_tda = MatrixXcd::Zero(L*L,L*L);
   for(int it1=0; it1<L*L; it1++)
   {
      for(int it2=0; it2<L*L; it2++)
         H_tda(it1,it2) = matrix_elem(mi(it1).first,mi(it1).second, mi(it2).first, mi(it2).second);
   }

   cout << H_tda << endl;
}
