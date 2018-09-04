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
   
   int sigma1, sigma2;
   cin >> sigma1 >> sigma2;
   sigma(0,2) = sigma1; sigma(1,2) = sigma2;

   MatrixXcd H_spa = construct_h0() - U_prime/2*matrixelement_sigmaz(sigma);
   VectorXd hf_eivals= Eigenvalues(H_spa);
   vector <double> hf_vector_eivals(hf_eivals.data(), hf_eivals.data()+hf_eivals.size()); 
   cout << "SPA:" << endl <<  H_spa << endl << endl << hf_eivals.transpose() << endl;// << Eigenvectors(H_spa) << endl << endl;

   U = Eigenvectors(H_spa);
   if(sigma1==1 && sigma2==1)
   {
      MatrixXcd temp_U = U;
      U.col(1) = U.col(2);
      U.col(2) = temp_U.col(1);
   }
   cout << U << endl << endl;

   MatrixXcd H_tda = MatrixXcd::Zero(L*L,L*L);
   for(int it1=0; it1<L*L; it1++)
   {
      for(int it2=0; it2<L*L; it2++)
         H_tda(it1,it2) = matrix_elem(mi(it1).first,mi(it1).second, mi(it2).first, mi(it2).second);
   }

   cout << "TDA \n" << H_tda << endl<< endl;
   VectorXd tda_eivals = Eigenvalues(H_tda);
   double E_HF = hf_eivals(0)+hf_eivals(1);
   tda_eivals.array() += E_HF;
   cout << tda_eivals.transpose() << endl;
   vector <double> tda_vector_eivals(tda_eivals.data(), tda_eivals.data()+tda_eivals.size());   
   tda_vector_eivals.push_back(E_HF);

   ofstream Zout("Z_1-1.dat");
   ofstream Fout("F_1-1.dat");

   for(double temperature = 10.00; temperature >= 0.001; temperature-=0.001)
   {
      double Z_tda=0.00; double Z_hf = 1.00;
      for(int it=0; it< tda_vector_eivals.size(); it++)
      {
         Z_tda += exp(-tda_vector_eivals.at(it)/temperature);
      }
      double mu = get_mu(temperature,hf_vector_eivals);
      for(int it=0; it< hf_eivals.size(); it++)
      {
         Z_hf *= (1+exp(-(hf_eivals(it)-mu)/temperature));
      }
      Zout << temperature << " " << Z_hf << " " << Z_tda << " " << mu << endl;
      Fout << temperature << " " << -temperature*log(Z_hf) << " " << -temperature*log(Z_tda) << endl;
   }
   Zout.close();
   Fout.close();
//    cout << endl << Eigenvectors(H_tda) << endl << endl;

}
