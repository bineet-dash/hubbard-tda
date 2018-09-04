#include "tda.hpp"
#include <cstring>

double t=1; 
double U_prime=4;
int L=4;
MatrixXd sigma;
MatrixXcd U;
long idum= -1;

int main(int argc, char* argv[])
{
   if(argc!=2) exit(1);

   sigma = MatrixXd::Zero(L,3);
   ofstream Zout;
   ofstream Fout;

   if(!strcmp(argv[1],"anti"))
   {
      for(int i=0; i<L; i++) sigma(i,2) = pow(-1,i);
      Zout.open("Z_1-1.dat");
      Fout.open("F_1-1.dat");
   }
   else if(!strcmp(argv[1],"ferro"))
   {
      for(int i=0; i<L; i++) sigma(i,2) = 1;//pow(-1,i);
      Zout.open("Z_11.dat");
      Fout.open("F_11.dat");
   }
   else 
   {
      cout << "invalid option" << endl;
      exit(1);
   }

   MatrixXcd H_spa = construct_h0() - U_prime/2*matrixelement_sigmaz(sigma);
   for(int it=0; it<H_spa.rows(); it++) H_spa(it,it) += ran0(&idum)*0.02-0.01;

   VectorXd hf_eivals= Eigenvalues(H_spa);
   vector <double> hf_vector_eivals(hf_eivals.data(), hf_eivals.data()+hf_eivals.size()); 
   // cout << "SPA:" << endl <<  H_spa << endl << endl;
   cout << "SPA eigenvalues: \n" << hf_eivals.transpose() << endl << endl;// << Eigenvectors(H_spa) << endl << endl;

   double E_HF = 0.0;
   for(int level=0; level<L; level++) E_HF += hf_eivals(level);

   U = Eigenvectors(H_spa);
   /* if(sigma(0,2)==1 && sigma(1,2)==1)
   {
      MatrixXcd temp_U = U;
      U.col(1) = U.col(2);
      U.col(2) = temp_U.col(1);
   } */
//    cout << U << endl << endl;

   MatrixXcd H_tda = MatrixXcd::Zero(L*L,L*L);
   for(int it1=0; it1<L*L; it1++)
   {
      for(int it2=0; it2<L*L; it2++)
         H_tda(it1,it2) = matrix_elem(mi(it1).first,mi(it1).second, mi(it2).first, mi(it2).second,E_HF);
   }

   // cout << "TDA \n" << H_tda << endl<< endl;
   VectorXd tda_eivals = Eigenvalues(H_tda);
   tda_eivals.array() += E_HF;
   vector <double> tda_vector_eivals(tda_eivals.data(), tda_eivals.data()+tda_eivals.size());   
   tda_vector_eivals.push_back(E_HF);
   sort(tda_vector_eivals.begin(),tda_vector_eivals.end());
   cout << "TDA eigenvalues" << endl;
   for(auto it=tda_vector_eivals.begin(); it!= tda_vector_eivals.end(); it++) cout << *it << " ";
   cout << endl;

   for(double temperature = 100.00; temperature >= 0.001; temperature-=0.01)
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
