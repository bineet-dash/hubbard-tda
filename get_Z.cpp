#include "tda.hpp"
#include <cstring>
#include <chrono>
#include <cstdlib>

double t=1; 
double U_prime=4;
int L=8;
double DELTA;
MatrixXd sigma;
MatrixXcd U;

using namespace std::chrono;

vector <pair<int,int>> select_excitations(vector <double> v, double delta)
{
  auto begin_it = v.begin(); auto end_it = v.end();
  double mean = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
  int lower_index = lower_bound(begin_it,end_it, mean-delta)-begin_it; //cout << endl << v[lower_index] << endl;
  int fermi_index = lower_bound(begin_it,end_it, mean)- begin_it; //cout << v[fermi_index] << endl << endl;

  vector <pair<int,int>> excitations;

  for(int index = lower_index; index < fermi_index; index++) 
  {
    int upper_index = upper_bound(begin_it,end_it, v[index]+delta)-begin_it;
    for(int i = fermi_index; i<upper_index; i++) excitations.push_back(make_pair(index,i));
  }
  return excitations;
}

void greens_sigma_generate(MatrixXd& suggested_sigma, int lattice_index, long & idum)
{
  if(ran0(&idum)<=0.5) suggested_sigma(lattice_index,2) *= -1;
}

double sinhx_by_x(double x, double beta)
{
   double res = sinh(beta*x/2)/x;
   if(!isinf(res) && !isnan(res)) return res;
   else return beta/2;
}

double Z_summand(MatrixXd sigma, double temperature)
{
   MatrixXcd H0 = construct_h0(); //will be used umpteen # of times
   MatrixXcd Id = MatrixXcd::Identity(H0.rows(),H0.cols());

   MatrixXcd H_spa = H0 - U_prime/2*matrixelement_sigmaz(sigma)+ U_prime/4*sigma.unaryExpr(&Sqr).sum()*Id;
   pair<MatrixXcd,vector<double>> spa_spectrum = stdEigenspectrum(H_spa);
   double E_HF = gs_energy(spa_spectrum.second); 
   U = spa_spectrum.first;
   vector <pair<int,int>> s = select_excitations(spa_spectrum.second,DELTA);
   MatrixXcd H_tda = construct_truncated_tda(s, E_HF);

   VectorXd Htda_eivals = Eigenvalues(H_tda); 

   double spa_part = 1.0;
   for(int i=0; i< (spa_spectrum.second).size(); i++) spa_part *= 1+exp(-spa_spectrum.second.at(i)/temperature);

   double num = 1.0;  double denom = 1.0;
   for(int i=0; i< (spa_spectrum.second).size(); i++)
   {
      for(int j=0; j< (spa_spectrum.second).size(); j++)
      {
         if(i==j) continue;
         num *= sinhx_by_x(spa_spectrum.second.at(i)-spa_spectrum.second.at(j),1/temperature);
      }
   }

   for(int i=0; i< Htda_eivals.size(); i++)
   {
      denom *= sinhx_by_x(Htda_eivals(i)-E_HF, 1/temperature);
   }

   double Z_sigma = exp(-U_prime*L/(temperature*4))*spa_part*num/denom;
   return Z_sigma;
}

double F_summand(MatrixXd sigma, double temperature)
{
   MatrixXcd H0 = construct_h0(); //will be used umpteen # of times
   MatrixXcd Id = MatrixXcd::Identity(H0.rows(),H0.cols());

   MatrixXcd H_spa = H0 - U_prime/2*matrixelement_sigmaz(sigma)+ U_prime/4*sigma.unaryExpr(&Sqr).sum()*Id;
   pair<MatrixXcd,vector<double>> spa_spectrum = stdEigenspectrum(H_spa);
   double E_HF = gs_energy(spa_spectrum.second); 
   U = spa_spectrum.first;
   vector <pair<int,int>> s = select_excitations(spa_spectrum.second,DELTA);
   MatrixXcd H_tda = construct_truncated_tda(s, E_HF);

   VectorXd Htda_eivals = Eigenvalues(H_tda); 

   double spa_part = 1.0;
   for(int i=0; i< (spa_spectrum.second).size(); i++) 
   {
      double proposed = 1+exp(-spa_spectrum.second.at(i)/temperature);
      if(!isinf(proposed)) spa_part += -temperature*proposed;
      else spa_part += spa_spectrum.second.at(i);
   }

   double num = 0.0;  
   double denom = -temperature*log(sinhx_by_x(E_HF, 1/temperature));

   for(int i=0; i< (spa_spectrum.second).size(); i++)
   {
      for(int j=0; j< (spa_spectrum.second).size(); j++)
      {
         if(i==j) continue;
         num += -temperature*log(sinhx_by_x(spa_spectrum.second.at(i)-spa_spectrum.second.at(j),1/temperature));
      }
   }
   for(int i=0; i< Htda_eivals.size(); i++)
   {
      denom += -temperature*log(sinhx_by_x(Htda_eivals(i), 1/temperature));
   }
   
   cout.precision(3);
   for(int i=0; i< (spa_spectrum.second).size(); i++) cout << spa_spectrum.second.at(i) << " "; cout << endl;
   cout << E_HF << " " << (Htda_eivals).transpose() << endl;
   cout << spa_part << " " << num << " " << denom << " " << num-denom << "\n\n";

   double F_sigma = 0.25*U_prime*L + spa_part + (num-denom);
   return F_sigma;
}


int main(int argc, char* argv[])
{
  if(argc!=4) {cerr << "Enter (1) lattice size, (2) U and (3) no of sweeps.\n"; exit(1);}
  L = atoi(argv[1]);
  U_prime = atof(argv[2]);
  DELTA = U_prime + 1;
  int no_sweeps = atoi(argv[3]);
  int N_therm = 0.5*no_sweeps;
  int N_meas = no_sweeps-N_therm;

  int initial_exp = -3;
  int final_exp = 0;
  double final_temp = 10*pow(10,final_exp);
  long idum = time(NULL);

  sigma = MatrixXd::Zero(L,3);
  sigma.col(2) = VectorXd::Constant(L,1);

  ofstream dataout("Z_pspa.dat");

  for(int j=final_exp; j>=initial_exp; j--)
  {
   for(double i=9; i>=1; i-=1)
   {
      double temperature = i*pow(10,j);
      double Z = 0.0, F=0.0;

      cout << "temperature =" << temperature << "\n===================\n";

      // for(int site=0; site <L; site++)
      // {
      //    for(int poss_field = -1; poss_field<=1; poss_field+=2)
      //    {
      //       sigma(site,2) = poss_field;
      //       Z += Z_summand(sigma, temperature);
      //       F += F_summand(sigma, temperature);
      //    }
      // }
      for(int site=0; site <L; site++) sigma(site,2) = pow(-1,i);
      

      cout << endl << endl;
            // exit(1);
      dataout << temperature << " " << Z << " " << F << endl;
   }
  }









}