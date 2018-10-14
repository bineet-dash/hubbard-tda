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

double log_sinhx_by_x(double x, double beta)
{
  if(x==0) return log(beta/2);
  else
  {
   double res = sinh(beta*x/2)/x;
   if(!isinf(res) && !isnan(res)) return log(res);
   if(isinf(sinh(beta*x/2))) return beta*abs(x)/2 - log(abs(x)); 
   if(isnan(res)) {cout << "Problem with " << beta << " " << x << endl << endl; exit(1); }
  }
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

  double Z = 0.0, F=0.0;


  for(int site=0; site <L; site++) sigma(site,2) = pow(-1,site);
  MatrixXcd H0 = construct_h0(); //will be used umpteen # of times
  MatrixXcd Id = MatrixXcd::Identity(H0.rows(),H0.cols());

  MatrixXcd H_spa = H0 - U_prime/2*matrixelement_sigmaz(sigma)+ U_prime/4*sigma.unaryExpr(&Sqr).sum()*Id;
  pair<MatrixXcd,vector<double>> spa_spectrum = stdEigenspectrum(H_spa);
  double E_HF = gs_energy(spa_spectrum.second); 
  U = spa_spectrum.first;
  vector <pair<int,int>> s = select_excitations(spa_spectrum.second,DELTA);
  MatrixXcd H_tda = construct_truncated_tda(s, E_HF);

  VectorXd Htda_eivals = Eigenvalues(H_tda); 
  cout << E_HF << " " <<  (Htda_eivals.array()+E_HF).transpose() << endl << endl;
  cout << Eigenvalues(H_spa).transpose() << endl << endl;
  
  ofstream fout("manual_variation_of_PSPA_kernel.dat");

  for(int j=final_exp; j>=initial_exp; j--)
  {
   for(double i=9; i>=1; i-=1)
   {
    double temperature = i*pow(10,j);
    cout << "T= " << temperature << "\n==========================\n";

    double num = 0.0;  
    for(int i=0; i< (spa_spectrum.second).size(); i++)
    {
      for(int j=0; j< (spa_spectrum.second).size(); j++)
      {
        if(i==j) continue;
        num += log_sinhx_by_x(spa_spectrum.second.at(i)-spa_spectrum.second.at(j),1/temperature);
        // num += log(sinh((spa_spectrum.second.at(i)-spa_spectrum.second.at(j))/(2*temperature)));

        cout << "(" << spa_spectrum.second.at(i)-spa_spectrum.second.at(j)<< ", " 
             << log_sinhx_by_x(spa_spectrum.second.at(i)-spa_spectrum.second.at(j),1/temperature) << ") -> ";
      }
    }
    cout << endl << endl;

    double denom = log_sinhx_by_x(E_HF, 1/temperature);
    for(int i=0; i< Htda_eivals.size(); i++)
    {
      denom += log_sinhx_by_x(Htda_eivals(i)+E_HF, 1/temperature);
      // denom += log(sinh((Htda_eivals(i)+E_HF)/(2*temperature)));
      cout << denom << " -> ";
    }
    
    cout << endl << endl << num << " " << denom << " " << num-denom << endl << endl << endl;

    fout << temperature << " " <<  num << " " << denom << " " << num-denom << endl;      
   }  
  }


      cout << endl << endl;
            // exit(1);

}