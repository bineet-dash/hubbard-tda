#include "tda.hpp"
#include <cstring>
#include <chrono>
#include <cstdlib>

double t=1; 
double U_prime=4;
int L=8;
double DELTA = 2.0;
MatrixXd sigma;
MatrixXcd U;

using namespace std::chrono;

inline double gs_energy(VectorXd hf_eivals) {return hf_eivals.block(0,0,hf_eivals.size()/2,1).sum();}
inline double gs_energy(vector<double> v) {return accumulate(v.begin(), v.begin()+v.size()/2, 0.00);}

vector <pair<int,int>> select_excitations(vector <double> v, double delta)
{
  auto begin_it = v.begin(); auto end_it = v.end();
  int lower_index = lower_bound(begin_it,end_it,-delta)-begin_it; //cout << endl << v[lower_index] << endl;
  int fermi_index = lower_bound(begin_it,end_it, 0.00)- begin_it; //cout << v[fermi_index] << endl << endl;

  vector <pair<int,int>> excitations;

  for(int index = lower_index; index < fermi_index; index++) 
  {
    int upper_index = upper_bound(begin_it,end_it, v[index]+delta)-begin_it;
    for(int i = fermi_index; i<=upper_index; i++) excitations.push_back(make_pair(index,i));
  }
  return excitations;
}

void greens_sigma_generate(MatrixXd& suggested_sigma, int lattice_index, long & idum)
{
  if(ran0(&idum)<=0.5) suggested_sigma(lattice_index,2) *= -1;
}

int main(int argc, char* argv[])
{
  if(argc!=3) {cerr << "Enter (1) lattice size, (2) no of sweeps.\n"; exit(1);}
  L = atoi(argv[1]);
  int no_sweeps = atoi(argv[2]);
  int N_therm = 0.5*no_sweeps;
  int N_meas = no_sweeps-N_therm;

  int initial_exp = -3;
  int final_exp = 0;
  double final_temp = 10*pow(10,final_exp);
  long idum = time(NULL);

  sigma = MatrixXd::Zero(L,3);
  sigma.col(2) = VectorXd::Constant(L,1);
  for(int i=0; i<L; i++)  greens_sigma_generate(sigma, i, idum);
  MatrixXd suggested_sigma = sigma;
  MatrixXcd H0 = construct_h0(); //will be used umpteen # of times

  MatrixXcd H_spa = H0 - U_prime/2*matrixelement_sigmaz(sigma);
  for(int it=0; it<H_spa.rows(); it++) H_spa(it,it) += ran0(&idum)*0.02-0.01;
  pair<MatrixXcd,vector<double>> spa_spectrum = stdEigenspectrum(H_spa);
  double E_HF = gs_energy(spa_spectrum.second); 
  U = spa_spectrum.first;
  vector <pair<int,int>> s = select_excitations(spa_spectrum.second,DELTA);
  MatrixXcd H_tda = construct_truncated_tda(s, E_HF);
  VectorXd Htda_eivals = Eigenvalues(H_spa);
  double free_energy = tda_free_energy(Htda_eivals,E_HF, final_temp);

  string filename, latticedata;
  latticedata = "_U="+to_string(int(U_prime))+"_size="+to_string(L)+"_sweeps="+to_string(no_sweeps)+"_delta="+to_string(int(DELTA));
  filename="truncated/m_length_trunctda_"+ current_time_str()+latticedata+".txt"; ofstream outfile_mlength;//(filename);
  filename="truncated/trunctda_results_"+current_time_str()+latticedata+".txt"; ofstream outfile_freeenergy;//(filename);
  // cout << "==============================\n"<< "filename is: " << filename << "\n========================\n";

    auto begin_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count();
    cout << begin_ms << endl << endl;


  for(int j=final_exp; j>=initial_exp; j--)
  {
    for(double i=10; i>=2; i-=1)
    {
      double temperature = i*pow(10,j);
      for(int sweep=0; sweep<N_therm; sweep++)
      {
        for(int lattice_index=0; lattice_index<L; lattice_index++)
        {
             cout << duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count()-begin_ms << endl;

          greens_sigma_generate(suggested_sigma,lattice_index, idum);
          MatrixXcd suggested_Hspa = H0-U_prime/2*matrixelement_sigmaz(suggested_sigma);
          for(int it=0; it<H_spa.rows(); it++) suggested_Hspa(it,it) += ran0(&idum)*0.02-0.01;
             cout << duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count()-begin_ms << endl;

          pair<MatrixXcd,vector<double>> suggested_spa_spectrum = stdEigenspectrum(suggested_Hspa);
             cout << duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count()-begin_ms << endl;
  
          double suggested_E_HF = gs_energy(suggested_spa_spectrum.second); 
          MatrixXcd original_U = U;
          U = suggested_spa_spectrum.first;
             cout << duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count()-begin_ms << endl;

          vector <pair<int,int>> s = select_excitations(suggested_spa_spectrum.second,DELTA);
              cout << "size=" << s.size() << endl;
              cout << duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count()-begin_ms << endl;

          MatrixXcd suggested_Htda = construct_truncated_tda(s, E_HF);
              cout << duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count()-begin_ms << endl;

          VectorXd suggested_Htda_eivals = Eigenvalues(suggested_Htda);
          double suggested_free_energy = tda_free_energy(suggested_Htda_eivals,suggested_E_HF, temperature);
              cout << duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count()-begin_ms << endl;

          double uniform_rv = ran0(&idum); double move_prob = exp((free_energy - suggested_free_energy)/temperature);
          if(uniform_rv <= move_prob)
          {
            free_energy = suggested_free_energy;
            sigma = suggested_sigma;
          }
          else
          {
            suggested_sigma=sigma;
            U = original_U;
          }
              cout << duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count()-begin_ms << endl;
              exit(1);

        }
        cout << "\r sweep = " << sweep << " done."; cout.flush();
      }

      double final_free_energy = 0.0;
      double magnetisation = 0.0;
      double S_pi = 0.0;

      for(int sweep= N_therm; sweep<no_sweeps; sweep++)
      {
        for(int lattice_index=0; lattice_index<L; lattice_index++)
        {
          greens_sigma_generate(suggested_sigma,lattice_index, idum);
          MatrixXcd suggested_Hspa = H0-U_prime/2*matrixelement_sigmaz(suggested_sigma);
          for(int it=0; it<H_spa.rows(); it++) suggested_Hspa(it,it) += ran0(&idum)*0.02-0.01;
          pair<MatrixXcd,vector<double>> suggested_spa_spectrum = stdEigenspectrum(suggested_Hspa);

          double suggested_E_HF = gs_energy(suggested_spa_spectrum.second); 
          MatrixXcd original_U = U;
          U = suggested_spa_spectrum.first;
          vector <pair<int,int>> s = select_excitations(suggested_spa_spectrum.second,DELTA);
          MatrixXcd suggested_Htda = construct_truncated_tda(s, E_HF);
          VectorXd suggested_Htda_eivals = Eigenvalues(suggested_Htda);
          double suggested_free_energy = tda_free_energy(suggested_Htda_eivals,suggested_E_HF, temperature);

          double uniform_rv = ran0(&idum); double move_prob = exp((free_energy - suggested_free_energy)/temperature);
          if(uniform_rv <= move_prob)
          {
            free_energy = suggested_free_energy;
            sigma = suggested_sigma;
          }
          else
          {
            suggested_sigma=sigma;
            U = original_U;
          }
        }

        final_free_energy += free_energy; 
        magnetisation += sigma.col(2).sum();

        double sq = 0.0;
        for(int i=0; i<L; i++)
        {
          for(int j=0; j<L; j++)
          {
            sq += sigma(i,2)*sigma(j,2)*pow(-1,i-j)/pow(L,2);
          }
        }
        S_pi += sq;
        cout << "\r sweep = " << sweep << " done."; cout.flush();
      }

      outfile_mlength << temperature <<  " " << sigma.col(2).transpose() << endl;
      outfile_freeenergy << temperature << " " << final_free_energy/double(N_meas) << " " << magnetisation/double(N_meas) << " " << S_pi/double(N_meas) << endl;

      cout << "\rtemperature = " << temperature << " done."; cout.flush();
    }
  }

  cout << endl;
  // end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
  // show_time(begin_ms, end_ms,"MC calculation");
  // spinarrangement_Mathematica_output(sigma,outfile_spinarr);
  // outfile_spinarr.close();

  // outfile_mcdetails.close();
  outfile_mlength.close();
  outfile_freeenergy.close();
  return 0;
}