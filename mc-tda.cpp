#include "tda.hpp"
#include <cstring>
#include <chrono>
#include <cstdlib>

double t=1;
double U_prime=2;
int L=4;
MatrixXd sigma;
MatrixXcd U;

using namespace std::chrono;

void greens_sigma_generate(MatrixXd& suggested_sigma, int lattice_index, long & idum)
{
  if(ran0(&idum)<=0.5) suggested_sigma(lattice_index,2) *= -1;
}

int main(int argc, char* argv[])
{
  if(argc!=4) {cerr << "Enter (1) lattice size, (2) U and (3) no of sweeps.\n"; exit(1);}
  L = atoi(argv[1]);
  U_prime = atof(argv[2]);
  int no_sweeps = atoi(argv[3]);
  int N_therm = 0.5*no_sweeps;
  int N_meas = no_sweeps-N_therm;

  int initial_exp = -3;
  int final_exp = 0;
  double final_temp = 10*pow(10,final_exp);
  milliseconds begin_ms, end_ms;
  long idum = time(NULL);

  sigma = MatrixXd::Zero(L,3);
  sigma.col(2) = VectorXd::Constant(L,1);
  for(int i=0; i<L; i++)  greens_sigma_generate(sigma, i, idum);
  MatrixXd suggested_sigma = sigma;
  MatrixXcd H0 = construct_h0();

  MatrixXcd H_spa = H0 - U_prime/2*matrixelement_sigmaz(sigma);
  for(int it=0; it<H_spa.rows(); it++) H_spa(it,it) += ran0(&idum)*0.02-0.01;

  pair<MatrixXcd,vector<double>> spa_spectrum = stdEigenspectrum(H_spa);
  double E_HF = gs_energy(spa_spectrum.second); 
  U = spa_spectrum.first;
  MatrixXcd H_tda = construct_tda(E_HF);
  VectorXd Htda_eivals = Eigenvalues(H_tda); 
  double free_energy = tda_free_energy(Htda_eivals,E_HF, final_temp);

  string filename, latticedata;
  latticedata = "_U="+to_string(int(U_prime))+"_size="+to_string(L)+"_sweeps="+to_string(no_sweeps);
  // filename="data/spin_arrangement"+current_time_str()+latticedata+".nb"; ofstream outfile_spinarr(filename);
  // spinarrangement_Mathematica_output(sigma,outfile_spinarr);
  filename="data/m_length_tda_"+ current_time_str()+latticedata+".txt"; ofstream outfile_mlength(filename);
  filename="data/tda_results_"+current_time_str()+latticedata+".txt"; ofstream outfile_freeenergy(filename);
  // filename="data/mcdetails"+current_time_str()+latticedata+".txt"; ofstream outfile_mcdetails(filename);
  cout << "==============================\n"<< "filename is: " << filename << "\n========================\n";

  for(int j=final_exp; j>=initial_exp; j--)
  {
    for(double i=10; i>=2; i-=1)
    {
      double temperature = i*pow(10,j);
      for(int sweep=0; sweep<N_therm; sweep++)
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
  
          MatrixXcd suggested_Htda = construct_tda(suggested_E_HF);
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
        // cout << "\r sweep = " << sweep << " done."; cout.flush();
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
  
          MatrixXcd suggested_Htda = construct_tda(suggested_E_HF);
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

      // cout << "\rtemperature = " << temperature << " done."; cout.flush();
    }
  }

  cout << endl;
  end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
  // show_time(begin_ms, end_ms,"MC calculation");
  // spinarrangement_Mathematica_output(sigma,outfile_spinarr);
  // outfile_spinarr.close();

  // outfile_mcdetails.close();
  outfile_mlength.close();
  outfile_freeenergy.close();
  return 0;
}
