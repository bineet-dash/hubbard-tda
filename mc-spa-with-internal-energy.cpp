#include "tda.hpp"
#include <cstring>
#include <chrono>
#include <cstdlib>

double t=1; 
double U_prime=2;
int L=8;
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
  MatrixXcd Id = MatrixXcd::Identity(H0.rows(), H0.cols());

  MatrixXcd H_spa = H0 - U_prime/2*matrixelement_sigmaz(sigma) + U_prime/4*sigma.unaryExpr(&Sqr).sum()*Id; 
  // for(int it=0; it<H_spa.rows(); it++) H_spa(it,it) += ran0(&idum)*0.02-0.01;
  double internal_energy = spa_internal_energy(H_spa, final_temp);

  string filename, latticedata;
  latticedata = "_U_"+to_string(int(U_prime))+"_size_"+to_string(L)+"_sweeps_"+to_string(no_sweeps);
  // filename="data/spin_arrangement"+current_time_str()+latticedata+".nb"; ofstream outfile_spinarr(filename);
  // spinarrangement_Mathematica_output(sigma,outfile_spinarr);
  filename="spa/m_length_spa_internal_energy"+ current_time_str()+latticedata+".dat"; ofstream outfile_mlength(filename);
  filename="spa/spa_results_internal_energy_"+current_time_str()+latticedata+".dat"; ofstream outfile_results(filename);
  // filename="data/mcdetails"+current_time_str()+latticedata+".txt"; ofstream outfile_mcdetails(filename);
  cout << "==============================\n"<< "filename is: " << filename << "\n========================\n";

  for(int j=final_exp; j>=initial_exp; j--)
  {
    for(double i=9; i>=1; i-=1)
    {
      double temperature = i*pow(10,j);
      for(int sweep=0; sweep<N_therm; sweep++)
      {
        for(int lattice_index=0; lattice_index<L; lattice_index++)
        {
          greens_sigma_generate(suggested_sigma,lattice_index, idum);
          MatrixXcd suggested_Hspa = H0-U_prime/2*matrixelement_sigmaz(suggested_sigma)+U_prime/4*sigma.unaryExpr(&Sqr).sum()*Id;
          // for(int it=0; it<H_spa.rows(); it++) suggested_Hspa(it,it) += ran0(&idum)*0.02-0.01;
          double suggested_internal_energy = spa_internal_energy(suggested_Hspa,temperature);

          double uniform_rv = ran0(&idum); double move_prob = exp((internal_energy - suggested_internal_energy)/temperature);
          if(uniform_rv <= move_prob)
          {
            internal_energy = suggested_internal_energy;
            sigma = suggested_sigma;
          }
          else
          {
            suggested_sigma=sigma;
          }
        }
        cout << "\r sweep = " << sweep << " done."; cout.flush();
      }

      double final_internal_energy = 0.0;
      double magnetisation = 0.0;
      double S_pi = 0.0;
      double internal_energy = 0.0;

      for(int sweep= N_therm; sweep<no_sweeps; sweep++)
      {
        for(int lattice_index=0; lattice_index<L; lattice_index++)
        {
          greens_sigma_generate(suggested_sigma,lattice_index, idum);
          MatrixXcd suggested_Hspa = H0-U_prime/2*matrixelement_sigmaz(suggested_sigma) + U_prime/4*sigma.unaryExpr(&Sqr).sum()*Id;
          // for(int it=0; it<H_spa.rows(); it++) suggested_Hspa(it,it) += ran0(&idum)*0.02-0.01;
          double suggested_internal_energy = spa_internal_energy(suggested_Hspa,temperature);

          double uniform_rv = ran0(&idum); double move_prob = exp((internal_energy - suggested_internal_energy)/temperature);
          if(uniform_rv <= move_prob)
          {
            internal_energy = suggested_internal_energy;
            sigma = suggested_sigma;
          }
          else
          {
            suggested_sigma=sigma;
          }
        }
        
        MatrixXcd H_spa_afterSweep = H0-U_prime/2*matrixelement_sigmaz(suggested_sigma) + U_prime/4*sigma.unaryExpr(&Sqr).sum()*Id;
        internal_energy += spa_internal_energy(H_spa_afterSweep,temperature)/L;
        final_internal_energy += internal_energy/L; 
        magnetisation += sigma.col(2).mean();

        double sq = 0.0;
        for(int i=0; i<L; i++)
        {
          for(int j=0; j<L; j++)
          {
            sq += sigma(i,2)*sigma(j,2)*pow(-1,i-j);
          }
        }
        S_pi += sq/pow(L,2);
        cout << "\r sweep = " << sweep << " done."; cout.flush();
      }

      outfile_mlength << temperature <<  " " << sigma.col(2).transpose() << endl;
      outfile_results << temperature << " " << final_internal_energy/N_meas << " " <<  internal_energy/N_meas
                      << " " << S_pi/N_meas << " " << magnetisation/N_meas  << endl;
                      
      cout << "\rtemperature = " << temperature << " done."; cout.flush();
    }
  }

  ofstream fout("final_eivals.dat");
  MatrixXcd final_Hspa = H0-U_prime/2*matrixelement_sigmaz(sigma)+U_prime/4*sigma.unaryExpr(&Sqr).sum()*Id;
  fout << Eigenvalues(final_Hspa) << endl;

  cout << endl;
  end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
  // show_time(begin_ms, end_ms,"MC calculation");
  // spinarrangement_Mathematica_output(sigma,outfile_spinarr);
  // outfile_spinarr.close();

  // outfile_mcdetails.close();
  outfile_mlength.close();
  outfile_results.close();
  return 0;
}
