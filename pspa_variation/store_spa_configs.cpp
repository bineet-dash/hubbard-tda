#include "../tda.hpp"
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

int get_unique_id(VectorXd v)
{
  int result = 0;
  for(int i=0; i<v.size(); i++)
  {
    if(v(i)==1) result += pow(2,i);
  }
  return result;
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
  double free_energy = spa_free_energy(H_spa, final_temp);

  string filename, latticedata;
  latticedata = "_U="+to_string(int(U_prime))+"_size="+to_string(L)+"_sweeps="+to_string(no_sweeps);
  filename="pspa_variation/configs/m_length_spa_"+latticedata+".dat"; ofstream outfile_mlength(filename);
  ofstream outfile_results;//(filename);
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
          double suggested_free_energy = spa_free_energy(suggested_Hspa,temperature);

          double uniform_rv = ran0(&idum); double move_prob = exp((free_energy - suggested_free_energy)/temperature);
          if(uniform_rv <= move_prob)
          {
            free_energy = suggested_free_energy;
            sigma = suggested_sigma;
          }
          else
          {
            suggested_sigma=sigma;
          }
        }
        cout << "\r sweep = " << sweep << " done."; cout.flush();
      }

      int config_counter = 1;
      for(int sweep= N_therm; sweep<no_sweeps; sweep++)
      {
        for(int lattice_index=0; lattice_index<L; lattice_index++)
        {
          greens_sigma_generate(suggested_sigma,lattice_index, idum);
          MatrixXcd suggested_Hspa = H0-U_prime/2*matrixelement_sigmaz(suggested_sigma) + U_prime/4*sigma.unaryExpr(&Sqr).sum()*Id;
          double suggested_free_energy = spa_free_energy(suggested_Hspa,temperature);

          double uniform_rv = ran0(&idum); double move_prob = exp((free_energy - suggested_free_energy)/temperature);
          if(uniform_rv <= move_prob)
          {
            free_energy = suggested_free_energy;
            sigma = suggested_sigma;
            if(sweep%5==0) 
            {
              outfile_mlength << temperature <<  " " << get_unique_id(sigma.col(2)) << " " << sigma.col(2).transpose() << endl;
              config_counter++; 
            }
          }
          else
          {
            suggested_sigma=sigma;
          }
        }
      
        cout << "\r sweep = " << sweep << " done."; cout.flush();
      }
                      
      cout << "\rtemperature = " << temperature << " done."; cout.flush();
    }
  }

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

