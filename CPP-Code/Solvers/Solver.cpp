#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream> 
#include <string>
#include <unistd.h>
#include <vector>
#include<sys/wait.h>
#include <chrono>  // for high_resolution_clock
#include <signal.h> // catch Ctrl+C command
using namespace std;

void my_handler(int s){
           printf("Caught signal %d\n",s);
           exit(1); 

}

class SolverTasks {
public:
	vector<float> p0_vals, a_vals, n_vals;
	// Contains all tasks that need to be done as {n, p0, a}
	vector<vector<float>> distr_total;
	// Contains the distributed tasks
	vector<vector<vector<float>>> distr_threads;
	
	SolverTasks(vector<float> p0_vals, vector<float> a_vals, float n_min, float n_max, float dn) {
		p0_vals = p0_vals;
		a_vals = a_vals;
		for (float n=n_min; n<=n_max; n+=dn) {
			n_vals.push_back(n);
		}
		// Fill contents of total distribution
		for (size_t i=0; i<p0_vals.size(); ++i) {
			for (size_t j=0; j<a_vals.size(); ++j) {
				for (size_t k=0; k<n_vals.size(); ++k) {
					distr_total.push_back({n_vals[k], p0_vals[i], a_vals[j]});
				}
			}
		}
	}
	
	void distribute_tasks(int threads) {
		for (int i=0; i<threads; ++i) {
			vector<vector<float>> subdistr;
			for (size_t j=0; j<distr_total.size(); ++j) {
				if (j*threads+i < distr_total.size()) {
					subdistr.push_back(distr_total[j*threads+i]);
				} else {
					break;
				}
			}
			distr_threads.push_back(subdistr);
		}
	}
};

class DiffEqSolver {
public:
	// Declare initial parameters given by EOS
	float gamma, A, p_init;
	
	// Define how the increment will be returned. We use the variable ODE
	float * increment(float r, float u, float p, float dr) {
		static float * res;
		return [&](float r, float u, float p, float dr){
			return [&](float du1, float dp1){
				return [&](float du2, float dp2) {
					return [&](float du3, float dp3) {
						return [&](float du4, float dp4) {
							static float res[2];
							res[0] = (du1+2*du2+2*du3+du4)/6;
							res[1] = (dp1+2*dp2+2*dp3+dp4)/6;
							return res;
						}(ODE(r+dr,u+du3,p+dp3)[0]*dr,ODE(r+dr,u+du3,p+dp3)[1]*dr);
					}(ODE(r+dr/2,u+du2/2,p+dp2/2)[0]*dr,ODE(r+dr/2,u+du2/2,p+dp2/2)[1]*dr);
				}(ODE(r+dr/2,u+du1/2,p+dp1/2)[0]*dr,ODE(r+dr/2,u+du1/2,p+dp1/2)[1]*dr);
			}(ODE(r,u,p)[0]*dr,ODE(r,u,p)[1]*dr);
		}(r, u, p, dr);
	}
	
	// Constructor function for the DiffEqSolver class
	DiffEqSolver(float g = 1.3333, float a = 1, float p0 = 1) {
		// Get values or default values for parameters for EOS
		gamma = g;
		A = a;
		p_init = p0;
	}
	
	// Define EOS needed later in TOV equation
	float eos(float p, float r) {
		float ret;
		if (p > 0) {
			return A*exp(1/gamma*log(p));
		} else {
			return 0;
		}
	}
	
	// ODE function
	float * ODE(float r, float u, float p){
		static float res[2];
		float rho = eos(p,r);
		res[0] = 4*3.1415*rho*r*r;
		if (r > 0) {
			res[1] = -(p+rho)*(4*3.1415*p*r*r+u)/(1-2*u/r);
		} else {
			res[1] = 0;
		}
		return res;
	}
	
	// Just for debugging
	float* solveODE() {
		float r = 0;
		float u = 0;
		float p = p_init;
		float dr = 0.0005;
		float r_max=10000.0;
// 		int steps = (r_max-r)/dr;
// 			float u_series[steps+1];
// 			float p_series[steps+1];
		float* res;
		static float ret[5];
// 			u_series[0] = u;
// 			p_series[0] = p;
		for (int n = 0; r < r_max; r += dr) {
			res = increment(r, u, p, dr);
			if ( p+res[1] <= 0.0 ) {
				break;
			};
			u = u+res[0];
			p = p+res[1];
// 				u_series[n+1]=u;
// 				p_series[n+1]=p;
		}
		ret[0]=r;
		ret[1]=u;
		ret[2]=p;
		ret[3] = isnan(u) + isnan(p) + r>r_max-2*dr;
		ret[4] = dr;
		return ret;
	}
	
	void writeZeroVal(float* ret, float n, float a, float p0, std::ofstream& myfile) {
		string output = "{'r'=" + std::to_string(ret[0]) + ",'u'=" + std::to_string(ret[1]) + ",'p'=" + std::to_string(ret[2]) + ",'n'=" + std::to_string(n) + ",'a'=" + std::to_string(a) + ",'p0'=" + std::to_string(p0) + ",'dr'=" + std::to_string(ret[4]) + ",'solved'=" + std::to_string(ret[3]<1) + "}\n";
		std::cout << output;
		myfile << output;
	}
	
	int calculateZeroVals() {
		// Define initial parameters
		float dn = 0.01;
		float n_min = 0.01;
		float n_max = 5;
		float* ret;
		
		// Write to a file
		std::ofstream myfile;
		myfile.open ("results.txt");
		
		int threads = 6;
		int processes[threads];
		
		vector<float> p0_vals = {0.005,0.01,0.1, 0.2, 0.4, 0.8, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192};
		vector<float> a_vals = {0.005,0.01,0.1, 0.2, 0.4, 0.8, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192};
		
		SolverTasks tasks(p0_vals, a_vals, n_min, n_max, dn);
		tasks.distribute_tasks(threads);
		
		int stat;
		pid_t pid[threads];
		for (int j=0; j<threads; j++) {
			if ((pid[j] = fork())==0) {
				vector<vector<float>> subdistr = tasks.distr_threads[j];
				for ( size_t h=0; h<subdistr.size(); ++h) {
					float n = subdistr[h][0];
					gamma = 1+1/n;
					p_init = subdistr[h][1];
					A = subdistr[h][2];
					ret = solveODE();
					writeZeroVal(ret, n, A , p_init, myfile);
				}
				exit(j);
			}
		}
		
		for (int i=0; i<threads; i++) {
			pid_t cpid = waitpid(pid[i], &stat, 0);
			if (WIFEXITED(stat)) {
				cout << "Subprocess "<< WEXITSTATUS(stat) << " terminated with PID=" << cpid << "\n";
			}
		}
		myfile.close();
		return 0;
	}
};

int main() {
	
	struct sigaction sigIntHandler;

	sigIntHandler.sa_handler = my_handler;
	sigemptyset(&sigIntHandler.sa_mask);
	sigIntHandler.sa_flags = 0;

	sigaction(SIGINT, &sigIntHandler, NULL);

	// Record start time
	auto start = std::chrono::high_resolution_clock::now();
	
	// Initialise class for Differential Equation solver
	DiffEqSolver DiffSolver;
	DiffSolver.calculateZeroVals();
	
	// Record end time
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout << "Elapsed time: " << elapsed.count() << " s\n";
}
