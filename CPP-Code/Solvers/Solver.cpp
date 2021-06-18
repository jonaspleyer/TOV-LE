#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream> 
#include <string>
#include <unistd.h>
#include <vector>
#include <sys/wait.h>
#include <chrono>  // for high_resolution_clock
#include <signal.h> // catch Ctrl+C command
#include <mongocxx/client.hpp> // Mongo database
#include <mongocxx/instance.hpp> // Mongo database
#include <bsoncxx/json.hpp>
#include <bsoncxx/builder/stream/helpers.hpp>
#include <bsoncxx/builder/stream/document.hpp>
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
	float gamma, A, p_init, dr, r_end;
	float dn = 0.25;
	float n_min = 0.01;
	float n_max = 5;
	// Define the values for which solutions should be generated
// 	vector<float> p0_vals = {0.001,0.005,0.01,0.1, 0.2, 0.4, 0.8, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192};
// 	vector<float> a_vals = {0.005,0.01,0.1, 0.2, 0.4, 0.8, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192};
	// For testing purposes
	vector<float> p0_vals = {0.1,1,8};
	vector<float> a_vals = {0.1,1,8};
	
	// Initialise the mongodb
	mongocxx::instance inst{};
	mongocxx::client conn{mongocxx::uri{}};
	mongocxx::database db = conn["TOV-LE"];
	mongocxx::collection coll = db["Exponents-polytropic-EOS-CPP2"];
	
	// Constructor function for the DiffEqSolver class
	DiffEqSolver(float g = 1.3333, float a = 1, float p0 = 1, float r_step = 0.0005, float r_end_pass = 10000.0) {
		// Get values or default values for parameters for EOS
		gamma = g;
		A = a;
		p_init = p0;
		dr = r_step;
		r_end = r_end_pass;
	}
	
	// Define EOS needed later in TOV equation
	float eos(float p, float r) {
		if (p > 0) {
			return A*pow(p, 1/gamma);
		} else {
			return 0;
		}
	}
	
	// Define how the increment will be returned. We use the variable ODE
	float * increment(float r, float u, float p, float dr) {
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
	
	// ODE function
	float * ODE(float r, float u, float p){
		static float res[2];
		float rho = eos(p,r);
		res[0] = 4*3.1415*rho*r*r;
		if (r > 0) {
			res[1] = -(p+rho)*(4*3.1415*p*r*r*r+u/r)/(1-2*u/r);
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
		float* res;
		static float ret[5];
		while (r < r_end) {
			res = increment(r, u, p, dr);
			if ( p+res[1] <= 0.0 ) {
				break;
			};
			u = u+res[0];
			p = p+res[1];
			r += dr;
		}
		ret[0] = r;
		ret[1] = u;
		ret[2] = p;
		ret[3] = (isnan(u) + isnan(p) + r>r_end-2*dr)==0;
		ret[4] = dr;
		ret[5] = r_end;
		return ret;
	}
	
	// Writes a new ZeroVal in the database
	void mongoWriteZeroVal(float* ret, float n, float a, float p0) {
		auto builder = bsoncxx::builder::stream::document{};
		bsoncxx::document::value doc_value = builder
			<< "n" << n
			<< "terms" << 0
			<< "A_init" << a
			<< "p_init" << p0
			<< "r_max" << ret[0]
			<< "equ_type" << "TOV"
			<< "succ" << ret[3]
			<< "details" << bsoncxx::builder::stream::open_document
				<< "r_end" << ret[5]
				<< "dr" << ret[4]
				<< "u" << ret[1]
				<< "p" << ret[2]
			<< bsoncxx::builder::stream::close_document
			<< bsoncxx::builder::stream::finalize;
		bsoncxx::stdx::optional<mongocxx::result::insert_one> result = coll.insert_one(std::move(doc_value));
		if (result) {
			cout << "Inserted new entry\n";
		} else {
			cout << "No new entry\n";
		}
	}
	
	// We want to check if the value we are calculating is already present before we calculate it.
	// Outputs 1 if no result is present, multiple ones are (deletes them all) or new result would be better.
	int checkAndDeleteZeroValIfPresent(float r_end, float dr, float n, float a, float p0, int delete_toggle=1) {
		int check=0;
		// Define the document for which to check in MongoDB
		auto builder = bsoncxx::builder::stream::document{};
		bsoncxx::document::value find_value = builder
			<< "n" << n
			<< "terms" << 0
			<< "A_init" << a
			<< "p_init" << p0
			<< "equ_type" << "TOV"
			<< bsoncxx::builder::stream::finalize;
		
		// Define a view document to actually use it
		bsoncxx::document::view find_view = find_value.view();
		
		// Find all documents in MongoDB that match a solution for these parameters
		mongocxx::cursor cursor = coll.find(find_view);
		
		// See how many elements are matching
		int cursor_len = std::distance(cursor.begin(), cursor.end());
		if (cursor_len>1) {
			// Delete all results if we find many.
			coll.delete_many(find_view);
			check=1;
		} else {
			// Otherwise check if the single result will be improved
			bsoncxx::stdx::optional<bsoncxx::document::value> doc_value = coll.find_one(find_view);
			if (doc_value) {
				bsoncxx::document::view doc_view = (*doc_value).view();
				if (delete_toggle*checkIfNewResultBetter(doc_view, r_end, dr, n, a, p0)>0) {
					cout << "Deleting result to calculate more precise one:\n";
					coll.delete_one(find_view);
					check=1;
				}
			} else {
				check=1;
			}
		}
		return check;
	}
	
	// Returns 1 if the new result would be better. Return 0 otherwise.
	int checkIfNewResultBetter(bsoncxx::document::view doc_view, float r_end, float dr, float n, float a, float p0) {
		// If any of the below conditions are true then we want to recalculate
		int check = 0;
		int check_end = 0;
// 		int check_succ = 0;
		int check_stepsize = 0;
		if (doc_view["r_max"].type() == bsoncxx::type::k_double) {
			float r_end_old = doc_view["details"]["r_end"].get_double();
			check_end = (doc_view["r_max"].get_double()>=r_end_old-dr && r_end_old<r_end);
		}
		if (doc_view["details"].type() == bsoncxx::type::k_document && doc_view["details"]["dr"].type() == bsoncxx::type::k_double) {
			check_stepsize = (doc_view["details"]["dr"].get_double() > dr);
		}
// 		if (doc_view["succ"].type() == bsoncxx::type::k_double) {
// 			check_succ=(doc_view["succ"].get_double()==0)*check_stepsize;
// 		}
		check = check_end + check_stepsize > 0;
		// Return 1 if the new result would improve. Return 0 otherwise
		return check;
	}
	
	// Print and delete if toggle is set.
	void printAll(int delete_toggle=0) {
		mongocxx::cursor cursor = coll.find({});
		int count=0;
		for(auto doc : cursor) {
			std::cout << bsoncxx::to_json(doc) << "\n";
			if (delete_toggle) {
				coll.delete_one(doc);
				cout << "Deleting entry\n";
			}
			count++;
		}
		cout << "Counted " << count << " results.\n";
	}
	
	// Actually calculate all the zero values we are interested in
	int calculateZeroVals() {
		float* ret;
		
		// Define number of threads
		int threads = 15;
		int stat;
		pid_t pid[threads];
		
		// Distribute the tasks with this struct
		SolverTasks tasks(p0_vals, a_vals, n_min, n_max, dn);
		tasks.distribute_tasks(threads);
		
		// Split the jobs over the different threads
		for (int j=0; j<threads; j++) {
			
			// If statements fork the process
			if ((pid[j] = fork())==0) {
				
				// This vector contains all jobs for this thread
				vector<vector<float>> subdistr = tasks.distr_threads[j];
				
				// Iterate through every job list per thread
				for ( size_t h=0; h<subdistr.size(); ++h) {
					// Define the variables for which we want to solve the ODE
					float n = subdistr[h][0];
					gamma = 1+1/n;
					p_init = subdistr[h][1];
					A = subdistr[h][2];
					
					// Solve the ODE if the check says so
					if (checkAndDeleteZeroValIfPresent(r_end, dr, n, A, p_init)) {
						ret = solveODE();
						mongoWriteZeroVal(ret, n, A , p_init);
					}
				}
				// When every result is calculated exit the respective threads
				exit(j);
			}
		}
		// Prints a message after every thread has terminated and waits to continue
		for (int i=0; i<threads; i++) {
			pid_t cpid = waitpid(pid[i], &stat, 0);
			if (WIFEXITED(stat)) {
				cout << "Subprocess "<< WEXITSTATUS(stat) << " terminated with PID=" << cpid << "\n";
			}
		}
		return 0;
	}
};

int main() {
	// Interrupt Handler to stop
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
// 	DiffSolver.printAll();
	
	// Record end time
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout << "Elapsed time: " << elapsed.count() << " s\n";
}
