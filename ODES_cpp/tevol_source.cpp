/**
 * @file    tevol_source.cpp
 * @brief   ODE system for age-structured COVID-19 model
 *
 * Source code. Compile with (requires gsl):
 * g++ -std=c++11 -O3 -o tevol_source ./tevol_source.cpp $(gsl-config --cflags) $(gsl-config --libs)
 *
 * Based on the model of:
 * https://science.sciencemag.org/content/sci/suppl/2020/04/28/science.abb8001.DC1/abb8001_Zhang_SM.pdf
 *
 * @author  LHD
 * @since   2020-05-12
 */

#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>
#include <sstream>

#include <boost/multi_array.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include "dyn.hpp"

using namespace std;

int main(int argc, const char *argv[]) {

    // Input arguments
    // Argument 1: Path to contact matrix file
    // Argument 2: Base transmission rate (beta)
		 
    // Model parameters
    double beta = atof(argv[2]); //base transmission rate
    double gamma = 1.0/5.1; //recovery rate
    const int dim = 14; //number of age groups

    // Inital condition parameter
    double epsilon = 1e-7;

    // Population structure
    vector<int> N(14);
    // Shanghai
    N[0] = 940622;
    N[1] = 770132;
    N[2] = 549734;
    N[3] = 975165;
    N[4] = 2519060;
    N[5] = 2686112;
    N[6] = 2410875;
    N[7] = 1897672;
    N[8] = 1684530;
    N[9] = 1751202;
    N[10] = 1753252;
    N[11] = 1871916;
    N[12] = 1285068;
    N[13] = 3101661;
    long int totalpop = accumulate(N.begin(),N.end(),0);


    // Susceptibility in age group i
    vector<double> sigma(14);
    sigma[0] = 0.34;
    sigma[1] = 0.34;
    sigma[2] = 0.34;
    sigma[3] = 1.0;
    sigma[4] = 1.0;
    sigma[5] = 1.0;
    sigma[6] = 1.0;
    sigma[7] = 1.0;
    sigma[8] = 1.0;
    sigma[9] = 1.0;
    sigma[10] = 1.0;
    sigma[11] = 1.0;
    sigma[12] = 1.0;
    sigma[13] = 1.44;

    // Setting initial conditions
    typedef boost::multi_array<double,2> mat_type;
    typedef mat_type::index index;
    mat_type y(boost::extents[3][dim]);
    fill(y.data(),y.data()+y.num_elements(),0.0);
    // y[i][x] = susceptible (i=0), infectious (i=1) or recovered (i=2) in age group x
    double lastR = 0.0;
    double lastI = 0.0;
	for(int x=0; x<=13; ++x) {
        y[0][x] = N[x]*(1.0-epsilon);
        y[1][x] = N[x]*(epsilon);
        lastI += y[1][x];
        y[2][x] = 0.0*N[x];
    }

    // Contact matrix
    vector< vector<double> > M(14, vector<double>(14, 0.0));
    ifstream matrix_file(argv[1]); double rate;
    for(int x=0; x<=13; ++x) {
        for(int y=0; y<=13; ++y) {
            matrix_file >> rate;
            M[x][y] = rate;
        }
    }

    // Structure of parameters for solver
    Sparam param = {beta, gamma, sigma, N, M, dim};

    // Integrator parameters
    double t = 0;
    double dt = 1e-4;
    double t_step = 1.0;
    const double eps_abs = 1e-8;
    const double eps_rel = 1e-10;

    // Define GSL odeiv parameters
    const gsl_odeiv_step_type * step_type = gsl_odeiv_step_rkf45;
    gsl_odeiv_step * step = gsl_odeiv_step_alloc (step_type, 3*dim);
    gsl_odeiv_control * control = gsl_odeiv_control_y_new (eps_abs,eps_rel);
    gsl_odeiv_evolve * evolve = gsl_odeiv_evolve_alloc (3*dim);
    gsl_odeiv_system sys = {dydt, NULL, 3*dim, &param};
	
	//Integration
    int status(GSL_SUCCESS);
    double diffI = 1.0;
    double diffR = 1.0;
    for (double t_target = t+t_step; diffR > 1e-14; t_target += t_step ) { //stop by difference
            
         // Integrate ODEs       
         while (t < t_target) {
            status = gsl_odeiv_evolve_apply (evolve,control,step,&sys,&t,t_target,&dt,y.data());
            if (status != GSL_SUCCESS) {
				cout << "SNAFU" << endl;
                break;
			}
        } // end while

        // Calculate and update R(t)
        double newR = 0.0;
        for(int x=0; x<=13; ++x) newR += y[2][x];
        //cout << t << " " << newR << "\n";
        diffR = abs(newR - lastR);
        lastR = newR;

	} //end while
    cout.flush();

    cout << beta << " " << lastR/(1.0*totalpop) << "\n";

    // Free memory
    gsl_odeiv_evolve_free(evolve);
    gsl_odeiv_control_free(control);
    gsl_odeiv_step_free(step);
    
    return 0;
}
