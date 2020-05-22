#ifndef DYN_HPP_INCLUDED
#define DYN_HPP_INCLUDED

#include <boost/multi_array.hpp>

/**
 * @file    dyn.hpp
 * @brief   ODE system for age-structured COVID-19 model
 *
 * @author  LHD
 * @since   2019-10-29
 */

struct Sparam {
    const double beta;
    const double gamma;
    const std::vector<double> sigma;
    const std::vector<int> N;
    const std::vector< std::vector<double> > M;
    const int dim;
}; // parameter structure

//********** function dydt definition **************************************************************
int dydt(double t, const double y[], double f[], void * param) {
// ODE system for intra-community growth process

    // Cast parameters
    Sparam& p = *static_cast<Sparam* >(param);

    // Create multi_array reference to y and f
    typedef boost::multi_array_ref<const double,2> CSTmatref_type;
    typedef boost::multi_array_ref<double,2> matref_type;
    typedef CSTmatref_type::index indexref;
    CSTmatref_type yref(y,boost::extents[3][p.dim]);
    matref_type fref(f,boost::extents[3][p.dim]);

    // Compute derivatives
    // y[i][x] = susceptible (i=0), infectious (i=1) or recovered (i=2) in age group x
    for(int i=0; i<p.dim; ++i) {
        fref[0][i] = 0.0;
        fref[1][i] = -p.gamma*yref[1][i];
        fref[2][i] = p.gamma*yref[1][i];
        for(int j=0; j<p.dim; ++j) {    
            fref[0][i] += -p.beta*p.M[i][j]*yref[1][j]*p.sigma[i]*yref[0][i]/p.N[j];
            fref[1][i] += +p.beta*p.M[i][j]*yref[1][j]*p.sigma[i]*yref[0][i]/p.N[j];
        }
    }

    return GSL_SUCCESS;

} //********** end function dydt definition ********************************************************

#endif // DYN_HPP_INCLUDED
