#include <cliques/helpers/math.h>
#include <cliques/tests/catch.hpp>
#include <cliques/helpers/io.h>
#include <vector>
#include <math.h>

TEST_CASE( "Math helpers behave correctly", "[Math]" ) {    
    //input matrix
    std::vector<double> matrix = {1, 1, 1,
                                  0, 1, 1, 
                                  0, 0, 1};
    double my_time = 1;
    int msize = 3;
    
    //compute matrix exponential
    std::vector<double> exp_matrix = clq::exp(matrix, my_time, msize);
    
    REQUIRE(exp_matrx.size() == 9);

    // exact analytical results
    double e = std::exp(1);
    std::vector<double> result = {e, e, 1.5*e,
                                  0, e, e,
                                  0, 0, e};
    
    //print result and compare
    for(double i : exp_matrix){ 
        REQUIRE(result[k] - i< 1e-12);
    }
}
