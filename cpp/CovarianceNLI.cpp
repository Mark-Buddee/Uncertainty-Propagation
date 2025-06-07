#include <dace/dace.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_multiroots.h>
#include "AstroCnvRef.h"
#include <fmt/printf.h>
#include <gsl/gsl_vector.h>
#include "RK78.h"
#include "constants.h"
#include <chrono>
#include "osctomean.h"
#include "convert.h"
#include "dynamics.h"
#include <eigen/Eigen/Dense>
#include "calcnonlinearindx.h"
#include "uncertainityConvert.h"
#include "mex.hpp"
#include "mexAdapter.hpp"

using namespace std;
using namespace matlab::data;

using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {
    public:
    void operator()(ArgumentList outputs, ArgumentList inputs) override {

        double cst[3];

        double mu = inputs[0][0];
        double J2 = inputs[0][1];
        double Re = inputs[0][2];
        cst[0] = mu; 
        cst[1] = J2; 
        cst[2] = Re;
        double tf = inputs[0][3]; // tf
   
        typedef DACE::DA state_type; // define state type as DA
        DACE::DA::init(2, 6);        // initialise DACE for X order in N variables
        DACE::DA::setEps(1e-40);
        DACE::AlgebraicVector<state_type> x0(6), x0t(6), fx(6), xesc(6);
        DACE::AlgebraicMatrix<double> CovCar(6,6);

        double statesGuess[6], sigmaCart[6], xval;

        ArrayFactory f;
        // Nonlinear indicies. 
        outputs[0] =f.createArray<double>({5});

        for (unsigned int j = 0; j < 6; j++){
            statesGuess[j] = inputs[1][j];
            sigmaCart[j] = inputs[2][j];
            CovCar.at(j,j) = sigmaCart[j]*sigmaCart[j];
        }
        
      
        for (unsigned int i = 0; i < 6; i++)
          {  x0[i] = statesGuess[i] + DACE::DA(i + 1);
          }

        // Cartesian
        xesc = axisTransform(CovCar ,x0);
        fx = RK78(6, xesc, 0.0, tf, cst, DynamicsCart);
        xval  = calcnonlinearindexLoads(fx, 3.0);
        outputs[0][0] = xval;

        //MEE
        x0t = eci2mee(x0, mu);
        xesc = axisTransform(convertCovariance( CovCar ,x0t),  x0t);
        fx = RK78(6, xesc, 0.0, tf, cst,DynamicsMEE);
        xval  =  calcnonlinearindexLoads(fx, 3.0);
        outputs[0][1] = xval;

        // Geq
        x0t =  RV2GEq(x0, cst);
        xesc = axisTransform(convertCovariance( CovCar ,x0t),  x0t);
        fx = RK78(6, xesc, 0.0, tf, cst,DynamicsEq);
         xval  =  calcnonlinearindexLoads(fx, 3.0);
        outputs[0][2] =xval;

        // Ceq
        x0t =  RV2CEq(x0, cst);
        xesc = axisTransform(convertCovariance( CovCar ,x0t),  x0t);
        fx = RK78(6, xesc, 0.0, tf, cst,DynamicsCeq);
         xval  =  calcnonlinearindexLoads(fx, 3.0);
        outputs[0][3]  =xval; 

        // Kepler
        x0t =  eci2kep2(x0, cst[0]); 
        xesc = axisTransform(convertCovariance( CovCar ,x0t),  x0t);
        fx = RK78(6, xesc, 0.0, tf, cst,DynamicsKep);
          xval  =  calcnonlinearindexLoads(fx, 3.0);
        outputs[0][4] = xval;

    }
  

    };