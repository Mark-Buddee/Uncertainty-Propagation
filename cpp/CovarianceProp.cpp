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
        DACE::AlgebraicVector<state_type> x0(6), x0t(6), fx(6), fxbar(6);
        DACE::AlgebraicMatrix<double> CovCar(6,6);

  
        //DACE::AlgebraicMatrix<double> ;
        double statesGuess[6], xarray[6], constantPart[6], dx0[6];
        Eigen::MatrixXd Mat(6,6);
        AlgebraicVector<state_type> EigenV(6,0.0);
        AlgebraicVector<state_type> Geqoe_Eigen_DA(6);
        AlgebraicVector<state_type> Geqoesc(6);
        DACE::AlgebraicMatrix<double> STM1(6,6), STM0(6,6);
        DACE::AlgebraicMatrix<state_type> STM(6,6);
        double xval;

        ArrayFactory f;

        // Cartesian
        outputs[0] =f.createArray<double>({6});
        outputs[1] =f.createArray<double>({6});
        outputs[2] =f.createArray<double>({1});

        // MEE 
        outputs[3] =f.createArray<double>({6});
        outputs[4] =f.createArray<double>({6});
        outputs[5] =f.createArray<double>({1});

         // Geq 
        outputs[6] =f.createArray<double>({6});
        outputs[7] =f.createArray<double>({6});
        outputs[8] =f.createArray<double>({1});

         // Ceq 
        outputs[9] =f.createArray<double>({6});
        outputs[10] =f.createArray<double>({6});
        outputs[11] =f.createArray<double>({1});

         // Kepler 
        outputs[12] =f.createArray<double>({6});
        outputs[13] =f.createArray<double>({6});
        outputs[14] =f.createArray<double>({1});

        CovCar.at(0,0) =  (1e3/6378.137e3)*(1e3/6378.137e3);
        CovCar.at(1,1) =  (1e3/6378.137e3)*(1e3/6378.137e3);
        CovCar.at(2,2) =  (1e3/6378.137e3)*(1e3/6378.137e3);
        CovCar.at(3,3) =  (1.0/6378.137e3*8.068110649226988e2)*(1.0/6378.137e3*8.068110649226988e2);
        CovCar.at(4,4) =  (1.0/6378.137e3*8.068110649226988e2)*(1.0/6378.137e3*8.068110649226988e2);
        CovCar.at(5,5) =  (1.0/6378.137e3*8.068110649226988e2)*(1.0/6378.137e3*8.068110649226988e2);

        for (unsigned int j = 0; j < 6; j++){
            statesGuess[j] = inputs[1][j];
        }

        
        vector <double>  dx(6), dx0s(6), fxe(6),fx0(6);
    
        for (unsigned int j = 0; j < 6; j++){
          
            dx[j] = inputs[2][j];
      
            dx0s[j] = 0.0;
        }

        int Jsum = 0;
        for (unsigned int i = 0; i < 6; i++)
        {
                x0[i] = statesGuess[i] + DACE::DA(i + 1);
         
        }
        
        Geqoesc = convertU(CovCar, x0);
        fx = RK78(6, Geqoesc, 0.0, tf, cst, DynamicsCart);
        xval  = calcnonlinearindexLoads(fx, 1.0);
        outputs[2][0] = xval;  

        fxe = fx.eval(dx);
        fx0 = fx.eval(dx0s);

        for (unsigned int i = 0; i < 6; i++){
            outputs[0][i] = fxe[i];
            outputs[1][i] = fx0[i]; 
                     
        }


        //convert to mee
        x0t = eci2mee(x0, mu);
        Geqoesc = convertU(CovCar, x0t);
        fx = RK78(6, Geqoesc, 0.0, tf, cst,DynamicsMEE);
        xval  = calcnonlinearindexLoads(fx, 1.0);
        outputs[5][0] = xval;
       

       // cout << fx.linear() << endl;
        // not right 
        fxe = fx.eval(dx);
        fx0 = fx.eval(dx0s);

        //cout << fx.linear() << endl;
        for (unsigned int i = 0; i < 6; i++){
            outputs[3][i] = fxe[i];
            outputs[4][i] = fx0[i];
          
        }


        // Geq
        x0t =  RV2GEq(x0, cst);
        Geqoesc = convertU(CovCar, x0t);
        fx = RK78(6, Geqoesc, 0.0, tf, cst,DynamicsEq);
        xval  = calcnonlinearindexLoads(fx, 1.0);
        outputs[8][0] = xval;

        fxe = fx.eval(dx);
        fx0 = fx.eval(dx0s);

    
        for (unsigned int i = 0; i < 6; i++){
            outputs[6][i] = fxe[i];
            outputs[7][i] = fx0[i];
          
        }

        // Ceq
        x0t =  RV2CEq(x0, cst);
        Geqoesc = convertU(CovCar, x0t);
        fx = RK78(6, Geqoesc, 0.0, tf, cst,DynamicsCeq);
        xval  = calcnonlinearindexLoads(fx, 1.0);
        outputs[11][0] = xval;

        fxe = fx.eval(dx);
        fx0 = fx.eval(dx0s);
       
        for (unsigned int i = 0; i < 6; i++){
            outputs[9][i] = fxe[i];
            outputs[10][i] = fx0[i];
          
        }

        // Kepler

        x0t =  eci2Kep(x0, cst[0]);
        Geqoesc = convertU(CovCar, x0t);
        fx = RK78(6,  Kep2eci(Geqoesc,cst[0]), 0.0, tf, cst,DynamicsCart);
         fx = eci2Kep(fx, cst[0]);
        xval  = calcnonlinearindexLoads(fx, 1.0);
        outputs[14][0] = xval;

        fxe = fx.eval(dx);
        fx0 = fx.eval(dx0s);
       
        for (unsigned int i = 0; i < 6; i++){
            outputs[12][i] = fxe[i];
            outputs[13][i] = fx0[i];
          
        }


    }
  

    };