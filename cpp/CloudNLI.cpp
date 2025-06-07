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

        cst[0] = inputs[0][0]; //mu
        cst[1] = inputs[0][1]; //J2
        cst[2] = inputs[0][2]; //Re
        double tf = inputs[0][3]; // tf

        typedef DACE::DA state_type; // define state type as DA
        DACE::DA::init(2, 6);        // initialise DACE for X order in N variables
        DACE::DA::setEps(1e-40);
        DACE::AlgebraicVector<state_type> x0(6), x0t(6), fx(6), fxbar(6);

        //DACE::AlgebraicMatrix<double> ;
        double statesGuess[6], xarray[6], constantPart[6],dx[6], dx0[6];
        Eigen::MatrixXd Mat(6,6);
        AlgebraicVector<state_type> EigenV(6,0.0);
        AlgebraicVector<state_type> Geqoe_Eigen_DA(6);
        AlgebraicVector<state_type> Geqoesc(6);
        DACE::AlgebraicMatrix<double> STM(6,6);

        double xval;

        vector <double> dxs(6), dx0s(6), fxe(6),fx0(6);
        ArrayFactory f;

        // Cartesian
        outputs[0] =f.createArray<double>({6});
        outputs[1] =f.createArray<double>({6,6});

       // MEE 
        outputs[2] =f.createArray<double>({6});
        outputs[3] =f.createArray<double>({6,6});

 // Geqoe 
        outputs[4] =f.createArray<double>({6});
        outputs[5] =f.createArray<double>({6,6});

 // Ceqoe 
        outputs[6] =f.createArray<double>({6});
        outputs[7] =f.createArray<double>({6,6});

 // Kepler 
        outputs[8] =f.createArray<double>({6});
        outputs[9] =f.createArray<double>({6,6});


        for (unsigned int j = 0; j < 6; j++){
            statesGuess[j] = inputs[1][j];
        }


        for (unsigned int i = 0; i < 6; i++)
        {
            if (i <3)
            {
                x0[i] = statesGuess[i] + DACE::DA(i + 1);
            }
            else
            {
                x0[i] = statesGuess[i] + DACE::DA(i + 1);
            }
        }

        

        // propagate x0, and output the result please.
        fx = RK78(6, x0, 0.0, tf, cst, DynamicsCart);
        
        fxe = cons(fx);
        STM  = fx.linear();

        for (unsigned int i = 0; i < 6; i++){
            outputs[0][i] = fxe[i];

            for (unsigned int j = 0; j < 6; j++){
                outputs[1][i][j] = cons(STM.at(i,j));
            }

        }

       //  MEE

       for (unsigned int j = 0; j < 6; j++){
            statesGuess[j] = inputs[2][j];
        }

        for (unsigned int i = 0; i < 6; i++)
        {
            if (i <3)
            {
                x0[i] = statesGuess[i] + DACE::DA(i + 1);
            }
            else
            {
                x0[i] = statesGuess[i] +DACE::DA(i + 1);
            }
        }


        //propagate
        fx = RK78(6, x0, 0.0, tf, cst,DynamicsMEE);
 
        fxe = cons(fx);
 

        STM  = fx.linear();
  
        for (unsigned int i = 0; i < 6; i++){
            outputs[2][i] = fxe[i];

            for (unsigned int j = 0; j < 6; j++){
                outputs[3][i][j] =STM.at(i,j);
            }

        }



             //  Geqoe

       for (unsigned int j = 0; j < 6; j++){
            statesGuess[j] = inputs[3][j];
        }

        for (unsigned int i = 0; i < 6; i++)
        {
            if (i <3)
            {
                x0[i] = statesGuess[i] + DACE::DA(i + 1);
            }
            else
            {
                x0[i] = statesGuess[i] +DACE::DA(i + 1);
            }
        }


        //propagate
        fx = RK78(6, x0, 0.0, tf, cst,DynamicsEq);
        fxe = cons(fx);
        STM  = fx.linear();
  
        for (unsigned int i = 0; i < 6; i++){
            outputs[4][i] = fxe[i];

            for (unsigned int j = 0; j < 6; j++){
                outputs[5][i][j] =STM.at(i,j);
            }

        }

                     //  Ceqoe

       for (unsigned int j = 0; j < 6; j++){
            statesGuess[j] = inputs[4][j];
        }

        for (unsigned int i = 0; i < 6; i++)
        {
            if (i <3)
            {
                x0[i] = statesGuess[i] + DACE::DA(i + 1);
            }
            else
            {
                x0[i] = statesGuess[i] +DACE::DA(i + 1);
            }
        }


        //propagate
        fx = RK78(6, x0, 0.0, tf, cst,DynamicsCeq);
        fxe = cons(fx);
        STM  = fx.linear();
  
        for (unsigned int i = 0; i < 6; i++){
            outputs[6][i] = fxe[i];

            for (unsigned int j = 0; j < 6; j++){
                outputs[7][i][j] =STM.at(i,j);
            }

        }

                             //  Kepler

       for (unsigned int j = 0; j < 6; j++){
            statesGuess[j] = inputs[5][j];
        }

        for (unsigned int i = 0; i < 6; i++)
        {
            if (i <3)
            {
                x0[i] = statesGuess[i] + DACE::DA(i + 1);
            }
            else
            {
                x0[i] = statesGuess[i] +DACE::DA(i + 1);
            }
        }


        //propagate
        fx = RK78(6, x0, 0.0, tf, cst,DynamicsKep);
        fxe = cons(fx);
        STM  = fx.linear();
  
        for (unsigned int i = 0; i < 6; i++){
            outputs[8][i] = fxe[i];

            for (unsigned int j = 0; j < 6; j++){
                outputs[9][i][j] =STM.at(i,j);
            }

        }

        


    }


};