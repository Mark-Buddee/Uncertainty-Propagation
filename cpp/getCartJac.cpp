#include <dace/dace.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_multiroots.h>
#include <fmt/printf.h>
#include <gsl/gsl_vector.h>
#include "RK78.h"
#include "constants.h"
#include <chrono>
#include "dynamics.h"
#include "osctomean.h"
#include "convert.h"


#include "mex.hpp"
#include "mexAdapter.hpp"


using namespace std;
using namespace matlab::data;
using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {
    public:
    void operator()(ArgumentList outputs, ArgumentList inputs) override {
        double cst[13];

        cst[3] = inputs[0][0]; //Tmax
        cst[4] = inputs[0][1]; //m0
        cst[5] = inputs[0][2]; //mu
        cst[6] = inputs[0][3]; //J2
        cst[7] = inputs[0][4]; //Re
        const int method = inputs[0][5];
        cst[9] = inputs[0][6]; //params.LU;
        cst[10] = inputs[0][7]; //params.Area;
        cst[11] = inputs[0][8]; // params.Cd;
        cst[12] = inputs[0][9];//params.MU;

        const int Nnodes = inputs[0][10];
        const double nu = inputs[0][11];

        typedef DACE::DA state_type; // define state type as DA
        DACE::DA::init(2, 9);        // initialise DACE for X order in N variables
        DACE::DA::setEps(1e-40);
        DACE::AlgebraicVector<state_type> x0(6), fx(6);
        DACE::AlgebraicMatrix<state_type> Jacobian(6, 6);
        DACE::DA derivative(1), result(1), b;
        DACE::AlgebraicMatrix<double> STMgamma(6, 7);
        double statesGuess[6], constantPart[6];
        double dt; 
        
        double t1 = inputs[2][0];
        double t2 = inputs[3][0];

        int Jsum = 0;
        for (unsigned int j = 0; j < 6; j++){
            statesGuess[j] = inputs[1][j];
        }

        // organize outputs
        ArrayFactory f;
        outputs[0] = f.createArray<double>({6,7});


        for (unsigned int i = 0; i < 6; i++)
        {
            x0[i] = statesGuess[i] + DACE::DA(i + 1);

        }

//         cout << x0<< endl;
//         cout << t1<< endl;
//         cout << t2<< endl;
//         
        fx = RK78(6, x0, t1, t2, cst, DynamicsCart);

        STMgamma = fx.linear();

        for (unsigned int i = 0; i < 6; i++)
            constantPart[i] = DACE::cons(fx[i]);



        for (unsigned int i = 0; i < 6; i++)
            constantPart[i] = DACE::cons(fx[i]);

        // output stm and constant term
        for (unsigned int j = 0; j < 6; j++) // Prints row of x
        {
            for (unsigned int i = 0; i < 6; i++)
                outputs[0][Jsum +j][i]= STMgamma.at(j, i);
            outputs[0][Jsum +j][6] = constantPart[j];

        }

        Jsum += 6;


    }

};
