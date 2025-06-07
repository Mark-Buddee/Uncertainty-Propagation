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
#include "convert.h"
#include "dynamics.h"
#include "osctomean.h"


#include "mex.hpp"
#include "mexAdapter.hpp"


using namespace std;
using namespace matlab::data;
using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {
    public:
    void operator()(ArgumentList outputs, ArgumentList inputs) override {

        double cst[17];

      
     
        cst[3] = inputs[0][0]; //Tmax
        cst[4] = inputs[0][1]; //m0
        cst[5] = inputs[0][2]; //mu
        cst[6] = inputs[0][3]; //J2
        cst[7] = inputs[0][4]; //Re
        cst[13] = inputs[0][12]; //Isp
        cst[14] = inputs[0][13]; // g0 
        cst[15] = inputs[0][15]; // drag bool

        const int nodecoord = inputs[0][16]; //coordinate type for the nodes 
       

        const int method = inputs[0][5];

        if (method < 2)
        {
            cst[8] = 1.0;
            cout << "Method 1 running" << endl;
        }
        else
        {
            cst[8] = 2.0;
            cout << "Method 2 (dv control) running " << endl;
        }

        double LU = inputs[0][6]; 
        double TU = inputs[0][14];

        cst[9] = inputs[0][6]; //params.LU;
        cst[10] = inputs[0][7]; //params.Area;
        cst[11] = inputs[0][8]; // params.Cd;
        cst[12] = inputs[0][9];//params.MU;

        const int Nnodes = inputs[0][10];
        const double nu = inputs[0][11];

        typedef DACE::DA state_type; // define state type as DA
        DACE::DA::init(2, 9);        // initialise DACE for X order in N variables
        DACE::DA::setEps(1e-40);
        DACE::AlgebraicVector<state_type> xTargetKep(3), x0(6), fx(6),  dv(3), dx(3),fxnew(6),fxendRV(6),fxendCC(6);
        DACE::AlgebraicMatrix<state_type> Jacobian(6, 9), xval(Nnodes, 9);
        DACE::DA derivative(1), result(1), b;
        DACE::AlgebraicMatrix<double> STMgamma(6, 9), STMxftransform(6,9), Jbar(6, 9);
        double statesGuess[Nnodes+ 1][9], xarray[9], constantPart[6], B, normJbar, time, maxval, eclipse;
        double pi = 4.0*atan(1.0);

        vector<double> zerostate(9);
        state_type r , v, Ws, a;
        vector<state_type> zerostateDA(9);
        double  dt[Nnodes],massguess[Nnodes];
        DACE::AlgebraicVector<double> weights(4);
        state_type Q; 
        DACE::AlgebraicVector<double> Qlinear(9);
        DACE::AlgebraicVector<state_type> Qlawfrac(3);

        for (unsigned int i = 0; i < Nnodes + 1; i++){
            for (unsigned int j = 0; j < 9; j++){
                statesGuess[i][j] = inputs[1][i][j];
            }
        }

        for (unsigned int i = 0; i < Nnodes; i++){
            dt[i] = inputs[2][i];
            massguess[i] = inputs[3][i];
        }

        const double targetcoord = inputs[4][0]; //1: kepler. 2: MEE with a
         
         

        // organize outputs
        ArrayFactory f;
        outputs[0] = f.createArray<double>({6*static_cast<unsigned long>(Nnodes),10});
        outputs[1] = f.createArray<double>({static_cast<unsigned long>(Nnodes), 9});
        outputs[2] = f.createArray<double>({6, 7});

        int Jsum = 0; 
        for (int k = 0; k < Nnodes; k++)

        {
            cst[16] = massguess[k];

            time += dt[k];
            
            for (unsigned int i = 0; i < 6; i++)
            {
                x0[i] = statesGuess[k][i] + DACE::DA(i + 1);

                if (i < 3)
                    cst[i] = statesGuess[k][i + 6];
            }



            // 1 = Cartesian, 2= MEE, 3 = GeqOE, 4 = CeqOE, 5 = Keplerian 

            if (nodecoord == 1){
                // do nothing 
            }
            else if (nodecoord == 2){
                x0 = mee2eci(x0, cst[5]);

            }
            else if (nodecoord == 3){
                // convert from GEQ to cartesian
                x0 = GEq2RV(x0, cst);
            }
            else if (nodecoord == 4){
                   x0 = CEq2RV(x0, cst);
            }
            else if (nodecoord == 5){
               // cout << "here" << endl;
               // cout << cons(x0) << endl; 
                 x0 = Kep2eci(x0, cst[5]);
                // cout << cons(x0) << endl ;
                
            }
            else{
                cout << "Error" << endl; 

            }

           
            
            fx = RK78(6, x0, 0.0, dt[k], cst, DynamicsCart);

             
        

            if (nodecoord == 1){
                // do nothing 
            }
            else if (nodecoord == 2){
                fx = eci2mee(fx, cst[5]);

            }
            else if (nodecoord == 3){
                // convert from GEQ to cartesian
                fx = RV2GEq(fx, cst);
            }
            else if (nodecoord == 4){
                   fx = RV2CEq(fx, cst);
            }
            else if (nodecoord == 5){

                 fx = eci2Kep(fx, cst[5]);
             //   cout << cons(fx)<< endl;
            }
              else{
                cout << "Error" << endl; 

            }


            STMgamma = fx.linear();
           

            if (k == Nnodes-1){

                // do cons(fx)+DA
                for (unsigned int i = 0; i < 6; i++)
                {
                    fxnew[i] = cons(fx[i])+ DACE::DA(i + 1);
                }
               

                // convert to mean keplerian
                //fxendRV  = GEq2RV(fxnew,cst);

                if (nodecoord == 1){
                    fxendRV = fxnew;
                }
                else if (nodecoord == 2){
                    fxendRV = mee2eci(fxnew, cst[5]);

                }
                else if (nodecoord == 3){
                    // convert from GEQ to cartesian
                    fxendRV = GEq2RV(fxnew, cst);
                }
                else if (nodecoord == 4){
                    fxendRV = CEq2RV(fxnew, cst);
                }
                else if (nodecoord == 5){
                    fxendRV = Kep2eci(fxnew, cst[5]);
                }

           


                
                fxendRV = oscCart2MeanCart(fxendRV,cst[5], cst[6], cst[7]);
               // cout << cons(fxendRV) << endl;

               

                if (targetcoord <1.5)

                {
                    fxendCC = cart2kep(fxendRV,cst[5]);

                    if( cons(fxendCC[3]) > 2*pi-0.0001){
                        cout << "RAAN > pi" << endl;
                        fxendCC[3]  =  fxendCC[3]  - 2*pi;
                    }

                }
                else {
                    cout << "converting to MEE with sma" << endl;
                    
                    fxendCC = cart2kep(fxendRV,cst[5]);

                     fxendCC = Kep2meesma(fxendCC);

                    fxendCC = eci2meewsma(fxendRV,cst[5]);

                }


 
                // output the maps.
                STMxftransform = fxendCC.linear();

                for (unsigned int j = 0; j < 6; j++) // Prints row of x
                {

                    for (unsigned int i = 0; i < 6; i++)
                    {
                        outputs[2][j][i] = STMxftransform.at(j, i);
                    }
                    outputs[2][j][6] = cons(fxendCC[j]);
                }
            }

            
            for (unsigned int i = 0; i < 6; i++)
                constantPart[i] = DACE::cons(fx[i]);


            // output stm and constant term
            for (unsigned int j = 0; j < 6; j++) // Prints row of x
            {
                for (unsigned int i = 0; i < 9; i++)
                    outputs[0][Jsum +j][i]= STMgamma.at(j, i);

                outputs[0][Jsum +j][9] = constantPart[j];
            }

            Jsum += 6;


            // calculate the jacobian and jbar norm
            normJbar = 0;
            for (unsigned int i = 0; i < 6; i++)
            {
                for (unsigned int j = 1; j < 10; j++)
                {
                    Jacobian.at(i, j - 1) = fx[i].deriv(j);
                    Jbar.at(i, j - 1) = fx[i].deriv(j).eval(zerostate);
                    normJbar += pow(Jbar.at(i, j - 1),2);
                }
            }

            // analyse the dependency on each variable individually
            for (int l = 0; l < 9; l++)
            {
                zerostateDA = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0};

                zerostateDA[l] = DACE::DA(l + 1);

                b = 0;
                for (unsigned int i = 0; i < 6; i++)
                {
                    for (unsigned int j = 1; j < 10; j++)
                    {
                        b += pow(Jacobian.at(i, j - 1).eval(zerostateDA).deriv(l + 1),2);
                    }
                }

                if (cons(b) == 0)
                {b = 1e-30;} 

                 xval.at(k, l) =  nu *sqrt(normJbar) / sqrt(b);
                
                outputs[1][k][l] = DACE::cons(xval.at(k, l));

            }


        }
    

    }
};