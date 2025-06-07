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
#include <Eigen/Dense>
#include "calcnonlinearindx.h"
#include "uncertainityConvert.h"
#include "mex.hpp"
#include "mexAdapter.hpp"

namespace mdata = matlab::data;
using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) override {
        double cst[3];
        cst[0] = inputs[0][0]; // mu
        cst[1] = inputs[0][1]; // J2
        cst[2] = inputs[0][2]; // Re
        double tf = inputs[0][3]; // tf

        typedef DACE::DA state_type;
        DACE::DA::init(2, 6);
        DACE::DA::setEps(1e-40);
        DACE::AlgebraicVector<state_type> x0(6), fx(6);

        double statesGuess[6];
        DACE::AlgebraicMatrix<double> STM(6, 6);
        vector<double> fxe(6);
        mdata::ArrayFactory f;

        outputs[0] = f.createArray<double>({6});
        outputs[1] = f.createArray<double>({6,6});

        for (unsigned int j = 0; j < 6; j++) {
            statesGuess[j] = inputs[1][j];
        }

        for (unsigned int i = 0; i < 6; i++) {
            x0[i] = statesGuess[i] + DACE::DA(i + 1);
        }

        fx = RK78(6, x0, 0.0, tf, cst, DynamicsCart);
        fxe = cons(fx);
        STM = fx.linear();

        for (unsigned int i = 0; i < 6; i++) {
            outputs[0][i] = fxe[i];
            for (unsigned int j = 0; j < 6; j++) {
                outputs[1][i][j] = cons(STM.at(i, j));
            }
        }
    }
};
