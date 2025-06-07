#ifndef dynamics_H
#define dynamics_H

template <typename T, typename U>
DACE::AlgebraicVector<T>  J2accelerationMEE(DACE::AlgebraicVector<T> x, U J2 , U Re, U mu, T r)
{

    T h = x[3];
    T k = x[4];
    T L = x[5];


    DACE::AlgebraicVector<T>  FJ2(3);

    FJ2[0] = -3*mu*J2*Re*Re/(2*r*r*r*r)*(1 - 12*(h*sin(L) - k*cos(L))*(h*sin(L) - k*cos(L))/ (1+ h*h + k*k)/(1+ h*h + k*k));
    FJ2[1] = -12*mu*J2*Re*Re/(r*r*r*r)*(h*sin(L) - k*cos(L))*(h*cos(L) + k*sin(L))/(1+ h*h + k*k)/(1+ h*h + k*k);
    FJ2[2] = -6*mu*J2*Re*Re/(r*r*r*r)*(1 - h*h - k*k)*(h*sin(L) - k*cos(L))/(1+ h*h + k*k)/(1+ h*h + k*k);

    return FJ2;


}


template <typename T>
T getAandB(DACE::AlgebraicVector<T> &x, T (&A)[6], T (&B)[6][3], const double mu)
{
  
  using DACE::cos;
  using DACE::pow;
  using DACE::sin;
  using DACE::sqrt;
  using std::cos;
  using std::pow;
  using std::sin;
  using std::sqrt;

  // unpack state vector
  T p = x[0];
  T f = x[1];
  T g = x[2];
  T h = x[3];
  T k = x[4];
  T L = x[5];

  // calculate intermediate variables that simply matrices
  T cosL = cos(L);
  T sinL = sin(L);
  T w = 1 + f * cosL + g * sinL;
  T w2 = w * w;
  T s2 = 1 + h * h + k * k;
  T sqrt1 = sqrt(p / mu);
  T sqrt2 = sqrt(mu * p);
  T pow1 = pow(2 * f * cosL + 2 * g * sinL + 2, 2);

  // fill A matrix
  A[0] = 0.0;
  A[1] = 0.0;
  A[2] = 0.0;
  A[3] = 0.0;
  A[4] = 0.0;
  A[5] = sqrt2 * (w / p)*(w / p);

  // fill B matrix
  B[0][0] = 0.0;
  B[0][1] = 2 * p / w * sqrt1;
  B[0][2] = 0.0;
  B[1][0] = sqrt1 * sinL;
  B[1][1] = sqrt1 * (1 / w) * ((w + 1) * cosL + f);
  B[1][2] = -sqrt1 * (g / w) * (h * sinL - k * cosL);
  B[2][0] = -sqrt1 * cosL;
  B[2][1] = sqrt1 * (1 / w) * ((w + 1) * sinL + g);
  B[2][2] = sqrt1 * (f / w) * (h * sinL - k * cosL);
  B[3][0] = 0.0;
  B[3][1] = 0.0;
  B[3][2] = sqrt1 * s2 * cosL / (2 * w);
  B[4][0] = 0.0;
  B[4][1] = 0.0;
  B[4][2] = sqrt1 * s2 * sinL / (2 * w);
  B[5][0] = 0.0;
  B[5][1] = 0.0;
  B[5][2] = sqrt1 * (1 / w) * (h * sinL - k * cosL);

  return sqrt(mu * x[0]) * (w / x[0]) * (w / x[0]);
}

template <typename U>
U Density_HP(U heightm)
{

  U h[50] = {100.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0, 260.0, 270.0, 280.0, 290.0, 300.0, 320.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 620.0, 640.0, 660.0, 680.0, 700.0, 720.0, 740.0, 760.0, 780.0, 800.0, 840.0, 880.0, 920.0, 960.0, 1000.0};

  U c_min[50] = {4.974e+05, 2.490e+04, 8.377e+03, 3.899e+03, 2.122e+03, 1.263e+03,
                 8.008e+02, 5.283e+02, 3.617e+02, 2.557e+02, 1.839e+02, 1.341e+02,
                 9.949e+01, 7.488e+01, 5.709e+01, 4.403e+01, 3.430e+01, 2.697e+01,
                 2.139e+01, 1.708e+01, 1.099e+01, 7.214e+00, 4.824e+00, 3.274e+00,
                 2.249e+00, 1.558e+00, 1.091e+00, 7.701e-01, 5.474e-01, 3.916e-01,
                 2.819e-01, 2.042e-01, 1.488e-01, 1.092e-01, 8.070e-02, 6.012e-02,
                 4.519e-02, 3.430e-02, 2.632e-02, 2.043e-02, 1.607e-02, 1.281e-02,
                 1.036e-02, 8.496e-03, 7.069e-03, 4.680e-03, 3.200e-03, 2.210e-03,
                 1.560e-03, 1.150e-03};


  heightm = heightm / 1000.0;

  U density = 0;

  if (heightm > 1000.0)
  {
    density = 1.150e-03;
  }
  else if (heightm < 100.0)
  {
    density = 4.974e+05;
  }
  else
  {

    for (int i = 0; i < 50; i++)
    {
      if (h[i] <= heightm && h[i + 1] >= heightm)
      {
        U diffh = heightm - h[i];
        U diffh2 = h[i + 1] - h[i];

        density = c_min[i] + (c_min[i + 1] - c_min[i]) * diffh / diffh2;
      }
    }
  }

  return density * 1e-12;
}

template <typename U>
U Smoothed_eclipse(U cs, U ct, DACE::AlgebraicVector<double> x_sun, DACE::AlgebraicVector<double> x_sat, U Rs, U Re)
{

  using DACE::cos;
  using DACE::pow;
  using DACE::sin;
  using DACE::sqrt;
  using std::cos;
  using std::pow;
  using std::sin;
  using std::sqrt;

  // cout << "Rs =" << Rs << endl;

  U a_SR = asin(Rs / x_sun.vnorm());
  U a_BR = asin(Re / x_sat.vnorm());

  U a_D = acos(x_sat.dot(x_sun) / x_sat.vnorm() / x_sun.vnorm());

  U gamma_L = 1 / (1 + exp(-cs * (a_D - ct * (a_SR + a_BR))));

  return gamma_L;
}

template <typename U>
DACE::AlgebraicVector<U> getSunPosVec(U JD, U LU)
{

  U pi = 4 * atan(1.0);

  U T_UT1 = (JD - 2451545.0) / 36525;

  U lambdaMSun = 280.460 + 36000.771 * T_UT1;

  U MSun = 357.5291092 + 35999.05034 * T_UT1;

  U lambdaEcliptic = lambdaMSun + 1.914666471 * sin(MSun * pi / 180.0) + 0.019994643 * sin(2 * MSun * pi / 180.0);

  U r_Sun = 1.000140612 - 0.016708617 * cos(MSun * pi / 180.0) - 0.000139589 * cos(2 * MSun * pi / 180.0);

  U epsilon = 23.439291 - 0.0130042 * T_UT1;

  DACE::AlgebraicVector<U> SunVec(3);
  SunVec[0] = r_Sun * cos(lambdaEcliptic * pi / 180.0);
  SunVec[1] = r_Sun * cos(epsilon * pi / 180.0) * sin(lambdaEcliptic * pi / 180.0);
  SunVec[2] = r_Sun * sin(epsilon * pi / 180.0) * sin(lambdaEcliptic * pi / 180.0);

  SunVec = 1.495978707e11 * SunVec / LU;

  // cout << SunVec << endl;
  return SunVec;
}

template <typename T, typename U>
double getEclipse(U t, DACE::AlgebraicVector<T> x, U cs, U ct, U LU, U Rs, U Re)
{

  DACE::AlgebraicVector<double> r_Earth2Sun(3);

  r_Earth2Sun = getSunPosVec(t + 2400000.5 + 60000, LU);

  // convert dace to double.
  DACE::AlgebraicVector<double> x_double(6);
  for (unsigned int i = 0; i < 6; i++)
  {
    x_double[i] = DACE::cons(x[i]);
  }

  DACE::AlgebraicVector<double> CART = x_double;

  DACE::AlgebraicVector<double> rSat2Earth(3);
  DACE::AlgebraicVector<double> rSat2Sun(3);

  for (unsigned int i = 0; i < 3; i++)
  {
    rSat2Earth[i] = -CART[i];
    rSat2Sun[i] = rSat2Earth[i] + r_Earth2Sun[i];
  }

  double numinus1 = Smoothed_eclipse(cs, ct, rSat2Sun, rSat2Earth, Rs, Re);

  return 1.0 - numinus1;
}

template <typename T, typename U>
DACE::AlgebraicVector<T> DynamicsCart(DACE::AlgebraicVector<T> x, U t, const double *cst)
{

    using DACE::cos;
    using DACE::pow;
    using DACE::sin;
    using DACE::sqrt;
    using std::cos;
    using std::pow;
    using std::sin;
    using std::sqrt;


    const double Tmax =  cst[3];
    const double Isp = cst[13];
    const double g0 = cst[14];
    const double mu = cst[5];
    const double J2 = cst[6];
    const double Re = cst[7];
    const double LU = cst[9];
    const double Area = cst[10];
    const double Cd = cst[11];
    const double m0 = cst[4];
    const double massunit = cst[12];
    const double mguess = cst[16];
    DACE::AlgebraicVector<T> dxdt(6);

    //cout << mguess << endl;

    T r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
    // cout << DACE::cons(term) << endl;

    if (DACE::cons(r2) <= 0)
    {
        cout << "So the iter.tex has a mistake, its giving all states as zero, likely after infeasible run" << endl;
    }

    T r = sqrt(r2);

    T r3 = r * r * r;

    for (unsigned int i = 0; i < 3; i++)
    {
        dxdt[i] = x[3 + i];
        dxdt[i + 3] = -mu * x[i] / (r3);
    }

    dxdt[3] = dxdt[3] - mu * x[0] / (r3)*1.5 * J2 * (Re / r) * (Re / r) * (1 - 5 * x[2] * x[2] / r2);
    dxdt[4] = dxdt[4] - mu * x[1] / (r3)*1.5 * J2 * (Re / r) * (Re / r) * (1 - 5 * x[2] * x[2] / r2);
    dxdt[5] = dxdt[5] - mu * x[2] / (r3)*1.5 * J2 * (Re / r) * (Re / r) * (3 - 5 * x[2] * x[2] / r2);

    //ADD DRAG
    if ( cst[15] > 0.0)

    {
        //cout << "the heck" << endl;
        double rho = Density_HP(DACE::cons(r*LU -Re*LU))*LU*LU*LU/massunit;

        T vmag = sqrt(x[3] * x[3] + x[4] * x[4] + x[5] * x[5]);
        for (unsigned int i = 0; i < 3; i++)
        {
            dxdt[3 + i] =  dxdt[3 + i] + 0.5 * rho * Cd * Area * vmag * vmag /mguess * -x[3 + i] / vmag;
        }
    }



    DACE::AlgebraicVector<T> u(3), uECI(3);
    for (unsigned int i = 0; i < 3; i++)
    {
        //cout << cst[3] << endl;
        u[i] =  cst[i] +  DACE::DA(i + 7); //
    }

    // convert to ECI
    uECI = RTN2ECI(x,u);

    for (unsigned int i = 0; i < 3; i++)
    {
        dxdt[3 + i] +=  uECI[i];
    }

    return dxdt;


}


template <typename T, typename U>
DACE::AlgebraicVector<T> DynamicsMEE(DACE::AlgebraicVector<T> x, U t, const double *cst)
{

    // initialise matrices to zero
    T A[6], B[6][3];



    const double Tmax =  cst[3];
    const double Isp = cst[13];
    const double g0 = cst[14];
    const double mu = cst[5];
    const double J2 = cst[6];
    const double Re = cst[7];
    const double LU = cst[9];
    const double Area = cst[10];
    const double Cd = cst[11];
    const double m0 = cst[4];
    const double massunit = cst[12];
    const double mguess = cst[16];


    DACE::AlgebraicVector<T>  FJ2(3);

    DACE::AlgebraicVector<T> dxdt(6);

    T r = x[0]/(1.0 + x[1]*cos(x[5]) + x[2]*sin(x[5]));

    FJ2 = J2accelerationMEE(x, J2 , Re, mu, r);

    getAandB(x, A, B,mu);


    // acceleration
    DACE::AlgebraicVector<T> u(3);
    for (unsigned int i = 0; i < 3; i++)
    {
        // cout << params.u0[i] << endl;
        u[i] = cst[i] + DACE::DA(i + 7); //
    }



    for (int i = 0; i < 6; i++)
    {
        dxdt[i] = A[i];

        for (int j = 0; j < 3; j++)
        {
            dxdt[i] = dxdt[i] + B[i][j] * FJ2[j] + B[i][j]*u[j];
        }


    }




    return dxdt;
}




#endif