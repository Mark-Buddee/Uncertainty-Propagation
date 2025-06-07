#ifndef dynamics_H
#define dynamics_H

template <typename T, typename U>
DACE::AlgebraicVector<T>  J2accelerationKep(DACE::AlgebraicVector<T> x, U J2 , U Re, U mu, T r)
{
    T a = x[0];
    T e = x[1];
    T inc = x[2];
    T RAAN = x[3];
    T AOP = x[4];
    T f = x[5];


    T theta = f + AOP; 


    DACE::AlgebraicVector<T>  FJ2(3);

    FJ2[0] = -3*mu*J2*Re*Re/(r*r*r*r)*(0.5 - 1.5*sin(inc)*sin(inc)*sin(theta)*sin(theta));
    FJ2[1] = -3*mu*J2*Re*Re/(r*r*r*r)*(sin(inc)*sin(inc)*sin(theta)*cos(theta));
    FJ2[2] = -3*mu*J2*Re*Re/(r*r*r*r)*(sin(inc)*sin(theta)*cos(inc));


    return FJ2; 


}


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

    const double mu = cst[0];
    const double J2 = cst[1];
    const double Re = cst[2];


  DACE::AlgebraicVector<T> dxdt(6);

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

    // cout << J2 << endl;
    //
    dxdt[3] = dxdt[3] - mu * x[0] / (r3)*1.5 * J2 * (Re / r) * (Re / r) * (1 - 5 * x[2] * x[2] / r2);
    dxdt[4] = dxdt[4] - mu * x[1] / (r3)*1.5 * J2 * (Re / r) * (Re / r) * (1 - 5 * x[2] * x[2] / r2);
    dxdt[5] = dxdt[5] - mu * x[2] / (r3)*1.5 * J2 * (Re / r) * (Re / r) * (3 - 5 * x[2] * x[2] / r2);

  return dxdt;
}

template <typename T, typename U>
DACE::AlgebraicVector<T> DynamicsCeq(DACE::AlgebraicVector<T> X, U t, const double *cst)
{

  U pi = 4.0 * atan(1.0);
    const double mu = cst[0];
    const double J2 = cst[1];
    const double Re = cst[2];


  DACE::AlgebraicVector<T> kep = RV2OE(GEq2RV(X, cst), mu);

  T inc = kep[2];
  T Raan = kep[3]; 
  T AOP = kep[4];
  T theta = kep[5];

  DACE::AlgebraicVector<T> dxdt(6);

  T nu = X[0];
  T p1 = X[1];
  T p2 = X[2];
  T varL = X[3];
  T q1 = X[4];
  T q2 = X[5];

  // auxiliary variables
  T g2 = p1 * p1 + p2 * p2;
  T a = pow((mu / nu / nu), 1.0 / 3.0);
  T c = pow((mu * mu / nu), 1.0 / 3.0) * sqrt(1.0 - p1 * p1 - p2 * p2);
  T rho = a * (1.0 - g2);
  T alpha = 1.0 / (1 + sqrt(1.0 - p1 * p1 - p2 * p2));

  // generalised Kepler equation
  U tolN = 1e-12;
  U diff = 0.5;
  // cout << "eps" << endl << eps << endl ;
  U K = 2.094395102393196;
  while (abs(diff) > tolN)
  {
    U K_old = K;

    T fx = K + p1 * cos(K) - p2 * sin(K) - varL;

    T fxd = 1.0 - p1 * sin(K) - p2 * cos(K);

    K = K_old - cons(fx / fxd);

    diff = K_old - K;
  }


  T r = a * (1 - p1 * sin(K) - p2 * cos(K));
  T r_dot = sqrt(mu * a) / r * (p2 * sin(K) - p1 * cos(K));

  T sL = a / r * (alpha * p1 * p2 * cos(K) + (1.0 - alpha * p2 * p2) * sin(K) - p1);
  T cL = a / r * (alpha * p1 * p2 * sin(K) + (1.0 - alpha * p1 * p1) * cos(K) - p2);

  T L = atan2(sL, cL);

  if (cons(L) < 0.0)
  {
    L = L + 2.0 * pi;
  }

  T Kq = 1.0 / (1.0 + q1 * q1 + q2 * q2);

  DACE::AlgebraicVector<T> ex(3), ey(3), ez(3);

  ex = {Kq * (1.0 - q1 * q1 + q2 * q2), Kq * 2.0 * q1 * q2, Kq * -2.0 * q1};
  ey = {Kq * (2.0 * q1 * q2), Kq * (1.0 + q1 * q1 - q2 * q2), Kq * 2.0 * q2};
  ez = DACE::cross(ex, ey);

  DACE::AlgebraicVector<T> er(3), ef(3), eh(3), R(3);
  er = ex * cos(L) + ey * sin(L);
  ef = ey * cos(L) - ex * sin(L);
  eh = ez;

  T x = r * cos(L);
  T y = r * sin(L);
  R = x * ex + y * ey;

  T phi = acos(R[2]/r);
  T UJ2 = 0.0;
  T Ut = 0.0; 

  T Fh = -1.5*J2*mu*Re*Re/(r*r*r*r)*sin(2*inc)*sin(AOP+theta);
  T Fr = -1.5*J2*mu*Re*Re/(r*r*r*r)*(1- 3*sin(inc)*sin(inc)*sin(AOP+theta)*sin(AOP+theta)); 
  T Pr = 0.0;
  T Pf = 0.0;

  // Perturbations
  // Variables that required perturbations
  T h = sqrt(c * c - 2.0 * r * r * UJ2); // angular momentum magnitude
  T omegh = q1 * cos(L) - q2 * sin(L); // angular velocity of equinoctial frame (h direction tilde)
  T eps_dot =  Ut + r_dot*Pr + h/r*Pf;   // total energy time derivative

  // Recurring parameters
  T A = (h - c) / r / r;
  T B =  r/h*omegh*Fh;
  T CAP = (2.0 * UJ2 -r*Fr) / c;
  T D = r / mu * eps_dot;

  // Non-dimensional variables
  T sigma = r / rho;
  T xi = 1.0 + sigma;

  dxdt[0] = -3 * pow((nu / mu / mu), (1.0 / 3.0)) * eps_dot;
  dxdt[1] = p2 * (A - B) + CAP * (r * r_dot / c * p1 + xi * p2 + sigma * cos(L)) + D * (sigma * p1 + xi * sin(L));
  dxdt[2] = p1 * (B - A) + CAP * (r * r_dot / c * p2 - xi * p1 - sigma * sin(L)) + D * (sigma * p2 + xi * cos(L));
  dxdt[3] = nu + A - B + D * r_dot * c / mu * alpha * xi + CAP * (1 / alpha + alpha * (1 - r / a));

  T Kqdot =r/2/h/Kq*Fh;
  dxdt[4] = Kqdot * sin(L);
  dxdt[5] = Kqdot * cos(L);

 
  return dxdt;


}


template <typename T, typename U>
DACE::AlgebraicVector<T> DynamicsMEE(DACE::AlgebraicVector<T> x, U t, const double *cst)
{

    // initialise matrices to zero
    T A[6], B[6][3];


    const double mu = cst[0];
    const double J2 = cst[1];
    const double Re = cst[2];

    DACE::AlgebraicVector<T>  FJ2(3);

    DACE::AlgebraicVector<T> dxdt(6);

    T r = x[0]/(1.0 + x[1]*cos(x[5]) + x[2]*sin(x[5]));

    FJ2 = J2accelerationMEE(x, J2 , Re, mu, r);

    getAandB(x, A, B,mu);


    for (int i = 0; i < 6; i++)
    {
        dxdt[i] = A[i];

        for (int j = 0; j < 3; j++)
        {
            dxdt[i] = dxdt[i] + B[i][j] * FJ2[j];
        }


    }

    return dxdt;
}





template <typename T, typename U>
    DACE::AlgebraicVector<T> DynamicsKep(DACE::AlgebraicVector<T> x, U t, const double *cst)
{

    const double mu = cst[0];
    const double J2 = cst[1];
    const double Re = cst[2];

    DACE::AlgebraicVector<T> FJ2(3); 

    T a = x[0];
    T e = x[1];
    T inc = x[2];
    T RAAN = x[3];
    T AOP = x[4];
    T theta = x[5];

    T p = a*(1-e*e);
    T r = p/(1 + e*cos(theta));

    T h = sqrt(mu*p);
    
    FJ2 = J2accelerationKep(x, J2 , Re, mu, r);

    T Fr = FJ2[0];
    T Ft = FJ2[1];
    T Fn = FJ2[2];

    T E = 2*atan(tan(theta/2)*sqrt((1-e)/(1+e)));

    DACE::AlgebraicVector<T> doedt(6);

    doedt[0] = 2*(a*a)/h*(e*sin(theta)*Fr + p*Ft/r);
    doedt[1] = sqrt(p/mu)*(Fr*sin(theta)+ Ft*(cos(E) + cos(theta)));
    doedt[2] = r*cos(theta+AOP)*Fn/h;
    doedt[3] = r*sin(theta+AOP)*Fn/(h*sin(inc));
    doedt[4] = -sqrt(p/mu)*(Fn*r/p*1.0/tan(inc)*sin(theta+AOP) + 1.0/e*(Fr*cos(theta) -Ft*(1+r/p)*sin(theta)));
    doedt[5] = h/r/r+ 1/(e*h)*(p*cos(theta)*Fr -(p+r)*sin(theta)*Fn);


        


    return doedt;

}



template <typename T, typename U>
DACE::AlgebraicVector<T> DynamicsEq(DACE::AlgebraicVector<T> X, U t, const double *cst)
{

  U pi = 4.0 * atan(1.0);
    const double mu = cst[0];
    const double J2 = cst[1];
    const double Re = cst[2];

  DACE::AlgebraicVector<T> kep = RV2OE(GEq2RV(X, cst), mu);

  T inc = kep[2];
  T Raan = kep[3]; 
  T AOP = kep[4];
  T theta = kep[5];

  DACE::AlgebraicVector<T> dxdt(6);

  T nu = X[0];
  T p1 = X[1];
  T p2 = X[2];
  T varL = X[3];
  T q1 = X[4];
  T q2 = X[5];

  // auxiliary variables
  T g2 = p1 * p1 + p2 * p2;
  T a = pow((mu / nu / nu), 1.0 / 3.0);
  T c = pow((mu * mu / nu), 1.0 / 3.0) * sqrt(1.0 - p1 * p1 - p2 * p2);
  T rho = a * (1.0 - g2);
  T alpha = 1.0 / (1 + sqrt(1.0 - p1 * p1 - p2 * p2));

  // generalised Kepler equation
  U tolN = 1e-13;
  U diff = 0.5;
  // cout << "eps" << endl << eps << endl ;
  U K = 2.094395102393196;
  while (abs(diff) > tolN)
  {
    U K_old = K;

    T fx = K + p1 * cos(K) - p2 * sin(K) - varL;

    T fxd = 1.0 - p1 * sin(K) - p2 * cos(K);

    K = K_old - cons(fx / fxd);

    diff = K_old - K;
  }


  T r = a * (1 - p1 * sin(K) - p2 * cos(K));
  T r_dot = sqrt(mu * a) / r * (p2 * sin(K) - p1 * cos(K));


  T sL = a / r * (alpha * p1 * p2 * cos(K) + (1.0 - alpha * p2 * p2) * sin(K) - p1);
  T cL = a / r * (alpha * p1 * p2 * sin(K) + (1.0 - alpha * p1 * p1) * cos(K) - p2);

  T L = atan2(sL, cL);

  if (cons(L) < 0.0)
  {
    L = L + 2.0 * pi;
  }

  T Kq = 1.0 / (1.0 + q1 * q1 + q2 * q2);

  DACE::AlgebraicVector<T> ex(3), ey(3), ez(3);

  ex = {Kq * (1.0 - q1 * q1 + q2 * q2), Kq * 2.0 * q1 * q2, Kq * -2.0 * q1};
  ey = {Kq * (2.0 * q1 * q2), Kq * (1.0 + q1 * q1 - q2 * q2), Kq * 2.0 * q2};
  ez = DACE::cross(ex, ey);

  DACE::AlgebraicVector<T> er(3), ef(3), eh(3), R(3);
  er = ex * cos(L) + ey * sin(L);
  ef = ey * cos(L) - ex * sin(L);
  eh = ez;

  T x = r * cos(L);
  T y = r * sin(L);
  R = x * ex + y * ey;

  T phi = acos(R[2]/r);
  T UJ2 = J2/2*mu/r*(Re/r)*(Re/r)*(3*cos(phi)*cos(phi)-1);



  T Ut = 0.0; //3*J2*mu*Re*Re/(2.0*r*r*r*r)*r_dot; 

  T Fh = -1.5*J2*mu*Re*Re/(r*r*r*r)*sin(2*inc)*sin(AOP+theta);
  T Fr = -1.5*J2*mu*Re*Re/(r*r*r*r)*(1- 3*sin(inc)*sin(inc)*sin(AOP+theta)*sin(AOP+theta)); 
  T Pr = 0.0;
  T Pf = 0.0;


  // Perturbations
  // Variables that required perturbations
  T h = sqrt(c * c - 2.0 * r * r * UJ2); // angular momentum magnitude
  T omegh = q1 * cos(L) - q2 * sin(L); // angular velocity of equinoctial frame (h direction tilde)
  T eps_dot =  Ut + r_dot*Pr + h/r*Pf;   // total energy time derivative

  // Recurring parameters
  T A = (h - c) / r / r;
  T B =  r/h*omegh*Fh;
  T CAP = (2.0 * UJ2 -r*Fr) / c;
  T D = r / mu * eps_dot;

  // Non-dimensional variables
  T sigma = r / rho;
  T xi = 1.0 + sigma;

  dxdt[0] = -3 * pow((nu / mu / mu), (1.0 / 3.0)) * eps_dot;
  dxdt[1] = p2 * (A - B) + CAP * (r * r_dot / c * p1 + xi * p2 + sigma * cos(L)) + D * (sigma * p1 + xi * sin(L));
  dxdt[2] = p1 * (B - A) + CAP * (r * r_dot / c * p2 - xi * p1 - sigma * sin(L)) + D * (sigma * p2 + xi * cos(L));
  dxdt[3] = nu + A - B + D * r_dot * c / mu * alpha * xi + CAP * (1 / alpha + alpha * (1 - r / a));

  T Kqdot =r/2/h/Kq*Fh;
  dxdt[4] = Kqdot * sin(L);
  dxdt[5] = Kqdot * cos(L);

 
  return dxdt;


}


#endif