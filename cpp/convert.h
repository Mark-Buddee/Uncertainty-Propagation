#ifndef convert_H
#define convert_H

template<typename T> DACE::AlgebraicVector<T> RV2OE(DACE::AlgebraicVector<T> rv, const double mu)
{

	DACE::AlgebraicVector<T> kep(6);
    double pi = 4 * atan(1.0);

	DACE::AlgebraicVector<T> R(3), V(3);
	for (int i = 0; i < 3; i++)
	{
		R[i] = rv[i];
		V[i] = rv[i + 3];
	}

    T r = vnorm(R);
    T v = vnorm(V);

    DACE::AlgebraicVector<T> H = cross(R, V);
    T h = vnorm(H);

    DACE::AlgebraicVector<T> E = -R/r -1/mu*cross(H,V);
    T e = vnorm(E);

    // periforcal frame 
    DACE::AlgebraicVector<T> u1(3), u2(3), u3(3), i(3), j(3), k(3);

    u1 = E/e; 
    u3 = H/h; 
    u2 = cross(u3,u1);

    // intertial frame
    i[0] = 1.0; 
    j[1]= 1.0;
    k[2] = 1.0;

    DACE::AlgebraicVector<T> uln = cross(k, u3); //node line 

    //angles from perifocal frame 
    T INC = acos(dot(u3,k));
    T RAAN = atan2(dot(uln,j),dot(uln,i));

    if (cons(RAAN) < 0.0)RAAN = RAAN + 2*pi;

    T AOP = atan2(dot(-u2,uln),dot(uln,u1));  // argument of perigee
    if (cons(AOP) < 0.0)AOP = AOP + 2*pi;

    T TA = atan2(dot(R,u2),dot(R,u1));    // true anomaly
    if (cons(TA) < 0.0)TA = TA + 2*pi;

    T u = 2*atan(sqrt((1-e)/(1+e))*tan(TA/2));    // eccentric anomaly

    T MA = u - e*sin(u); // mean anomaly

    T a = r/(2.0 - v*v*r/mu);

    kep[0] = a; 
    kep[1] = e; 
    kep[2] = INC;
    kep[3] = RAAN;
    kep[4] = AOP;
    kep[5] = TA;

    return kep; 

}



// template <typename T>
// DACE::AlgebraicVector<T> CEq2GEq(DACE::AlgebraicVector<T> XC, const double *cst){
// 
//     double pi = 4 * atan(1.0);
//     const double mu = cst[0];
//     const double J2 = cst[1];
//     const double Re = cst[2];
//     T nuC = XC[0];
//     T p1C = XC[1];
//     T p2C = XC[2];
//     T varLC = XC[3];
//     T q1C = XC[4];
//     T q2C = XC[5];
// 
//     T epsk = -pow((mu*nuC, 2.0/3.0))/2.0; 
//     T eps = epsk + UJ2;  
// 
//     T UJ2 = J2/2*mu/r*(Re/r)*(Re/r)*(3*cos(phi)*cos(phi)-1);
// 
// 
// 
// 
// 
//     // GEqOE array
//     X[0]= nu; 
//     X[1]= p1C;
//     X[2]= p2C; 
//     X[3]= varL; 
//     X[4]= q1C; 
//     X[5] = q2C; 
// 
// 
// 
// 
//     return X;
//     
// }


template <typename T>
DACE::AlgebraicVector<T> RV2CEq(DACE::AlgebraicVector<T> RV, const double *cst){

    double pi = 4 * atan(1.0);
    const double mu = cst[0];
    const double J2 = cst[1];
    const double Re = cst[2];
    
    DACE::AlgebraicVector<T> kep = RV2OE(RV, mu);

   // cout << kep << endl;
    T inc = kep[2];
    T Raan = kep[3]; 
    
    DACE::AlgebraicVector<T> R(3), V(3), X(6);
    
    R[0] = RV[0];
    R[1] = RV[1];
    R[2] = RV[2];
    
    V[0] = RV[3];
    V[1] = RV[4];
    V[2] = RV[5];
    T r = vnorm(R);
    T v = vnorm(V);
    T h = vnorm(cross(R,V));
    T r_dot = dot(R,V)/r;

    T phi = acos(R[2]/r); 
    T UJ2 = 0.0;

    T epsk = 0.5*v*v - mu/r; 
    T eps = epsk + UJ2;  

    T nu = 1.0/mu*pow((-2.0*eps),1.5);       
    
    // Calculate q in order to obtain equinoctial frame
    T q1 = tan(inc/2.0)*sin(Raan);  // (5)
    T q2 = tan(inc/2.0)*cos(Raan);  // (6)

    T Kq = 1.0 / (1.0 + q1 * q1 + q2 * q2);

    DACE::AlgebraicVector<T> ex(3), ey(3), er(3);

    ex = {Kq * (1 - q1 * q1 + q2 * q2), Kq * 2 * q1 * q2, Kq * -2 * q1};
    ey = {Kq * (2 * q1 * q2), Kq * (1 + q1 * q1 - q2 * q2), Kq * 2 * q2};
    er = R/r;

    T cL = dot(er,ex);
    T sL = dot(er, ey);

    T L = atan2(sL,cL);

    if (cons(L) < 0)
    {
        L = L + 2 * pi;
    }

    // Calculate rho,c and a
    T Ueff = h*h/2/r/r + UJ2;   // effective potential
    T c = sqrt(2*r*r*Ueff);
    T rho = c*c/mu;
    T a = -mu/2/eps;

    T p1 = (rho/r-1)*sin(L) - c*r_dot/mu*cos(L);  // (2)
    T p2 = (rho/r-1)*cos(L) + c*r_dot/mu*sin(L);  // (3)
    
    // Obtain K for varL
    T w = sqrt(mu/a);
    T sK = (mu+c*w-r*r_dot*r_dot)*sin(L) - r_dot*(c+w*r)*cos(L);
    T cK = (mu+c*w-r*r_dot*r_dot)*cos(L) + r_dot*(c+w*r)*sin(L);
    T K = atan2(sK,cK);

    if (cons(K) < 0.0)    K = K + 2*pi;

    // Obtain generalized mean longitud varL
    T varL = K + 1/(mu+c*w)*(cK*p1-sK*p2);  // (4)
    
    // GEqOE array
    X[0]= nu; 
    X[1]= p1;
    X[2]= p2; 
    X[3]= varL; 
    X[4]= q1; 
    X[5] = q2; 
    return X;
}


template <typename T>
DACE::AlgebraicVector<T> RV2GEq(DACE::AlgebraicVector<T> RV, const double *cst){

    double pi = 4 * atan(1.0);
    const double mu = cst[0];
    const double J2 = cst[1];
    const double Re = cst[2];
    
    DACE::AlgebraicVector<T> kep = RV2OE(RV, mu);

   // cout << kep << endl;
    T inc = kep[2];
    T Raan = kep[3]; 
    
    DACE::AlgebraicVector<T> R(3), V(3), X(6);
    
    R[0] = RV[0];
    R[1] = RV[1];
    R[2] = RV[2];
    
    V[0] = RV[3];
    V[1] = RV[4];
    V[2] = RV[5];
    T r = vnorm(R);
    T v = vnorm(V);
    T h = vnorm(cross(R,V));
    T r_dot = dot(R,V)/r;

    T phi = acos(R[2]/r); 
    T UJ2 = J2/2*mu/r*(Re/r)*(Re/r)*(3*cos(phi)*cos(phi)-1);

    T epsk = 0.5*v*v - mu/r; 
    T eps = epsk + UJ2;  

    T nu = 1.0/mu*pow((-2.0*eps),1.5);       
    
    // Calculate q in order to obtain equinoctial frame
    T q1 = tan(inc/2.0)*sin(Raan);  // (5)
    T q2 = tan(inc/2.0)*cos(Raan);  // (6)

    T Kq = 1.0 / (1.0 + q1 * q1 + q2 * q2);

    DACE::AlgebraicVector<T> ex(3), ey(3), er(3);

    ex = {Kq * (1 - q1 * q1 + q2 * q2), Kq * 2 * q1 * q2, Kq * -2 * q1};
    ey = {Kq * (2 * q1 * q2), Kq * (1 + q1 * q1 - q2 * q2), Kq * 2 * q2};
    er = R/r;

    T cL = dot(er,ex);
    T sL = dot(er, ey);

    T L = atan2(sL,cL);

    if (cons(L) < 0)
    {
        L = L + 2 * pi;
    }

    
    // Calculate rho,c and a
    T Ueff = h*h/2/r/r + UJ2;   // effective potential
    T c = sqrt(2*r*r*Ueff);
    T rho = c*c/mu;
    T a = -mu/2/eps;

    T p1 = (rho/r-1)*sin(L) - c*r_dot/mu*cos(L);  // (2)
    T p2 = (rho/r-1)*cos(L) + c*r_dot/mu*sin(L);  // (3)
    
    // Obtain K for varL
    T w = sqrt(mu/a);
    T sK = (mu+c*w-r*r_dot*r_dot)*sin(L) - r_dot*(c+w*r)*cos(L);
    T cK = (mu+c*w-r*r_dot*r_dot)*cos(L) + r_dot*(c+w*r)*sin(L);
    T K = atan2(sK,cK);

    if (cons(K) < 0.0)    K = K + 2*pi;

    // Obtain generalized mean longitud varL
    T varL = K + 1/(mu+c*w)*(cK*p1-sK*p2);  // (4)
    
    // GEqOE array
    X[0]= nu; 
    X[1]= p1;
    X[2]= p2; 
    X[3]= varL; 
    X[4]= q1; 
    X[5] = q2; 




    return X;
}
template <typename T>
DACE::AlgebraicVector<T> GEq2RV(DACE::AlgebraicVector<T> X, const double *cst)
{

    //cout << X << endl;
    double pi = 4 * atan(1.0);
    const double mu = cst[0];
    const double J2 = cst[1];
    const double Re = cst[2];

    DACE::AlgebraicVector<T> posvel(6);

    T nu = X[0];
    T p1 = X[1];
    T p2 = X[2];
    T varL = X[3];
    T q1 = X[4];
    T q2 = X[5];

    // auxiliary variables
    T a = pow((mu / nu / nu), 1.0 / 3.0);

    T c = pow((mu * mu / nu), 1.0 / 3.0) * sqrt(1.0 - p1 * p1 - p2 * p2);

    T alpha = 1.0 / (1.0 + sqrt(1.0 - p1 * p1 - p2 * p2));

    // generalised Kepler equation
    double tolN = 1e-12;
    T diff = 0.5;
    // cout << "eps" << endl << eps << endl ;
    T K = (varL);
    while (cons(diff) > tolN)
    {
        T K_old = K;

        T fx = K + p1 * cos(K) - p2 * sin(K) - varL;

        T fxd = 1.0 - p1 * sin(K) - p2 * cos(K);

        K = K_old - (fx / fxd);

        diff = K_old - K;
    }

    
    T r = a * (1.0 - p1 * sin(K) - p2 * cos(K));
  
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

    DACE::AlgebraicVector<T> er(3), ef(3), eh(3), R(3), V(3);
    er = ex * cos(L) + ey * sin(L);
    ef = ey * cos(L) - ex * sin(L);

    R = r * er;

    T phi = acos(R[2]/r); 
    T UJ2 = J2/2*mu/r*(Re/r)*(Re/r)*(3*cos(phi)*cos(phi)-1);

    T h = sqrt(c * c - 2.0 * r * r * UJ2); // angular momentum magnitude

    V = r_dot*er + h/r * ef;

    posvel[0] = R[0];
    posvel[1] = R[1];
    posvel[2] = R[2];
    posvel[3] = V[0];
    posvel[4] = V[1];
    posvel[5] = V[2];
    return posvel;
}

template <typename T>
DACE::AlgebraicVector<T> CEq2RV(DACE::AlgebraicVector<T> X, const double *cst)
{
    double pi = 4 * atan(1.0);
        const double mu = cst[0];
    const double J2 = cst[1];
    const double Re = cst[2];

    DACE::AlgebraicVector<T> posvel(6);

    T nu = X[0];
    T p1 = X[1];
    T p2 = X[2];
    T varL = X[3];
    T q1 = X[4];
    T q2 = X[5];

    // auxiliary variables
    T a = pow((mu / nu / nu), 1.0 / 3.0);
    T c = pow((mu * mu / nu), 1.0 / 3.0) * sqrt(1.0 - p1 * p1 - p2 * p2);
    T alpha = 1.0 / (1.0 + sqrt(1.0 - p1 * p1 - p2 * p2));

    // generalised Kepler equation
    double tolN = 1e-14;
    T diff = 0.5;
    // cout << "eps" << endl << eps << endl ;
    T K = (varL);
    while (cons(diff) > tolN)
    {
        T K_old = K;

        T fx = K + p1 * cos(K) - p2 * sin(K) - varL;

        T fxd = 1.0 - p1 * sin(K) - p2 * cos(K);

        K = K_old - (fx / fxd);

        diff = K_old - K;
    }

    T r = a * (1.0 - p1 * sin(K) - p2 * cos(K));

  
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

    DACE::AlgebraicVector<T> er(3), ef(3), eh(3), R(3), V(3);
    er = ex * cos(L) + ey * sin(L);
    ef = ey * cos(L) - ex * sin(L);

    R = r * er;

    T phi = acos(R[2]/r); 
    T UJ2 = 0.0;

    T h = sqrt(c * c - 2.0 * r * r * UJ2); // angular momentum magnitude

    V = r_dot*er + h/r * ef;

    posvel[0] = R[0];
    posvel[1] = R[1];
    posvel[2] = R[2];
    posvel[3] = V[0];
    posvel[4] = V[1];
    posvel[5] = V[2];
    return posvel;
}


#endif