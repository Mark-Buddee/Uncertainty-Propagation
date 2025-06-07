#ifndef convert_H
#define convert_H
#include <math.h> 

 DACE::DA constrainAngle( DACE::DA x){

     double pi = 4 * atan(1.0);

      while (cons(x) > pi) {
        x -= 2 * pi;
    }
    while (cons(x) < -pi) {
        x += 2 * pi;
    }
    return x;

}

template <typename T> DACE::AlgebraicVector<T> eci2Kep(DACE::AlgebraicVector<T> x, const double mu) {
  /*member function to convert ECI state vector into keplerian classical element
  !< AlgebraicVector of ECI reference frame  {x, y, z, dx, dy, dz}
  !> return AlgebraicVector of keplerian element res = {a, e , i, RA, PA, TA}
     RA: rigth ascension of ascending node; PA: argument of periapsis; TA: true anomaly; i: orbital inclination*/
    using DACE::cons;
    //const double mu=398600.0;

    DACE::AlgebraicVector<T> pos(3), vel(3);

    pos[0] = x[0];
    pos[1] = x[1];
    pos[2] = x[2];

    vel[0] = x[3];
    vel[1] = x[4];
    vel[2] = x[5];

    T rm = vnorm(pos);
    T vm = vnorm(vel);

    DACE::AlgebraicVector<T> H(3), ec(3), Hver(3);

    H = pos.cross(vel);
    T Hm = vnorm(H);
    Hver = H/Hm;

    ec = -pos/rm + vel.cross(H)/mu;

    T a = 1.0/(2.0/rm - vm*vm/mu);

    T h = H[0]/(1.0+H[2]);
    T k = -H[1]/(1.0+H[2]);

    T in = atan2(sqrt(H[0]*H[0] + H[1]*H[1]), H[2] );//orbit inclination

    DACE::AlgebraicVector<double> ek(3); ek[0]= 0.0;ek[1]= 0.0;ek[2]= 1.0;
    DACE::AlgebraicVector<T> N = ek.cross(Hver);

    T RA = atan2(N[1], N[0] ); //right ascension
    if (abs(cons(RA)) < 1e-9) RA = RA - cons(RA);
    RA = RA - floor(cons(RA)/(2.0*M_PI))*2.0*M_PI;

    T PA = acos(N.dot(ec)/(vnorm(N)*vnorm(ec)) ); //periapsis argument
    if ( cons(ec)[2] < 0 ) { PA = 2.0*M_PI - PA; }

    T TA = acos(ec.dot(pos)/(vnorm(ec)*vnorm(pos)) ); //true anomaly
    if (cons(pos.dot(vel) ) < 0 ) {TA = 2.0*M_PI - TA; }


    RA = constrainAngle( RA);
    PA = constrainAngle( PA);
    TA = constrainAngle( TA);

    DACE::AlgebraicVector<T> res(6);
    res[0] = a;
    res[1] = vnorm(ec);
    res[2] = in;
    res[3] = RA;
    res[4] = PA;
    res[5] = TA;

    return res;
  }



template <typename T> DACE::AlgebraicVector<T> Kep2eci(DACE::AlgebraicVector<T> v, const double mu) {
  /*member function to convert keplerian  classical element into Earth-Centred inertial reference frame element
  !< keplerian element v = {a, e , i, RA, PA, TA}
     RA: rigth ascension of ascending node; PA: argument of periapsis; TA: true anomaly; i:orbital inclination
  !> return AlgebraicVector of ECI reference frame res = {x, y, z, dx, dy, dz}*/
    //const double mu = 398600.0;

    T p = v[0]*(1.0 - v[1]*v[1] );
    DACE::AlgebraicVector<T> rm(3), vm(3); // position and velocity in perifocal refererence frame
    rm[0] = p*cos(v[5])/(1.0 + v[1]*cos(v[5]) );
    rm[1] = p*sin(v[5])/(1.0 + v[1]*cos(v[5]) );
    rm[2] = 0.0;
    vm[0] = -sin(v[5])*sqrt(mu/p);
    vm[1] = (v[1] + cos(v[5]))*sqrt(mu/p);
    vm[2] = 0.0;

    T cRA = cos(v[3]);  T sRA = sin(v[3]);
    T cPA = cos(v[4]);  T sPA = sin(v[4]);
    T ci = cos(v[2]);  T si = sin(v[2]);

    T RR[3][3]; // rotational matrix from perifocal to eci reference frame
    RR[0][0] = cRA*cPA-sRA*ci*sPA;  RR[0][1] = -cRA*sPA-sRA*ci*cPA; RR[0][2] = sRA*si;
    RR[1][0] = sRA*cPA+cRA*ci*sPA;  RR[1][1] = -sRA*sPA+cRA*ci*cPA; RR[1][2] = -cRA*si;
    RR[2][0] = si*sPA;              RR[2][1] = si*cPA;              RR[2][2] = ci;

    DACE::AlgebraicVector<T> rr(3),vv(3);
    for(unsigned int i=0;i<3;i++){
        rr[i]=0.0;
        vv[i]=0.0;
        for(unsigned int j=0;j<3;j++){
            rr[i]=rr[i]+RR[i][j]*rm[j];
            vv[i]=vv[i]+RR[i][j]*vm[j];
        }
    }

    DACE::AlgebraicVector<T> res(6);
    res[0]=rr[0]; res[1]=rr[1]; res[2]=rr[2];
    res[3]=vv[0]; res[4]=vv[1]; res[5]=vv[2];

    return res;
}

template <typename T> DACE::AlgebraicVector<T> eci2meewsma(DACE::AlgebraicVector<T> x, const double mu) {
  /*function to convert Earth-Centred Inertial into Modified Equinoctial Elements elements
  !> AlgebraicVector ECI element {x, y , z, dx, dy, dz}
  !< return AlgebraicVector MEE element res = {p, f, g, h, k, L}*/
    //const double mu=398600.0;

    DACE::AlgebraicVector<T> pos(3), vel(3);

    pos[0] = x[0];
    pos[1] = x[1];
    pos[2] = x[2];

    vel[0] = x[3];
    vel[1] = x[4];
    vel[2] = x[5];

    T rm = vnorm(pos);

    DACE::AlgebraicVector<T> ef(3), eg(3), ew(3),er(3),ev(3), ec(3), H(3);
    H = pos.cross(vel);
    ew = H/vnorm(H);

    T p = sqr(vnorm(H))/mu;

    er = pos/rm;
    ev = (rm*vel - pos.dot(vel)*pos/rm)/vnorm(H);

    ec = -er + vel.cross(H)/mu;

    T k = ew[0]/(1.0 + ew[2]);
    T h = -ew[1]/(1.0 + ew[2]);

    T dum=1/(1+DACE::sqr(k)+DACE::sqr(h));


    T AA[3][3];
    AA[0][0] = dum*(1-DACE::sqr(k)+DACE::sqr(h)); AA[0][1] = dum*2*k*h;                         AA[0][2] = dum*2*k;
    AA[1][0] = dum*2*k*h;                         AA[1][1] = dum*(1+DACE::sqr(k)-DACE::sqr(h)); AA[1][2] = dum*(-2*h);
    AA[2][0] = dum*(-2*k);                        AA[2][1] = dum*2*h;                           AA[2][2] = dum*(1-DACE::sqr(k)-sqr(h));

    ef[0] = AA[0][0]; ef[1] = AA[1][0]; ef[2] = AA[2][0];
    eg[0] = AA[0][1]; eg[1] = AA[1][1]; eg[2] = AA[2][1];

    T f = ec.dot(ef);
    T g = ec.dot(eg);

    T cosL = er[0] + ev[1];
    T sinL = er[1] - ev[0];

    T L = atan(sinL/cosL);
    if ( DACE::cons(cosL) < 0.0 ) { L += M_PI; }
    if ( DACE::cons(L) < 0.0 ) { L += 2.0*M_PI; } //only to have alpha \in [0,2pi]

    //if ( cons(sinL) < 0.0 ) { L = 2.0*M_PI - L; }

    DACE::AlgebraicVector<T> res(6);

    T eccen = sqrt(f*f + g*g);
    T a = p/(1-eccen*eccen);

    res[0] = a;
    res[1] = f;
    res[2] = g;
    res[3] = h;
    res[4] = k;
    res[5] = L;
    //cout << cons(cos(L)) << " " << cons(sin(L))<< endl;

    return res;

}

template<typename T> DACE::AlgebraicVector<T> RTN2ECI(DACE::AlgebraicVector<T> rv, DACE::AlgebraicVector<T> dvRTN)

{

    DACE::AlgebraicVector<T> r(3), v(3), rhat(3), Nhat(3), rcrossv(3), that(3), dvECI(3); 

    DACE::AlgebraicMatrix<T> RTNmatrix(3,3); 


    for (unsigned int i = 0; i < 3; i++){
        r[i] = rv[i];
        v[i] = rv[i+3];
            
    }

    T rnorm = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);

    rhat = r/rnorm; 

    rcrossv = r.cross(v); 

    T rcrossvnorm = sqrt(rcrossv[0]*rcrossv[0] + rcrossv[1]*rcrossv[1] + rcrossv[2]*rcrossv[2]);

    Nhat = rcrossv/rcrossvnorm;

    that = Nhat.cross(rhat);
    //that = rhat.cross(Nhat);

    for(unsigned int j = 0; j < 3; j++){
        RTNmatrix.at(0,j) = rhat[j];
        RTNmatrix.at(1,j) = that[j]; 
        RTNmatrix.at(2,j) = Nhat[j];
    }
//     cout << cons(rhat) << endl;
//     cout << cons(RTNmatrix) << endl;


    dvECI = DACE::transpose(RTNmatrix)*dvRTN;

    return dvECI;

}


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
    T RAAN = atan(dot(uln,j)/dot(uln,i));
    if ( DACE::cons(dot(uln,i)) < 0.0 ) { RAAN += M_PI; }
    if ( DACE::cons(RAAN) < 0.0 ) { RAAN += 2.0*M_PI; } //only to have alpha \in [0,2pi]



    if (cons(RAAN) < 0.0)RAAN = RAAN + 2*pi;

  
    T AOP = atan(dot(-u2,uln)/dot(uln,u1));
    if ( DACE::cons(dot(uln,u1)) < 0.0 ) { AOP += M_PI; }
    if ( DACE::cons(AOP) < 0.0 ) { AOP += 2.0*M_PI; } //only to have alpha \in [0,2pi]

    if (cons(AOP) < 0.0)AOP = AOP + 2*pi;

    T TA = atan(dot(R,u2)/dot(R,u1));
    if ( DACE::cons(dot(R,u1)) < 0.0 ) { TA += M_PI; }
    if ( DACE::cons(TA) < 0.0 ) { TA += 2.0*M_PI; } //only to have alpha \in [0,2pi]

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



template <typename T> DACE::AlgebraicVector<T> Kep2meesma(DACE::AlgebraicVector<T> v) {
  /*member function to convert keplerian  classical element into MEE reference frame element
  !> keplerian element v = {a, e , i, RA, PA, TA}
     OMEGA: rigth ascension of ascending node; omega: argument of periapsis; M: mean anomaly;i:orbital inclination
  !< return AlgebraicVector of MEE reference frame res = {p, f, g, h, k, L}*/

    T p = v[0]*(1.0 - DACE::sqr(v[1]) );
    T f = v[1]*cos(v[3] + v[4]);
    T g = v[1]*sin(v[3] + v[4]);
    T h = tan(v[2]/2.0)*cos(v[3]);
    T k = tan(v[2]/2.0)*sin(v[3]);

    T L = v[5] + v[3] + v[4];

    DACE::AlgebraicVector<T> res(6);
    res[0] = v[0];
    res[1] = f;
    res[2] = g;
    res[3] = h;
    res[4] = k;
    res[5] = L;

    return res;

}


template <typename T> DACE::AlgebraicVector<T> eci2mee(DACE::AlgebraicVector<T> x, const double mu) {
  /*function to convert Earth-Centred Inertial into Modified Equinoctial Elements elements
  !> AlgebraicVector ECI element {x, y , z, dx, dy, dz}
  !< return AlgebraicVector MEE element res = {p, f, g, h, k, L}*/
    //const double mu=398600.0;

    DACE::AlgebraicVector<T> pos(3), vel(3);

    pos[0] = x[0];
    pos[1] = x[1];
    pos[2] = x[2];

    vel[0] = x[3];
    vel[1] = x[4];
    vel[2] = x[5];

    T rm = vnorm(pos);

    DACE::AlgebraicVector<T> ef(3), eg(3), ew(3),er(3),ev(3), ec(3), H(3);
    H = pos.cross(vel);
    ew = H/vnorm(H);

    T p = sqr(vnorm(H))/mu;

    er = pos/rm;
    ev = (rm*vel - pos.dot(vel)*pos/rm)/vnorm(H);

    ec = -er + vel.cross(H)/mu;

    T k = ew[0]/(1.0 + ew[2]);
    T h = -ew[1]/(1.0 + ew[2]);

    T dum=1/(1+DACE::sqr(k)+DACE::sqr(h));

    T AA[3][3];
    AA[0][0] = dum*(1-DACE::sqr(k)+DACE::sqr(h)); AA[0][1] = dum*2*k*h;                         AA[0][2] = dum*2*k;
    AA[1][0] = dum*2*k*h;                         AA[1][1] = dum*(1+DACE::sqr(k)-DACE::sqr(h)); AA[1][2] = dum*(-2*h);
    AA[2][0] = dum*(-2*k);                        AA[2][1] = dum*2*h;                           AA[2][2] = dum*(1-DACE::sqr(k)-sqr(h));

    ef[0] = AA[0][0]; ef[1] = AA[1][0]; ef[2] = AA[2][0];
    eg[0] = AA[0][1]; eg[1] = AA[1][1]; eg[2] = AA[2][1];

    T f = ec.dot(ef);
    T g = ec.dot(eg);

    T cosL = er[0] + ev[1];
    T sinL = er[1] - ev[0];

    
    T L = atan(sinL/cosL);
    if ( DACE::cons(cosL) < 0.0 ) { L += M_PI; }
    if ( DACE::cons(L) < 0.0 ) { L += 2.0*M_PI; } //only to have alpha \in [0,2pi]

    //if ( cons(sinL) < 0.0 ) { L = 2.0*M_PI - L; }

    DACE::AlgebraicVector<T> res(6);
    res[0] = p;
    res[1] = f;
    res[2] = g;
    res[3] = h;
    res[4] = k;
    res[5] = L;
    //cout << cons(cos(L)) << " " << cons(sin(L))<< endl;

    return res;

}

template <typename T> DACE::AlgebraicVector<T> mee2eci(DACE::AlgebraicVector<T> mee, const double mu) {
    /*function to convert Modified Equinoctial Elements into Earth-Centred Inertial elements
    !> mee = {p
       f=|e|*cos(RA+PA)
       g=|e|*sin(RA+PA)
       h=tan(i/2)*cos(RA)
       k=tan(i/2)*sin(RA)
       L=TA+RA+PA}
       RA: rigth ascension of ascending node; PA: argument of periapsis; TA: true anomaly; i:orbital inclination
   !< return AlgebraicVector ECI element res = {x, y, z, dx, dy, dz}
      */
    //const double mu=398600.0;  //km^3/s^2

    T p=mee[0];
    T f=mee[1];
    T g=mee[2];
    T h=mee[3];
    T k=mee[4];
    T L=mee[5];
    //cout << cons(cos(L)) << " " << cons(sin(L))<< endl;

    T e = sqrt(f*f + g*g);

    T H  = sqrt(p*mu);
    T rm = p/(1.0 + f*cos(L) + g*sin(L) );

    T dum=1/(1+DACE::sqr(h)+DACE::sqr(k));

    DACE::AlgebraicVector<T> r(3),v(3);
    r[0] = rm*(cos(L) + (h*h - k*k)*cos(L) + 2.0*h*k*sin(L) )*dum;
    r[1] = rm*(sin(L) - (h*h - k*k)*sin(L) + 2.0*h*k*cos(L) )*dum;
    r[2] = 2.0*rm*(h*sin(L) - k*cos(L) )*dum;
    v[0] = -sqrt(mu/p)*(sin(L) + (h*h - k*k)*sin(L) - 2.0*h*k*cos(L) + g - 2.0*f*h*k + (h*h - k*k)*g )*dum;
    v[1] = -sqrt(mu/p)*(-cos(L) + (h*h - k*k)*cos(L) + 2.0*h*k*sin(L) - f + 2.0*g*h*k + (h*h - k*k)*f )*dum;
    v[2] = 2.0*sqrt(mu/p)*(h*cos(L) + k*sin(L) + f*h + k*g )*dum;

    DACE::AlgebraicVector<T> res(6);
    res[0]=r[0]; res[1]=r[1]; res[2]=r[2];
    res[3]=v[0]; res[4]=v[1]; res[5]=v[2];

    return res;

}



template <typename T>
DACE::AlgebraicVector<T> RV2GEq(DACE::AlgebraicVector<T> RV, const double *cst){

    double pi = 4 * atan(1.0);
    const double mu = cst[5];
    const double J2 = cst[6];
    const double Re = cst[7];
    
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

    
    T L = atan(sL/cL);
    if ( DACE::cons(cL) < 0.0 ) { L += M_PI; }
    if ( DACE::cons(L) < 0.0 ) { L += 2.0*M_PI; } //only to have alpha \in [0,2pi]


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
 

    T K = atan(sK/cK);
    if ( DACE::cons(cK) < 0.0 ) { K += M_PI; }
    if ( DACE::cons(K) < 0.0 ) { K += 2.0*M_PI; } //only to have alpha \in [0,2pi]



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
    const double mu = cst[5];
    const double J2 = cst[6];
    const double Re = cst[7];

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

   
     T L = atan(sL/cL);
    if ( DACE::cons(cL) < 0.0 ) { L += M_PI; }
    if ( DACE::cons(L) < 0.0 ) { L += 2.0*M_PI; } //only to have alpha \in [0,2pi]


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

// template <typename T>
// DACE::AlgebraicVector<T> CEq2RV(DACE::AlgebraicVector<T> X, const double *cst)
// {
//     double pi = 4 * atan(1.0);
//     const double mu = cst[5];
//     const double J2 = cst[6];
//     const double Re = cst[7];
// 
//     DACE::AlgebraicVector<T> posvel(6);
// 
//     T nu = X[0];
//     T p1 = X[1];
//     T p2 = X[2];
//     T varL = X[3];
//     T q1 = X[4];
//     T q2 = X[5];
// 
//     // auxiliary variables
//     T a = pow((mu / nu / nu), 1.0 / 3.0);
//     T c = pow((mu * mu / nu), 1.0 / 3.0) * sqrt(1.0 - p1 * p1 - p2 * p2);
//     T alpha = 1.0 / (1.0 + sqrt(1.0 - p1 * p1 - p2 * p2));
// 
//     // generalised Kepler equation
//     double tolN = 1e-14;
//     T diff = 0.5;
//     // cout << "eps" << endl << eps << endl ;
//     T K = (varL);
//     while (cons(diff) > tolN)
//     {
//         T K_old = K;
// 
//         T fx = K + p1 * cos(K) - p2 * sin(K) - varL;
// 
//         T fxd = 1.0 - p1 * sin(K) - p2 * cos(K);
// 
//         K = K_old - (fx / fxd);
// 
//         diff = K_old - K;
//     }
// 
//     T r = a * (1.0 - p1 * sin(K) - p2 * cos(K));
// 
//   
//     T r_dot = sqrt(mu * a) / r * (p2 * sin(K) - p1 * cos(K));
// 
//     T sL = a / r * (alpha * p1 * p2 * cos(K) + (1.0 - alpha * p2 * p2) * sin(K) - p1);
//     T cL = a / r * (alpha * p1 * p2 * sin(K) + (1.0 - alpha * p1 * p1) * cos(K) - p2);
// 
// 
//     T L = atan(sL/cL);
//     if ( DACE::cons(cL) < 0.0 ) { L += M_PI; }
//     if ( DACE::cons(L) < 0.0 ) { L += 2.0*M_PI; } //only to have alpha \in [0,2pi]
// 
// 
//     if (cons(L) < 0.0)
//     {
//         L = L + 2.0 * pi;
//     }
// 
//     T Kq = 1.0 / (1.0 + q1 * q1 + q2 * q2);
// 
//     DACE::AlgebraicVector<T> ex(3), ey(3), ez(3);
// 
//     ex = {Kq * (1.0 - q1 * q1 + q2 * q2), Kq * 2.0 * q1 * q2, Kq * -2.0 * q1};
//     ey = {Kq * (2.0 * q1 * q2), Kq * (1.0 + q1 * q1 - q2 * q2), Kq * 2.0 * q2};
// 
//     DACE::AlgebraicVector<T> er(3), ef(3), eh(3), R(3), V(3);
//     er = ex * cos(L) + ey * sin(L);
//     ef = ey * cos(L) - ex * sin(L);
// 
//     R = r * er;
// 
//     T phi = acos(R[2]/r); 
//     T UJ2 = 0.0;
// 
//     T h = sqrt(c * c - 2.0 * r * r * UJ2); // angular momentum magnitude
// 
//     V = r_dot*er + h/r * ef;
// 
//     posvel[0] = R[0];
//     posvel[1] = R[1];
//     posvel[2] = R[2];
//     posvel[3] = V[0];
//     posvel[4] = V[1];
//     posvel[5] = V[2];
//     return posvel;
// }

// template <typename T>
// DACE::AlgebraicVector<T> RV2CEq(DACE::AlgebraicVector<T> RV, const double *cst){
// 
//     double pi = 4 * atan(1.0);
//     const double mu = cst[5];
//     const double J2 = cst[6];
//     const double Re = cst[7];
//     
//     DACE::AlgebraicVector<T> kep = RV2OE(RV, mu);
// 
//    // cout << kep << endl;
//     T inc = kep[2];
//     T Raan = kep[3]; 
//     
//     DACE::AlgebraicVector<T> R(3), V(3), X(6);
//     
//     R[0] = RV[0];
//     R[1] = RV[1];
//     R[2] = RV[2];
//     
//     V[0] = RV[3];
//     V[1] = RV[4];
//     V[2] = RV[5];
//     T r = vnorm(R);
//     T v = vnorm(V);
//     T h = vnorm(cross(R,V));
//     T r_dot = dot(R,V)/r;
// 
//     T phi = acos(R[2]/r); 
//     T UJ2 = 0.0;
// 
//     T epsk = 0.5*v*v - mu/r; 
//     T eps = epsk + UJ2;  
// 
//     T nu = 1.0/mu*pow((-2.0*eps),1.5);       
//     
//     // Calculate q in order to obtain equinoctial frame
//     T q1 = tan(inc/2.0)*sin(Raan);  // (5)
//     T q2 = tan(inc/2.0)*cos(Raan);  // (6)
// 
//     T Kq = 1.0 / (1.0 + q1 * q1 + q2 * q2);
// 
//     DACE::AlgebraicVector<T> ex(3), ey(3), er(3);
// 
//     ex = {Kq * (1 - q1 * q1 + q2 * q2), Kq * 2 * q1 * q2, Kq * -2 * q1};
//     ey = {Kq * (2 * q1 * q2), Kq * (1 + q1 * q1 - q2 * q2), Kq * 2 * q2};
//     er = R/r;
// 
//     T cL = dot(er,ex);
//     T sL = dot(er, ey);
// 
//     T L = atan(sL/cL);
//     if ( DACE::cons(cL) < 0.0 ) { L += M_PI; }
//     if ( DACE::cons(L) < 0.0 ) { L += 2.0*M_PI; } //only to have alpha \in [0,2pi]
// 
// 
//     if (cons(L) < 0)
//     {
//         L = L + 2 * pi;
//     }
// 
//     // Calculate rho,c and a
//     T Ueff = h*h/2/r/r + UJ2;   // effective potential
//     T c = sqrt(2*r*r*Ueff);
//     T rho = c*c/mu;
//     T a = -mu/2/eps;
// 
//     T p1 = (rho/r-1)*sin(L) - c*r_dot/mu*cos(L);  // (2)
//     T p2 = (rho/r-1)*cos(L) + c*r_dot/mu*sin(L);  // (3)
//     
//     // Obtain K for varL
//     T w = sqrt(mu/a);
//     T sK = (mu+c*w-r*r_dot*r_dot)*sin(L) - r_dot*(c+w*r)*cos(L);
//     T cK = (mu+c*w-r*r_dot*r_dot)*cos(L) + r_dot*(c+w*r)*sin(L);
// 
//     T K = atan(sK/cK);
//     if ( DACE::cons(cK) < 0.0 ) { K += M_PI; }
//     if ( DACE::cons(K) < 0.0 ) { K += 2.0*M_PI; } //only to have alpha \in [0,2pi]
// 
// 
//     if (cons(K) < 0.0)    K = K + 2*pi;
// 
//     // Obtain generalized mean longitud varL
//     T varL = K + 1/(mu+c*w)*(cK*p1-sK*p2);  // (4)
//     
//     // GEqOE array
//     X[0]= nu; 
//     X[1]= p1;
//     X[2]= p2; 
//     X[3]= varL; 
//     X[4]= q1; 
//     X[5] = q2; 
//     return X;
// }
// 


template <typename T>
DACE::AlgebraicVector<T> RV2CEq(DACE::AlgebraicVector<T> RV, const double *cst){

    double pi = 4 * atan(1.0);
    const double mu = cst[5];
    const double J2 = cst[6];
    const double Re = cst[7];
    
    
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
DACE::AlgebraicVector<T> CEq2RV(DACE::AlgebraicVector<T> X, const double *cst)
{
    double pi = 4 * atan(1.0);
   const double mu = cst[5];
    const double J2 = cst[6];
    const double Re = cst[7];
    

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