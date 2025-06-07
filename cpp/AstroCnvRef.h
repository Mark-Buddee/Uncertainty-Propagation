/********************************************************************************************/
//  AstroCnvRef.h                                                                            /
//  Astrodynamic function for the Conversion of reference frame from ane the each other      /
//                                                                                           /
//                                                                                           /
//  Created by Daniele Antonio Santeramo on 28/11/16.                                        /
//                                                                                           /
/********************************************************************************************/

#ifndef ASTROCNVREF_H_INCLUDED_
#define ASTROCNVREF_H_INCLUDED_

#include <dace/dace.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>


using namespace std;
void multiply(std::vector<double>& v, double scalar) {
    for (int i = 0; i < v.size(); i++) {
        v[i] *= scalar;
    }
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

template <typename T> DACE::AlgebraicVector<T> eci2kep2(DACE::AlgebraicVector<T> x, const double mu) {

    T v2s = x[3]*x[3] + x[4]*x[4] + x[5]*x[5]; 
    T r = sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
    T Ws = 0.5*v2s - mu/r; 
    T a = -mu/2.0/Ws; 

    DACE::AlgebraicVector<T> L(3), pos(3), vel(3),xkep(6); 

    for (unsigned int i = 0; i < 3; i++){
        pos[i] = x[i];
        vel[i] = x[i+3];
    }

    L[0] = x[1]*x[5] - x[2]*x[4];
    L[1] = x[2]*x[3] - x[0]*x[5]; 
    L[2] = x[0]*x[4] -x[1]*x[3]; 

    T L2s = L[0]*L[0] + L[1]*L[1] + L[2]*L[2];

    T p = L2s/mu; 

    T e = sqrt((sqrt(1.0 - p/a))*(sqrt(1.0 - p/a)));

    a = sqrt(a*a); 
 

    T cosE = (1.0 - r/a)/e; 
    T sinE = (pos[0]*vel[0] + pos[1]*vel[1] + pos[2]*vel[2])/(e*sqrt(mu*a));

    T E = atan2(sinE, cosE);

    T sinI = sqrt( L[0]*L[0] + L[1]*L[1])/sqrt(L2s);
    T cosI = L[2]/sqrt(L2s); 

    T I = atan2(sinI, cosI);

    T sinw = ((vel[0]*L[1] - vel[1]*L[0])/mu - pos[2]/r)/(e*sinI); 
    T cosw = ((sqrt(L2s)*vel[2])/mu - (L[0]*pos[1] - L[1]*pos[0])/(sqrt(L2s)*r))/(e*sinI);
    T RAAN = atan2(sinw, cosw); 

    T cosO = -L[1]/(sqrt(L2s)*sinI);
    T sinO = L[0]/(sqrt(L2s)*sinI); 
    T AOP = atan2(sinO, cosO); 

    T P = 2*M_PI*sqrt(a*a*a/mu);

    T tau = -(E -e*sin(E))/sqrt(mu*pow(a,-3.0));

    T True_an = 2.0*atan(sqrt((1.0+e)/(1.0-e))*tan(E/2.0));

    xkep[0] = a; 
    xkep[1] = e; 
    xkep[2] = I; 
    xkep[3] = RAAN; 
    xkep[4] = AOP; 
    xkep[5] = True_an; 

   // cout << cons(xkep) << endl;
    return xkep; 

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

    T L = atan2(sinL, cosL);
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

template <typename T> DACE::AlgebraicVector<T> Kep2mee(DACE::AlgebraicVector<T> v) {
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
    res[0] = p;
    res[1] = f;
    res[2] = g;
    res[3] = h;
    res[4] = k;
    res[5] = L;

    return res;

}

template <typename T> DACE::AlgebraicVector<T> mee2Kep(DACE::AlgebraicVector<T> mee) {
  /*member function to convert MEE reference frame element into keplerian  classical element
  !> AlgebraicVector of MEE reference frame mee = {p, f, g, h, k, L}
  !< return AlgebraicVector of keplerian element res = {a, e , i, RA, PA, TA}
    RA: rigth ascension of ascending node; PA: argument of periapsis; TA: true anomaly; i:orbital inclination*/

  using DACE::cons;

  T p=mee[0];
  T f=mee[1];
  T g=mee[2];
  T h=mee[3];
  T k=mee[4];
  T L=mee[5];

  T e = sqrt(DACE::sqr(f) + DACE::sqr(g) );
  T a = p/(1.0 - DACE::sqr(e) );
  T i = 2.0*atan(sqrt(DACE::sqr(h) + DACE::sqr(k)) );
  i = i - floor(cons(i)/(2.0*M_PI))*2.0*M_PI;

  T RA = atan2(k, h); //rigth ascension of ascending node
  RA = RA - floor(cons(RA)/(2.0*M_PI))*2.0*M_PI;

  T RAPA = atan2(g, f); // sum of rigth ascension and periapsis argument
  RAPA = RAPA - floor(cons(RAPA)/(2.0*M_PI))*2.0*M_PI;

  T PA = RAPA - RA; // periapsis argument
  T TA = L - RAPA; // true anomaly

  DACE::AlgebraicVector<T> res(6);
  res[0] = a;
  res[1] = e;
  res[2] = i;
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

    DACE::AlgebraicVector<T> res(6);
    res[0] = a;
    res[1] = vnorm(ec);
    res[2] = in;
    res[3] = RA;
    res[4] = PA;
    res[5] = TA;


    return res;
  }


template<typename T> DACE::AlgebraicVector<T> in2orb( DACE::AlgebraicVector<T> rrin, T i, T Om, T l)
{

    // convert classica orbital elements to modified equinoctial elements


    DACE::AlgebraicVector<T> row1(3), row2(3), row3(3);
    DACE::AlgebraicVector<T> rrorb(3);

    row1[0] = (cos(l)*cos(Om)-sin(l)*cos(i)*sin(Om));
    row1[1] = (cos(l)*sin(Om)+sin(l)*cos(i)*cos(Om));
    row1[2] = (sin(l)*sin(i));
    row2[0] = (-sin(l)*cos(Om)-cos(l)*cos(i)*sin(Om));
    row2[1] = (-sin(l)*sin(Om)+cos(l)*cos(i)*cos(Om));
    row2[2] = (cos(l)*sin(i));
    row3[0] = (sin(i)*sin(Om));
    row3[1] = (-sin(i)*cos(Om));
    row3[2] = (cos(i));

    rrorb[0] = dot(row1,rrin) ;
    rrorb[1] = dot(row2,rrin) ;
    rrorb[2] = dot(row3,rrin) ;

    return rrorb;

}

template<typename T>
DACE::AlgebraicVector<T> in2orb(const DACE::AlgebraicVector<T>& in, const DACE::AlgebraicVector<T>& r)
{
    DACE::AlgebraicVector<T> rr(3), vv(3);
    for (int i=0; i<3; i++) {
        rr[i] = r[i];
        vv[i] = r[i+3];
    }

    /* Convert vector from inertial to orbital reference frame. */
    DACE::AlgebraicVector<T> hh = rr.cross(vv);
    DACE::AlgebraicVector<T> tt = hh.cross(rr);

    const T rrnorm = rr.vnorm();
    const T hhnorm = hh.vnorm();
    const T ttnorm = tt.vnorm();

    DACE::AlgebraicVector<T> rru = rr / rrnorm;
    DACE::AlgebraicVector<T> hhu = hh / hhnorm;
    DACE::AlgebraicVector<T> ttu = tt / ttnorm;

    DACE::AlgebraicMatrix<T> rotmat(3, 3);
    rotmat.setrow(0, rru);
    rotmat.setrow(1, ttu);
    rotmat.setrow(2, hhu);

    DACE::AlgebraicVector<T> orb = rotmat*in;

    return orb;
}


#endif //ASTROCNVREF_H_INCLUDED_
