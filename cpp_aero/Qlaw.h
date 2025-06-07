#ifndef Qlaw_H
#define Qlaw_H

template<typename T,typename U> T QLawQuotient( DACE::AlgebraicVector<U> Weights, DACE::AlgebraicVector<T> xTarget,  DACE::AlgebraicVector<T> Xcurrent, U rpmin, U k , U f, U mu)
{
        T sma      = Xcurrent[0];
        T e        = Xcurrent[1];
        T ape      = Xcurrent[4];
    
        T sma_t 	 = xTarget[0];


        U artificialScaling = sqrt(3); 
        U normweights = sqrt(Weights[1]*Weights[1] + Weights[2]*Weights[2] + Weights[3]*Weights[3]); 
    
        U Wp       = Weights[0];
        U Wa       = artificialScaling * Weights[1]/normweights;
        U We       = artificialScaling * Weights[2]/normweights;
        U Winc     = artificialScaling * Weights[3]/normweights;

        // Relative distance
        DACE::AlgebraicVector<T> doe(3);
        
        for (unsigned int i = 0; i < 3; i++)
        {
            doe[i] = Xcurrent[i] - xTarget[i];
        }
      

        // Penalty
        T P        = exp( k * ( 1 - sma * (1 - e) / rpmin ) );

        // put h,p,r,b in terms of orbital elements
        T p        = sma*(1-e*e);
        T h        = sqrt(mu*p);

        // Denominators
        T da_xx    = 2 * f * sqrt(sma*sma*sma * (1 + e) / ( mu * (1 - e)));
        T de_xx    = 2 * f * p / h;
        T di_xx    = f*p / ( h*(sqrt( 1 - e*e*sin(ape)*sin(ape) ) - e*abs(cos(ape))));

        T s_sma = pow((1 + pow(( (sma - sma_t) / (3 * sma_t)),4) ),0.5);

        T Q = (1+ Wp*P)*(Wa*s_sma*(doe[0]/da_xx)*(doe[0]/da_xx)
                       + We*(doe[1]/de_xx)*(doe[1]/de_xx) + Winc*(doe[2]/di_xx)*(doe[2]/di_xx)); 

        return Q;

}



#endif
