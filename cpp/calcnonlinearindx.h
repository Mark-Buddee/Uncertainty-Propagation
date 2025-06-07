#ifndef calcnonlinearindx_H
#define calcnonlinearindx_H


template <typename T, typename U>
U calcnonlinearindexLoads(DACE::AlgebraicVector<T> fx,  U nSigma)
{
    // nSigma = number of standard deviations. 

    U xval;
    DACE::AlgebraicMatrix<T> STM(6, 6);
    DACE::AlgebraicMatrix<U> Jbar(6, 6), Bmat(6,6);
    U b;

    // calculate the STM and jbar norm
    U normJbar = 0;

    for (unsigned int i = 0; i < 6; i++)
    {
        for (unsigned int j = 1; j < 7; j++)
        {
            STM.at(i, j - 1) = fx[i].deriv(j);
            Jbar.at(i,j-1) = cons(STM.at(i, j - 1));
            normJbar += pow(Jbar.at(i, j - 1),2);
        }
    }

    // analyse the dependency on each variable individually
    for (unsigned int i = 0; i < 6; i++)
    {
        for (unsigned int j = 1; j < 7; j++)
        {

            for (int l = 0; l < 6; l++)
            {
                Bmat.at(i,j-1) += DACE::cons(pow((STM.at(i, j - 1)).deriv(l + 1),2)); 

            }

        }
    }

    b = 0;
    for (unsigned int i = 0; i < 6; i++)
    {
        for (unsigned int j = 0; j < 6; j++)
        {
            b += Bmat.at(i,j);          
        }
    }

    xval = nSigma*sqrt(b)/sqrt(normJbar);

    return xval;
}


#endif