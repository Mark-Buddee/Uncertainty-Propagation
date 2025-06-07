#ifndef uncertainityConvert_H
#define uncertainityConvert_H


template<typename T, typename U> DACE::AlgebraicMatrix<U>  convertCovariance( DACE::AlgebraicMatrix<U> Cov0 ,DACE::AlgebraicVector<T> x0t)
{

    DACE::AlgebraicMatrix<U> J(6, 6), covariance(6,6);

    J = x0t.linear();

    // get the covariance of the transformation.
    covariance = J*Cov0*J.transpose();

    return covariance;

}



template<typename T, typename U> DACE::AlgebraicVector<T> axisTransform(DACE::AlgebraicMatrix<U> CovarianceMatrix, DACE::AlgebraicVector<T> x0t)
{

    Eigen::MatrixXd Mat(6,6);
    DACE::AlgebraicVector<T> EigenV(6,0.0);
    DACE::AlgebraicVector<T> x_Eigen_DA(6);
    DACE::AlgebraicVector<T> xsc(6);



    for (unsigned int i=0; i<6; i++) {
        for (unsigned int j=0; j<6; j++) {
            Mat(i,j)= CovarianceMatrix.at(i,j);
        }
    }

 

    Eigen::EigenSolver<Eigen::MatrixXd> es(Mat, true);
    Eigen::VectorXcd EigVal = es.eigenvalues();
    Eigen::MatrixXcd EigVec = es.eigenvectors();
    //cout << EigVal << endl;


    for (int i=0; i<6;i++) {
        EigenV[i] = sqrt(EigVal[i].real())*DACE::DA(i+1);
    }

    for (int i=0; i<6; i++) {
        for (int j=0; j<6; j++) {
            x_Eigen_DA[i] = x_Eigen_DA[i] + EigVec(i,j).real()*EigenV[j];
        }
    }


    xsc =  x0t.cons() + x_Eigen_DA;

    return xsc;

}





template<typename T, typename U> DACE::AlgebraicVector<T> convertU(DACE::AlgebraicMatrix<U> Cov0, DACE::AlgebraicVector<T> x0t)
{

    DACE::AlgebraicMatrix<U> J(6, 6), covariance(6,6);
    Eigen::MatrixXd Mat(6,6);
    DACE::AlgebraicVector<T> EigenV(6,0.0);
    DACE::AlgebraicVector<T> x_Eigen_DA(6);
    DACE::AlgebraicVector<T> xsc(6);

    J = x0t.linear();
    // get the covariance of the transformation.
    covariance = J*Cov0*J.transpose();

    for (unsigned int i=0; i<6; i++) {
        for (unsigned int j=0; j<6; j++) {
            Mat(i,j)= covariance.at(i,j);
        }
    }

    Eigen::EigenSolver<Eigen::MatrixXd> es(Mat, true);
    Eigen::VectorXcd EigVal = es.eigenvalues();
    Eigen::MatrixXcd EigVec = es.eigenvectors();
    //cout << EigVal << endl;

    for (int i=0; i<6;i++) {
        EigenV[i] = 3*sqrt(EigVal[i].real())*DACE::DA(i+1);
    }

    for (int i=0; i<6; i++) {
        for (int j=0; j<6; j++) {
            x_Eigen_DA[i] = x_Eigen_DA[i] + EigVec(i,j).real()*EigenV[j];
        }
    }

    xsc =  x0t.cons() + x_Eigen_DA;

    return xsc; 

}


#endif