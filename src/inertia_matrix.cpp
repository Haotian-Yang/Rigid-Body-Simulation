#include <inertia_matrix.h>
#include <cassert>
#include <iostream>

//compute inertia matrix and volume by integrating on surfaces
void inertia_matrix(Eigen::Matrix3d &I, Eigen::Vector3d & center, double &mass, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, double density) {


    double one, X, Y, Z, XX, YY, ZZ, XY, ZX, YZ = 0;
    for(int i=0; i<F.rows(); i++)
    {
        Eigen::Vector3d X0, X1, X2, N;
        X0 = V.row(F(i, 0));
        X1 = V.row(F(i, 1));
        X2 = V.row(F(i, 2));

        // The normal of the current triangle
        N = (X1 - X0).cross(X2 - X0);

        one += 1.0/6.0*N.x() * (X0.x() + X1.x() + X2.x());
        X +=  1.0/12.0*N.x() * (X0.x()*X0.x() + X1.x()*X1.x() + X2.x()*X2.x() + X1.x()*X2.x() + X0.x()*X1.x() + X0.x()*X2.x());
        Y +=  1.0/12.0*N.y() * (X0.y()*X0.y() + X1.y()*X1.y() + X2.y()*X2.y() + X1.y()*X2.y() + X0.y()*X1.y() + X0.y()*X2.y());
        Z +=  1.0/12.0*N.z() * (X0.z()*X0.z() + X1.z()*X1.z() + X2.z()*X2.z() + X1.z()*X2.z() + X0.z()*X1.z() + X0.z()*X2.z());
        XX += 1.0/20.0*N.x() * (X0.x()*X0.x()*X0.x() + X0.x()*X0.x()*(X1.x() + X2.x()) + (X1.x() + X2.x())*(X1.x()*X1.x() + X2.x()*X2.x()) + X0.x()*(X1.x()*X1.x() + X1.x()*X2.x() +X2.x()*X2.x()));
        YY += 1.0/20.0*N.y() * (X0.y()*X0.y()*X0.y() + X0.y()*X0.y()*(X1.y() + X2.y()) + (X1.y() + X2.y())*(X1.y()*X1.y() + X2.y()*X2.y()) + X0.y()*(X1.y()*X1.y() + X1.y()*X2.y() +X2.y()*X2.y()));
        ZZ += 1.0/20.0*N.z() * (X0.z()*X0.z()*X0.z() + X0.z()*X0.z()*(X1.z() + X2.z()) + (X1.z() + X2.z())*(X1.z()*X1.z() + X2.z()*X2.z()) + X0.z()*(X1.z()*X1.z() + X1.z()*X2.z() +X2.z()*X2.z()));
        XY += 1.0/60.0*N.x() *(X0.y()*(X1.x()*X1.x() + X1.x()*X2.x() + X2.x()*X2.x()) + X1.y()*(3*X1.x()*X1.x() + 2*X1.x()*X2.x() + X2.x()*X2.x()) + (X1.x()*X1.x() + 2*X1.x()*X2.x() + 3*X2.x()*X2.x())*X2.y() + X0.x()*X0.x()*(3*X0.y() + X1.y() + X2.y()) + X0.x()*(2*X0.y()*(X1.x() + X2.x()) + X1.x()*(2*X1.y() + X2.y())+ X2.x()*(X1.y() + 2*X2.y())));
        YZ += 1.0/60.0*N.y() *(X0.z()*(X1.y()*X1.y() + X1.y()*X2.y() + X2.y()*X2.y()) + X1.z()*(3*X1.y()*X1.y() + 2*X1.y()*X2.y() + X2.y()*X2.y()) + (X1.y()*X1.y() + 2*X1.y()*X2.y() + 3*X2.y()*X2.y())*X2.z() + X0.y()*X0.y()*(3*X0.z() + X1.z() + X2.z()) + X0.y()*(2*X0.z()*(X1.y() + X2.y()) + X1.y()*(2*X1.z() + X2.z())+ X2.y()*(X1.z() + 2*X2.z())));
        ZX += 1.0/60.0*N.z()*(X0.z()*X0.z()*(X1.x() + X2.x()) +X1.z()*X1.z()*(3*X1.x()+X2.x()) + 2*X1.z()*(X1.x()+X2.x())*X2.z() + (X1.x()+3*X2.x())* X2.z()*X2.z() + X0.x()*(3*X0.z()*X0.z() + X1.z()*X1.z() + X1.z()*X2.z() + X2.z()*X2.z() + 2*X0.z()*(X1.z()+X2.z())) + X0.z()*(X1.x()*(2*X1.z()+X2.z()) + X2.z()*(X1.z() + 2*X2.z())) );
    }

    mass = one;
    center << X/mass, Y/mass, Z/mass;
    I.setZero();

    I(0,0) = YY + ZZ;
    I(0,1) = -XY;
    I(0,2) = -ZX;
    I(1,0) = I(0,1);
    I(1,1) = XX + ZZ;
    I(1,2) = -YZ;
    I(2,0) = I(0,2);
    I(2,1) = I(1,2);
    I(2,2) = XX + YY;
    I(0,0) -= mass * (center.y() * center.y() + center.z() * center.z() );
    I(0,1) += mass * center.x() * center.y();
    I(0,2) += mass * center.z() * center.x();
    I(1,0) = I(0,1);
    I(1,1) -= mass * ( center.z() * center.z() + center.x() * center.x() );
    I(1,2) += mass * center.y() * center.z();
    I(2,0) = I(0,2);
    I(2,1) = I(1,2);
    I(2,2) -= mass * ( center.x() * center.x() + center.y() * center.y() );

}