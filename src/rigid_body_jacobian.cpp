#include <rigid_body_jacobian.h>
#include <iostream>

void rigid_body_jacobian(Eigen::Matrix36d &J, 
                         Eigen::Ref<const Eigen::Matrix3d> R, Eigen::Ref<const Eigen::Vector3d> p, 
                         Eigen::Ref<const Eigen::Vector3d> X) {

    Eigen::Matrix3d X_skew;
   // X_skew << 0.0, -(X(2)-p(2)), X(1)-p(1),
   //           X(2)-p(2), 0.0, -(X(0)-p(0)),
   //           -(X(1)-p(1)), X(0)-p(0), 0.0;

    X_skew << 0.0, -X(2), X(1),
            X(2), 0.0, -X(0),
            -X(1), X(0), 0.0;

    J.setZero();
    J.block(0, 0, 3, 3) = R * X_skew.transpose() * R.transpose();
    J.block(0, 3, 3, 3) = Eigen::Matrix3d::Identity();
    //std::cout << "J\n" << J << std::endl;
   
}

