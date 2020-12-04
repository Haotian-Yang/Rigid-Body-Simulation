#include <rodrigues.h>
#include <cmath>

void rodrigues(Eigen::Matrix3d &R, Eigen::Ref<const Eigen::Vector3d> omega) {
    R.setZero();

    double theta = omega.norm();
    Eigen::Vector3d omega_norm = omega/omega.norm();


    Eigen::Matrix3d omega_skew;
    omega_skew << 0.0, -omega_norm(2), omega_norm(1),
                  omega_norm(2), 0.0, -omega_norm(0),
                  -omega_norm(1), omega_norm(0), 0.0;


    // exp([w]theta)
    R = Eigen::Matrix3d::Identity() +
        sin(theta)*omega_skew +
        (1-cos(theta))*omega_skew*omega_skew;


}