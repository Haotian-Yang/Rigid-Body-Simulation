#include <Eigen/Dense>
#include <EigenTypes.h>

//  R - rotation matrix 
//  omega - angular velocity vector. Note: This is wrong. We need to the pass omega x dt into this function
void rodrigues(Eigen::Matrix3d &R, Eigen::Ref<const Eigen::Vector3d> omega);