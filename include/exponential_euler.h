#include <Eigen/Dense>
#include <EigenTypes.h>
#include <rodrigues.h>
#include <iostream>
#include <inertia_matrix.h>
#include <rigid_to_world.h>
#include <unsupported/Eigen/MatrixFunctions>


//Input:
//  q - 12n vector where n is the number of rigid bodies. Each rigid body is stored as 12 doubles. 
//      The first 9 doubles are the columns of a 3x3 rotation matrix and the final 3 doubles are the world space position of the object's center of mass.
//  qdot - 6n vector of generalied velocities. The first 3 doubles of each body are the world space angular velocity and 
//         the second 3 are the world space linear velocity.
//  dt - the integration time step
//  masses - a vector to mass matrices for each rigid body
//  forces - a 6n vector of generalized forces for n rigid bodies. The first 3 doubles of each rigid body are the torques acting on the object
//           while the second 3 doubles are the linear forces.
//Output:
//  q - updated generalized coordinates 
//  qdot - updated generalized velocities 
inline void exponential_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, 
                            std::vector<Eigen::Matrix66d> &masses, Eigen::Ref<const Eigen::VectorXd> forces) {

    Eigen::Matrix3d Rt, RIR;
    Rt.setZero();
    RIR.setZero();
    Rt = Eigen::Map<const Eigen::Matrix3d>(q.segment<9>(0).data());


    RIR = Rt * masses[0].block(0, 0, 3, 3) * Rt.transpose();
    double m = masses[0].block(3, 3, 3, 3)(0, 0);

    Eigen::Vector3d pdot, pt, omega, torques_ext, f_ext;
    pt = q.segment(9, 3);
    omega = qdot.segment(0, 3);
    pdot = qdot.segment(3, 3);



    torques_ext = forces.segment(0, 3);
    f_ext = forces.segment(3, 3);

    //Velocity update
    Eigen::Vector3d omega_new, pdot_new;
    omega_new << omega - (RIR.inverse()*(dt*omega.cross(RIR * omega) + dt*torques_ext));
    pdot_new << pdot + dt*f_ext/m;
    qdot.segment(0, 3) = omega_new;
    qdot.segment(3, 3) = pdot_new;

    //Position update;
    Eigen::Matrix3d omega_skew;
    omega_skew << 0.0, -omega(2), omega(1),
            omega(2), 0.0, -omega(0),
            -omega(1), omega(0), 0.0;

    Eigen::Matrix3d R_new;
    rodrigues(R_new, omega*dt);
    R_new = R_new * Rt;

    Eigen::Vector3d p_new = pt + dt*pdot;
    q.segment(0, 3) = R_new.block(0, 0, 3, 1);
    q.segment(3, 3) = R_new.block(0, 1, 3, 1);
    q.segment(6, 3) = R_new.block(0, 2, 3, 1);
    q.segment(9, 3) = p_new;




}