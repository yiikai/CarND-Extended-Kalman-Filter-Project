#include "tools.h"
#include "kalman_filter.h"
#include <iostream>
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.
const float PI2 = 2 * M_PI;
KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_*x_;
  P_ = F_*P_*F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
	cout<<"Update"<<endl;
	Eigen::VectorXd y = z - H_*x_;
	Eigen::MatrixXd S = H_*P_*H_.transpose() + R_;
	Eigen::MatrixXd K = P_*H_.transpose()*S.inverse();
	x_ = x_ + K * y;
	Eigen::MatrixXd I = Eigen::MatrixXd::Identity(x_.size(),x_.size());
	P_ = (I - K * H_) * P_;	 
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
	float px = x_[0];
	float py = x_[1];
	float vx = x_[2];
	float vy = x_[3];
	//calc pred_z
	float rho = sqrt(px*px + py*py);
	float phi = atan2(py,px);
  if(rho < EPSILON)
    rho = EPSILON;
	float rho_dot = (px*vx + py*vy) / sqrt(px*px + py*py);
	Eigen::VectorXd z_pred(3);
	z_pred << rho, phi, rho_dot;
	
	Eigen::VectorXd y = z - z_pred;
  if(y(1) > M_PI)
  {
    y(1) -= PI2;
  }

  if(y(1) < -M_PI)
  {
    y(1) += PI2;
  }

	Eigen::MatrixXd S = H_*P_*H_.transpose() + R_;
	Eigen::MatrixXd K = P_*H_.transpose()*S.inverse();
	x_ = x_ + K * y;
	Eigen::MatrixXd I = Eigen::MatrixXd::Identity(x_.size(),x_.size());
	P_ = (I - K * H_) * P_;	 

}
