#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

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
   *  predict the state
   */
	x_ = F_ * x_;
	P_ = F_ * P_ * F_.transpose() + Q_;

 
}	

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * update the state by using Kalman Filter equations
   */
	VectorXd zPred = (H_ *  x_);
	MatrixXd y = z - zPred;
	
	UpdateShared(y);

}


void KalmanFilter::UpdateEKF(const VectorXd &z) {
	/**
	 * TODO: update the state by using Extended Kalman Filter equations
	 */

	 // recover state parameters
	float px = x_(0);
	float py = x_(1);
	float vx = x_(2);
	float vy = x_(3);
	// in case Px is 0 make it to 0.0001 to not impact the results too much
	/*if (fabs(px) < 0.0001) {
		px = 0.0001;
	}*/
	float rho = sqrt(px * px + py * py);
	// in case rho is 0 make it to 0.0001 to not impact the results too much
	if (fabs(rho) < 0.0001) {
		rho = 0.0001;
	}
	float phi = atan2(py, px);
	float rhoDot = (px * vx + py * vy) / rho;
	VectorXd zPred(3);
	zPred << rho, phi, rhoDot;
	VectorXd y = z - zPred;
	// make sure the angle is between -pi to pi
	if (y(1) > M_PI) {
		y(1) -= 2 * M_PI;
	}

	if (y(1) < -M_PI) {
		y(1) += 2 * M_PI;
	}
	UpdateShared(y);
}
void KalmanFilter::UpdateShared(const VectorXd& y) {
	/*
		 duplicated portion for KF and EKF
	*/
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd K = P_ * Ht *S.inverse();
	// new state
	x_ = x_ + (K * y);
	// in case x_ is 4X4 so check x_' size every time when update and get coresponed I
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);

	P_ = (I - K * H_) *P_;


	}
