#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
   */
  H_laser_ << 1, 0, 0, 0,
             0,1,0,0;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    std::cout<<"initialized....."<<std::endl;
    MeasurementPackage::SensorType type = measurement_pack.sensor_type_;
    if( MeasurementPackage::LASER == type)
    {
      //Lidar no need non-linear change
      float px = measurement_pack.raw_measurements_[0];
      float py = measurement_pack.raw_measurements_[1];
      ekf_.x_ = Eigen::VectorXd(4);
      ekf_.x_ << px,py,0,0;
    }
    else
    {
      //radar need non-linear change
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      //float rho_dot = measurement_pack.raw_measurements_[2];
      float px = rho*cos(phi);
      float py = rho*sin(phi);
      ekf_.x_ = Eigen::VectorXd(4);
      ekf_.x_ << px,py,0,0;
    }
    if(fabs(ekf_.x_(0)) < EPSILON and fabs(ekf_.x_(1)) < EPSILON )
    {
      ekf_.x_(0) = EPSILON;
      ekf_.x_(1) = EPSILON;
    }
    ekf_.P_ =  Eigen::MatrixXd(4,4);
    ekf_.P_ << 1,0,0,0,
               0,1,0,0,
               0,0,1000,0,
               0,0,0,1000;  
	  previous_timestamp_ = measurement_pack.timestamp_;	
    // done initializing, no need to predict or update
    is_initialized_ = true;
    std::cout<<"initialized..... over"<<std::endl;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  std::cout<<"prediction....."<<std::endl;
  float dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  ekf_.F_ = Eigen::MatrixXd(4,4);
  ekf_.F_ << 1,0,dt,0,
             0,1,0,dt,
             0,0,1,0,
             0,0,0,1;
  //define process noise covariance matrix Q.
  int noise_ax = 9;
  int noise_ay = 9;
  int dt2 = dt*dt;
  int dt3 = dt2*dt;
  int dt4 = dt3*dt;
  ekf_.Q_ = Eigen::MatrixXd(4,4);
  ekf_.Q_ << dt4/4*noise_ax, 0 , dt3/2*noise_ax, 0,
             0,dt4/4*noise_ay, 0, dt3/2*noise_ay,
             dt3/2*noise_ax, 0 , dt2*noise_ax, 0,
             0, dt3/2*noise_ay, 0, dt2*noise_ay;


  ekf_.Predict();

  std::cout<<"prediction.....over"<<std::endl;
  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
  std::cout<<"update....."<<std::endl;
  if(MeasurementPackage::RADAR == measurement_pack.sensor_type_)
  {
    //jacobin matrix for non-linear calc
    Hj_ =  tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  }
  else
  {
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  std::cout<<"update.....over"<<std::endl;
  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
