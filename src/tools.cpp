#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

//const float EPS = 0.0001;
//const float EPS2 = 0.0000001;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
    VectorXd rmse(4);
	rmse << 0,0,0,0;
	if(estimations.size() == 0
	  || estimations.size() != ground_truth.size()){
	    cout << "error input" << endl;
	    return rmse;
	}
    VectorXd squared_residuals = rmse;
    VectorXd c = rmse;

	//accumulate squared residuals
    unsigned int i;
	for(i=0; i < estimations.size(); ++i){
        // ... your code here
        c = estimations[i] - ground_truth[i];
        c = c.array()*c.array();
		squared_residuals = squared_residuals + c;
	}
    rmse = squared_residuals / estimations.size();
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
	MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//TODO: YOUR CODE HERE

	//check division by zero
  //  if (fabs(px) < EPSILON and fabs(py) < EPSILON){
    //    px = EPSILON;
     //   py = EPSILON;
   // }

	//compute the Jacobian matrix
	float square_xy = px*px+py*py;
	float sqrt_xy = sqrt(square_xy);
	float sqrt32_xy = square_xy * sqrt_xy;

	Hj << px/sqrt_xy, py/sqrt_xy, 0, 0,
	      -py/square_xy, px/square_xy, 0, 0,
	      py*(vx*py-vy*px)/sqrt32_xy, px*(vy*px-vx*py)/sqrt32_xy, px/sqrt_xy,py/sqrt_xy;

	return Hj;
}
