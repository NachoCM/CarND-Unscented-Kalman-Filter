#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <numeric>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 10, 0, 0,
        0, 0, 0, 10, 0,
        0, 0, 0, 0, 10;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.25;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  

  // state dimension
  n_x_ = 5;

  // spreading parameter
  lambda_ = 3 - n_x_;

  n_aug_ = 7; 

  n_sigma_points_ = 2 * n_aug_ + 1;

  Xsig_pred_ = MatrixXd(n_x_, n_sigma_points_);

  weights_ = VectorXd(n_sigma_points_);
  weights_.fill(1/(2*(lambda_+n_aug_)));
  weights_(0)=lambda_/(lambda_+n_aug_);

}

UKF::~UKF() {}

double UKF::normalizeAngle(double x){
    x = fmod(x + M_PI,2*M_PI);
    if (x < 0)
        x += 2*M_PI;
    return x - M_PI;
}


void UKF::ShowNIS(const std::string type,const std::vector<double> &v, const double threshold){
    int over_threshold=0;
    for (auto& n : v)
      if (n>threshold) over_threshold++;
    double pct_over=100.0*over_threshold/v.size();
    double mean=std::accumulate(v.begin(), v.end(), 0.0)/v.size();
    cout<<type<<" NIS \t"<<v.back()<<" mean: "<<mean<< " over threshold: "<<pct_over<<"% \n";

}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(const MeasurementPackage meas_package) {
  
  if (!is_initialized_) {
    cout << "UKF: " << endl;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
        float range=meas_package.raw_measurements_[0];
        float bearing=meas_package.raw_measurements_[1];
        float radial_velocity=meas_package.raw_measurements_[2];
        float px=cos(bearing)*range;
        float py=-sin(bearing)*range;
        x_ << px, py, 0,0,0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
        float px=meas_package.raw_measurements_[0];
        float py=meas_package.raw_measurements_[1];
        x_ << px, py, 0, 0, 0;
    }
    

    // done initializing, no need to predict or update
    previous_timestamp_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;    //dt - expressed in seconds
  previous_timestamp_ = meas_package.timestamp_;
  Prediction(dt);
  
  /*****************************************************************************
   *  Update
   ****************************************************************************/
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    UpdateRadar(meas_package);
      
  } else {
    // LIDAR update
    UpdateLidar(meas_package);
  }

  // print the output
  //cout << "x_ = " << x_ << endl;
  //cout << "P_ = " << P_ << endl;
}


void UKF::GenerateSigmaPoints(const double delta_t) {

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_,n_sigma_points_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug.tail(2) << 0, 0;
  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_)=P_;
  P_aug.block(n_x_,n_x_,2,2)<< pow(std_a_,2), 0, 0, pow(std_yawdd_,2);
  //create square root matrix
  MatrixXd A=P_aug.llt().matrixL();
  //create augmented sigma points
  Xsig_aug.col(0)=x_aug;
  for (int i=0;i<n_aug_;i++){
      Xsig_aug.col(i+1)=x_aug+sqrt(lambda_+n_aug_)*A.col(i);
      Xsig_aug.col(n_aug_+i+1)=x_aug-sqrt(lambda_+n_aug_)*A.col(i);
  }

  VectorXd raw_pred = VectorXd(5);
  VectorXd noise_pred = VectorXd(5);
  //predict sigma points
  for(int i=0;i<n_sigma_points_;i++){
      double px=Xsig_aug(0,i);
      double py=Xsig_aug(1,i);
      double v=Xsig_aug(2,i);
      double sigma=Xsig_aug(3,i);
      double sigma_dot=Xsig_aug(4,i);
      double noise_longitudinal_acc=Xsig_aug(5,i);
      double noise_yaw_acc=Xsig_aug(6,i);
      
      if(sigma_dot==0){
        raw_pred<< v*cos(sigma)*delta_t,
                 v*sin(sigma)*delta_t,
                 0,
                 0,
                 0;  
      } else {
        raw_pred<< (sin(sigma+sigma_dot*delta_t)-sin(sigma))*v/sigma_dot,
                 (-cos(sigma+sigma_dot*delta_t)+cos(sigma))*v/sigma_dot,
                 0,
                 sigma_dot*delta_t,
                 0;
      }
      noise_pred<<pow(delta_t,2)*cos(sigma)*noise_longitudinal_acc/2,
                  pow(delta_t,2)*sin(sigma)*noise_longitudinal_acc/2,
                  delta_t*noise_longitudinal_acc,
                  pow(delta_t,2)*noise_yaw_acc/2,
                  delta_t*noise_yaw_acc;  

      Xsig_pred_.col(i)=Xsig_aug.col(i).head(5) + raw_pred + noise_pred;

  }
}



/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(const double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  GenerateSigmaPoints(delta_t);
  

  //predict state mean
  for (int i=0;i<n_x_;i++){
      x_(i)=weights_.dot(Xsig_pred_.row(i));
  }
  
  //predict state covariance matrix
  
  P_.fill(0.0);
  for (int i = 0; i < n_sigma_points_; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3)=normalizeAngle(x_diff(3));

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(const MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  //measurement dimension, lidar can measure px and py
  int n_z = 2;

  //measurement vector
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_(0),meas_package.raw_measurements_(1);
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_sigma_points_);
  Zsig = Xsig_pred_.topRows(2);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  z_pred.fill(0);
  //calculate mean predicted measurement
  for (int i=0;i<n_sigma_points_;i++){
      z_pred+=Zsig.col(i)*weights_(i);
  }
  //calculate measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  MatrixXd diff;
  S.fill(0);
  for (int i=0;i<n_sigma_points_;i++){
      diff=Zsig.col(i)-z_pred;
      S+= weights_(i)*diff*diff.transpose();
  }
  VectorXd R_vec = VectorXd(2);
  R_vec << pow(std_laspx_,2), pow(std_laspy_,2);
  S+= R_vec.asDiagonal();

  //calculate cross correlation matrix Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < n_sigma_points_; i++) {

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3)=normalizeAngle(x_diff(3));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }


  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;


  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  //calculate NIS
  float NIS = z_diff.transpose() * S.inverse() * z_diff;
  //NIS for 5 degrees of freedom should be under 11.070 95% of the time
  LIDAR_NIS.push_back(NIS);
  ShowNIS("LIDAR",LIDAR_NIS,11.070);
  
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(const MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  /*****************************************************************************
   *  Transform prediction into measurement space
   ****************************************************************************/
  //measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //measurement vector
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_(0),meas_package.raw_measurements_(1),meas_package.raw_measurements_(2);
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_sigma_points_);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  z_pred.fill(0);
  //transform sigma points into measurement space
  for (int i=0;i<n_sigma_points_;i++){
      double px=Xsig_pred_(0,i);
      double py=Xsig_pred_(1,i);
      double v=Xsig_pred_(2,i);
      double sigma=Xsig_pred_(3,i);
      double sigma_dot=Xsig_pred_(4,i);
      double root_square_sum = pow(pow(px,2)+pow(py,2), 0.5);
      if (root_square_sum<1.0e-6){
            std::cout<<"root square sum too close to 0"<<'\n';
            root_square_sum = 1.0e-6;
        }
      Zsig.col(i)<<root_square_sum,
          atan2(py,px),
          (px*cos(sigma)*v+py*sin(sigma)*v)/root_square_sum;
      //calculate mean predicted measurement
      z_pred+=Zsig.col(i)*weights_(i);
  }
  //calculate measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  MatrixXd diff;
  S.fill(0);
  for (int i=0;i<n_sigma_points_;i++){
      diff=Zsig.col(i)-z_pred;
      diff(1)=normalizeAngle(diff(1));
      S+= weights_(i)*diff*diff.transpose();
  }
  VectorXd R_vec = VectorXd(3);
  R_vec << pow(std_radr_,2), pow(std_radphi_,2), pow(std_radrd_,2);
  S+= R_vec.asDiagonal();

  //calculate cross correlation matrix Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < n_sigma_points_; i++) {

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    z_diff(1)=normalizeAngle(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3)=normalizeAngle(x_diff(3));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  z_diff(1)=normalizeAngle(z_diff(1));

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  //calculate NIS
  double NIS = z_diff.transpose() * S.inverse() * z_diff;
  //NIS for 3 degrees of freedom should be under 7.815 95% of the time.
  RADAR_NIS.push_back(NIS);
  ShowNIS("RADAR",RADAR_NIS,7.815);

}
