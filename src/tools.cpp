#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
    rmse << 0,0,0,0;
    
    //  the estimation vector size should not be zero
    if (estimations.size()==0){
        cerr<<"Empty estimations vector on CalculateRMSE";
        return rmse;
    }
    //  the estimation vector size should equal ground truth vector size
    if (estimations.size()!=ground_truth.size()){
        cerr<<"Size of estimation and ground truth vectors should be equal";
        return rmse;
    }
    VectorXd dif(4);
    //accumulate squared residuals
    for(int i=0; i < estimations.size(); ++i){
        dif=estimations[i]-ground_truth[i];
        //rmse=rmse.array()+pow(dif.array(),2);
        rmse=rmse.array()+ dif.array() * dif.array();
    }
    
    //calculate the mean
    rmse=rmse/estimations.size();
    
    //calculate the squared root
    rmse=rmse.array().sqrt();
    
    //return the result
    return rmse;
}

