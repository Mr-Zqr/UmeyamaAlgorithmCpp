/* A C++ library implements Umeyama's algorithm using Eigen. 
 * "Umeyama S. Least-squares estimation of transformation parameters between two point patterns[J]. 
 * IEEE Transactions on Pattern Analysis & Machine Intelligence, 1991, 13(04): 376-380."
 * Qingrui Zhao 2023.4.21
 * qingray.zhao@gmail.com
 */

#include "UmeyamaAlignment.h"
#include <ostream>

int UmeyamaAlignment::align(std::vector<Eigen::Vector3d> & point_set_x, std::vector<Eigen::Vector3d> & point_set_y, Eigen::Matrix3d & R, Eigen::Vector3d & t, double & c)
{
    if(point_set_x.size() != point_set_y.size())
    {
        std::cout << "[Umeyama Alignment] The number of points in two point sets must be equal!" << std::endl;
        throw;
    }
    else if(point_set_x.size() < 3)
    {
        std::cout << "[Umeyama Alignment] The number of points in two point sets must be greater than 3!" << std::endl;
        throw;
    }
    else 
    {
    int point_size = point_set_x.size();

    Eigen::Vector3d miu_x = calc_miu(point_set_x, point_size);
    Eigen::Vector3d miu_y = calc_miu(point_set_y, point_size);

    Eigen::Matrix3d H = calc_H(point_set_x, point_set_y, miu_x, miu_y, point_size);
    Eigen::Matrix3d DS_;

    double sigma_x = calc_sigma(point_set_x, miu_x, point_size);

    R = calc_R(H, DS_);
    t = miu_y - R * miu_x;  
    c = calc_c(sigma_x, DS_);
    return 1;
    }
}

Eigen::Vector3d UmeyamaAlignment::calc_miu(std::vector<Eigen::Vector3d> & origin, int point_size)
{
    Eigen::Vector3d avg = Eigen::Vector3d::Zero();
    for(int i = 0; i < point_size; i++)
    {
        avg += origin[i];
    }
    avg /= point_size;
    return avg;
}

Eigen::Matrix3d UmeyamaAlignment::calc_H(std::vector<Eigen::Vector3d> & x, std::vector<Eigen::Vector3d> & y, Eigen::Vector3d miu_x, Eigen::Vector3d miu_y, int point_size)
{
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(3, 3);
    for(int i = 0; i < point_size; i++)
    {
        H = H + (y[i] - miu_y) * (x[i] - miu_x).transpose();
    }
    H = H / point_size;
    return H;
}

Eigen::Matrix3d UmeyamaAlignment::calc_R(Eigen::Matrix3d & origin, Eigen::Matrix3d & DS, const float er)
{
// 进行svd分解
   Eigen::JacobiSVD<Eigen::MatrixXd> svd_holder(origin,
                                                Eigen::ComputeThinU |
                                                Eigen::ComputeThinV);

    // 构建SVD分解结果
    Eigen::MatrixXd U = svd_holder.matrixU();
    Eigen::MatrixXd V = svd_holder.matrixV();
    Eigen::MatrixXd D = svd_holder.singularValues();

    // 构建D矩阵
    Eigen::MatrixXd Diag(V.cols(), U.cols());
    Diag.setZero();

    Eigen::MatrixXd S = Eigen::MatrixXd::Identity(3, 3);
    for (unsigned int i = 0; i < D.size(); ++i) 
    {
        if (D(i, 0) > er) {
            Diag(i, i) = D(i, 0);
        } else {
            Diag(i, i) = 0;
        }
    }
    if(origin.determinant() < 0)
    {
        S(2, 2) = -S(2, 2);
    }
    DS = Diag*S;
    return U * V.transpose();
}

double UmeyamaAlignment::calc_sigma(std::vector<Eigen::Vector3d> & x, Eigen::Vector3d miu_x, int point_size)
{
    double sigma = 0;
    for(int i = 0; i < point_size; i++)
    {
        sigma = sigma + (x[i] - miu_x).norm() * (x[i] - miu_x).norm();
    }
    sigma = sigma / point_size;
    return sigma;
}

double UmeyamaAlignment::calc_c(double sigma_x, Eigen::Matrix3d & DS)
{
    double c = 0;
    c = DS.trace() / sigma_x;
    return c;
}