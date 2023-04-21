/* A C++ library implements Umeyama's algorithm using Eigen. 
 * "Umeyama S. Least-squares estimation of transformation parameters between two point patterns[J]. 
 * IEEE Transactions on Pattern Analysis & Machine Intelligence, 1991, 13(04): 376-380."
 * Qingrui Zhao 2023.4.21
 * qingray.zhao@gmail.com
 */

#include "UmeyamaAlignment.h"

void UmeyamaAlignment::align(std::vector<Eigen::Vector3d> & point_set_x, std::vector<Eigen::Vector3d> & point_set_y, Eigen::Matrix3d & R, Eigen::Vector3d & t, double & c)
{
    if(point_set_x.size() != point_set_y.size())
    {
        std::cout << "The number of points in two point sets are not equal!" << std::endl;
        return;
    }

    int point_size = point_set_x.size();

    Eigen::Vector3d miu_x = calc_miu(point_set_x, point_size);
    Eigen::Vector3d miu_y = calc_miu(point_set_y, point_size);

    Eigen::Matrix3d H = calc_H(point_set_x, point_set_y, miu_x, miu_y, point_size);
    Eigen::Matrix3d DS_;

    double delta_x = calc_delta(point_set_x, miu_x, point_size);

    R = calc_R(H, DS_);
    t = miu_y - R * miu_x;  
    c = calc_c(delta_x, DS_);
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
        H = H + (x[i] - miu_x) * (y[i] - miu_y).transpose();
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

    // 构建S矩阵
    Eigen::MatrixXd S(V.cols(), U.cols());
    S.setZero();

    for (unsigned int i = 0; i < D.size(); ++i) {

        if (D(i, 0) > er) {
            S(i, i) = D(i, 0);
        } else {
            S(i, i) = 0;
        }
    }
    return V * U.transpose();
}

double UmeyamaAlignment::calc_delta(std::vector<Eigen::Vector3d> & x, Eigen::Vector3d miu_x, int point_size)
{
    double delta = 0;
    for(int i = 0; i < point_size; i++)
    {
        delta += (x[i] - miu_x).norm();
    }
    delta = delta / point_size;
    return delta;
}

double UmeyamaAlignment::calc_c(double delta_x, Eigen::Matrix3d & DS)
{
    double c = 0;
    c = DS.trace() / delta_x;
    return c;
}