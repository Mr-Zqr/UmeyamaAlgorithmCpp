# include <iostream>
# include <Eigen/Dense>
# include <Eigen/SVD>
# include <Eigen/Core>
# include <vector>

class UmeyamaAlignment
{
    private:
        static Eigen::Matrix3d calc_R(Eigen::Matrix3d & origin, Eigen::Matrix3d & DS, const float er=0);
        static Eigen::Vector3d calc_miu(std::vector<Eigen::Vector3d> & origin, int point_size);
        static Eigen::Matrix3d calc_H(std::vector<Eigen::Vector3d> & x, std::vector<Eigen::Vector3d> & y, Eigen::Vector3d miu_x, Eigen::Vector3d miu_y, int point_size);
        static double calc_c(double delta_x, Eigen::Matrix3d & DS);
        static double calc_delta(std::vector<Eigen::Vector3d> & x, Eigen::Vector3d miu_x, int point_size);
    public:
    UmeyamaAlignment();
    ~UmeyamaAlignment();
    static void align(std::vector<Eigen::Vector3d> & point_set_x, std::vector<Eigen::Vector3d> & point_set_y, Eigen::Matrix3d & R, Eigen::Vector3d & t, double & c);
};
