#include "UmeyamaAlignment.h"
#include <Eigen/Dense>
#include <iostream>

int main()
{
    std::vector<Eigen::Vector3d> point_set_x;
    std::vector<Eigen::Vector3d> point_set_y;
    Eigen::Matrix3d R;
    Eigen::Vector3d t;
    double c;
    point_set_x.push_back(Eigen::Vector3d(0.721046620579811, 0.522495305777102, 0.993704624120852));
    point_set_x.push_back(Eigen::Vector3d(0.218676632399634, 0.105798273250228, 0.109697464523194));
    point_set_x.push_back(Eigen::Vector3d(0.063591370975106, 0.404579995857626, 0.448372912066495));
    point_set_x.push_back(Eigen::Vector3d(0.365816176838171, 0.763504640848813, 0.627896379614169));

    point_set_y.push_back(Eigen::Vector3d(-0.518632480327867, 0.977166142386389, 1.308468839104955));
    point_set_y.push_back(Eigen::Vector3d(0.523348616013710, 0.838081633002886, 0.988450754807947));
    point_set_y.push_back(Eigen::Vector3d(0.338759281547603, 1.021038526625798, 1.389046021281762));
    point_set_y.push_back(Eigen::Vector3d(-0.065716737541079, 0.811449645673728, 1.600864482292104));

    UmeyamaAlignment::align(point_set_x, point_set_y, R, t, c);

    std::cout << "R: " << std::endl << R << std::endl;
    std::cout << "t: " << std::endl << t << std::endl;
    std::cout << "c: " << c << std::endl;
}   