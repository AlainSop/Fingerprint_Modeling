#include "Image.h"
#include "GeometricalWarps.h"
#include "Simulation.h"
#include "ShowManyImages.h"
#include "SaveManyImages.h"
int main() {

    //opening images
    Image CLEAN("data_in/clean_finger.png");
    Image Im_1_N("data_in/clean_finger.png"); // geometrical transform with nearest-neighbours interpolation
    Image Im_2_B("data_in/clean_finger.png"); // geometrical transform with bilinear interpolation
    Image Im_3_BC("data_in/clean_finger.png"); // geometrical transform with bicubic interpolation
    Image Im_4_N_stability("data_in/clean_finger.png"); // stability of nearest-neighbours after 100 rotations
    Image Im_5_B_stability("data_in/clean_finger.png"); // stability of bilinear after 100 rotations
    Image Im_6_BC_stability("data_in/clean_finger.png"); // stability of bicubic after 100 rotations
    Image Im_7_elasticity("data_in/clean_finger.png"); // simulation of elasticity

    cv::Point2d center = Im_1_N.get_center_coordinates();
    
    Im_1_N.transform("N", center, cv::Point2d(0.9, 1.6), cv::Point2d(1, 1), 32, cv::Point2d(0, 0));
    Im_2_B.transform("B", center, cv::Point2d(0.9, 1.6), cv::Point2d(1, 1), 32, cv::Point2d(0, 0));
    Im_3_BC.transform("BC", center, cv::Point2d(0.9, 1.6), cv::Point2d(1, 1), 32, cv::Point2d(0, 0));

    for (int i = 0; i < 100; i++) {
        Im_4_N_stability.transform("N", center, cv::Point2d(1, 1), cv::Point2d(1, 1), 32, cv::Point2d(0, 0));
        Im_5_B_stability.transform("B", center, cv::Point2d(1, 1), cv::Point2d(1, 1), 32, cv::Point2d(0, 0));
        Im_6_BC_stability.transform("BC", center, cv::Point2d(1, 1), cv::Point2d(1, 1), 32, cv::Point2d(0, 0));
    }

    Im_7_elasticity.warp_transform("B", cv::Point2d(275, 75));

    std::cout << "Fig 1 : Original clean fingerprint" << std::endl;
    std::cout << "Fig 2 : Geometrical transform with nearest-neighbours interpolation" << std::endl;
    std::cout << "Fig 3 : Geometrical transform with bilinear interpolation" << std::endl;
    std::cout << "Fig 4 : Geometrical transform with bicubic interpolation" << std::endl;
    std::cout << "Fig 5 : Stability of nearest-neighbours interpolation after 100 rotations" << std::endl;
    std::cout << "Fig 6 : Stability of bilinear interpolation after 100 rotations" << std::endl;
    std::cout << "Fig 7 : Stability of bicubic interpolation after 100 rotations" << std::endl;
    std::cout << "Fig 8 : Simulation of elasticity of a finger" << std::endl;

    SaveManyImages("data_out/GeometricalWarps.png", 8, CLEAN.im, Im_1_N.im, Im_2_B.im, Im_3_BC.im, Im_4_N_stability.im, Im_5_B_stability.im, Im_6_BC_stability.im, Im_7_elasticity.im);


    return 0;
}
