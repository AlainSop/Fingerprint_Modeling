#include "Image.h"
#include "Simulation.h"
#include "Restoration.h"
#include "LinearFiltering.h"
#include "ShowManyImages.h"
#include "SaveManyImages.h"
int main() {

    //opening images
    Image CLEAN("data_in/clean_finger.png");
    Image WEAK("data_in/weak_finger.png");
    Image Im_1("data_in/clean_finger.png");
    Image Im_2("data_in/clean_finger.png");
    Image Im_3("data_in/clean_finger.png");

    cv::Point2d center = Im_1.get_center_coordinates();
    cv::Point2d pt1(0, 0), pt2(Im_1.im.cols - 1, Im_1.im.rows - 1);

    //creating pointers to functions
    Fonction iso_gauss = isotropic_gaussian;
    Fonction aniso_log = anisotropic_logistic;

    Im_1.weak_finger(center, pt1, pt2, iso_gauss, 0.0005);

    Im_2.weak_finger(center, pt1, pt2, aniso_log, 105, 0.5, 5);

    double m1 = 10, m2 = 100, eps1 = 0.2, eps2 = 0.2, eps3 = 0.1, eps4 = 0.1;
    Im_3.weak_finger_deformations(center + cv::Point2d(20, 0), pt1, pt2, 105, 0.5, 5, m1, eps1, eps2, eps1, eps2);
    Im_3.weak_finger_deformations(center + cv::Point2d(20, 0), pt1, pt2, 105, 0.5, 5, m2, eps3, eps4, eps3, eps4);

    std::cout << "Fig 1 : Original weak pressure fingerprint" << std::endl;
    std::cout << "Fig 2 : Simulation of a weak pressure fingerprint with isotropic Gaussian coefficient function" << std::endl;
    std::cout << "Fig 3 : Simulation of a weak pressure fingerprint with anisotropic logistic coefficient function" << std::endl;
    std::cout << "Fig 4 : Simulation of a weak pressure fingerprint with anisotropic logistic coefficient function and random deformations at boundaries" << std::endl;

    SaveManyImages("data_out/SimulationPressure.png", 4, WEAK.im, Im_1.im, Im_2.im, Im_3.im);


    return 0;
}