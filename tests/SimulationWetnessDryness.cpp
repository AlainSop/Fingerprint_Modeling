#include "Image.h"
#include "Simulation.h"
#include "Restoration.h"
#include "LinearFiltering.h"
#include "ShowManyImages.h"
#include "SaveManyImages.h"
int main() {

    //opening images
    Image DRY("data_in/dry_finger.png");
    Image MOIST("data_in/moist_finger.png");
    Image Im_1("data_in/clean_finger.png");
    Image Im_2("data_in/clean_finger.png");

    Im_1.binarize();
    Im_2.binarize();

    Im_1.im = Im_1.bin_erosion("4");
    Im_2.im = Im_2.bin_dilatation("12");


    Image Im_3("data_in/clean_finger.png");
    Image Im_4("data_in/clean_finger.png");

    Im_3.im = Im_3.gs_erosion("8");
    Im_4.im = Im_4.gs_dilatation("12");

    Image Im_5("data_in/clean_finger.png");
    Image Im_6("data_in/clean_finger.png");
    double k = Im_6.find_threshold_parameter();

    cv::Point2d center = Im_5.get_center_pressure();

    Im_5.im = Im_5.gs_faded_erosion("8", center, 100, 0.01, 1.4);

    Im_6.im = Im_6.gs_faded_dilatation("4", center, 60, 0.06, 1.3);
    Im_6.im = Im_6.gs_faded_dilatation("4", center, 55, 0.15, 0.15);
    Im_6.im = Im_6.gs_faded_dilatation("2w", center, 35, 0.1, 0.2);


    cv::Mat net_gauss = cv::Mat::zeros(5, 5, CV_64FC1);
    std::vector<double> values{ 1, 4, 6, 4, 1, 4, 16, 24, 16, 4, 6, 24, -476, 24, 6, 4, 16, 24, 16, 4, 1, 4, 6, 4, 1 };
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            net_gauss.at<double>(i, j) = -1.0 * values.back() / 256.0;
            values.pop_back();
        }
    }
    Image nett(net_gauss.clone());

    Im_6 = convolutionDft(Im_6, nett);
    Im_6.contrast(0.1, k, k, 0.9);
    Im_6.make_dry();

    std::cout << "Fig 1 : Original moist finger" << std::endl;
    std::cout << "Fig 2 : Original dry finger" << std::endl;
    std::cout << "Fig 3 : Eroded binarized fingerprint" << std::endl;
    std::cout << "Fig 4 : Dilated binarized fingerprint" << std::endl;
    std::cout << "Fig 5 : Eroded fingerprint" << std::endl;
    std::cout << "Fig 6 : Dilated fingerprint" << std::endl;
    std::cout << "Fig 7 : Simulation of a moist finger" << std::endl;
    std::cout << "Fig 8 : Simulation of a dry finger" << std::endl;

    SaveManyImages("data_out/SimulationWetnessDryness.png", 8, MOIST.im, DRY.im, Im_1.im, Im_2.im, Im_3.im, Im_4.im, Im_5.im, Im_6.im);



    return 0;
}