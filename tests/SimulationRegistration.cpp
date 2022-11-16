#include "ShowManyImages.h"
#include "RegistrationOptimization.h"
#include "SaveManyImages.h"

int main() {

    //opening the image
    Image Im_1("data_in/clean_finger.png");
    Image Im_2("data_in/txy_finger.png");
    Image Im_3("data_in/txy_finger.png");
    Image Im_4("data_in/txy_finger.png");

    Im_3.transform("B", Im_3.get_center_pressure(), cv::Point2d(1,1), cv::Point2d(1,1), 0, cv::Point2d(-21,17));

    Im_4 = absolute_error(Im_1,Im_3);


    std::cout << "Fig 1 : Clean finger" << std::endl;
    std::cout << "Fig 2 : Translated finger" << std::endl;
    std::cout << "Fig 3 : Warped translated finger" << std::endl;
    std::cout << "Fig 4 : Absolute error" << std::endl;

    SaveManyImages("data_out/SimulationRegistration.png", 4, Im_1.im, Im_2.im, Im_3.im, Im_4.im);

    return 0;
}