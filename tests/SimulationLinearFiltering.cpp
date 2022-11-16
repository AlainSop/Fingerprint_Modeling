#include "LinearFiltering.h"
#include "ShowManyImages.h"
#include "SaveManyImages.h"

int main(){
    Image img("data_in/clean_finger.png");
    Image blur_filter(cv::Mat::zeros(9,9,CV_64FC1));
    blur_filter.set_rectangle_value(0,4,8,4,1);
    blur_filter.im /= cv::norm(blur_filter.im,cv::NORM_L1);
    Image convmat = convolutionMatrix(img,blur_filter);
    Image convdft = convolutionDft(img,blur_filter);
    Image convani = convolutionMatrixAni(img, blur_filter);

    std::cout << "Fig 1 : Original clean finger" << std::endl;
    std::cout << "Fig 2 : Vertical blur filter" << std::endl;
    std::cout << "Fig 3 : Blurred image using matrix scalar product" << std::endl;
    std::cout << "Fig 4 : Blurred image using DFT" << std::endl;
    std::cout << "Fig 5 : Anisotropic blurred image" << std::endl;

    SaveManyImages("data_out/SimulationLinearFiltering.png", 5, img.im, blur_filter.im,convmat.im, convdft.im,convani.im);
    return 0;
}
