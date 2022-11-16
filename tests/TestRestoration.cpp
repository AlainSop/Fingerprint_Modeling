#include "Image.h"
#include "Simulation.h"
#include "Restoration.h"
#include "LinearFiltering.h"
#include "ShowManyImages.h"
#include "SaveManyImages.h"
int main() {

    Image img("data_in/weak_finger.png");
    DicPatches dic(img, 3025, 9);
    cv::Mat matdic = dic.printDic();
    Image mask(img.ring());

    Image im_test = restore2(img,mask,dic);
    std::cout << "Fig 1 : Original weak finger" << std::endl;
    std::cout << "Fig 2 : Mask applied on the image" << std::endl;
    std::cout << "Fig 3 : Restored image" << std::endl;
    SaveManyImages("data_out/TestRestoration.png",3,img.im,mask.im,im_test.im);
    return 0;
}


