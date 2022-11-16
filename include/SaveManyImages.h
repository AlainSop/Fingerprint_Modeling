#include <opencv2/opencv.hpp>

#include <stdio.h>
#include <stdarg.h>

#include "Image.h"

/**
 * Save multiple images (opencv matrices) in a single Image, up to 12 images
 * @param title The title of the Image
 * @param nArgs The number of images to be recorded (if you add more than nArgs Image only nArgs images are still recorded)
 * @param ... The images (cv::Mat)
 */

void SaveManyImages(std::string title, int nArgs, ...) {
    int size;
    int i;
    int m, n;
    int x, y;

    // w - Maximum number of images in a row
    // h - Maximum number of images in a column
    int w, h;

    // scale - How much we have to resize the Image
    float scale;
    int max;

    // If the number of arguments is lesser than 0 or greater than 12
    // return without displaying
    if (nArgs <= 0) {
        printf("Number of arguments too small....\n");
        return;
    }
    else if (nArgs > 14) {
        printf("Number of arguments too large, can only handle maximally 12 images at a time ...\n");
        return;
    }
    // Determine the size of the Image,
    // and the number of rows/cols
    // from number of arguments
    else if (nArgs == 1) {
        w = h = 1;
        size = 300;
    }
    else if (nArgs == 2) {
        w = 2; h = 1;
        size = 300;
    }
    else if (nArgs == 3 || nArgs == 4) {
        w = 2; h = 2;
        size = 300;
    }
    else if (nArgs == 5 || nArgs == 6) {
        w = 3; h = 2;
        size = 200;
    }
    else if (nArgs == 7 || nArgs == 8) {
        w = 4; h = 2;
        size = 200;
    }
    else {
        w = 4; h = 3;
        size = 150;
    }

    // Create a new 3 channel Image
    cv::Mat DispImage = cv::Mat::zeros(cv::Size(100 + size * w, 60 + size * h), CV_64FC1);

    // Used to get the arguments passed
    va_list args;
    va_start(args, nArgs);

    // Loop for nArgs number of arguments
    for (i = 0, m = 20, n = 20; i < nArgs; i++, m += (20 + size)) {
        // Get the Pointer to the IplImage
        cv::Mat img = va_arg(args, cv::Mat);

        // Check whether it is NULL or not
        // If it is NULL, release the Image, and return
        if (img.empty()) {
            printf("Invalid arguments");
            return;
        }

        // Find the width and height of the Image
        x = img.cols;
        y = img.rows;

        // Find whether height or width is greater in order to resize the Image
        max = (x > y) ? x : y;

        // Find the scaling factor to resize the Image
        scale = (float)((float)max / size);

        // Used to Align the images
        if (i % w == 0 && m != 20) {
            m = 20;
            n += 20 + size;
        }

        // Set the Image ROI to display the current Image
        // Resize the input Image and copy the it to the Single Big Image
        cv::Rect ROI(m, n, (int)(x / scale), (int)(y / scale));
        cv::Mat temp; resize(img, temp, cv::Size(ROI.width, ROI.height));
        temp.copyTo(DispImage(ROI));
    }

    // Create a new window, and show the Single Big Image
    cv::namedWindow(title, 1);
    Image DispIm("data_in/clean_finger.png");
    DispIm.im = DispImage;
    DispIm.save(title);

    // End the number of arguments
    va_end(args);
}
