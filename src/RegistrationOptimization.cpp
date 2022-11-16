#include "RegistrationOptimization.h"

double Image::average_value() {
    double S = 0;
    for (int i = 0; i < im.rows; i++) {
        for (int j = 0; j < im.cols; j++) {
            S += im.at<double>(i, j);
        }
    }
    return S / (im.rows * im.cols);
}

double loss_function_translation(cv::Point2d t, Image image0, Image image1) {
    Image im_temp(cv::Mat::zeros(image1.im.rows, image1.im.cols, CV_64FC1));
    im_temp.im = image1.im.clone();
    im_temp.transform("B", im_temp.get_center_pressure(), cv::Point2d(1, 1), cv::Point2d(1, 1), 0, t);
    double S = 0;
    for (int i = 0; i < image0.im.rows; i++) {
        for (int j = 0; j < image0.im.cols; j++) {
            S += pow(image0.im.at<double>(i, j) - im_temp.im.at<double>(i, j), 2);
        }
    }
    return S;

}


double loss_function(cv::Point2d t, Image image0, Image image1, double f_barre) {
    Image im_temp(cv::Mat::zeros(image1.im.rows, image1.im.cols, CV_64FC1));
    im_temp.im = image1.im.clone();
    im_temp.transform("B", im_temp.get_center_pressure(), cv::Point2d(1, 1), cv::Point2d(1, 1), 0, t);
    double S1 = 0;
    double S2 = 0;
    double S3 = 0;
    double g_barre = im_temp.average_value();
    for (int i = 0; i < image0.im.rows; i++) {
        for (int j = 0; j < image0.im.cols; j++) {
            S1 += (image0.im.at<double>(i, j) - f_barre) * (im_temp.im.at<double>(i, j) - g_barre);
            S2 += pow(image0.im.at<double>(i, j) - f_barre, 2);
            S3 += pow(im_temp.im.at<double>(i, j) - g_barre, 2);
        }
    }
    if ((S2 == 0) || (S3 == 0)) {
        return 0;
    }
    else {
        return abs(S1 / (sqrt(S2) * sqrt(S3)));
    }
}


cv::Point2d best_translation(Image image0, Image image1, int a) {
    std::map<double, cv::Point2d> Tab; // map of all the translation parameters ordered by their loss function value
    cv::Point2d translation_temp = cv::Point2d(-image0.im.rows, -image0.im.cols);
    double f_barre = image0.average_value();
    for (int i = 0; i < 2 * image0.im.rows - 1; i++) {
        for (int j = 0; j < 2 * image0.im.cols - 1; j++) {
            if ((i + j) % 5 == 0) { // decrease of the complexity by selecting certains translations
                translation_temp.x += i;
                translation_temp.y += j;
                if (a == 0) {
                    Tab[loss_function_translation(translation_temp, image0, image1)] = translation_temp;
                }
                else if (a == 1) {
                    Tab[loss_function(translation_temp, image0, image1, f_barre)] = translation_temp;
                }
                translation_temp = cv::Point2d(-image0.im.rows, -image0.im.cols);
            }
        }
    }
    if (a == 0) {
        return Tab.begin()->second; // return the translation parameter of the lowest value of the first loss function
    }
    else {
        return Tab.end()->second; // return the translation parameter of the highest value of the second loss function
    }
}

Image absolute_error(Image image0, Image image1) {
    Image im_temp(cv::Mat::zeros(image0.im.rows, image0.im.cols, CV_64FC1));
    im_temp.im = image0.im.clone();
    for (int i = 0; i < image0.im.rows; i++) {
        for (int j = 0; j < image0.im.cols; j++) {
            im_temp.im.at<double>(i, j) = abs(image0.im.at<double>(i, j) - image1.im.at<double>(i, j));
        }
    }
    return im_temp;
}

