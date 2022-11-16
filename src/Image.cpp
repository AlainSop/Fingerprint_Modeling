#include "Image.h"

void Image::save(std::string name) {
    cv::Mat S = cv::Mat::zeros(cv::Size(im.cols, im.rows), cv::IMREAD_GRAYSCALE);
    S.convertTo(S, CV_64FC1);
    for (int i = 0; i < S.rows; i++){
        for (int j = 0; j < S.cols; j++) {
            S.at<double>(i, j) = get_value(i, j) * 255.0;
        }
    }
    imwrite(name, S);
}

///////////////////// GETTER /////////////////////
double Image::get_value(int x, int y) {
    return im.at<double>(x, y);
}


std::pair<double, double> Image::min_max() {
    double mini, maxi;
    cv::minMaxLoc(im, &mini, &maxi);

    return std::pair<double, double>(mini, maxi);
}
///////////////////// SETTERS /////////////////////
void Image::set_value(int x, int y, double value) {
    im.at<double>(x, y) = value;
}

void Image::set_rectangle_value(int x1, int y1, int x2, int y2, double value) {
    cv::Rect r = cv::Rect(cv::Point2d(y1, x1), cv::Point2d(y2+1, x2+1));	//defines the rectangle
    rectangle(im, r, value, cv::FILLED);									//sets all pixels in the rectangle to a specific value
}

void Image::black(int x, int y) {
    set_value(x, y, 0);
}

void Image::white(int x, int y) {
    set_value(x, y, 1);
}


///////////////////// FULL IMAGE MODIFIERS /////////////////////
void Image::negatif() {
    for (int i = 0; i < im.rows; i++) {
        for (int j = 0; j < im.cols; j++) {
            im.at<double>(i, j) = 1 - im.at<double>(i, j);
        }
    }
}


Image Image::symmetry_y()const {
    cv::Mat S = cv::Mat::zeros(im.cols, im.cols, CV_64FC1);
    for (int i = 0; i < im.cols; i++) {
        S.at<double>(im.cols - 1 - i, i) = 1;
    }
    return Image(im * S);
}


Image Image::symmetry_x_y()const {
    return Image(im.t());
}


Image Image::symmetry_x()const {
    cv::Mat S = cv::Mat::zeros(im.rows, im.rows, CV_64FC1);
    for (int i = 0; i < im.rows; i++) {
        S.at<double>(im.rows - 1 - i, i) = 1;
    }
    return Image(S * im);
}

