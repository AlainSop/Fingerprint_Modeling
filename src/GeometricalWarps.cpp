#include "GeometricalWarps.h"


cv::Point2d
similarityR2(int i, int j, cv::Point2d r, cv::Point2d r_theta, double deg, cv::Point2d spot, cv::Point2d t) {
    double theta = -deg * 3.141592 / 180; // from degree to radian
    cv::Point2d P;               // the point returned in R�

    //Apply the transformation as seen in the report
    P.x = r_theta.x * (cos(theta) * r.x * (i - spot.x) + sin(theta) * r.y * (j - spot.y)) + spot.x + t.x;
    P.y = r_theta.y * (-sin(theta) * r.x * (i - spot.x) + cos(theta) * r.y * (j - spot.y)) + spot.y + t.y;

    return P;
}


cv::Point2d
similarity_inverseR2(int i, int j, cv::Point2d r, cv::Point2d r_theta, double deg, cv::Point2d spot, cv::Point2d t) {
    return similarityR2(i, j, cv::Point2d(1 / r_theta.x, 1 / r_theta.y), cv::Point2d(1 / r.x, 1 / r.y), -deg, t + spot,
                        -t);
}


double interpolation_neighbors(cv::Point2d pt, cv::Mat im) {
    cv::Point ref; // the nearest neighbour
    if (pt.x - floor(pt.x) <= ceil(pt.x) - pt.x) {
        ref.x = floor(pt.x);
    } else {
        ref.x = ceil(pt.x);
    }
    if (pt.y - floor(pt.y) <= ceil(pt.y) - pt.y) {
        ref.y = floor(pt.y);
    } else {
        ref.y = ceil(pt.y);
    }
    return im.at<double>(ref.x, ref.y); // approximation of the value of the nearest neighbour
}


double interpolation_linear(double x, double x1, double x2, double fx1, double fx2) {
    if (x1 == x2) { // if we already know the value of x
        return fx1;
    } else {
        return fx1 + (fx2 - fx1) / (x2 - x1) * (x - x1);
    }
}


double interpolation_bilinear(cv::Point2d pt, cv::Mat im) {
    cv::Point2d Q11(floor(pt.y), floor(pt.x)), Q12(floor(pt.y), ceil(pt.x)), Q21(ceil(pt.y), floor(pt.x)), Q22(
            ceil(pt.y), ceil(pt.x)); // the neighbours
    // we perform two linear interpolations w.r.t. the x coordinate
    double fp1 = interpolation_linear(pt.x, Q11.y, Q12.y, im.at<double>(Q11), im.at<double>(Q12));
    double fp2 = interpolation_linear(pt.x, Q21.y, Q22.y, im.at<double>(Q21), im.at<double>(Q22));
    // we finaly do a linear interpolation w.r.t. the y coordinate to find an approximation of f(pt)
    return interpolation_linear(pt.y, Q11.x, Q22.x, fp1, fp2);
}


double h00(double t) {
    return (1 + 2 * t) * pow((1 - t), 2);
}

double h10(double t) {
    return t * pow(1 - t, 2);
}

double h01(double t) {
    return pow(t, 2) * (3 - 2 * t);
}

double h11(double t) {
    return pow(t, 2) * (t - 1);
}


double first_derivative(double fx0, double fx2) {
    return (fx2 - fx0) / 2;
}


double interpolation_cubic(double x, double x0, double x1, double x2, double x3, double fx0, double fx1, double fx2,
                           double fx3) {
    double t = x - x1;
    return h00(t) * fx1 + h10(t) * first_derivative(fx0, fx2) + h01(t) * fx2 + h11(t) * first_derivative(fx1, fx3);
}


double interpolation_bicubic(cv::Point2d pt, cv::Mat im) {
    cv::Point2d Q11(floor(pt.y), floor(pt.x)), Q21(ceil(pt.y), floor(pt.x)), Q12(floor(pt.y), ceil(pt.x)), Q22(
            ceil(pt.y), ceil(pt.x));
    cv::Point2d Q00(Q11.x - 1, Q11.y - 1), Q10(Q11.x, Q11.y - 1), Q20(Q11.x + 1, Q11.y - 1), Q30(Q11.x + 2, Q11.y - 1);
    cv::Point2d Q01(Q11.x - 1, Q11.y), Q02(Q11.x - 1, Q11.y + 1), Q03(Q11.x - 1, Q11.y + 2);
    cv::Point2d Q13(Q11.x, Q11.y + 2), Q23(Q11.x + 1, Q11.y + 2), Q33(Q11.x + 2, Q11.y + 2);
    cv::Point2d Q31(Q11.x + 2, Q11.y), Q32(Q11.x + 2, Q11.y + 1);
    double f1 = interpolation_cubic(pt.x, Q00.y, Q01.y, Q02.y, Q03.y, im.at<double>(Q00), im.at<double>(Q01),
                                    im.at<double>(Q02), im.at<double>(Q03));
    double f2 = interpolation_cubic(pt.x, Q10.y, Q11.y, Q12.y, Q13.y, im.at<double>(Q10), im.at<double>(Q11),
                                    im.at<double>(Q12), im.at<double>(Q13));
    double f3 = interpolation_cubic(pt.x, Q20.y, Q21.y, Q22.y, Q23.y, im.at<double>(Q20), im.at<double>(Q21),
                                    im.at<double>(Q22), im.at<double>(Q23));
    double f4 = interpolation_cubic(pt.x, Q30.y, Q31.y, Q32.y, Q33.y, im.at<double>(Q30), im.at<double>(Q31),
                                    im.at<double>(Q32), im.at<double>(Q33));
    return interpolation_cubic(pt.y, Q00.x, Q10.x, Q20.x, Q30.x, f1, f2, f3, f4);
}


void Image::transform(std::string interpolation, cv::Point2d spot, cv::Point2d r, cv::Point2d r_theta, double deg,
                      cv::Point2d t) {
    cv::Mat im_temp = cv::Mat::zeros(im.rows, im.cols, CV_64FC1); // matrix associated to the result Image
    cv::Point2d sim_temp;
    for (int i = 0; i < im.rows; i++) { // we iterate on every pixel coordinate of the result Image
        for (int j = 0; j < im.cols; j++) {
            sim_temp = similarity_inverseR2(i, j, r, r_theta, deg, spot,
                                            t); // we find the coordinates sim_temps associated to the transformed pixel (i, j)
            if ((sim_temp.x <= im.rows - 1) and (0 <= sim_temp.x) and (sim_temp.y <= im.cols - 1) and
                (0 <= sim_temp.y)) {
                if (interpolation == "neighbors" or interpolation == "N") { // nearest-neighbours interpolation
                    im_temp.at<double>(i, j) = interpolation_neighbors(sim_temp, im);
                } else if (interpolation == "bilinear" or interpolation == "B") { // bilinear interpolation
                    im_temp.at<double>(i, j) = interpolation_bilinear(sim_temp, im);
                } else if (interpolation == "bicubic" or interpolation == "BC") { // bicubic interpolation
                    if (sim_temp.x > im.rows - 2 or sim_temp.x < 1 or sim_temp.y > im.cols - 2 or sim_temp.y < 1) {
                        im_temp.at<double>(i, j) = interpolation_bilinear(sim_temp, im);
                    } else {
                        im_temp.at<double>(i, j) = interpolation_bicubic(sim_temp, im);
                    }
                } else { // otherwise we only take the upper left neighbour
                    im_temp.at<double>(i, j) = im.at<double>(floor(sim_temp.y), floor(sim_temp.x));
                }
            } else {
                im_temp.at<double>(i, j) = 1;
            }
        }
    }
    im = im_temp;
}


void
Image::inverse_transform(std::string interpolation, cv::Point2d spot, cv::Point2d r, cv::Point2d r_theta, double deg,
                         cv::Point2d t) {
    transform(interpolation, spot + t, cv::Point2d(1 / r_theta.x, 1 / r_theta.y), cv::Point2d(1 / r.x, 1 / r.y), -deg,
              -t);
}


double euclidian_dist(cv::Point center, cv::Point pt) {
    return sqrt(pow(center.x - pt.x, 2) + pow(center.y - pt.y, 2));
}


double gaussian(double r, double alpha) {
    return exp(-alpha * r * r);
}


void Image::warp_transform(std::string interpolation, cv::Point2d spot) {
    cv::Mat im_temp = cv::Mat::zeros(im.rows, im.cols, CV_64FC1); // matrix associated to the result Image
    cv::Point2d sim_temp;
    for (int i = 0; i < im.rows; i++) { // we iterate on every pixel coordinate of the result Image
        for (int j = 0; j < im.cols; j++) {
            double theta = -15 * gaussian(euclidian_dist(spot, cv::Point2d(i, j)), 0.001);
            double rx = 1.2 - 0.4 * gaussian(euclidian_dist(spot, cv::Point2d(i, j)), 0.0005);
            double ry = 1.1 - 0.4 * gaussian(euclidian_dist(spot, cv::Point2d(i, j)), 0.0001);
            double tx = 10 - 30 * gaussian(abs(i - spot.x), 0.0001);
            double ty = 0;
            sim_temp = similarity_inverseR2(i, j, cv::Point2d(rx, ry), cv::Point2d(1, 1), theta, spot, cv::Point2d(ty,
                                                                                                                   tx)); // we find the coordinates sim_temps associated to the transformed pixel (i, j)
            if ((sim_temp.x <= im.rows - 1) and (0 <= sim_temp.x) and (sim_temp.y <= im.cols - 1) and
                (0 <= sim_temp.y)) {
                if (interpolation == "neighbors" or interpolation == "N") { // nearest-neighbours interpolation
                    im_temp.at<double>(i, j) = interpolation_neighbors(sim_temp, im);
                } else if (interpolation == "bilinear" or interpolation == "B") { // bilinear interpolation
                    im_temp.at<double>(i, j) = interpolation_bilinear(sim_temp, im);
                } else if (interpolation == "bicubic" or interpolation == "BC") { // bicubic interpolation
                    if (sim_temp.x > im.rows - 2 or sim_temp.x < 1 or sim_temp.y > im.cols - 2 or sim_temp.y < 1) {
                        im_temp.at<double>(i, j) = interpolation_bilinear(sim_temp, im);
                    } else {
                        im_temp.at<double>(i, j) = interpolation_bicubic(sim_temp, im);
                    }
                } else { // otherwise we only take the upper left neighbour
                    im_temp.at<double>(i, j) = im.at<double>(floor(sim_temp.y), floor(sim_temp.x));
                }
            } else {
                im_temp.at<double>(i, j) = 1;
            }
        }
    }
    im = im_temp;
}

