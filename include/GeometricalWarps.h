#ifndef GEOMETRICAL_WARPS
#define GEOMETRICAL_WARPS
#include "Image.h"

/**
 * Function that takes as parameters two coordinates (i,j) and perform a geometrical transformation on it
 * @param i : x-coordinate of the pixel
 * @param j : y-coordinate of the pixel
 * @param r : scaling factor
 * @param r_theta : scaling factor after rotation (used mostly for the inverse transform)
 * @param deg : angle of rotation in degree
 * @param spot : center of rotation and of scaling
 * @param t : translation
 * @return the transformed coordinates
 */
cv::Point2d similarityR2(int i, int j, cv::Point2d r, cv::Point2d r_theta, double deg, cv::Point2d spot, cv::Point2d t);


/**
 * Function that apply the inverse transformation parametrized by (r, r_theta, deg, spot, t)
 * @param i : x-coordinate of the pixel
 * @param j : y-coordinate of the pixel
 * @param r : scaling factor
 * @param r_theta : scaling factor after rotation (used mostly for the inverse transform)
 * @param deg : angle of rotation in degree
 * @param spot : center of rotation and of scaling
 * @param t : translation
 * @return the transformed coordinates
 */
cv::Point2d similarity_inverseR2(int i, int j, cv::Point2d r, cv::Point2d r_theta, double deg, cv::Point2d spot, cv::Point2d t);


/**
 * Interpolation of the value of the point of coordinates pt by the neareast neighbor
 * @param pt : coordinates of the point whose value we want to interpolate
 * @param im : the Image matrix
 * @return the interpolated value at pt
 */
double interpolation_neighbors(cv::Point2d pt, cv::Mat im);


/**
 * Linear interpolation of a function between two point x1 and x2 with values fx1 and fx2
 * @param x : coordinate of the point whose value we want to interpolate
 * @param x1 : a point
 * @param x2 : a point
 * @param fx1 : the value the function evaluated at point x1
 * @param fx2 : the value of the function evaluated at point x2
 * @return the interpolated value at point x
 */
double interpolation_linear(double x, double x1, double x2, double fx1, double fx2);


/**
 * Bilinear interpolation of the value of the point pt according to its 4 neighbours in the matrix im
 * @param pt : coordinated of the point whose value we want to interpolate 
 * @param im : the Image matrix
 * @return the interpolated value at point pt
 */
double interpolation_bilinear(cv::Point2d pt, cv::Mat im);


/**
 * Hermite basis for polynomials of degree 3
 * @param t : point where we evaluate the function
 * @return the value of h00 at point t
 */
double h00(double t);
/**
 * Hermite basis for polynomials of degree 3
 * @param t : point where we evaluate the function
 * @return the value of h10 at point t
 */
double h10(double t);
/**
 * Hermite basis for polynomials of degree 3
 * @param t : point where we evaluate the function
 * @return the value of h01 at point t
 */
double h01(double t);
/**
 * Hermite basis for polynomials of degree 3
 * @param t : point where we evaluate the function
 * @return the value of h11 at point t
 */
double h11(double t);


/**
 * Approximation of the first derivative with unit step
 * @param fx0 : value of the function at point x0
 * @param fx2 : value of the function at point x1
 * @return the approximation of the first derivative of the function at point x
 */
double first_derivative(double fx0, double fx2);


/**
 * Cubic interpolation (in 1D)
 * @param x : the coordinate of the point whose value we want to interpolate
 * @param x0 : a point
 * @param x1 : a point
 * @param x2 : a point
 * @param x3 : a point
 * @param fx0 : value of the function at x0
 * @param fx1 : value of the function at x1
 * @param fx2 : value of the function at x2
 * @param fx3 : value of the function at x3
 * @return the interpolated value of the function at point x
 */
double interpolation_cubic(double x, double x0, double x1, double x2, double x3, double fx0, double fx1, double fx2, double fx3);


/**
 * Bicubic interpolation (in 2D)
 * @param pt : the coordinates of the point whose value we want to interpolate
 * @param im : the Image matrix
 * @return the interpolated value of the function at point pt
 */
double interpolation_bicubic(cv::Point2d pt, cv::Mat im);



/**
 * Returns the euclidian distance between two points
 * @param center : point with respect to which we take the distance
 * @param pt : point at which we compute the distance
 * @return the distance between point center and pt
 */
double euclidian_dist(cv::Point center, cv::Point pt);


/**
 * Gaussian function with rate alpha
 * @param r : variable at which we compute the function
 * @param alpha : rate of decrease
 * @return the value of the Gaussian function at point r
 */
double gaussian(double r, double alpha);
#endif