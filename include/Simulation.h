#ifndef SIMULATION
#define SIMULATION

#define M_PI           3.14159265358979323846  /* pi */
#include <random>
#include "Image.h"


/**
 * Computes the angle between the point pt and the y axis passing through center returns theta between 0 and 2pi
 * @param pt : point whose angle we need to compute
 * @param center : center of rotation
 * @return the angle between the point pt and the y axis passing through center
 */
double angle(cv::Point2d pt, cv::Point2d center);


/**
 * Returns the linear interpolation of the function at point x where we know the function at x0 and x1
 * @param x : coordinate of the point whose value we want to interpolate
 * @param x0 : a point
 * @param x1 : a point
 * @param fx0 : value of the function at point x0
 * @param fx1 : value of the function at point x1
 * @return the interpolated value of the function at point x
 */
double linear(double x, double x0, double x1, double fx0, double fx1);
/**
 * Generates a random variable (theta_k, A_theta_k, B_theta_k) where (A_theta_k, B_theta_k) are the scaling of the elipse at angle theta_k
 * @param theta_k : previous angle generated
 * @param m : mean number of fluctuation points
 * @param eps1 : minimum bound of x-variations
 * @param eps2 : maximum bound of x-variations
 * @param eps3 : minimum bound of y_variations
 * @param eps4 : maximum bound of y_variations
 * @return a random variable (theta_k, A_theta_k, B_theta_k)
 */
cv::Point3d generate_random_fluctuation(double theta_k, double m, double eps1, double eps2, double eps3, double eps4);

/**
 * Generates all random variables (theta_k, A_theta_k, B_theta_k) and put them into a vector of 3-dimensional points
 * @param m : mean number of fluctuation points
 * @param eps1 : minimum bound of x-variations
 * @param eps2 : maximum bound of x-variations
 * @param eps3 : minimum bound of y_variations
 * @param eps4 : maximum bound of y_variations
 * @return a vector of 3-dimensional points of the form (theta_k, A_theta_k, B_theta_k)
 */
std::vector<cv::Point3d> generate_random_boundary(double m, double eps1, double eps2, double eps3, double eps4);
/**
 * Returns the value of the scaling factors A and B at angle theta based on the random fluctuations we simulated
 * @param theta : angle between 0 and 2pi
 * @param fluctuations : vector of 3-dimensional points of the form (theta_k, A_theta_k, B_theta_k)
 * @return the value of A and B at angle theta
 */

cv::Point2d simulate_A_B(double theta, std::vector<cv::Point3d> fluctuations);

///////////////////// DISTANCES /////////////////////

/**
 * Returns the euclidian distance between two points
 * @param center : first point (center of pressure or center of coordinates)
 * @param pt : second point
 * @return the distance between the two points
 */
double euclidian_dist(cv::Point2d center, cv::Point2d pt);

/**
 * Returns the distance between two points, with different metrics for x-axis and y-axis
 * @param center : first point (center of pressure or center of coordinates)
 * @param pt : second point
 * @param scale_y : unit of the y axis metric
 * @return the distance between the two points with a different metric
 */
double distorted_distance(cv::Point2d center, cv::Point2d pt, double scale_y);



/**
 * Returns the distance between two points, with scaled metric for x_axis and y_axis
 * @param center : center of the Image
 * @param pt : point at which we want to compute the distance
 * @param scale_x : x scaling factor of the ellipse
 * @param scale_y : y scaling factor of the ellipse
 * @return the distance between center and pt scaled by scale_x and scale_y
 */
double elliptic_distance(cv::Point2d center, cv::Point2d pt, double scale_x, double scale_y);

///////////////////// ISOTROPIC FUNCTIONS /////////////////////
/**
 * Encodes the isotropic Gaussian coefficient function
 * @param center : center of the Image
 * @param pt : point at which we want to compute the function
 * @param alpha : coefficient
 * @param beta : unused parameter (needed to create a general function pointing)
 * @return the value of the function evaluated at point pt
 */
double isotropic_gaussian(cv::Point2d center, cv::Point2d pt, double alpha, double beta = 0);
/**
 * Encodes the isotropic inverted polynomial coefficient function
 * @param center : center of the Image
 * @param pt : point at which we want to compute the function
 * @param alpha : degree of the polynomial
 * @param beta : unused parameter (needed to create a general function pointing)
 * @return the value of the function evaluated at point pt
 */
double isotropic_inverted_poly(cv::Point2d center, cv::Point2d pt, double alpha, double beta = 0);
/**
 * Encodes the isotropic logistic coefficient function
 * @param center : center of the Image
 * @param pt : point at which we want to compute the difunctionstance
 * @param alpha : rate of decrease
 * @param beta : shift
 * @return the value of the function evaluated at point pt
 */
double isotropic_logistic(cv::Point2d center, cv::Point2d pt, double alpha, double beta);

///////////////////// ANISOTROPIC FUNCTIONS /////////////////////
/**
 * Encodes the anisotropic logistic coefficient function, with the y-axis scale being 1.6 times the x-axis one (best parameter for the simulation  of a weak pressure)
 * @param center : center of the Image
 * @param pt : point at which we want to compute the function
 * @param alpha : rate of decrease
 * @param beta : shift
 * @return the value of the function evaluated at point pt
 */
double anisotropic_logistic(cv::Point2d center, cv::Point2d pt, double alpha, double beta);
/**
 * Encodes the an isotropic logistic coefficient function, with the y-axis scale as a parameter
 * @param center : center of the Image
 * @param pt : point at which we want to compute the function
 * @param alpha : rate of decrease
 * @param beta : shift
 * @param scale : metric in y-axis = scale * metric in x-axis
 * @return the value of the function evaluated at point pt
 */
double anisotropic_logistic_v2(cv::Point2d center, cv::Point2d pt, double alpha, double beta, double scale);

///////////////////// ANISOTROPIC RANDOM FUNCTIONS /////////////////////
/**
 * Computes the value of the logistic function with random fluctuations on the boundary
 * @param center : center of pressure
 * @param pt : point at which we want to compute the function
 * @param fluctuations : vector of 3-dimensional points of the form (theta_k, A_theta_k, B_theta_k)
 * @param alpha : rate of decrease
 * @param beta : shift
 * @return the value of the function at point pt
 */
double random_anisotropic_logistic(cv::Point2d center, cv::Point2d pt, std::vector<cv::Point3d> fluctuations, double alpha, double beta);

/**
 * Computes the new coordinates of a pixel after its rotation at a certain angle
 * @param center : center of pressure
 * @param pt : point to rotate
 * @param deg : angle of rotation
 * @return the rotated point
 */
cv::Point2d rotation(cv::Point2d center, cv::Point2d pt, double deg);

#endif
