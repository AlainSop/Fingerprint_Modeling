#ifndef REGISTRATION_OPTIMIZATION
#define REGISTRATION_OPTIMIZATION
#include "Simulation.h"
#include "GeometricalWarps.h"
#include <fstream>
#include <iostream> 




/**
 * Function that renders the first loss function value for a given translation t and two images
 * @param t : vector of translation
 * @param image0 : first Image to compare
 * @param image1 : second Image to compare
 * @return the value of the loss function
 */
double loss_function_translation(cv::Point2d t, Image image0, Image image1);


/**
 * Function that renders the second loss function value for a given translation t, two images and the average value of the first Image
 * @param t : vector of translation
 * @param image0 : first Image to compare
 * @param image1 : second Image to compare
 * @param f_barre : average value of the first Image
 * @return the value of the loss function
 */
double loss_function(cv::Point2d t, Image image0, Image image1, double f_barre);

//
/**
 * Warning : this function has got a very long execution time (between 10 and 20 minutes). 
 *  Function with a greedy algorithm that renders the best parameters of translation with the first or second loss function
 * @param image0 : first Image to compare
 * @param image1 : second Image to compare
 * @param a : if a = 0 it uses the first lost function, otherwise it uses the second one
 * @return the best translation vector to compare the two images
 */
cv::Point2d best_translation(Image image0, Image image1, int a = 0);


/**
 * Function that renders the absolute error Image given two images
 * @param image0 : first Image to compare
 * @param image1 : second Image to compare
 * @return the Image of the absolute error
 */
Image absolute_error(Image image0, Image image1);

#endif
