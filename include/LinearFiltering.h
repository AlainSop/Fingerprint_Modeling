#ifndef LINEAR_FILTERING
#define LINEAR_FILTERING
#include <cstring>
#include <opencv2/opencv.hpp>

#include "Restoration.h"

/**
 * @brief Function using the classical definition of convolution to compute the filtering of an Image
 * @param img The Image to be filtered
 * @param kernel The filter
 * @return The convolved Image
 */
Image convolutionMatrix(const  Image &img,const Image &kernel);

/**
 * @brief Return the filtered Image using a kernel changing following the distance to the center of pressure
 * @param img The Image to be filtered
 * @param kernel The filter
 * @return The filtered Image
 */
Image convolutionMatrixAni(const  Image &img, const Image &kernel);

/**
 * Returns the dft non printable (two dimensions Real part and imaginary part) of an Image
 * @param img Image whose dft is returned
 * @return The dft as an Image
 */
Image dft(const Image & img);

/**
 * Compute the normalized_dft which can be printed to see effect of filters
 * @param I Image
 * @return Normalized dft of I (printable)
 */
Image normalized_dft(const Image & I);

/**
 * Compute the inverse dft of an Image
 * @param img Image whose idft is performed
 * @return The idft of an Image
 */
Image idft(const Image & img);

/**
 * Filtering of an Image using the dft and idft defined earlier
 * @param img The Image to filter
 * @param kernel The filter
 * @return The filtered Image
 */
Image convolutionDft(const Image & img, Image & kernel);

#endif