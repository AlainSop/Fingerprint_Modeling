#ifndef IMAGE
#define IMAGE

#include <stdlib.h>
#include <stdio.h>
#include <utility> //std::pair
#include <opencv2/opencv.hpp>
#include <iostream>
#include <cstring>
#include <cmath>


/**
 * allows us to define a pointer to a function (for coefficient functions)
 */
typedef double(*Fonction)(cv::Point2d, cv::Point2d, double, double);


///////////////////// MAIN CLASS IMAGE /////////////////////
/**
 * This class describes an Image
 * @param im The matrix from opencv used to contain the pixels
 */
class Image{
public:
    cv::Mat im;

    //constructors
    /**
     * Constructor
     * @param path : path of the Image that we choose for the instance of the class
     */
    Image(std::string path) {
        im = cv::imread(path, cv::IMREAD_GRAYSCALE);
        im.convertTo(im, CV_64FC1);//conversion en Mat_<double>
        im *= 1.0 / 255.0;
    };
    /**
     * Copy constructor
     * @param img : Image that we want to clone
     */
    Image(const Image& img) :im(img.im.clone()) {}
    /**
     * Constructor
     * @param M : matrix that represents an Image
     */
    Image(cv::Mat M) :im(M.clone()) {}
    /**
     * saves an Image under a given name
     * @param name : name that we want for the saved Image, can be a path
     */
    void save(std::string name);


    //getter and setter
    /**
     * gets a pixel color
     * @param x : abscissa of the pixel 
     * @param y : ordinate of the pixel
     * @return the intensity of pixel (x,y)
     */
    double get_value(int x, int y);
    /**
     * sets a pixel intesity to a specific value
     * @param x : abscissa of the pixel 
     * @param y : ordinate of the pixel
     * @param value : color that we want for the pixel
     */
    void set_value(int x, int y, double value);
    /**
     *draws a filled rectangle of a specific color on the Image
     * @param x1 : abscissa of the top left pixel of the rectangle
     * @param y1 : ordinate of the top left pixel of the rectangle
     * @param x2 : abscissa of the bottom right pixel of the rectangle
     * @param y2 : ordinate of the bottom right pixel of the rectangle
     * @param value : intensity that we want for the rectangle
     */
    void set_rectangle_value(int x1, int y1, int x2, int y2, double value);
    /**
     * puts the Image in negative
     */
    void negatif();

    
    // -------------- starter 1 ---------------------
    /**
     * gives the minimum and the maximum intensity of the Image
     * @return a pair of double containing (minimum, maximum) intensities of the Image
     */
    std::pair<double, double> min_max();
    /**
     * sets a pixel to black
     * @param x : abscissa of the pixel 
     * @param y : ordinate of the pixel
     */
    void black(int x, int y);
    /**
     * sets a pixel to white
     * @param x : abscissa of the pixel 
     * @param y : ordinate of the pixel
     */
    void white(int x, int y);
    /**
     * symmetrizes the Image with respect to a vertical centered axis
     * @return the vertically symmetrized Image
     */
    Image symmetry_y() const;
    /**
     * symmetrizes the Image with respect to the diagonal from top left to bottom right
     * @return the diagonally symmetrized Image
     */
    Image symmetry_x_y() const;
    /**
     *symmetrizes the Image with respect to an horizontal centered axis
     * @return the horizontally symmetrized Image
     */
    Image symmetry_x() const;


    // -------------- Main_course 1 ---------------------
    /**
     * gets the central pixel of the Image in terms of coordinates
     * @return the center pixel if there are many possible centers, it gives the top left one)
     */
    cv::Point2d get_center_coordinates();
    /**
     * Obtains the center of pressure of a fingerprint
     * @return the approximative center of pressure of the fingerprint
     */
    cv::Point2d get_center_pressure()const;
    /**
     * gets the location of the highest pixel of the Image
     * @return the darkest pixel of the Image
     */
    cv::Point2d get_loc_highest_pressure();
    /**
     * simulates a fingerprint with weak pressure by applying a coefficient function onto the Image
     * @param center : center of the Image (can be in terms of coordinates or in terms of pressure)
     * @param pt1 : top left pixel of the rectangle where we want to apply the simulation
     * @param pt2 : bottom right pixel of the rectangle where we want to apply the simulation
     * @param f : isotropic or anisotropic coefficient function
     * @param alpha : parameter for the coefficient function
     * @param beta : parameter for the coefficient function
     * @param deg : angle of rotation for the simulation
     */
    void weak_finger(cv::Point2d center, cv::Point2d pt1, cv::Point2d pt2, Fonction f, double alpha, double beta = 0, double deg = 0);
    /**
     * simulates a fingerprint with weak pressure with random deformations at boundaries
     * @param center : center of the Image (can be in terms of coordinates or in terms of pressure)
     * @param pt1 : top left pixel of the rectangle where we want to apply the simulation
     * @param pt2 : bottom right pixel of the rectangle where we want to apply the simulation
     * @param alpha : parameter for the coefficient function
     * @param beta : parameter for the coefficient function
     * @param deg : angle of rotation for the simulation
     * @param m : mean number of fluctuation points
     * @param eps1 : minimum bound of x-variations
	 * @param eps2 : maximum bound of x-variations
	 * @param eps3 : minimum bound of y_variations
	 * @param eps4 : maximum bound of y_variations
     */
    void weak_finger_deformations(cv::Point2d center, cv::Point2d pt1, cv::Point2d pt2, double alpha, double beta, double deg, double m, double eps1, double eps2, double eps3, double eps4);
    /**
     *
     * @param i : abscissa of the pixel
     * @param j : ordinate of the pixel
     * @param k : if 0, it considers 8 neighbors and if 1, it considers 24 neighbors
     * @return
     */
    std::vector <cv::Point2d> neighbours(int i, int j, int k);
    /**
     * gets intensities of neighbors of a pixel
     * @param neighborhood : size of neighborhood (can be "2h", "2w", "4", "8", "12", "24") 
     * @param i : abscissa of the pixel 
     * @param j : ordinate of the pixel 
     * @return a vector of double containing intensities of neighbors of pixel (i,j)
     */
    std::vector <double> neighbors(std::string neighborhood, int i, int j);
    /**
     * gets intensities of the 2 vertical neighbors of a pixel
     * @param i : abscissa of the pixel 
     * @param j : ordinate of the pixel 
     * @return a vector of double containing intensities of the 2 vertical neighbors of pixel (i,j)
     */
    std::vector <double> voisins2_height(int i, int j);
    /**
     * gets intensities of the 2 horizontal neighbors of a pixel
     * @param i : abscissa of the pixel 
     * @param j : ordinate of the pixel 
     * @return a vector of double containing intensities of the 2 horizontal neighbors of pixel (i,j)
     */
    std::vector <double> voisins2_width(int i, int j);
    /**
     * gets intensities of the 4 closest neighbors of a pixel
     * @param i : abscissa of the pixel 
     * @param j : ordinate of the pixel 
     * @return a vector of double containing intensities of the 4 closest neighbors of pixel (i,j)
     */
    std::vector <double> voisins4(int i, int j);
    /**
     * gets intensities of the 8 closest neighbors of a pixel (3x3 square)
     * @param i : abscissa of the pixel 
     * @param j : ordinate of the pixel 
     * @return a vector of double containing intensities of the 8 closest neighbors of pixel (i,j)
     */
    std::vector <double> voisins8(int i, int j);
    /**
     * gets intensities of the 12 closest neighbors of a pixel (diamond-shaped)
     * @param i : abscissa of the pixel 
     * @param j : ordinate of the pixel 
     * @return a vector of double containing intensities of the 12 closest neighbors of pixel (i,j)
     */
    std::vector <double> voisins12(int i, int j);
    /**
     * gets intensities of the 24 closest neighbors of a pixel (5x5 square) 
     * @param i : abscissa of the pixel 
     * @param j : ordinate of the pixel 
     * @return a vector of double containing intensities of the 24 closest neighbors of pixel (i,j)
     */
    std::vector <double> voisins24(int i, int j);
    /**
    * Method to obtain a ring mask around the fingerprint
    * @return the ring mask
    */
    cv::Mat ring();
    /**
    * Method to obtain an approximation of the boundary around the fingerprint
    * @return the boundary
    */
    cv::Mat boundary();
    /**
     * draws a white cross on the center of pressure of a fingerprint
     */
    void draw_center_pressure();

    
    // -------------- Starter 2 ---------------------
    /**
    * Geometrical transformation of an Image using 2d interpolation
    * @param interpolation : which interpolation method we use between nearest-neighbors ("N"), bilinear ("B") or bicubic ("BC")
    * @param spot : the center of rotation and of scaling
    * @param r : anisotropic scaling factor
    * @param r_theta : scaling after rotation (used mainly of the inverse transform)
    * @param deg : angle of rotation in degree
    * @param t : translation
    */
    void transform(std::string interpolation, cv::Point2d spot, cv::Point2d r, cv::Point2d r_theta, double deg, cv::Point2d t);

    /**
    * Geometrical inverse transformation of an Image using 2d interpolation
    * @param interpolation : which interpolation method we use between nearest-neighbors ("N"), bilinear ("B") or bicubic ("BC")
    * @param spot : the center of rotation and of scaling
    * @param r : anisotropic scaling factor
    * @param r_theta : scaling after rotation (used mainly of the inverse transform)
    * @param deg : angle of rotation in degree
    * @param t : translation
    */
    void inverse_transform(std::string interpolation, cv::Point2d spot, cv::Point2d r, cv::Point2d r_theta, double deg, cv::Point2d t);

    /**
    * Perform a transformation of the Image simulating the elastic property of skin when there is a pressure point
    * @param interpolation : which interpolation method we use between nearest-neighbors ("N"), bilinear ("B") or bicubic ("BC")
    * @param spot : center of pressure
    */
    void warp_transform(std::string interpolation, cv::Point2d spot);

    // -------------- starter 4 ---------------------
    /**
     * computation of the proportions of intensities of pixels, needed for the threshold selection method
     * @return a vector containing proportions of intensities in the Image
     */
    std::vector<double> get_P();
    /**
     * computation of the squared sigmas, needed for the threshold selection method
     * @param k : temporary value of the threshold parameter
     * @param P : vector of proportions of pixels intensity
     * @return the value of sigma(k, P)^2
     */
    double sigma2_B(double k, std::vector<double> P);
    /**
     * algorithm that finds the optimal threshold parameter for binarization
     * @return the value of the threshold for the Image
     */
    double find_threshold_parameter();
    /**
     * binarizes an Image given the threshold parameter
     * @param k : threshold parameter (if not given, it is found automatically by our algorithm)
     */
    void binarize(double k = -1);
    /**
     * contrasts an Image (p < "a" -> black, "a" < p < "b" -> darker, "c" < p < "d" -> brighter, "d" <p -> white)
     * @param a : lower bound of dark pixels
     * @param b : upper bound of dark pixels
     * @param c : lower bound of bright pixels
     * @param d : upper bound of bright pixels
     */
    void contrast(double a, double b, double c, double d);
    /**
     * checks if the structural element fits
     * @param neighborhood : represents the structural element
     * @param i : abscissa of a pixel
     * @param j : ordinate of a pixel
     * @return True if the structural element fits, i.e. if all the pixels of the neighborhood have value 1
     */
    bool fits(std::string neighborhood, int i, int j);
    /**
     * checks if the structural element fits
     * @param neighborhood : represents the structural element
     * @param i : abscissa of a pixel
     * @param j : ordinate of a pixel
     * @return True if the structural element hits, i.e. if at least one pixel of the neighborhood has value 1
     */
    bool hits(std::string neighborhood, int i, int j);
    /**
     * gets the matrix of eroded Image
     * @param neighborhood : represents the structural element
     * @return the matrix corresponding to the eroded Image (binary)
     */
    cv::Mat bin_erosion(std::string neighborhood);
    /**
     * gets the matrix of dilated Image
     * @param neighborhood : represents the structural element
     * @return the matrix corresponding to the dilated Image (binary)
     */
    cv::Mat bin_dilatation(std::string neighborhood);
    /**
     * gets the matrix of opened Image
     * @param neighborhood : represents the structural element
     * @return the matrix corresponding to the opened Image (binary)
     */
    cv::Mat bin_opening(std::string neighborhood);
    /**
     * gets the matrix of closed Image
     * @param neighborhood : represents the structural element
     * @return the matrix corresponding to the closed Image (binary)
     */
    cv::Mat bin_closing(std::string neighborhood);
    /**
     * gets the matrix of eroded Image
     * @param neighborhood : represents the structural element
     * @return the matrix corresponding to the eroded Image (grayscale)
     */
    cv::Mat gs_erosion(std::string neighborhood);
    /**
     * gets the matrix of dilated Image
     * @param neighborhood : represents the structural element
     * @return the matrix corresponding to the dilated Image (grayscale)
     */
    cv::Mat gs_dilatation(std::string neighborhood);
    /**
     * gets the matrix of opened Image
     * @param neighborhood : represents the structural element
     * @return the matrix corresponding to the opened Image (grayscale)
     */
    cv::Mat gs_opening(std::string neighborhood);
    /**
     * gets the matrix of closed Image
     * @param neighborhood : represents the structural element
     * @return the matrix corresponding to the closed Image (grayscale)
     */
    cv::Mat gs_closing(std::string neighborhood);
    /**
     * performs a faded erosion on an Image (erosion that fades anisotropically)
     * @param neighborhood : represents the structural element
     * @param center : center of the Image (can be in terms of coordinates or in terms of pressure)
     * @param alpha : parameter for the coefficient function
     * @param beta : parameter for the coefficient function
     * @param scale : scale of the metric for the coefficient function
     * @return matrix of the Image after faded-erosion
     */
    cv::Mat gs_faded_erosion(std::string neighborhood, cv::Point2d center, double alpha, double beta, double scale);
    /**
     * performs a faded dilatation on an Image (dilatation that fades anisotropically)
     * @param neighborhood : represents the structural element
     * @param center : center of the Image (can be in terms of coordinates or in terms of pressure)
     * @param alpha : parameter for the coefficient function
     * @param beta : parameter for the coefficient function
     * @param scale : scale of the metric for the coefficient function
     * @return matrix of the Image after faded-dilatation
     */
    cv::Mat gs_faded_dilatation(std::string neighborhood, cv::Point2d center, double alpha, double beta, double scale);
    /**
     *
     * @param i : abscissa of the center of the black dot that we want to add
     * @param j : ordinate of the center of the black dot that we want to add
     * @param S : matrix where we save and add black dots
     * @return matrix S with added black dot
     */
    cv::Mat add_black_dot(int i, int j, cv::Mat S);
    /**
     * simulates the dryness of a finger by adding black spots randomly on black furrows
     */
    void make_dry();
    /**
     * An overloading of < using matrices intensities.
     * @param img1 The left Image
     * @param img2 The right Image
     * @return True or false
     */
    friend bool operator<(const Image& img1, const Image& img2) {
    return norm(img1.im) < norm(img2.im);
    }


    // ---------------starter 5 ----------------------
    /**
     * Gets the average value of the Image
     * @return mean intensity of all pixels
     */
    double average_value();
};

#endif