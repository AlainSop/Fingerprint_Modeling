#ifndef RESTORATION
#define RESTORATION

#include <cstring>
#include <opencv2/opencv.hpp>
#include "Simulation.h"
#include <cmath>

/**
 * Create a patch around the point of coordinates (x,y) with a given size from a given Image
 * @param x Row coordinate
 * @param y Column coordinate
 * @param size Size of the patch
 * @param in_image Image where the patch is cropped from
 * @return A newly created patch
 */
Image createPatch(int x, int y, int size, const Image &in_image);


/**
 * Squares the Image with a given border size on the longest side
 * @param in_image The Image to square
 * @param border The size of the border
 * @param color 0 for a constant white border, anything else for an expanding
 * @return The squared Image
 */
Image squareImage(const Image &in_image, int border, double color);


/**
 * Add border to the entry Image with a given border size
 * @param in_image The Image to border
 * @param border The size of the border
 * @param border_type 0 for a constant white border, anything else for an expanding
 * @return The bordered Image
 */
Image borderImage(const Image &in_image, int border, int border_type);


/**
 * Computes the distance (intensity difference) between two patches, takes into account only the pixel outside of the mask
 * @param patch1 First patch
 * @param patch2 Second patch
 * @param mask Local mask
 * @return The distance between the two patches
 */
double euclideanDistance(const Image &patch1, const Image &patch2, const Image& mask);



/**
 * This class describes a dictionary of patches used to perform the restoration of an Image
 * @param patches std::set of patches of whatever size ordered by intensity
 */
class DicPatches{
public:

    std::set<Image> patches;
    /**
     * Creates and fill the dictionary from an Image given a number of patches and the size of the patches.
     * @param in_image The Image whose patch are cropped from
     * @param nb_patches The number of patches
     * @param patch_size The patch size
     */
    DicPatches(const Image &in_image, long unsigned int nb_patches, int patch_size) {
        srand(0);
        int rows = in_image.im.rows;
        int cols = in_image.im.cols;
        while(patches.size()!=nb_patches) {
            int y = (rand() % (cols - patch_size)) + patch_size/2;
            int x = (rand() % (rows - patch_size)) + patch_size/2;
            Image patch = createPatch(x,y,patch_size,in_image);
            patches.insert(patch);
        }
    }

    /**
     * A method that returns a matrix of the patch in the dic
     */
    cv::Mat printDic();
};


/**
 * Restore the value of a pixel
 * @param dic The patch dictionary
 * @param patch The patch around our pixel
 * @param mask The local mask around our pixel
 * @return The value of the restored pixel
 */
double restorePixelValue(const DicPatches &dic, const Image& patch, const Image& mask);

/**
 * Restore an Image following a priority algorithm
 * @param in_image The Image to be restored
 * @param mask The mask associated to the missing pixel
 * @param dic A dictionary of patches
 * @return Restored Image
 */
Image restore2(Image &in_image, Image mask, DicPatches &dic);

/**
 * Computes the Image gradient around a point for patches of size 3x3
 * @param patch A patch around the point
 * @return The horizontal, vertical gradient and its orientation
 */
std::pair<std::pair<double,double>,double> computegradient(Image &patch);

/**
 * Restore an Image using the chosen route
 * @param in_image The Image to restore
 * @param mask The mask applied to the Image
 * @param dic The patches dictionary
 * @param parcours The type of route chosen, 1 for borders, 2 for external spiral, 3 internal spiral, 4 linear route, default is 4
 * @return The restored Image
 */
Image restoreImage(Image &in_image, Image mask, DicPatches &dic, int parcours);

#endif