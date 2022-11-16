#include "Simulation.h"

//gets the proportion of every intensity in the image
std::vector<double> Image::get_P() {
    std::vector<double> P(256, 0);
    //for every pixel in the image
    for (int i = 0; i < im.rows; i++) {
        for (int j = 0; j < im.cols; j++) {
            //we add 1/(number of pixels in the image) to the counter of the intensity
            P[round(255.0*get_value(i, j))] += 1.0 / (im.rows*im.cols);
        }
    }

    return P;
}

//computation of the sigma_B squared, according to the definition of the given document for the automatic threshold selection method
double Image::sigma2_B(double k, std::vector<double> P) {

    double mu_T = 0, w_k = 0, mu_k = 0;

    for (int i = 1; i <= 256; i++) {
        mu_T += i * P[i-1];
        if (i <= k) {
            mu_k += i * P[i-1];
            w_k += P[i-1];
        }
    }

    double coef1 = mu_T * w_k - mu_k;
    double coef2 = w_k * (1 - w_k);
    return coef1 * coef1 / coef2;
}

//computation of the threshold parameter, according to the definition of the given document for the automatic threshold selection method
double Image::find_threshold_parameter() {

    std::vector<double> list;
    for (double k = 1; k <= 255; k += 1)
        list.push_back(k);

    std::vector<double> P = get_P();

    double valmax = sigma2_B(1, P);
    double argmax = 0;
    for (double k : list) {
        if (sigma2_B(k, P) > valmax){
            argmax = k;
            valmax = sigma2_B(k, P);
        }
    }
    return argmax / 255.0;
}

void Image::binarize(double k) {
    if (k == -1)
        k = find_threshold_parameter();
    //for every pixel in the Image,
    for (int i = 0; i < im.rows; i++) {
        for (int j = 0; j < im.cols; j++) {
            //pixels of intensity <= k are set to black
            if (get_value(i, j) <= k)
                black(i, j);
                //other pixels are set to white
            else
                white(i, j);
        }
    }
}

void Image::contrast(double sub1, double sup1, double sub2, double sup2) {
    double val;
    for (int i = 0; i < im.rows; i++) {
        for (int j = 0; j < im.cols; j++) {
            val = get_value(i, j);
            if (val >= sup2)
                white(i, j);
            if (val > sub2 && val <= sup2)
                set_value(i, j, (1 + val) / 2);

            if (val <= sub1)
                black(i, j);
            if (val >= sub1 && val <= sup1)
                set_value(i, j, val / 2);

        }
    }
}

std::vector<double> Image::voisins2_width(int i, int j) {
    std::vector<double> V;
    if (j != 0) { V.push_back(get_value(i, j - 1)); }
    if (j != im.cols - 1) { V.push_back(get_value(i, j + 1)); }
    return V;
}

std::vector<double> Image::voisins2_height(int i, int j) {
    std::vector<double> V;
    if (i != 0) { V.push_back(get_value(i - 1, j)); }
    if (i != im.rows - 1) { V.push_back(get_value(i + 1, j)); }
    return V;
}

std::vector<double> Image::voisins4(int i, int j) {
    std::vector<double> V;
    if (i != 0) { V.push_back(get_value(i - 1, j)); }
    if (i != im.rows - 1) { V.push_back(get_value(i + 1, j)); }
    if (j != 0) { V.push_back(get_value(i, j - 1)); }
    if (j != im.cols - 1) { V.push_back(get_value(i, j + 1)); }
    return V;
}

std::vector<double> Image::voisins8(int i, int j) {
    std::vector<double> V;
    for (int k = -1; k <= 1; k++) {
        for (int n = -1; n <= 1; n++) {
            if ((i+k >= 0) && (i+k <= im.rows-1) && (j+n >= 0) && (j+n <= im.cols-1)){
                if (k!=0 || n!=0){
                    V.push_back(get_value(i + k, j + n));
                }
            }
        }
    }
    return V;
}
std::vector<double> Image::voisins12(int i, int j) {
    std::vector<double> V = voisins8(i, j);
    if (i >= 2) { V.push_back(get_value(i - 2, j)); }
    if (i <= im.rows - 3) { V.push_back(get_value(i + 2, j)); }
    if (j >= 2) { V.push_back(get_value(i, j - 2)); }
    if (j <= im.cols - 3) { V.push_back(get_value(i, j + 2)); }
    return V;
}

std::vector<double> Image::voisins24(int i, int j) {

    std::vector<double> V;
    for (int k = -2; k <= 2; k++) {
        for (int n = -2; n <= 2; n++) {
            if ((i + k >= 0) && (i + k <= im.rows - 1) && (j + n >= 0) && (j + n <= im.cols - 1)) {
                if (k != 0 || n != 0) {
                    V.push_back(get_value(i + k, j + n));
                }
            }
        }
    }
    return V;
}


std::vector<double> Image::neighbors(std::string neighborhood, int i, int j) {
    if (neighborhood == "2w") { return voisins2_width(i, j); }
    else if (neighborhood == "2h") { return voisins2_height(i, j); }
    else if (neighborhood == "4") { return voisins4(i, j); }
    else if (neighborhood == "8"){ return voisins8(i, j); }
    else if (neighborhood == "12") { return voisins12(i, j); }
    else if (neighborhood == "24") { return voisins24(i, j); }
    else{
        std::cout << "Invalid neighborhood. Exiting..." << std::endl;
        exit(1);
    }
}


bool Image::fits(std::string neighborhood, int i, int j){
    std::vector<double> V = neighbors(neighborhood, i, j);
    return (* min_element(V.begin(), V.end()) == 1); //min_elem = 1 <=> it fits
}

bool Image::hits(std::string neighborhood, int i, int j) {
    std::vector<double> V = neighbors(neighborhood, i, j);
    return (* max_element(V.begin(), V.end()) == 1); //max_elem = 1 <=> it hits
}

cv::Mat Image::bin_erosion(std::string neighborhood) {
    //creates a matrix that has the same size as the Image, and sets pixel intensity to be grayscale and of type double (CV_64FC1)
    cv::Mat S = cv::Mat::zeros(cv::Size(im.cols, im.rows), cv::IMREAD_GRAYSCALE);
    S.convertTo(S, CV_64FC1);
    //loop that goes through every pixel of the Image
    for (int i = 0; i < S.rows; i++) {
        for (int j = 0; j < S.cols; j++) {
            S.at<double>(i, j) = fits(neighborhood, i, j); //sets the current pixel to 1 if the structural element fits
        }
    }
    return S;
}

cv::Mat Image::bin_dilatation(std::string neighborhood) {
    //creates a matrix that has the same size as the Image, and sets pixel intensity to be grayscale and of type double (CV_64FC1)
    cv::Mat S = cv::Mat::zeros(cv::Size(im.cols, im.rows), cv::IMREAD_GRAYSCALE);
    S.convertTo(S, CV_64FC1);
    //loop that goes through every pixel of the Image
    for (int i = 0; i < S.rows; i++) {
        for (int j = 0; j < S.cols; j++) {
            S.at<double>(i, j) = hits(neighborhood, i, j); //sets the current pixel to 1 if the structural element hits
        }
    }
    return S;
}

cv::Mat Image::bin_opening(std::string neighborhood) {
    Image IM(bin_erosion(neighborhood));    //erodes
    return IM.bin_dilatation(neighborhood); //dilates
}

cv::Mat Image::bin_closing(std::string neighborhood) {
    Image IM(bin_dilatation(neighborhood)); //dilates
    return IM.bin_erosion(neighborhood); 	//erodes
}


cv::Mat Image::gs_dilatation(std::string neighborhood) {
    //creates a matrix that has the same size as the Image, and sets pixel intensity to be grayscale and of type double (CV_64FC1)
    cv::Mat S = cv::Mat::zeros(cv::Size(im.cols, im.rows), cv::IMREAD_GRAYSCALE);
    S.convertTo(S, CV_64FC1);
    std::vector<double> V;
    //loop that goes through every pixel of the Image
    for (int i = 0; i < S.rows; i++) {
        for (int j = 0; j < S.cols; j++) {
            V = neighbors(neighborhood, i, j);      					// V contains all intensity values of pixel in the neighborhood of (i,j)
            S.at<double>(i, j) = *max_element(V.begin(), V.end());		// sets the value of new pixel to the maximum neighbor's intensity
        }
    }
    return S;
}

cv::Mat Image::gs_erosion(std::string neighborhood) {
    //creates a matrix that has the same size as the Image, and sets pixel intensity to be grayscale and of type double (CV_64FC1)
    cv::Mat S = cv::Mat::zeros(cv::Size(im.cols, im.rows), cv::IMREAD_GRAYSCALE);
    S.convertTo(S, CV_64FC1);
    std::vector<double> V;
    //loop that goes through every pixel of the Image
    for (int i = 0; i < S.rows; i++) {
        for (int j = 0; j < S.cols; j++) {
            V = neighbors(neighborhood, i, j);      					// V contains all intensity values of pixel in the neighborhood of (i,j)
            S.at<double>(i, j) = *min_element(V.begin(), V.end());		// sets the value of new pixel to the minimum neighbor's intensity
        }
    }
    return S;
}

cv::Mat Image::gs_opening(std::string neighborhood) {
    Image IM(gs_erosion(neighborhood));		//erodes
    return IM.gs_dilatation(neighborhood);	//dilates
}

cv::Mat Image::gs_closing(std::string neighborhood) {
    Image IM(gs_dilatation(neighborhood));	//dilates
    return IM.gs_erosion(neighborhood);		//erodes
}


cv::Mat Image::gs_faded_erosion(std::string neighborhood, cv::Point2d center, double alpha, double beta, double scale) {
    //creates a matrix that has the same size as the Image, and sets pixel intensity to be grayscale and of type double (CV_64FC1)
    cv::Mat S = cv::Mat::zeros(cv::Size(im.cols, im.rows), cv::IMREAD_GRAYSCALE);
    S.convertTo(S, CV_64FC1);
    std::vector<double> V; double coef;
    //loop that goes through every pixel of the Image
    for (int i = 0; i < S.rows; i++) {
        for (int j = 0; j < S.cols; j++) {
            coef = anisotropic_logistic_v2(center, cv::Point(i, j), alpha, beta, scale);	    //computes the coefficient of anisotropic function for the fading
            V = neighbors(neighborhood, i, j);      	// V contains all intensity values of pixel in the neighborhood of (i,j)

            // sets the value of new pixel to be in between the minimum neighbor's intensity and the original intensity (depending on the value of the coefficient function)
            S.at<double>(i, j) = coef * (*min_element(V.begin(), V.end())) + (1 - coef) * get_value(i, j);
        }
    }
    return S;
}

cv::Mat Image::gs_faded_dilatation(std::string neighborhood, cv::Point2d center, double alpha, double beta, double scale) {
    //creates a matrix that has the same size as the Image, and sets pixel intensity to be grayscale and of type double (CV_64FC1)
    cv::Mat S = cv::Mat::zeros(cv::Size(im.cols, im.rows), cv::IMREAD_GRAYSCALE);
    S.convertTo(S, CV_64FC1);
    std::vector<double> V; double coef;
    //loop that goes through every pixel of the Image
    for (int i = 0; i < S.rows; i++) {
        for (int j = 0; j < S.cols; j++) {
            coef = anisotropic_logistic_v2(center, cv::Point(i, j), alpha, beta, scale);	    //computes the coefficient of anisotropic function for the fading
            V = neighbors(neighborhood, i, j);      	// V contains all intensity values of pixel in the neighborhood of (i,j)

            // sets the value of new pixel to be in between the minimum neighbor's intensity and the original intensity (depending on the value of the coefficient function)
            S.at<double>(i, j) = (1 - coef) * (*max_element(V.begin(), V.end())) + coef * get_value(i, j);
        }
    }
    return S;
}


cv::Mat Image::add_black_dot(int i, int j, cv::Mat S) {
    //n and k go through a 4x4 square
    for (int k = -1; k < 3; k++) {
        for (int n = -1; n < 3; n++) {
            //corners of the square will not change (shape of the circle)
            if (k + n != -2 && k + n != 4 && (k != 2 || n != -1) && (k != -1 || n != 2)){
                //checks that the pixel is in the Image
                if ((i + k >= 0) && (i + k <= im.rows - 1) && (j + n >= 0) && (j + n <= im.cols - 1))
                    S.at<double>(i + k, j + n) = 0; //sets the pixel to black
            }
        }
    }
    return S;
}

void Image::make_dry(){
    cv::Point2d center = get_center_pressure();// gets the center of pressure of the Image

    //creates a matrix that has the same size as the Image, and sets pixel intensity to be grayscale and of type double (CV_64FC1)
    cv::Mat S = cv::Mat::zeros(cv::Size(im.cols, im.rows), cv::IMREAD_GRAYSCALE);
    S.convertTo(S, CV_64FC1);
    S = im.clone();

    //for every pixel in the Image,
    for (int i = 0; i < im.rows - 1; i++) {
        for (int j = 0; j < im.cols - 1; j++) {
            //if the pixel is black enough,
            if (get_value(i, j) < 0.35) {
                //if the pixel is far enough from the center of pressure,
                if (distorted_distance(center, cv::Point2d(i, j), 0.7) > 60) {
                    //with a probability of 93%,
                    if (rand() % 101 > 95) {
                        add_black_dot(i, j, S); //we add a black dot located at the current pixel to the temporary matrix
                    }
                }
            }
        }
    }
    im = S;
}
