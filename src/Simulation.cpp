#include "Simulation.h"

double angle(cv::Point2d pt, cv::Point2d center) {
    // we translate the point by the center to compute the angle with respect to the Oy axis
    cv::Point2d new_pt = pt - center;

    // to avoid dividing by 0
    if (new_pt.x == 0) {
        if (new_pt.y > 0) {
            return 3 * M_PI / 2;
        }
        if (new_pt.y < 0) {
            return M_PI / 2;
        }
        else {
            return 0;
        }
    }

        // if we are in the right half of the unit circle
    else if (new_pt.x > 0) {
        return fmod(-atan(new_pt.y / new_pt.x) + 2. * M_PI, 2. * M_PI);
    }

        // if we are in the left half of the unit circle
    else {
        return fmod(-atan(new_pt.y / new_pt.x) + M_PI, 2. * M_PI);
    }
}

double linear(double x, double x0, double x1, double fx0, double fx1) {
    if (x1 == x0) {
        return fx0;
    }
    return fx0 + (fx1 - fx0) / (x1 - x0) * (x - x0);
}

cv::Point3d generate_random_fluctuation(double theta_k, double m, double eps1, double eps2, double eps3, double eps4) {
    // random number generators
    std::random_device mch;
    std::default_random_engine generator(mch());
    std::exponential_distribution<double> exp_dist(m / (2. * M_PI));
    std::uniform_real_distribution<double> unif_dist_A(1. - eps1, 1. + eps2);
    std::uniform_real_distribution<double> unif_dist_B(1.6 - eps3, 1.6 + eps4);

    theta_k += exp_dist(generator);
    double A_theta_k = unif_dist_A(generator);
    double B_theta_k = unif_dist_B(generator);
    return cv::Point3d(theta_k, A_theta_k, B_theta_k);
}

std::vector<cv::Point3d> generate_random_boundary(double m, double eps1, double eps2, double eps3, double eps4) {

    std::vector<cv::Point3d> fluctuations;

    cv::Point3d pt = generate_random_fluctuation(0, m, eps1, eps2, eps3, eps4);

    // we generate a random sample in [0, 2pi]
    while (pt.x <= 2. * M_PI) {
        fluctuations.push_back(pt);
        pt = generate_random_fluctuation(pt.x, m, eps1, eps2, eps3, eps4);
    }

    return fluctuations;
}

cv::Point2d simulate_A_B(double theta, std::vector<cv::Point3d> fluctuations) {
    double theta_0, theta_1, A_theta_0, A_theta_1, B_theta_0, B_theta_1;
    long unsigned int i = 0;

    // constant fluctuation
    if (fluctuations.size() == 0) {
        return cv::Point2d(1., 1.6);
    }

    // we find the interval [theta_k, theta_k+1] where theta is
    while (theta > fluctuations[i].x) {
        i++;
        if (i == fluctuations.size()) {
            break;
        }
    }

    // periodic interpolation at left
    if (i == 0) {
        theta_0 = fluctuations.back().x - 2 * M_PI;
        theta_1 = fluctuations[0].x;
        A_theta_0 = fluctuations.back().y;
        A_theta_1 = fluctuations[0].y;
        B_theta_0 = fluctuations.back().z;
        B_theta_1 = fluctuations[0].z;
    }
        // periodic interpolation at right
    else if (i == fluctuations.size()) {
        theta_0 = fluctuations[i - 1].x;
        theta_1 = fluctuations[0].x + 2 * M_PI;
        A_theta_0 = fluctuations[i - 1].y;
        A_theta_1 = fluctuations[0].y;
        B_theta_0 = fluctuations[i - 1].z;
        B_theta_1 = fluctuations[0].z;
    }
        // interpolation elsewhere
    else {
        theta_0 = fluctuations[i - 1].x;
        theta_1 = fluctuations[i].x;
        A_theta_0 = fluctuations[i - 1].y;
        A_theta_1 = fluctuations[i].y;
        B_theta_0 = fluctuations[i - 1].z;
        B_theta_1 = fluctuations[i].z;
    }

    //cout << theta_0 << " | " << theta_1 << endl;
    return cv::Point2d(linear(theta, theta_0, theta_1, A_theta_0, A_theta_1), linear(theta, theta_0, theta_1, B_theta_0, B_theta_1));
}


///////////////////// DISTANCES /////////////////////
double euclidian_dist(cv::Point2d center, cv::Point2d pt) {
    return sqrt(pow(center.x - pt.x, 2) + pow(center.y - pt.y, 2));
}

double distorted_distance(cv::Point2d center, cv::Point2d pt, double scale_y) {
    return sqrt(pow(center.x - pt.x, 2) + pow(scale_y * (center.y - pt.y), 2));
}

double elliptic_distance(cv::Point2d center, cv::Point2d pt, double scale_x, double scale_y) {
    return sqrt(pow(scale_x * (center.x - pt.x), 2) + pow(scale_y * (center.y - pt.y), 2));
}

///////////////////// ISOTROPIC FUNCTIONS /////////////////////
double isotropic_gaussian(cv::Point2d center, cv::Point2d pt, double alpha, double beta ) {
    double r = euclidian_dist(center, pt);
    return exp(-alpha * r * r);
}

double isotropic_inverted_poly(cv::Point2d center, cv::Point2d pt, double alpha, double beta) {
    double r = euclidian_dist(center, pt);
    return 1 / (pow(r, alpha) + 1);
}

double isotropic_logistic(cv::Point2d center, cv::Point2d pt, double alpha, double beta) {
    double r = euclidian_dist(center, pt);
    return (1 + exp(-alpha * beta)) / (1 + exp(beta * (r - alpha)));
}

///////////////////// ANISOTROPIC FUNCTIONS /////////////////////
double anisotropic_logistic(cv::Point2d center, cv::Point2d pt, double alpha, double beta) {
    double r = distorted_distance(center, pt, 1.6);
    return (1 + exp(-alpha * beta)) / (1 + exp(beta * (r - alpha)));
}

double anisotropic_logistic_v2(cv::Point2d center, cv::Point2d pt, double alpha, double beta, double scale) {
    double r = distorted_distance(center, pt, scale);
    return (1 + exp(-alpha * beta)) / (1 + exp(beta * (r - alpha)));
}

///////////////////// ANISOTROPIC RANDOM FUNCTIONS /////////////////////

double random_anisotropic_logistic(cv::Point2d center, cv::Point2d pt, std::vector<cv::Point3d> fluctuations, double alpha, double beta) {
    double theta = angle(pt, center);
    cv::Point2d A_B = simulate_A_B(theta, fluctuations);
    double r = elliptic_distance(center, pt, A_B.x, A_B.y);
    return (1 + exp(-alpha * beta)) / (1 + exp(beta * (r - alpha)));
}

///////////////////// CENTER GETTERS /////////////////////


cv::Point2d Image::get_center_coordinates() {
    return cv::Point2d((im.rows + 1) / 2, (im.cols + 1) / 2);
}


cv::Point2d Image::get_loc_highest_pressure() {
    double mini, maxi;
    cv::Point minloc, maxloc;
    minMaxLoc(im, &mini, &maxi, &minloc, &maxloc);
    return (cv::Point2d)minloc;
}


cv::Point2d Image::get_center_pressure()const{
    cv::Point2d spot;
    double min_intensity = 1;
    double average_intensity = 0;
    for (int i = 10; i < im.rows-10; i += 10){
        for(int j = 10; j < im.cols-10; j += 10){
            for(int k = -10; k < 11; k++){
                for(int l = -10; l < 11; l++){
                    average_intensity += im.at<double>(i + k, j + l);
                }
            }
            average_intensity /= 441;
            if (average_intensity < min_intensity){
                min_intensity = average_intensity;
                spot.x = i;
                spot.y = j;
            }
            average_intensity = 0;

        }
    }
    return spot;
}


void Image::draw_center_pressure() {
    cv::Point2d spot = get_center_pressure();
    for (int i = 0; i < im.rows; i++) {
        set_value(i, spot.y - 1, 1);
        set_value(i, spot.y, 1);
        set_value(i, spot.y + 1, 1);
    }
    for (int j = 0; j < im.cols; j++) {
        set_value(spot.x - 1, j, 1);
        set_value(spot.x, j, 1);
        set_value(spot.x + 1, j, 1);
    }
}



cv::Point2d rotation(cv::Point2d center, cv::Point2d pt, double deg) { // deg : angle of rotation in degree
    double theta = deg * 3.141592 / 180;
    return cv::Point2d(cos(theta) * (pt.x - center.x) - sin(theta) * (pt.y - center.y) + center.x, sin(theta) * (pt.x - center.x) + cos(theta) * (pt.y - center.y) + center.y);
}

///////////////////// SIMULATION OF WEAK_FINGER /////////////////////
void Image::weak_finger(cv::Point2d center, cv::Point2d pt1, cv::Point2d pt2, Fonction f, double alpha, double beta, double deg) {

    cv::Point2d pt;
    double coef;

    negatif();

    for (int i = pt1.y; i <= pt2.y; i++) {
        for (int j = pt1.x; j <= pt2.x; j++) {

            pt = rotation(center, cv::Point2d(i, j), deg);
            coef = f(center, pt, alpha, beta);
            set_value(i, j, coef * get_value(i, j));
        }
    }
    negatif();
}


void Image::weak_finger_deformations(cv::Point2d center, cv::Point2d pt1, cv::Point2d pt2, double alpha, double beta, double deg, double m, double eps1, double eps2, double eps3, double eps4) {

    cv::Point2d pt;
    double coef;

    negatif();
    std::vector<cv::Point3d> fluctuations = generate_random_boundary(m, eps1, eps2, eps3, eps4); // we generate the random fluctuations
    for (int i = pt1.y; i <= pt2.y; i++) {
        for (int j = pt1.x; j <= pt2.x; j++) {
            pt = rotation(center, cv::Point2d(i, j), deg);
            coef = random_anisotropic_logistic(center, pt, fluctuations, alpha, beta); // we apply the random anisotropic logistic function
            set_value(i, j, coef * get_value(i, j));
        }
    }
    negatif();
}
