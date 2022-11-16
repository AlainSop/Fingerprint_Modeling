#include "LinearFiltering.h"

Image convolutionMatrix(const  Image &img,const Image &kernel){
    int border = kernel.im.rows/2;
    int rows = img.im.rows;
    int cols = img.im.cols;
    Image kernel_copy = kernel.symmetry_x().symmetry_y();
    Image bordered_img = borderImage(img, border,1);
    Image return_img(cv::Mat::zeros(img.im.rows,img.im.cols,CV_64FC1));
    for(int row = 0; row<rows;row++){
        for(int col = 0; col<cols;col++){
            cv::Mat temp = bordered_img.im(cv::Range(row, row+kernel_copy.im.rows), cv::Range(col, col+kernel_copy.im.cols)).clone();
            return_img.im.at<double>(row,col) = cv::sum(temp.mul(kernel_copy.im))[0];
        }
    }
    return return_img;
}

Image convolutionMatrixAni(const  Image &img, const Image &kernel){
    int border = kernel.im.rows/2;
    int rows = img.im.rows;
    int cols = img.im.cols;
    cv::Point2d pt= img.get_center_pressure();
    Image kernel_cp = kernel.symmetry_x().symmetry_y();
    Image bordered_img = borderImage(img, border,1);
    Image return_img(cv::Mat::zeros(img.im.rows,img.im.cols,CV_64FC1));
    for(int row = 0; row<rows;row++){
        for(int col = 0; col<cols;col++){
            cv::Point2d pt_local(row,col);
            Image kernel_copy(kernel_cp.im*(3-2*anisotropic_logistic(pt,pt_local,100,0.05)));
            cv::Mat temp = bordered_img.im(cv::Range(row, row+kernel_copy.im.rows), cv::Range(col, col+kernel_copy.im.cols)).clone();
            return_img.im.at<double>(row,col) = cv::sum(temp.mul(kernel_copy.im))[0];
        }
    }
    return return_img;
}

Image dft(const Image & img){
//    cv::Mat planes[] = {cv::Mat_<double>(img.im), cv::Mat::zeros(img.im.size(), CV_64FC1)};
//    cv::Mat complexI;    //Complex plane to contain the DFT coefficients {[0]-Real,[1]-Img}
//    merge(planes, 2, complexI);
//    cv::dft(complexI, complexI);  // Applying DFT
    cv::Mat padded;                            //expand input Image to optimal size
    int m = cv::getOptimalDFTSize( img.im.rows );
    int n = cv::getOptimalDFTSize( img.im.cols ); // on the border add zero values
    cv::copyMakeBorder(img.im, padded, 0, m - img.im.rows, 0, n - img.im.cols, cv::BORDER_CONSTANT, cv::Scalar::all(0));
    cv::Mat planes[] = {cv::Mat_<double>(padded), cv::Mat::zeros(padded.size(), CV_64FC1)};
    cv::Mat complexI;
    cv::merge(planes, 2, complexI);         // Add to the expanded another plane with zeros
    cv::dft(complexI, complexI);
    return Image(complexI);
}


Image normalized_dft(const Image & I){
    cv::Mat padded;                            //expand input Image to optimal size
    int m = cv::getOptimalDFTSize( I.im.rows );
    int n = cv::getOptimalDFTSize( I.im.cols ); // on the border add zero values
    cv::copyMakeBorder(I.im, padded, 0, m - I.im.rows, 0, n - I.im.cols, cv::BORDER_CONSTANT, cv::Scalar::all(0));
    cv::Mat planes[] = {cv::Mat_<float>(padded), cv::Mat::zeros(padded.size(), CV_32F)};
    cv::Mat complexI;
    cv::merge(planes, 2, complexI);         // Add to the expanded another plane with zeros
    cv::dft(complexI, complexI);            // this way the result may fit in the source cv::Matrix
    // compute the magnitude and switch to logarithmic scale
    // => log(1 + sqrt(Re(DFT(I))^2 + Im(DFT(I))^2))
    cv::split(complexI, planes);                   // planes[0] = Re(DFT(I), planes[1] = Im(DFT(I))
    cv::magnitude(planes[0], planes[1], planes[0]);// planes[0] = magnitude
    cv::Mat magI = planes[0];
    magI += cv::Scalar::all(1);                    // switch to logarithmic scale
    log(magI, magI);
    // crop the spectrum, if it has an odd number of rows or columns
    magI = magI(cv::Rect(0, 0, magI.cols & -1, magI.rows & -1));
    // rearcv::Range the quadrants of Fourier Image  so that the origin is at the Image center
    int cx = magI.cols/2;
    int cy = magI.rows/2;
    cv::Mat q0(magI, cv::Rect(0, 0, cx, cy));   // Top-Left - Create a ROI per quadrant
    cv::Mat q1(magI, cv::Rect(cx, 0, cx, cy));  // Top-Right
    cv::Mat q2(magI, cv::Rect(0, cy, cx, cy));  // Bottom-Left
    cv::Mat q3(magI, cv::Rect(cx, cy, cx, cy)); // Bottom-Right
    cv::Mat tmp;                           // swap quadrants (Top-Left with Bottom-Right)
    q0.copyTo(tmp);
    q3.copyTo(q0);
    tmp.copyTo(q3);
    q1.copyTo(tmp);                    // swap quadrant (Top-Right with Bottom-Left)
    q2.copyTo(q1);
    tmp.copyTo(q2);
    cv::normalize(magI, magI, 0, 1, cv::NORM_MINMAX); // Transform the matrix with float values into a
    // viewable Image form (float between values 0 and 1).
    return Image(magI);
}

Image idft(const Image & img){
    cv::Mat padded;                            //expand input Image to optimal size
    int m = cv::getOptimalDFTSize( img.im.rows );
    int n = cv::getOptimalDFTSize( img.im.cols ); // on the border add zero values
    cv::copyMakeBorder(img.im, padded, 0, m - img.im.rows, 0, n - img.im.cols, cv::BORDER_CONSTANT, cv::Scalar::all(0));
    cv::Mat invDFT, invDFTcvt;
    cv::idft(padded.clone(), invDFT, cv::DFT_SCALE | cv::DFT_REAL_OUTPUT ); // Applying IDFT
    invDFT.convertTo(invDFTcvt, CV_64FC1);
    return Image(invDFTcvt);
}

Image convolutionDft(const Image & img, Image & kernel){
    int border = kernel.im.rows/2;
    int rows = img.im.rows;
    int cols = img.im.cols;
//    Image kernel_copy = kernel.symmetry_x().symmetry_y();
    Image bordered_img = borderImage(img, border,1);
    Image return_img(bordered_img);
    for(int row = border; row<rows+border;row++){
        for(int col = border; col<cols+border;col++){
            cv::Mat temp = bordered_img.im(cv::Range(row-border, row+border+1), cv::Range(col-border, col+border+1)).clone();
            Image fft(temp);
            fft = dft(fft);
            Image kernel_fft(kernel);
            kernel_fft = dft(kernel_fft);
            Image fft_prod(fft.im.mul(kernel_fft.im));
            Image convo = idft(fft_prod);
            convo.im.copyTo(return_img.im(cv::Range(row-border, row+border+1), cv::Range(col-border, col+border+1)));
        }
    }
    Image return_img2(return_img.im(cv::Range(border, rows+border), cv::Range(border, cols+border)));
//    std::cout <<"taille Image retournée"<<return_img2.im.size()<< std::endl;
//    std::cout <<"taille Image en entrée"<<img.im.size()<< std::endl;
    return return_img2;
}


