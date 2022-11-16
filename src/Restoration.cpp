#include "Restoration.h"

Image squareImage(const Image &in_image, int border, double color){
    int rows = in_image.im.rows;
    int cols = in_image.im.cols;
    int side = std::max(rows,cols)+2*border;
    int vertical_border = (side-rows)/2 ;
    int horizontal_border = (side-cols)/2;
    cv::Mat squared_matrix;
    cv::copyMakeBorder(in_image.im,squared_matrix,vertical_border,vertical_border,horizontal_border,horizontal_border,cv::BORDER_CONSTANT,color);
    return Image(squared_matrix);
}

Image borderImage(const Image &in_image, int border, int border_type){
    cv::Mat bordered_matrix;
    if(border_type)cv::copyMakeBorder(in_image.im,bordered_matrix,border,border,border,border,cv::BORDER_REPLICATE);
    else cv::copyMakeBorder(in_image.im,bordered_matrix,border,border,border,border,cv::BORDER_CONSTANT,1);
    return Image(bordered_matrix);
}

Image createPatch(int x, int y, int size, const Image &in_image){
    cv::Mat patch_matrix = cv::Mat::ones(size, size, CV_64FC1);
    Image new_im = borderImage(in_image, size/2,0);
    new_im.im(cv::Range(x , x + size ), cv::Range(y , y + size)).copyTo(
            patch_matrix);
    return Image(patch_matrix);
}


cv::Mat DicPatches::printDic(){
    int size = (int) sqrt(patches.size());
    int new_size = 11;
    cv::Mat display_matrix = cv::Mat::zeros(new_size*size,new_size*size,CV_64FC1);
    int row = 0;
    int col = 0;
    for(auto & patch : patches){
        squareImage(patch,1,1).im.copyTo(display_matrix(cv::Rect(col,row,new_size,new_size)));
        row += new_size;
        if (row==new_size*size){
            col += new_size;
            row=0;
        }
    }
    return display_matrix;
}

double euclideanDistance(const Image &patch1, const Image &patch2, const Image& mask){
    cv::Mat mult = cv::Mat::ones(mask.im.size(),CV_64FC1);
    return cv::norm(patch1.im.mul(mult-mask.im),patch2.im.mul(mult-mask.im),cv::NORM_L2);
}

double restorePixelValue(const DicPatches &dic, const Image& patch, const Image& mask){
    double min_dist = 1000000000;
    const Image * closest_patch_ptr = &(*dic.patches.begin());
    for(const auto & patch_from_dic : dic.patches){
        double dist = euclideanDistance(patch, patch_from_dic, mask);
        if(dist<min_dist){
            min_dist = dist;
            closest_patch_ptr = &patch_from_dic;
        }
    }
    return closest_patch_ptr->im.at<double>(patch.im.rows/2,patch.im.cols/2);
}

Image restorePatchValue(const DicPatches &dic, const Image& patch, const Image& mask){
    double min_dist = 1000000000;
    const Image * closest_patch_ptr = &(*dic.patches.begin());
    for(const auto & patch_from_dic : dic.patches){
        double dist = euclideanDistance(patch, patch_from_dic, mask);
        if(dist<min_dist){
            min_dist = dist;
            closest_patch_ptr = &patch_from_dic;
        }
    }
    return *closest_patch_ptr;
}

Image restore2(Image &in_image, Image mask, DicPatches &dic){
    Image new_im(in_image);
    int rows = mask.im.rows;
    int cols = mask.im.cols;

    for (int row = 0; row < rows; row++) {
        for (int col = 0; col < cols; col++) {
            if(mask.get_value(row,col) ==1)new_im.set_value(row,col,1);
        }
    }
    for(int i=0;i<40;i++) {//Arbitrary condition to restore a ring around the fingerprint
        Image new_mask(new_im.boundary());
        Image priority(cv::Mat::zeros(rows,cols,CV_64FC1));
        if (cv::norm(mask.im,cv::NORM_L1)==0){
            break;
        }
        for (int row = 0; row < rows; row++) {
            for (int col = 0; col < cols; col++) {
                if (new_mask.get_value(row, col) >0 &&  mask.get_value(row, col) >0) {
                    Image patch = createPatch(row, col, 9, new_im);
                    std::pair<std::pair<double, double>, double> grad = computegradient(patch);
                    priority.set_value(row, col, abs(grad.first.first*cos(grad.second*M_PI/180) + grad.first.second*sin(grad.second*M_PI/180)));
                }
            }
        }
        while(cv::norm(priority.im,cv::NORM_L1)>0){
            double min, max;
            cv::Point min_loc, max_loc;
            cv::minMaxLoc(priority.im,&min,&max,&min_loc,&max_loc);
            Image small_mask = createPatch(max_loc.y,max_loc.x,9,mask);
            Image small_patch = createPatch(max_loc.y,max_loc.x,9,new_im);
            double value = restorePixelValue(dic, small_patch, small_mask);
            new_im.set_value(max_loc.y,max_loc.x,value);
            priority.set_value(max_loc.y,max_loc.x,0);
            mask.set_value(max_loc.y,max_loc.x,0);
        }
    }
    return new_im;
}


std::pair<std::pair<double,double>,double> computegradient(Image &patch){
    Image filter_vertical(cv::Mat::zeros(9,9,CV_64FC1));
    Image filter_horizontal(cv::Mat::zeros(9,9,CV_64FC1));
    for(int row =-4;row<5;row++){
        for(int col=-4;col<5;col++){
            filter_horizontal.set_value(row+4,col+4,col/(pow(abs(row),abs(row))+pow(abs(col),abs(col))));
            filter_vertical.set_value(row+4,col+4,row/(pow(abs(row),abs(row))+pow(abs(col),abs(col))));
        }
    }
    double gy = cv::sum(filter_vertical.im.mul(patch.im.clone()))[0];
    double gx = cv::sum(filter_horizontal.im.mul(patch.im.clone()))[0];
    double angle = atan2(gy,gx)*180/M_PI;
    std::pair<double,double> ret1(gx,gy);
    std::pair<std::pair<double,double>,double> ret(ret1,angle);
    return ret;
}

Image restoreImage(Image &in_image, Image mask, DicPatches &dic, int parcours){
    int rows = mask.im.rows;
    int cols = mask.im.cols;
    for(int row =0;row<rows;row++){
        for(int col=0;col<cols;col++){
            if(mask.im.at<double>(row,col)==1) in_image.set_value(row,col,1);
        }
    }

    Image new_im(squareImage(in_image,10,1));
    Image new_mask(squareImage(mask, 10,0));
    switch(parcours){
        case 1:// Parcours bord horizontal supérieur, puis inférieur puis vertical gauche puis vertical droite
        {
            int k = 0;
            while(norm(new_mask.im)>0){
                for(int col=k;col<cols-k;col++){
                    if(new_mask.im.at<double>(k,col)==1){
                        Image patch = createPatch(k,col,9,new_im);
                        new_im.set_value(k,col, restorePixelValue(dic, patch, new_mask.im(cv::Range(k-4,k+5),cv::Range(col-4,col+5))) );
                        new_mask.set_value(k,col,0);
                    }
                    if(new_mask.im.at<double>(rows-1-k,col)==1){
                        Image patch = createPatch(rows-1-k,col,9,new_im);
                        new_im.set_value(rows-1-k,col,restorePixelValue(dic, patch, new_mask.im(cv::Range(rows-1-k-4,rows-1-k+5),cv::Range(col-4,col+5))));
                        new_mask.set_value(rows-1-k,col,0);
                    }
                }

                for(int row=k+1;row<rows-k-1;row++){
                    if(new_mask.im.at<double>(row,k)==1){
                        Image patch = createPatch(row,k,9,new_im);
                        new_im.set_value(row,k, restorePixelValue(dic, patch, new_mask.im(cv::Range(row-4,row+5), cv::Range(k-4,k+5))));
                        new_mask.set_value(row,k,0);
                    }
                    if(new_mask.im.at<double>(row,cols-k-1)==1){
                        Image patch = createPatch(row,cols-k-1,9,new_im);
                        new_im.set_value(row,cols-k-1, restorePixelValue(dic, patch, new_mask.im(cv::Range(row-4,row+5), cv::Range(cols-k-1-4,cols-k-1+5))));
                        new_mask.set_value(row,cols-k-1,0) ;
                    }
                }
                k++;
            }
            break;
        }

        case 2: // parcours en spirale depuis le point en haut à gauche de l'Image jusqu'au centre
        {
            int l =0;
            while(norm(new_mask.im)>0) {
                for (int col = l; col < cols - l; col++) {
                    if (new_mask.im.at<double>(l, col)==1) {
                        Image patch = createPatch(l, col, 9, new_im);
                        new_im.set_value(l, col, restorePixelValue(dic, patch,
                                                                   new_mask.im(cv::Range(l - 4, l + 5), cv::Range(col - 4, col + 5))));
                        new_mask.set_value(l, col,0);
                    }
                }
                for (int row = l + 1; row < rows - l; row++) {
                    if (new_mask.im.at<double>(row, cols - 1 - l)==1) {
                        Image patch = createPatch(row, cols - 1 - l, 9, new_im);
                        new_im.set_value(row, cols - 1 - l, restorePixelValue(dic, patch, new_mask.im(cv::Range(row - 4, row + 5),
                                                                                                      cv::Range(cols - 1 -
                                                                                                                l - 4,
                                                                                                                cols - 1 -
                                                                                                                l + 5))));
                        new_mask.set_value(row, cols - 1 - l, 0);
                    }
                }
                for (int col = l + 1; col < cols - 1; col++) {
                    if (new_mask.im.at<double>(rows - 1 - l, cols - 1 - col)==1) {
                        Image patch = createPatch(rows - 1 - l, cols - 1 - col, 9, new_im);
                        new_im.set_value(rows - 1 - l, cols - 1 - col, restorePixelValue(dic, patch,
                                                                                         new_mask.im(cv::Range(
                                                                                                 rows - 1 -
                                                                                                 l - 4,
                                                                                                 rows - 1 -
                                                                                                 l + 5),
                                                                                                     cv::Range(cols -
                                                                                                               1 - col - 4,
                                                                                                               cols -
                                                                                                               1 - col + 5))));
                        new_mask.set_value(rows - 1 - l, cols - 1 - col, 0);
                    }
                }
                for (int row = l + 1; row < rows - 1; row++) {
                    if (new_mask.im.at<double>(rows - 1 - row, l)==1) {
                        Image patch = createPatch(rows - 1 - row, l, 9, new_im);
                        new_im.set_value(rows - 1 - row, l, restorePixelValue(dic, patch,
                                                                              new_mask.im(cv::Range(rows - 1 - row - 4,
                                                                                                    rows - 1 - row + 5),
                                                                                          cv::Range(l - 4, l + 5))));
                        new_mask.set_value(rows - 1 - row, l, 0);
                    }
                }
                l++;
            }
            break;}


        case 3: //Parcours en spirale depuis le milieu de l'Image
        {
            int half_row = (new_im.im.rows - 1) / 2;
            int half_col = (new_im.im.cols - 1) / 2;
            int row = half_row;
            int col = half_col;
            int d = 1;
            int m = 1;
            while (norm(new_mask.im) > 0) {

                while ( 2*(row - half_row) * d < m) {
                    if (new_mask.im.at<double>(row, col)==1) {
                        Image patch = createPatch(row, col, 9, new_im);
                        new_im.set_value(row, col, restorePixelValue(dic, patch,
                                                                     new_mask.im(cv::Range(row - 4, row + 5), cv::Range(col - 4, col + 5))));
                        new_mask.set_value(row, col,0);
                    }
                    row += d;
                }
                while (2*(col - half_col) * d < m) {
                    if (new_mask.im.at<double>(row, col)==1) {
                        Image patch = createPatch(row, col, 9, new_im);
                        new_im.set_value(row, col, restorePixelValue(dic, patch,
                                                                     new_mask.im(cv::Range(row - 4, row + 5), cv::Range(col - 4, col + 5))));
                        new_mask.set_value(row, col,0);
                    }
                    col += d;
                }

                d *= -1;
                m += 1;
            }
            break;
        }
        case 4: //Parcours linéaire depuis en haut à gauche
        {

            for (int row = 4; row < rows-4; row++) {
                for (int col = 4; col < cols-4; col++) {
                    if (new_mask.im.at<double>(row, col) == 1) {
                        Image patch = createPatch(row, col, 9, new_im);
                        new_im.set_value(row, col, restorePixelValue(dic, patch,
                                                                     new_mask.im(cv::Range(row - 4, row + 5), cv::Range(col - 4, col + 5)).clone()));
                        new_mask.set_value(row, col, 0);
                    }
                }
            }
            break;
        }
        default:
            for(int row =0;row<new_im.im.rows;row++) {
                for (int col = 0; col < new_im.im.cols; col++){
                    if(new_mask.im.at<double>(row, col)==1){
                        Image patch = createPatch(row, col,9, new_im);
                        new_im.set_value(row, col, restorePixelValue(dic, patch, new_mask.im(cv::Range(row-4,row+5),cv::Range(col-4,col+5)).clone()));
                        new_mask.set_value(row, col, 0);
                    }
                }
            }
            break;
    }
    return new_im;
}

cv::Mat Image::ring(){
    cv::Mat S1 = cv::Mat::zeros(cv::Size(im.cols, im.rows), cv::IMREAD_GRAYSCALE);
    cv::Mat S2 = cv::Mat::zeros(cv::Size(im.cols, im.rows), cv::IMREAD_GRAYSCALE);
    S1.convertTo(S1, CV_64FC1);
    S2.convertTo(S2, CV_64FC1);
    Image image0 = im.clone();
    Image image1(S1);
    Image image2(S2);

    image0.binarize(0.6);
    image1.im = image0.bin_erosion("24");
    for (int i = 0; i < 6; i++){ // composition of the erosion filter function
        image1.im = image1.bin_erosion("24");
    }

    image2.im = image0.bin_erosion("24");
    image2.im = image2.bin_erosion("24");
    for (int i = 0; i < 6; i++){ // composition of the dilatation filter function
        image2.im = image2.bin_dilatation("24");
    }

    for (int i = 0; i != im.rows; i++){
        for (int j = 0; j != im.cols; j++){
            if (image1.im.at<double>(i,j) != image2.im.at<double>(i,j)){
                image1.im.at<double>(i,j) = 1;
            }
            else{
                image1.im.at<double>(i,j) = 0;
            }
        }
    }

    return image1.im;
}

cv::Mat Image::boundary(){
    cv::Mat S1 = cv::Mat::zeros(cv::Size(im.cols, im.rows), cv::IMREAD_GRAYSCALE);
    cv::Mat S2 = cv::Mat::zeros(cv::Size(im.cols, im.rows), cv::IMREAD_GRAYSCALE);
    S1.convertTo(S1, CV_64FC1);
    S2.convertTo(S2, CV_64FC1);
    Image image0 = im.clone();
    Image image1(S1);
    Image image2(S2);

    image0.binarize(0.6);
    image1.im = image0.bin_erosion("24");
    for (int i = 0; i < 1; i++){ // composition of the erosion filter function
        image1.im = image1.bin_erosion("24");
    }

    image2.im = image0.bin_erosion("24");

    for (int i = 0; i < 1; i++){ // composition of the dilatation filter function
        image2.im = image2.bin_dilatation("24");
    }

    for (int i = 0; i != im.rows; i++){
        for (int j = 0; j != im.cols; j++){
            if (image1.im.at<double>(i,j) != image2.im.at<double>(i,j)){
                image1.im.at<double>(i,j) = 1;
            }
            else{
                image1.im.at<double>(i,j) = 0;
            }
        }
    }

    return image1.im;
}