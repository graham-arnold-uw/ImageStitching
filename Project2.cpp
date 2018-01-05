#include "mainwindow.h"
#include "math.h"
#include "ui_mainwindow.h"
#include <QtGui>
#include "Matrix.h"

/*******************************************************************************
    The following are helper routines with code already written.
    The routines you'll need to write for the assignment are below.
*******************************************************************************/

/*******************************************************************************
Draw detected Harris corners
    cornerPts - corner points
    numCornerPts - number of corner points
    imageDisplay - image used for drawing

    Draws a red cross on top of detected corners
*******************************************************************************/
void MainWindow::DrawCornerPoints(CIntPt *cornerPts, int numCornerPts, QImage &imageDisplay)
{
   int i;
   int r, c, rd, cd;
   int w = imageDisplay.width();
   int h = imageDisplay.height();

   for(i=0;i<numCornerPts;i++)
   {
       c = (int) cornerPts[i].m_X;
       r = (int) cornerPts[i].m_Y;

       for(rd=-2;rd<=2;rd++)
           if(r+rd >= 0 && r+rd < h && c >= 0 && c < w)
               imageDisplay.setPixel(c, r + rd, qRgb(255, 0, 0));

       for(cd=-2;cd<=2;cd++)
           if(r >= 0 && r < h && c + cd >= 0 && c + cd < w)
               imageDisplay.setPixel(c + cd, r, qRgb(255, 0, 0));
   }
}

/*******************************************************************************
Compute corner point descriptors
    image - input image
    cornerPts - array of corner points
    numCornerPts - number of corner points

    If the descriptor cannot be computed, i.e. it's too close to the boundary of
    the image, its descriptor length will be set to 0.

    I've implemented a very simple 8 dimensional descriptor.  Feel free to
    improve upon this.
*******************************************************************************/
void MainWindow::ComputeDescriptors(QImage image, CIntPt *cornerPts, int numCornerPts)
{
    int r, c, cd, rd, i, j;
    int w = image.width();
    int h = image.height();
    double *buffer = new double [w*h];
    QRgb pixel;

    // Descriptor parameters
    double sigma = 2.0;
    int rad = 4;

    // Computer descriptors from green channel
    for(r=0;r<h;r++)
       for(c=0;c<w;c++)
        {
            pixel = image.pixel(c, r);
            buffer[r*w + c] = (double) qGreen(pixel);
        }

    // Blur
    GaussianBlurImage(buffer, w, h, sigma);

    // Compute the desciptor from the difference between the point sampled at its center
    // and eight points sampled around it.
    for(i=0;i<numCornerPts;i++)
    {
        int c = (int) cornerPts[i].m_X;
        int r = (int) cornerPts[i].m_Y;

        if(c >= rad && c < w - rad && r >= rad && r < h - rad)
        {
            double centerValue = buffer[(r)*w + c];
            int j = 0;

            for(rd=-1;rd<=1;rd++)
                for(cd=-1;cd<=1;cd++)
                    if(rd != 0 || cd != 0)
                {
                    cornerPts[i].m_Desc[j] = buffer[(r + rd*rad)*w + c + cd*rad] - centerValue;
                    j++;
                }

            cornerPts[i].m_DescSize = DESC_SIZE;
        }
        else
        {
            cornerPts[i].m_DescSize = 0;
        }
    }

    delete [] buffer;
}

/*******************************************************************************
Draw matches between images
    matches - matching points
    numMatches - number of matching points
    image1Display - image to draw matches
    image2Display - image to draw matches

    Draws a green line between matches
*******************************************************************************/
void MainWindow::DrawMatches(CMatches *matches, int numMatches, QImage &image1Display, QImage &image2Display)
{
    int i;
    // Show matches on image
    QPainter painter;
    painter.begin(&image1Display);
    QColor green(0, 250, 0);
    QColor red(250, 0, 0);

    for(i=0;i<numMatches;i++)
    {
        painter.setPen(green);
        painter.drawLine((int) matches[i].m_X1, (int) matches[i].m_Y1, (int) matches[i].m_X2, (int) matches[i].m_Y2);
        painter.setPen(red);
        painter.drawEllipse((int) matches[i].m_X1-1, (int) matches[i].m_Y1-1, 3, 3);
    }

    QPainter painter2;
    painter2.begin(&image2Display);
    painter2.setPen(green);

    for(i=0;i<numMatches;i++)
    {
        painter2.setPen(green);
        painter2.drawLine((int) matches[i].m_X1, (int) matches[i].m_Y1, (int) matches[i].m_X2, (int) matches[i].m_Y2);
        painter2.setPen(red);
        painter2.drawEllipse((int) matches[i].m_X2-1, (int) matches[i].m_Y2-1, 3, 3);
    }

}


/*******************************************************************************
Given a set of matches computes the "best fitting" homography
    matches - matching points
    numMatches - number of matching points
    h - returned homography
    isForward - direction of the projection (true = image1 -> image2, false = image2 -> image1)
*******************************************************************************/
bool MainWindow::ComputeHomography(CMatches *matches, int numMatches, double h[3][3], bool isForward)
{
    int error;
    int nEq=numMatches*2;

    dmat M=newdmat(0,nEq,0,7,&error);
    dmat a=newdmat(0,7,0,0,&error);
    dmat b=newdmat(0,nEq,0,0,&error);

    double x0, y0, x1, y1;

    for (int i=0;i<nEq/2;i++)
    {
        if(isForward == false)
        {
            x0 = matches[i].m_X1;
            y0 = matches[i].m_Y1;
            x1 = matches[i].m_X2;
            y1 = matches[i].m_Y2;
        }
        else
        {
            x0 = matches[i].m_X2;
            y0 = matches[i].m_Y2;
            x1 = matches[i].m_X1;
            y1 = matches[i].m_Y1;
        }


        //Eq 1 for corrpoint
        M.el[i*2][0]=x1;
        M.el[i*2][1]=y1;
        M.el[i*2][2]=1;
        M.el[i*2][3]=0;
        M.el[i*2][4]=0;
        M.el[i*2][5]=0;
        M.el[i*2][6]=(x1*x0*-1);
        M.el[i*2][7]=(y1*x0*-1);

        b.el[i*2][0]=x0;
        //Eq 2 for corrpoint
        M.el[i*2+1][0]=0;
        M.el[i*2+1][1]=0;
        M.el[i*2+1][2]=0;
        M.el[i*2+1][3]=x1;
        M.el[i*2+1][4]=y1;
        M.el[i*2+1][5]=1;
        M.el[i*2+1][6]=(x1*y0*-1);
        M.el[i*2+1][7]=(y1*y0*-1);

        b.el[i*2+1][0]=y0;

    }
    int ret=solve_system (M,a,b);
    if (ret!=0)
    {
        freemat(M);
        freemat(a);
        freemat(b);

        return false;
    }
    else
    {
        h[0][0]= a.el[0][0];
        h[0][1]= a.el[1][0];
        h[0][2]= a.el[2][0];

        h[1][0]= a.el[3][0];
        h[1][1]= a.el[4][0];
        h[1][2]= a.el[5][0];

        h[2][0]= a.el[6][0];
        h[2][1]= a.el[7][0];
        h[2][2]= 1;
    }

    freemat(M);
    freemat(a);
    freemat(b);

    return true;
}

/*******************************************************************************
*******************************************************************************
*******************************************************************************

    HW 1 Routines

*******************************************************************************
*******************************************************************************
*******************************************************************************/
// Normalize the values of the kernel to sum-to-one
void NormalizeKernel(double *kernel, int kernelWidth, int kernelHeight)
{
    double denom = 0.000001; int i;
    for(i=0; i<kernelWidth*kernelHeight; i++)
        denom += kernel[i];
    for(i=0; i<kernelWidth*kernelHeight; i++)
        kernel[i] /= denom;
}

/*********************************************************************************/
// Convolve the image with the kernel
void MainWindow::Convolution(double* image, double *kernel, int width, int height, int kernelWidth, int kernelHeight, bool add)
/*
 * image: input image in matrix form of size (w*h)*3 having double values
 * kernel: 1-D array of kernel values
 * kernelWidth: width of the kernel
 * kernelHeight: height of the kernel
 * add: a boolean variable (taking values true or false)
*/
{
    //create copy of image for convolution with zero padding
    int bufferWidth = width+2*(kernelWidth/2);
    int bufferHeight = height+2*(kernelHeight/2);
    double* buffer = new double [bufferWidth*bufferHeight];

    for (int r = 0; r < height; r++)
        for (int c = 0; c < width; c++) {
            buffer[(r+(kernelHeight/2))*bufferWidth+(c+(kernelWidth/2))] = image[r*width+c];
        }

    //looping through image pixels for convolution on each pixel
    for (int r = 0; r < height; r++) {
        for (int c = 0; c < width; c++) {
           double rgb;
           rgb = 0.0;
           //for each kernel value
           for(int rd = 0; rd < kernelHeight; rd++) {
               for(int cd = 0; cd < kernelWidth; cd++)
               {
                   double weight = kernel[rd*kernelWidth + cd];

                   rgb += weight*buffer[(r+rd)*bufferWidth + (c+cd)];
               }
            }


           if (add) {
               rgb += 128.0;
           }

           image[r*width+c] = rgb;
        }
    }

    delete[] buffer;

}

/*******************************************************************************
*******************************************************************************
*******************************************************************************/


/*******************************************************************************
Blur a single channel floating point image with a Gaussian.
    image - input and output image
    w - image width
    h - image height
    sigma - standard deviation of Gaussian

    This code should be very similar to the code you wrote for assignment 1.
*******************************************************************************/
void MainWindow::GaussianBlurImage(double *image, int width, int height, double sigma)
{

    int radius = (int)(ceil(3*sigma));
    int size = 2*radius + 1; // This is the size of the kernel
    // Compute the kernel to convolve with the image
    double *kernel = new double [size*size];

    for(int rd=-radius;rd<=radius;rd++)
        for(int cd=-radius;cd<=radius;cd++)
        {
            double denom = 1/(2*M_PI*sigma);
            kernel[(rd + radius)*size + cd + radius] = denom * exp(-1*(rd*rd+cd*cd)/(2*sigma*sigma));
        }


    NormalizeKernel(kernel, size, size);

    MainWindow::Convolution(image, kernel, width, height, size, size, false);


}

/****************************************************/
/*******************************************************************************
Detect Harris corners.
    image - input image
    sigma - standard deviation of Gaussian used to blur corner detector
    thres - Threshold for detecting corners
    cornerPts - returned corner points
    numCornerPts - number of corner points returned
    imageDisplay - image returned to display (for debugging)
*******************************************************************************/
void MainWindow::HarrisCornerDetector(QImage image, double sigma, double thres, CIntPt **cornerPts, int &numCornerPts, QImage &imageDisplay)
{
    int r, c;
    int w = image.width();
    int h = image.height();
    double *bufferx = new double [w*h];
    double *buffery = new double [w*h];
    double *cornerImage = new double [w*h];
    double *Ixx = new double [w*h];
    double *Iyy = new double [w*h];
    double *Ixy = new double [w*h];
    QRgb pixel;

    numCornerPts = 0;

    // Compute the corner response using just the green channel
    for(int r=0; r<h; r++)
       for(int c = 0 ;c < w; c++)
        {
            pixel = image.pixel(c, r);

            bufferx[r*w + c] = (double) qGreen(pixel);
            buffery[r*w + c] = (double) qGreen(pixel);
            cornerImage[r*w + c] = (double) qGreen(pixel);
            Ixx[r*w + c] = 0.0;
            Iyy[r*w + c] = 0.0;
            Ixy[r*w + c] = 0.0;
        }

    // Compute the First derivative of an image along the horizontal direction
    double *xkernel = new double [3];
    xkernel[0] = -1;
    xkernel[1] = 0;
    xkernel[2] = 1;
    MainWindow::Convolution(bufferx, xkernel, w, h, 3, 1, false);

    // Compute the First derivative of an image along the vertical direction
    double *ykernel = new double [3];
    ykernel[0] = -1;
    ykernel[1] = 0;
    ykernel[2] = 1;
    MainWindow::Convolution(buffery, ykernel, w, h, 1, 3, false);

    // Compute 3 Harris Matrix Images
    for(int r=0;r<h;r++)
       for(int c=0;c<w;c++)
        {
            Ixx[r*w + c] = bufferx[r*w + c] * bufferx[r*w + c];
            Iyy[r*w + c] = buffery[r*w + c] * buffery[r*w + c];
            Ixy[r*w + c] = bufferx[r*w + c] * buffery[r*w + c];
        }

    //Smooth each Harris Matrix Image
    MainWindow::GaussianBlurImage(Ixx, w, h, sigma);
    MainWindow::GaussianBlurImage(Iyy, w, h, sigma);
    MainWindow::GaussianBlurImage(Ixy, w, h, sigma);

    // Compute matrix and threshold at each pixel
    for(int r=1; r < h - 1; r++)
       for(int c = 0;c < w; c++)
        {
            //calculate R function
            double R = 0.0;
            double detH = Ixx[r*w + c]*Iyy[r*w + c] - Ixy[r*w + c]*Ixy[r*w + c];
            double traceH = Ixx[r*w + c] + Iyy[r*w + c];
            if (traceH > 0.0) {
                R = detH/traceH;
            }

            cornerImage[r*w + c] = R;

            //scale R
            R *= 0.1;
            R += 50.0;

            if (R > 255.0) {
                R = 255.0;
            }

        }

    //threshold R and perform
    //no maximum suppression
    for(r=1;r<h-1;r++)
           for(c=1;c<w-1;c++)
            {
                int rd, cd;
                bool greaterFound = false;
                double currR = cornerImage[r*w + c];

                if(currR > thres)
                {
                    for(rd=-1;rd<=1;rd++)
                        for(cd=-1;cd<=1;cd++)
                        {
                            if(currR < cornerImage[(r+rd)*w + c + cd])
                                greaterFound = true;
                        }

                    if(!greaterFound)
                    {
                        numCornerPts++;
                        cornerImage[r*w + c] = 255.0;
                    }
                }
            }

    //file up corner points data structure
    *cornerPts = new CIntPt [numCornerPts];
    int index = 0;

    for(r=1; r<h-1;r++)
        for(c=1;c<w-1;c++)
        {
            if (cornerImage[r*w+c] == 255.0) {
                (*cornerPts)[index].m_X = (double) c;
                (*cornerPts)[index].m_Y = (double) r;
                index += 1;
            }
        }

    // Once you uknow the number of corner points allocate an array as follows:
    // *cornerPts = new CIntPt [numCornerPts];
    // Access the values using: (*cornerPts)[i].m_X = 5.0;
    //
    // The position of the corner point is (m_X, m_Y)
    // The descriptor of the corner point is stored in m_Desc
    // The length of the descriptor is m_DescSize, if m_DescSize = 0, then it is not valid.

    // Once you are done finding the corner points, display them on the image
    DrawCornerPoints(*cornerPts, numCornerPts, imageDisplay);

    delete [] bufferx;
    delete [] buffery;
    delete [] Ixx;
    delete [] Iyy;
    delete [] Ixy;
    delete [] cornerImage;
}




/*******************************************************************************
Find matching corner points between images.
    image1 - first input image
    cornerPts1 - corner points corresponding to image 1
    numCornerPts1 - number of corner points in image 1
    image2 - second input image
    cornerPts2 - corner points corresponding to image 2
    numCornerPts2 - number of corner points in image 2
    matches - set of matching points to be returned
    numMatches - number of matching points returned
    image1Display - image used to display matches
    image2Display - image used to display matches
*******************************************************************************/
void MainWindow::MatchCornerPoints(QImage image1, CIntPt *cornerPts1, int numCornerPts1,
                             QImage image2, CIntPt *cornerPts2, int numCornerPts2,
                             CMatches **matches, int &numMatches, QImage &image1Display, QImage &image2Display)
{


    // Compute the descriptors for each corner point.
    // You can access the descriptor for each corner point using cornerPts1[i].m_Desc[j].
    // If cornerPts1[i].m_DescSize = 0, it was not able to compute a descriptor for that point
    ComputeDescriptors(image1, cornerPts1, numCornerPts1);
    ComputeDescriptors(image2, cornerPts2, numCornerPts2);

    numMatches = numCornerPts1;
    *matches = new CMatches [numMatches];
    numMatches = 0;

    //find point with min distance descriptor
    //add to matches data structure
    for (int i = 0; i < numCornerPts1; i++) {
        if(cornerPts1[i].m_DescSize > 0) {
            double minDist = 9999999999.0;
            int minDistIndex = 0;
            for (int j = 0; j < numCornerPts2; j++)
            {
                if (cornerPts2[j].m_DescSize > 0) {
                    double distSum = 0.0;
                    //minDistIndex = j;
                    for(int k = 0; k < DESC_SIZE; k++)
                    {
                        double currDist = abs(cornerPts1[i].m_Desc[k] - cornerPts2[j].m_Desc[k]);
                        distSum += currDist;
                    }
                    if ( distSum < minDist)
                    {
                        minDist = distSum;
                        minDistIndex = j;
                    }
                }

            }
            (*matches)[numMatches].m_X1 = cornerPts1[i].m_X;
            (*matches)[numMatches].m_Y1 = cornerPts1[i].m_Y;
            (*matches)[numMatches].m_X2 = cornerPts2[minDistIndex].m_X;
            (*matches)[numMatches].m_Y2 = cornerPts2[minDistIndex].m_Y;
            numMatches++;
        }

    }


    // Once you uknow the number of matches allocate an array as follows:
    // *matches = new CMatches [numMatches];
    //
    // The position of the corner point in iamge 1 is (m_X1, m_Y1)
    // The position of the corner point in image 2 is (m_X2, m_Y2)

    // Draw the matches
    DrawMatches(*matches, numMatches, image1Display, image2Display);
}


/*******************************************************************************
Project a point (x1, y1) using the homography transformation h
    (x1, y1) - input point
    (x2, y2) - returned point
    h - input homography used to project point
*******************************************************************************/
void MainWindow::Project(double x1, double y1, double &x2, double &y2, double h[3][3])
{
    double a = h[0][0]*x1 + h[0][1]*y1 + h[0][2];
    double b = h[1][0]*x1 + h[1][1]*y1 + h[1][2];
    double c = h[2][0]*x1 + h[2][1]*y1 + h[2][2];

    x2 = a / c;
    y2 = b / c;
}

/*******************************************************************************
Count the number of inliers given a homography.  This is a helper function for RANSAC.
    h - input homography used to project points (image1 -> image2
    matches - array of matching points
    numMatches - number of matchs in the array
    inlierThreshold - maximum distance between points that are considered to be inliers

    Returns the total number of inliers.
*******************************************************************************/
int MainWindow::ComputeInlierCount(double h[3][3], CMatches *matches, int numMatches, double inlierThreshold)
{
    int count = 0;

    for(int i = 0; i < numMatches; i++)
    {
        double x, y;

        Project(matches[i].m_X1, matches[i].m_Y1, x, y, h);

        double distance = (x - matches[i].m_X2)*(x - matches[i].m_X2) +
                          (y - matches[i].m_Y2)*(y - matches[i].m_Y2);

        if(distance < inlierThreshold*inlierThreshold)
            count++;

    }

    return count;
}

/***************************************************************************/
//helper function that performs same task as computeInliers
//but it also stores them in the provided data stucture
int MainWindow::StoreInliers(double h[3][3], CMatches *matches, CMatches *inlierPoints, int numMatches, double inlierThreshold)
{
    int ct = 0;

       for(int i=0;i<numMatches;i++)
       {
           double x1, y1;

           Project(matches[i].m_X1, matches[i].m_Y1, x1, y1, h);

           double distance = (x1 - matches[i].m_X2)*(x1 - matches[i].m_X2) +
                             (y1 - matches[i].m_Y2)*(y1 - matches[i].m_Y2);

           if(distance < inlierThreshold*inlierThreshold)
           {
               inlierPoints[ct].m_X1 = matches[i].m_X1;
               inlierPoints[ct].m_Y1 = matches[i].m_Y1;
               inlierPoints[ct].m_X2 = matches[i].m_X2;
               inlierPoints[ct].m_Y2 = matches[i].m_Y2;

               ct++;
           }

       }
       return ct;
}


/*******************************************************************************
Compute homography transformation between images using RANSAC.
    matches - set of matching points between images
    numMatches - number of matching points
    numIterations - number of iterations to run RANSAC
    inlierThreshold - maximum distance between points that are considered to be inliers
    hom - returned homography transformation (image1 -> image2)
    homInv - returned inverse homography transformation (image2 -> image1)
    image1Display - image used to display matches
    image2Display - image used to display matches
*******************************************************************************/
void MainWindow::RANSAC(CMatches *matches, int numMatches, int numIterations, double inlierThreshold,
                        double hom[3][3], double homInv[3][3], QImage &image1Display, QImage &image2Display)
{
    int maxInliers = 0;
    CMatches matchPairs[4];
    double currHom[3][3];

    //check 4 pairs of random points
    //find max inlier count
    for (int i = 0; i < numIterations; i++) {
        for (int j = 0; j < 4; j++) {
            int randIndex = rand()%numMatches;
            matchPairs[j].m_X1 = matches[randIndex].m_X1;
            matchPairs[j].m_Y1 = matches[randIndex].m_Y1;
            matchPairs[j].m_X2 = matches[randIndex].m_X2;
            matchPairs[j].m_Y2 = matches[randIndex].m_Y2;
        }

        bool randHom = ComputeHomography(matchPairs, 4, currHom, true);


        if (randHom == true) {
            int inlierCount = ComputeInlierCount(currHom, matches, numMatches, inlierThreshold);
            if (inlierCount > maxInliers) {
                maxInliers = inlierCount;

                for (int k = 0; k < 3; k++) {
                    for (int l = 0; l < 3; l++) {
                        hom[k][l] = currHom[k][l];
                    }
                }
            }
        }
    }

    CMatches *inlierPoints = new CMatches[maxInliers];
    int numInliers = 0;

    //recompute inliers and homography of all ponits
    numInliers = StoreInliers(hom, matches, inlierPoints, numMatches, inlierThreshold);

    ComputeHomography(inlierPoints, numInliers, hom, true);
    ComputeHomography(inlierPoints, numInliers, homInv, false);



    // After you're done computing the inliers, display the corresponding matches.
    DrawMatches(inlierPoints, numInliers, image1Display, image2Display);

    delete []inlierPoints;

}

/*******************************************************************************
Bilinearly interpolate image (helper function for Stitch)
    image - input image
    (x, y) - location to interpolate
    rgb - returned color values

    You can just copy code from previous assignment.
*******************************************************************************/
bool MainWindow::BilinearInterpolation(QImage *image, double x, double y, double rgb[3])
{
    int x1 = (int)floor(x);
    int y1 = (int)floor(y);
    int x2 = (int)ceil(x+0.000001);
    int y2 = (int)ceil(y+0.000001);
    int imageWidth = image->width();
    int imageHeight = image->height();
    QRgb pixel;

    double rgb11[3] = {0.0, 0.0, 0.0};
    double rgb12[3] = {0.0, 0.0, 0.0};
    double rgb21[3] = {0.0, 0.0, 0.0};
    double rgb22[3] = {0.0, 0.0, 0.0};


    if(0 <= x1 && x1 < imageWidth && 0 <= y1 && y1 < imageHeight) {
        pixel = image->pixel(x1, y1);
        rgb11[0] = (double) qRed(pixel);
        rgb11[1] = (double) qGreen(pixel);
        rgb11[2] = (double) qBlue(pixel);
    } else {
        return false;
    }

    if(0 <= x2 && x2 < imageWidth && 0 <= y1 && y1 < imageHeight) {
        pixel = image->pixel(x2,y1);
        rgb21[0] = (double) qRed(pixel);
        rgb21[1] = (double) qGreen(pixel);
        rgb21[2] = (double) qBlue(pixel);
    } else {
        return false;
    }

    if(0 <= x1 && x1 < imageWidth && 0 <= y2 && y2 < imageHeight) {
        pixel = image->pixel(x1,y2);
        rgb12[0] = (double) qRed(pixel);
        rgb12[1] = (double) qGreen(pixel);
        rgb12[2] = (double) qBlue(pixel);
    } else {
        return false;
    }

    if(0 <= x2 && x2 < imageWidth && 0 <= y2 && y2 < imageHeight) {
        pixel = image->pixel(x2,y2);
        rgb22[0] = (double) qRed(pixel);
        rgb22[1] = (double) qGreen(pixel);
        rgb22[2] = (double) qBlue(pixel);
    } else {
        return false;
    }

    double denom = (1 / ((x2-x1)*(y2-y1)));

    for (int i = 0; i < 3; i++) {
        rgb[i] = denom * ((rgb11[i]*(x2-x)*(y2-y)) +
                          (rgb21[i]*(x-x1)*(y2-y)) +
                          (rgb12[i]*(x2-x)*(y-y1)) +
                          (rgb22[i]*(x-x1)*(y-y1)));
    }

    return true;
}

/***********************************
 * helper function to round pixel positions */
void MainWindow::RoundPoints(double &x, double &y) {
    x = floor(x+0.5);
    y = floor(y+0.5);
}

/**********************************
 * helper function to calculate stitched image dimensions */
void MainWindow::calcDim(int cornerX, int cornerY, int &column1, int &row1, int &column2, int &row2) {
    if (column1 > cornerX) {
        column1 = cornerX;
    }
    if (column2 < cornerX) {
        column2 = cornerX;
    }
    if (row1 > cornerY) {
        row1 = cornerY;
    }
    if (row2 < cornerY) {
        row2 = cornerY;
    }
}


/*******************************************************************************
Stitch together two images using the homography transformation
    image1 - first input image
    image2 - second input image
    hom - homography transformation (image1 -> image2)
    homInv - inverse homography transformation (image2 -> image1)
    stitchedImage - returned stitched image
*******************************************************************************/
void MainWindow::Stitch(QImage image1, QImage image2, double hom[3][3], double homInv[3][3], QImage &stitchedImage)
{
    // Width and height of stitchedImage
    int ws = 0;
    int hs = 0;

    //grab width/ height data
    int w1 = image1.width();
    int h1 = image1.height();
    int w2 = image2.width();
    int h2 = image2.height();

    double topLeft_x;
    double topLeft_y;
    double topRight_x;
    double topRight_y;
    double bottomLeft_x;
    double bottomLeft_y;
    double bottomRight_x;
    double bottomRight_y;

    //project corners of image 2 onto image 1
    Project(0,0,topLeft_x, topLeft_y, homInv);
    Project(w2-1,0,topRight_x, topRight_y, homInv);
    Project(0,h2-1,bottomLeft_x, bottomLeft_y, homInv);
    Project(w2-1,h2-1,bottomRight_x, topRight_y, homInv);

    MainWindow::RoundPoints(topLeft_x, topLeft_y);
    MainWindow::RoundPoints(topRight_x, topRight_y);
    MainWindow::RoundPoints(bottomLeft_x, bottomLeft_y);
    MainWindow::RoundPoints(bottomRight_x, topRight_y);


    int column1 = 0;
    int row1 = 0;
    int column2 = w1;
    int row2 = h1;

    MainWindow::calcDim((int)topLeft_x, (int)topLeft_y, column1, row1, column2, row2);
    MainWindow::calcDim((int)topRight_x, (int)topRight_y, column1, row1, column2, row2);
    MainWindow::calcDim((int)bottomLeft_x,(int) bottomLeft_y, column1, row1, column2, row2);
    MainWindow::calcDim((int)bottomRight_x, (int)topRight_y, column1, row1, column2, row2);

    ws = column2 - column1;
    hs = row2 - row1;

   stitchedImage = QImage(ws, hs, QImage::Format_RGB32);
   stitchedImage.fill(qRgb(0,0,0));

   for(int r=0;r<h1;r++)
       for(int c=0;c<w1;c++)
       {
           QRgb pixel = image1.pixel(c, r);
           stitchedImage.setPixel(c - column1, r - row1, pixel);
       }

    for(int r=0;r<hs;r++)
        for(int c=0;c<ws;c++)
        {
            double x2Proj;
            double y2Proj;
            //double rgb[3];
            double currColumn = (double) c+column1;
            double currRow = (double) r+row1;

            Project(currColumn, currRow, x2Proj, y2Proj, hom);

            if (0 <= x2Proj && x2Proj < w2 && 0 <= y2Proj && y2Proj < h2)
            {
              double rgb[3];
              BilinearInterpolation(&image2, x2Proj, y2Proj, rgb);

              stitchedImage.setPixel(c, r, qRgb((int)(floor(rgb[0]+0.5)),
                                                (int)(floor(rgb[1]+0.5)),
                                                (int)(floor(rgb[2]+0.5))));
            }
        }
}

