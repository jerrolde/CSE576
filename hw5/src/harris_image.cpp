#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cassert>

#include "image.h"
//#include "matrix.h"

using namespace std;


// Create a feature descriptor for an index in an image.
// const Image& im: source image.
// int x,y: coordinates for the pixel we want to describe.
// returns: Descriptor for that index.
Descriptor describe_index(const Image& im, int x, int y, int w)
  {
  Descriptor d;
  d.p={(double)x,(double)y};
  d.data.reserve(w*w*im.c);
  
  // If you want you can experiment with other descriptors
  // This subtracts the central value from neighbors
  // to compensate some for exposure/lighting changes.
  for(int c=0;c<im.c;c++)
    {
    float cval = im.clamped_pixel(x,y,c);
    for(int dx=-w/2;dx<=w/2;dx++)for(int dy=-w/2;dy<=w/2;dy++)
      d.data.push_back(im.clamped_pixel(x+dx,y+dy,c)-cval);
    }
  return d;
  }

// Marks the spot of a point in an image.
// Image& im: image to mark.
// Point p: spot to mark in the image.
void mark_spot(Image& im, const Point& p)
  {
  int x = p.x;
  int y = p.y;
  
  for(int i = -9; i < 10; ++i)
    {
    im.set_pixel(x+i, y, 0, 1);
    im.set_pixel(x, y+i, 0, 1);
    im.set_pixel(x+i, y, 1, 0);
    im.set_pixel(x, y+i, 1, 0);
    im.set_pixel(x+i, y, 2, 1);
    im.set_pixel(x, y+i, 2, 1);
    }
  }

// Marks corners denoted by an array of descriptors.
// Image& im: image to mark.
// vector<Descriptor> d: corners in the image.
Image mark_corners(const Image& im, const vector<Descriptor>& d)
  {
  Image im2=im;
  for(auto&e1:d)mark_spot(im2,e1.p);
  return im2;
  }


// HW5 1.1b
// Creates a 1d Gaussian filter.
// float sigma: standard deviation of Gaussian.
// returns: single row Image of the filter.
Image make_1d_gaussian(float sigma)
  {
  
  Image ret;
  int w;
  w = 6 * sigma;
  if (!(w % 2))
    w++;

  ret = Image(w, 1, 1);

  float pixel;
  int x2;

  for (int i = 0; i < w; i++)
  {
      x2 = w / 2 - i;
      pixel = 1.0 / (sqrt(2 * M_PI) * (sigma));
      pixel = pixel * exp((-pow(x2, 2)) / (2 * powf(sigma, 2)));
      ret.set_pixel(i, 0, 0, pixel);
  }
  
  // Image lin(1,1); // set to proper dimension
  // lin.data[0]=1;
  
  
  return ret;
  }

// HW5 1.1b
// Smooths an image using separable Gaussian filter.
// const Image& im: image to smooth.
// float sigma: std dev. for Gaussian.
// returns: smoothed Image.
Image smooth_image(const Image& im, float sigma)
  {
  // TODO: use two convolutions with 1d gaussian filter.
  // Hint: to make the filter from vertical to horizontal or vice versa
  // use "swap(filter.h,filter.w)"
  
  Image filter = make_1d_gaussian(sigma);
  Image ret = convolve_image(im, filter, 1);
  filter.h = filter.w;
  filter.w = 1;
  ret = convolve_image(ret, filter, 1);
  return ret;
  }


// HW5 1.1
// Calculate the structure matrix of an image.
// const Image& im im: the input image.
// float sigma: std dev. to use for weighted sum.
// returns: structure matrix. 1st channel is Ix^2, 2nd channel is Iy^2,
//          third channel is IxIy.
Image structure_matrix(const Image& im2, float sigma)
  {
  TIME(1);
  // only grayscale or rgb
  assert((im2.c==1 || im2.c==3) && "only grayscale or rgb supported");
  Image im;
  // convert to grayscale if necessary
  if(im2.c==1)im=im2;
  else im=rgb_to_grayscale(im2);
  
  Image S(im.w, im.h, 3);

  Image Ix, Iy;
  Ix = convolve_image(im, make_gx_filter(), 0);
  Iy = convolve_image(im, make_gy_filter(), 0);

  for (int i = 0; i < im.w; i++)
  {
    for (int j = 0; j < im.h; j++)
    {
      S(i, j, 0) = Ix(i, j, 0) * Ix(i, j, 0);
      S(i, j, 1) = Iy(i, j, 0) * Iy(i, j, 0);
      S(i, j, 2) = Ix(i, j, 0) * Iy(i, j, 0);
    }
  }

  return smooth_image(S, sigma);
  }


// HW5 1.2
// Estimate the cornerness of each pixel given a structure matrix S.
// const Image& im S: structure matrix for an image.
// returns: a response map of cornerness calculations.
Image cornerness_response(const Image& S, int method)
  {
  Image R(S.w, S.h);
  // TODO: fill in R, "cornerness" for each pixel using the structure matrix.
  // We'll use formulation det(S) - alpha * trace(S)^2, alpha = .06.
  // E(S) = det(S) / trace(S)
  if (method == 1)
  {
    // method = 1 (exact 2nd eigenvalue) is optional
    NOT_IMPLEMENTED();
  }
  else
  {
    //method = 0 (det(S)/tr(S))

    float det, tr;
    for (int i = 0; i < S.w; i++)
    {
      for (int j = 0; j < S.h; j++)
      {
        det = (S(i, j, 0) * S(i, j, 1)) - (S(i, j, 2) * S(i, j, 2));
        tr = S(i, j, 0) + S(i, j, 1);
        R(i, j, 0) = det / tr;
        //R(i, j, 0) = det - 0.06 * tr * tr; // comment says I should use this?
      }
    }
  }
  
  return R;
  }


// HW5 1.3
// Perform non-max supression on an image of feature responses.
// const Image& im: 1-channel image of feature responses.
// int w: distance to look for larger responses.
// returns: Image with only local-maxima responses within w pixels.
Image nms_image(const Image& im, int w)
  {
  //TIME(1);
  Image r=im;
  // TODO: perform NMS on the response map.
  // for every pixel in the image:
  //     for neighbors within w:
  //         if neighbor response greater than pixel response:
  //             set response to be very low
  
  for (int i = 0; i < im.w; i++)
  {
    for (int j = 0; j < im.h; j++)
    {
      float val = im(i, j, 0);

      //check 2w+1 window around val
      for (int ii = -w; ii <= w; ii++)
      {
        for (int jj = -w; jj <= w; jj++)
        {
          // do nothing when out of bounds
          if ((i + ii >= 0) && (j + jj >= 0) && (i + ii < im.w) && (j + jj < im.h))
          {
            float neighbor = im(i + ii, j + jj, 0);
            if (neighbor > val) r(i, j, 0) = -9999;
          }
          
        }
      }
    }
  }
  
  return r;
  }


// HW5 1.4
// Perform corner detection and extract features from the corners.
// const Image& im: input image.
// const Image& nms: nms image
// float thresh: threshold for cornerness.
// returns: vector of descriptors of the corners in the image.
vector<Descriptor> detect_corners(const Image& im, const Image& nms, float thresh, int window)
  {
  vector<Descriptor> d;
  //TODO: count number of responses over threshold (corners)
  //TODO: and fill in vector<Descriptor> with descriptors of corners, use describe_index.
  
  for (int i = 0; i < im.w; i++)
  {
    for (int j = 0; j < im.h; j++)
    {
      if (nms(i, j, 0) >= thresh)
      {
        d.push_back(describe_index(im, i, j, window));
      }
    }
  }
  
  return d;
  }


// Perform harris corner detection and extract features from the corners.
// const Image& im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
// returns: vector of descriptors of the corners in the image.
vector<Descriptor> harris_corner_detector(const Image& im, float sigma, float thresh, int window, int nms, int corner_method)
  {
  // Calculate structure matrix
  Image S = structure_matrix(im, sigma);
  
  // Estimate cornerness
  Image R = cornerness_response(S,corner_method);
  
  // Run NMS on the responses
  Image Rnms = nms_image(R, nms);
  
  return detect_corners(im, Rnms, thresh, window);
  }

// Find and draw corners on an image.
// Image& im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
Image detect_and_draw_corners(const Image& im, float sigma, float thresh, int window, int nms, int corner_method)
  {
  vector<Descriptor> d = harris_corner_detector(im, sigma, thresh, window, nms, corner_method);
  return mark_corners(im, d);
  }
