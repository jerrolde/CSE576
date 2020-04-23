#include <cmath>
#include "image.h"

using namespace std;

// HW1 #1
// float x,y: inexact coordinates
// int c: channel
// returns the nearest neighbor to pixel (x,y,c)
float Image::pixel_nearest(float x, float y, int c) const
  {
  // Since you are inside class Image you can
  // use the member function pixel(a,b,c)

  int x_int, y_int;
  x_int = (int) lroundf(x);
  y_int = (int) lroundf(y);
  
  return clamped_pixel(x_int, y_int, c);
  }

// HW1 #1
// float x,y: inexact coordinates
// int c: channel
// returns the bilinearly interpolated pixel (x,y,c)
float Image::pixel_bilinear(float x, float y, int c) const
  {
  // Since you are inside class Image you can
  // use the member function pixel(a,b,c)
  
  // referencing http://supercomputingblog.com/graphics/coding-bilinear-interpolation/
  int x1, x2, y1, y2;
  float q12, q22, q11, q21;
  float r2, r1, p;

  x1 = floor(x);
  x2 = ceil(x);
  y1 = ceil(y);
  y2 = floor(y);

  q12 = clamped_pixel(x1, y2, c);
  q22 = clamped_pixel(x2, y2, c);
  q11 = clamped_pixel(x1, y1, c);
  q21 = clamped_pixel(x2, y1, c);
  r1 = ((x2 - x) / (x2 - x1)) * q11 + ((x - x1) / (x2 - x1)) * q21;
  r2 = ((x2 - x) / (x2 - x1)) * q12 + ((x - x1) / (x2 - x1)) * q22;
  p = ((y2 - y) / (y2 - y1)) * r1 + ((y - y1) / (y2 - y1)) * r2;
  
  return p;
  }

// HW1 #1
// int w,h: size of new image
// const Image& im: input image
// return new Image of size (w,h,im.c)
Image nearest_resize(const Image& im, int w, int h)
  {
  Image ret(w,h,im.c);
  
  float ratio_x, ratio_y;
  ratio_x = (float) im.w / (float) w;
  ratio_y = (float) im.h / (float) h;

  int x, y, c;
  float x_new, y_new, pixel;

  for (x = 0; x < w; x++)
  {
    for (y = 0; y < h; y++)
    {
      for (c = 0; c < im.c; c++)
      {
        x_new = -0.5 + ratio_x * (x + 0.5);
        y_new = -0.5 + ratio_y * (y + 0.5);
        pixel = im.pixel_nearest(x_new, y_new, c);
        ret.set_pixel(x, y, c, pixel);
      }
    }
  }
  
  return ret;
  }


// HW1 #1
// int w,h: size of new image
// const Image& im: input image
// return new Image of size (w,h,im.c)
Image bilinear_resize(const Image& im, int w, int h)
  {
  Image ret(w,h,im.c);

  float ratio_x, ratio_y;
  ratio_x = (float) im.w / (float) w;
  ratio_y = (float) im.h / (float) h;

  printf("im.w = %i, im.h = %i, w = %i, h = %i, ratio_x = %f, ratio_y = %f", im.w, im.h, w, h, ratio_x, ratio_y);
 
  int x, y, c;
  float x_new, y_new, pixel;

  for (x = 0; x < w; x++)
  {
    for (y = 0; y < h; y++)
    {
      for (c = 0; c < im.c; c++)
      {
        x_new = -0.5 + ratio_x * (x + 0.5);
        y_new = -0.5 + ratio_y * (y + 0.5);
        pixel = im.pixel_bilinear(x_new, y_new, c);
        ret.set_pixel(x, y, c, pixel);
      }
    }
  }
  
  return ret;
  }


