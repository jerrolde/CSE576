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
    int x2 = ceil(x);
    int x1 = floor(x);
    int y1 = ceil(y);
    int y2 = floor(y);

    float dx1 = x - x1;
    float dx2 = x2 - x;
    float dy1 = y - y2;
    float dy2 = y1 - y;

    float areaUL = dx1 * dy1;
    float areaUR = dx2 * dy1;
    float areaLL = dx1 * dy2;
    float areaLR = dx2 * dy2;

    float UL = clamped_pixel(x1, y2, c);
    float UR = clamped_pixel(x2, y2, c);
    float LL = clamped_pixel(x1, y1, c);
    float LR = clamped_pixel(x2, y1, c);

    if (x2 == x1 && y1 == y2) {
      return UL;
    } else if (x2 == x1) {
      return dy1 * LL + dy2 * UL;
    } else if (y1 == y2) {
      return dx1 * UL + dx2 * UR;
    } else {
      return UL * areaLR + UR * areaLL + LL * areaUR + LR * areaUL;
    }

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


