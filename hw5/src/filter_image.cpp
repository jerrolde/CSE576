#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"

#define M_PI 3.14159265358979323846

// HW1 #2.1
// Image& im: image to L1-normalize
void l1_normalize(Image &im)
{

  double norm;
  float new_pixel;
  int x, y, c;

  for (c = 0; c < im.c; c++)
  {
    norm = 0;
    for (x = 0; x < im.w; x++)
    {
      for (y = 0; y < im.h; y++)
      {
        norm += im.clamped_pixel(x, y, c);
      }
    }

    for (x = 0; x < im.w; x++)
    {
      for (y = 0; y < im.h; y++)
      {
        new_pixel = im.clamped_pixel(x, y, c) / norm;
        im.set_pixel(x, y, c, new_pixel);
      }
    }
  }
}

// HW1 #2.1
// int w: size of filter
// returns the filter Image of size WxW
Image make_box_filter(int w)
{
  assert(w % 2); // w needs to be odd

  Image im(w, w, 1);

  int x, y;
  float val;

  val = 1.0 / (w * w);

  for (int x = 0; x < w; x++)
  {
    for (int y = 0; y < w; y++)
    {
      im.set_pixel(x, y, 0, val);
    }
  }

  return im;
}

// HW1 #2.2
// const Image&im: input image
// const Image& filter: filter to convolve with
// bool preserve: whether to preserve number of channels
// returns the convolved image
Image convolve_image(const Image &im, const Image &filter, bool preserve)
{
  assert(filter.c == 1);
  Image ret;
  // This is the case when we need to use the function clamped_pixel(x,y,c).
  // Otherwise you'll have to manually check whether the filter goes out of bounds

  // TODO: Make sure you set the sizes of ret properly. Use ret=Image(w,h,c) to reset ret
  // TODO: Do the convolution operator

  ret = Image(im.w, im.h, im.c);

  int x, y, c, i, j;
  float old_pixel, new_pixel;
  for (c = 0; c < im.c; c++)
  {
    for (x = 0; x < im.w; x++)
    {
      for (y = 0; y < im.h; y++)
      {
        for (i = 0 - (filter.w / 2); i <= (filter.w / 2); i++)
        {
          for (j = 0 - (filter.h / 2); j <= (filter.h / 2); j++)
          {
            old_pixel = im.clamped_pixel(x + i, y + j, c);
            new_pixel = old_pixel * filter.clamped_pixel(i + filter.w / 2, j + filter.h / 2, 0);
            ret.set_pixel(x, y, c, ret.clamped_pixel(x, y, c) + new_pixel);
          }
        }
      }
    }
  }

  if (!preserve)
  {
    Image new_ret;
    new_ret = Image(ret.w, ret.h, 1);
    for (c = 0; c < ret.c; c++)
    {
      for (x = 0; x < ret.w; x++)
      {
        for (y = 0; y < ret.h; y++)
        {
          new_ret.set_pixel(x, y, 0, new_ret.clamped_pixel(x, y, 0) + ret.clamped_pixel(x, y, c));
        }
      }
    }

    return new_ret;
  }

  return ret;
}

// HW1 #2.3
// returns basic 3x3 high-pass filter
Image make_highpass_filter()
{

  Image filter;

  filter = Image(3, 3, 1);
  filter.set_pixel(1, 0, 0, -1.0);
  filter.set_pixel(0, 1, 0, -1.0);
  filter.set_pixel(2, 1, 0, -1.0);
  filter.set_pixel(1, 2, 0, -1.0);
  filter.set_pixel(1, 1, 0, 4.0);

  return filter;
}

// HW1 #2.3
// returns basic 3x3 sharpen filter
Image make_sharpen_filter()
{

  Image filter;

  filter = Image(3, 3, 1);
  filter.set_pixel(1, 0, 0, -1.0);
  filter.set_pixel(0, 1, 0, -1.0);
  filter.set_pixel(2, 1, 0, -1.0);
  filter.set_pixel(1, 2, 0, -1.0);
  filter.set_pixel(1, 1, 0, 5.0);

  return filter;
}

// HW1 #2.3
// returns basic 3x3 emboss filter
Image make_emboss_filter()
{
  Image filter;

  filter = Image(3, 3, 1);
  filter.set_pixel(0, 0, 0, -2.0);
  filter.set_pixel(1, 0, 0, -1.0);
  filter.set_pixel(0, 1, 0, -1.0);
  filter.set_pixel(2, 1, 0, 1.0);
  filter.set_pixel(1, 2, 0, 1.0);
  filter.set_pixel(1, 1, 0, 1.0);
  filter.set_pixel(2, 2, 0, 2.0);

  return filter;
}

// HW1 #2.4
// float sigma: sigma for the gaussian filter
// returns basic gaussian filter
Image make_gaussian_filter(float sigma)
{

  Image ret;
  int w;
  w = 6 * sigma;
  if (!(w % 2))
    w++;

  ret = Image(w, w, 1);

  int x, y, x2, y2;
  float pixel, exponent;

  for (x = 0; x < w; x++)
  {
    x2 = x - w / 2;
    for (y = 0; y < w; y++)
    {
      y2 = y - w / 2;
      exponent = -1.0 * ((pow(x2, 2) + pow(y2, 2)) / (2.0 * pow(sigma, 2)));
      pixel = (1.0 / (2.0 * M_PI * pow(sigma, 2))) * exp(exponent);
      ret.set_pixel(x, y, 0, pixel);
    }
  }

  l1_normalize(ret);

  return ret;
}

// HW1 #3
// const Image& a: input image
// const Image& b: input image
// returns their sum
Image add_image(const Image &a, const Image &b)
{
  assert(a.w == b.w && a.h == b.h && a.c == b.c); // assure images are the same size

  Image ret;
  ret = Image(a.w, a.h, a.c);

  float dat_a, dat_b;
  int x, y, c;
  for (c = 0; c < a.c; c++)
  {
    for (x = 0; x < a.w; x++)
    {
      for (y = 0; y < a.h; y++)
      {
        dat_a = a.clamped_pixel(x, y, c);
        dat_b = b.clamped_pixel(x, y, c);
        ret.set_pixel(x, y, c, dat_a + dat_b);
      }
    }
  }

  return ret;
}

// HW1 #3
// const Image& a: input image
// const Image& b: input image
// returns their difference res=a-b
Image sub_image(const Image &a, const Image &b)
{
  assert(a.w == b.w && a.h == b.h && a.c == b.c); // assure images are the same size

  Image ret;
  ret = Image(a.w, a.h, a.c);

  float dat_a, dat_b;
  int x, y, c;
  for (c = 0; c < a.c; c++)
  {
    for (x = 0; x < a.w; x++)
    {
      for (y = 0; y < a.h; y++)
      {
        dat_a = a.clamped_pixel(x, y, c);
        dat_b = b.clamped_pixel(x, y, c);
        ret.set_pixel(x, y, c, dat_a - dat_b);
      }
    }
  }

  return ret;
}

// HW1 #4.1
// returns basic GX filter
Image make_gx_filter()
{
  Image filter;

  filter = Image(3, 3, 1);
  filter.set_pixel(0, 0, 0, -1.0);
  filter.set_pixel(2, 0, 0, 1.0);
  filter.set_pixel(0, 1, 0, -2.0);
  filter.set_pixel(2, 1, 0, 2.0);
  filter.set_pixel(0, 2, 0, -1.0);
  filter.set_pixel(2, 2, 0, 1.0);

  return filter;
}

// HW1 #4.1
// returns basic GY filter
Image make_gy_filter()
{
  Image filter;

  filter = Image(3, 3, 1);
  filter.set_pixel(0, 0, 0, -1.0);
  filter.set_pixel(1, 0, 0, -2.0);
  filter.set_pixel(2, 0, 0, -1.0);
  filter.set_pixel(0, 2, 0, 1.0);
  filter.set_pixel(1, 2, 0, 2.0);
  filter.set_pixel(2, 2, 0, 1.0);

  return filter;
}

// HW1 #4.2
// Image& im: input image
void feature_normalize(Image &im)
{
  assert(im.w * im.h); // assure we have non-empty image

  int c, x, y;
  float max, min, range, pixel;

  min = __FLT_MAX__;
  max = __FLT_MIN__;
  for (c = 0; c < im.c; c++)
  {
    for (x = 0; x < im.w; x++)
    {
      for (y = 0; y < im.h; y++)
      {
        pixel = im.clamped_pixel(x, y, c);
        if (pixel > max)
          max = pixel;
        if (pixel < min)
          min = pixel;
      }
    }
  }

  range = max - min;
  if (range == 0)
  {
    im.clear();
  }
  else
  {
    for (c = 0; c < im.c; c++)
    {
      for (x = 0; x < im.w; x++)
      {
        for (y = 0; y < im.h; y++)
        {
          pixel = im.clamped_pixel(x, y, c);
          im.set_pixel(x, y, c, (pixel - min) / range);
        }
      }
    }
  }
}

// Normalizes features across all channels
void feature_normalize_total(Image &im)
{
  assert(im.w * im.h * im.c); // assure we have non-empty image

  int nc = im.c;
  im.c = 1;
  im.w *= nc;

  feature_normalize(im);

  im.w /= nc;
  im.c = nc;
}

// HW1 #4.3
// Image& im: input image
// return a pair of images of the same size
pair<Image, Image> sobel_image(const Image &im)
{
  Image gx, gy, g, theta;

  gx = convolve_image(im, make_gx_filter(), false);
  gy = convolve_image(im, make_gy_filter(), false);

  g = Image(gx.w, gx.h, gx.c);
  theta = Image(gx.w, gx.h, gx.c);

  int c, x, y;
  float p_gx, p_gy, mag, ang;
  for (c = 0; c < gx.c; c++)
  {
    for (x = 0; x < gx.w; x++)
    {
      for (y = 0; y < gx.h; y++)
      {
        p_gx = gx.clamped_pixel(x, y, c);
        p_gy = gy.clamped_pixel(x, y, c);
        mag = pow(pow(p_gx, 2) + pow(p_gy, 2), 0.5);
        ang = atan2f(p_gy, p_gx);

        g.set_pixel(x, y, c, mag);
        theta.set_pixel(x, y, c, ang);
      }
    }
  }

  return {g, theta};
}

// HW1 #4.4
// const Image& im: input image
// returns the colorized Sobel image of the same size
Image colorize_sobel(const Image &im)
{

  Image f = make_gaussian_filter(4);
  Image blur = convolve_image(im, f, true);
  blur.clamp();

  pair<Image,Image> res = sobel_image(blur);
  Image mag = res.first;
  Image theta = res.second;

  feature_normalize(mag);
  int c, x, y;
  float pixel;
  for (c = 0; c < theta.c; c++)
  {
    for (x = 0; x < theta.w; x++)
    {
      for (y = 0; y < theta.h; y++)
      {
        pixel = theta.clamped_pixel(x, y, c);
        theta.set_pixel(x, y, c, (pixel / (2 * M_PI)) + 0.5);
      }
    }
  }

  Image ret;
  ret = Image(im.w, im.h, 3);
  ret.set_channel(0, theta);
  ret.set_channel(1, mag);
  ret.set_channel(2, mag);

  hsv_to_rgb(ret);
  

  return ret;
}

// HW1 #4.5
// const Image& im: input image
// float sigma1,sigma2: the two sigmas for bilateral filter
// returns the result of applying bilateral filtering to im
Image bilateral_filter(const Image &im, float sigma1, float sigma2)
{
  Image bf = im;

  // TODO: Your bilateral code
  NOT_IMPLEMENTED();

  return bf;
}

// HELPER MEMBER FXNS

void Image::feature_normalize(void) { ::feature_normalize(*this); }
void Image::feature_normalize_total(void) { ::feature_normalize_total(*this); }
void Image::l1_normalize(void) { ::l1_normalize(*this); }

Image operator-(const Image &a, const Image &b) { return sub_image(a, b); }
Image operator+(const Image &a, const Image &b) { return add_image(a, b); }
