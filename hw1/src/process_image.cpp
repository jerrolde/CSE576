#include <cstdio>
#include <cstring>
#include <cassert>
#include <cmath>

#include "image.h"

using namespace std;

// HW0 #3
// const Image& im: input image
// return the corresponding grayscale image
Image rgb_to_grayscale(const Image &im)
{
  assert(im.c == 3);         // only accept RGB images
  Image gray(im.w, im.h, 1); // create a new grayscale image (note: 1 channel)

  for (int i = 0; i < im.w; i++)
  {
    for (int j = 0; j < im.h; j++)
    {
      float r, g, b, y;
      r = get_clamped_pixel(im, i, j, 0);
      g = get_clamped_pixel(im, i, j, 1);
      b = get_clamped_pixel(im, i, j, 2);
      y = 0.299 * r + 0.587 * g + 0.114 * b;
      set_pixel(gray, i, j, 0, y);
    }
  }

  return gray;
}

// Example function that changes the color of a grayscale image
Image grayscale_to_rgb(const Image &im, float r, float g, float b)
{
  assert(im.c == 1);
  Image rgb(im.w, im.h, 3);

  for (int q2 = 0; q2 < im.h; q2++)
    for (int q1 = 0; q1 < im.w; q1++)
    {
      rgb(q1, q2, 0) = r * im(q1, q2);
      rgb(q1, q2, 1) = g * im(q1, q2);
      rgb(q1, q2, 2) = b * im(q1, q2);
    }

  return rgb;
}

// HW0 #4
// Image& im: input image to be modified in-place
// int c: which channel to shift
// float v: how much to shift
void shift_image(Image &im, int c, float v)
{
  assert(c >= 0 && c < im.c); // needs to be a valid channel

  for (int i = 0; i < im.w; i++)
  {
    for (int j = 0; j < im.h; j++)
    {
      float old_pixel, new_pixel;
      old_pixel = get_clamped_pixel(im, i, j, c);
      new_pixel = old_pixel + v;
      set_pixel(im, i, j, c, new_pixel);
    }
  }
}

// HW0 #8
// Image& im: input image to be modified in-place
// int c: which channel to scale
// float v: how much to scale
void scale_image(Image &im, int c, float v)
{
  assert(c >= 0 && c < im.c); // needs to be a valid channel

  // TODO: scale all the pixels at the specified channel

  NOT_IMPLEMENTED();
}

// HW0 #5
// Image& im: input image to be modified in-place
void clamp_image(Image &im)
{
  for (int i = 0; i < im.w; i++)
  {
    for (int j = 0; j < im.h; j++)
    {
      for (int k = 0; k < im.c; k++)
      {
        float old_pixel, new_pixel;
        old_pixel = get_clamped_pixel(im, i, j, k);
        if (old_pixel < 0)
          new_pixel = 0;
        else if (old_pixel > 1)
          new_pixel = 1;
        else
          new_pixel = old_pixel;

        set_pixel(im, i, j, k, new_pixel);
      }
    }
  }
}

// These might be handy
float max(float a, float b, float c)
{
  return max({a, b, c});
}

float min(float a, float b, float c)
{
  return min({a, b, c});
}

// HW0 #6
// Image& im: input image to be modified in-place
void rgb_to_hsv(Image &im)
{
  assert(im.c == 3 && "only works for 3-channels images");

  for (int i = 0; i < im.w; i++)
  {
    for (int j = 0; j < im.h; j++)
    {
      float r, g, b, h, s, v;
      float m, c, hp;
      r = get_clamped_pixel(im, i, j, 0);
      g = get_clamped_pixel(im, i, j, 1);
      b = get_clamped_pixel(im, i, j, 2);

      v = max(r, g, b);
      m = min(r, g, b);
      c = v - m;
      s = (v != 0) ? c / v : 0;

      if (c == 0)
        hp = 0; //undefined
      else if (v == r)
        hp = (g - b) / c;
      else if (v == g)
        hp = (b - r) / c + 2;
      else if (v == b)
        hp = (r - g) / c + 4;
      else
        hp = 0; //shouldn't ever happen

      h = (hp < 0) ? hp / 6 + 1 : hp / 6;

      set_pixel(im, i, j, 0, h);
      set_pixel(im, i, j, 1, s);
      set_pixel(im, i, j, 2, v);
    }
  }
}

// HW0 #7
// Image& im: input image to be modified in-place
void hsv_to_rgb(Image &im)
{
  assert(im.c == 3 && "only works for 3-channels images");

  for (int i = 0; i < im.w; i++)
  {
    for (int j = 0; j < im.h; j++)
    {
      float r, g, b, h, s, v;
      float c, x, m;
      h = get_clamped_pixel(im, i, j, 0);
      s = get_clamped_pixel(im, i, j, 1);
      v = get_clamped_pixel(im, i, j, 2);

      c = v * s;
      x = c * (1 - abs(fmod(6.0 * h, 2.0) - 1));
      m = v - c;

      if (h < (1.0 / 6.0))
      {
        r = c;
        g = x;
        b = 0;
      }
      else if (h < (2.0 / 6.0))
      {
        r = x;
        g = c;
        b = 0;
      }
      else if (h < (3.0 / 6.0))
      {
        r = 0;
        g = c;
        b = x;        
      }
      else if (h < (4.0 / 6.0))
      {
        r = 0;
        g = x;
        b = c;        
      }
      else if (h < (5.0 / 6.0))
      {
        r = x;
        g = 0;
        b = c;       
      }
      else
      {
        r = c;
        g = 0;
        b = x;        
      }

      r += m;
      g += m;
      b += m;
      

      set_pixel(im, i, j, 0, r);
      set_pixel(im, i, j, 1, g);
      set_pixel(im, i, j, 2, b);
    }
  }
}

// HW0 #9
// Image& im: input image to be modified in-place
void rgb_to_lch(Image &im)
{
  assert(im.c == 3 && "only works for 3-channels images");

  // TODO: Convert all pixels from RGB format to LCH format

  NOT_IMPLEMENTED();
}

// HW0 #9
// Image& im: input image to be modified in-place
void lch_to_rgb(Image &im)
{
  assert(im.c == 3 && "only works for 3-channels images");

  // TODO: Convert all pixels from LCH format to RGB format

  NOT_IMPLEMENTED();
}

// Implementation of member functions
void Image::clamp(void) { clamp_image(*this); }
void Image::shift(int c, float v) { shift_image(*this, c, v); }
void Image::scale(int c, float v) { scale_image(*this, c, v); }

void Image::HSVtoRGB(void) { hsv_to_rgb(*this); }
void Image::RGBtoHSV(void) { rgb_to_hsv(*this); }
void Image::LCHtoRGB(void) { lch_to_rgb(*this); }
void Image::RGBtoLCH(void) { rgb_to_lch(*this); }
