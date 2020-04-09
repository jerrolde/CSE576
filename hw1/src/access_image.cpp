#include "image.h"

// HW0 #1
// const Image& im: input image
// int x,y: pixel coordinates
// int ch: channel of interest
// returns the 0-based location of the pixel value in the data array
int pixel_address(const Image& im, int x, int y, int ch)
  {
  int y_offset, ch_offset;

  // no offset corresponding to columns 
  y_offset = im.w; // each row is offset by width
  ch_offset = im.w * im.h; // each channel is offset by width * height
  return x + (y_offset * y) + (ch_offset * ch);
  }

// HW0 #1
// const Image& im: input image
// int x,y,ch: pixel coordinates and channel of interest
// returns the value of the clamped pixel at channel ch
float get_clamped_pixel(const Image& im, int x, int y, int ch)
  {
  int x_clamped, y_clamped, ch_clamped;
  if (x < 0) x_clamped = 0;
  else if (x > im.w) x_clamped = im.w - 1;
  else x_clamped = x;

  if (y < 0) y_clamped = 0;
  else if (y > im.h) y_clamped = im.h - 1;
  else y_clamped = y;
  
  if (ch < 0) ch_clamped = 0;
  else if (ch > im.c) ch_clamped = im.c - 1;
  else ch_clamped = ch;
    
  return im.data[pixel_address(im, x_clamped, y_clamped, ch_clamped)];
  }


// HW0 #1
// Image& im: input image
// int x,y,ch: pixel coordinates and channel of interest
void set_pixel(Image& im, int x, int y, int c, float value)
  {
  int x_valid, y_valid, ch_valid;
  if (x < 0) x_valid = 0;
  else if (x > im.w) x_valid = 0;
  else x_valid = 1;

  if (y < 0) y_valid = 0;
  else if (y > im.h) y_valid = 0;
  else y_valid = 1;
  
  if (c < 0) ch_valid = 0;
  else if (c > im.c) ch_valid = 0;
  else ch_valid = 1;

  if (x_valid & y_valid & ch_valid) im.data[pixel_address(im, x, y, c)] = value;
  }



// HW0 #2
// Copies an image
// Image& to: destination image
// const Image& from: source image
void copy_image(Image& to, const Image& from)
  {
  // allocating data for the new image
  to.data=(float*)calloc(from.w*from.h*from.c,sizeof(float));
  to.c=from.c;

  to.w = from.w;
  to.h = from.h;
  to.data = (float *) memcpy((void *)to.data, (void *)from.data, from.w*from.h*from.c*sizeof(float));  
  }
