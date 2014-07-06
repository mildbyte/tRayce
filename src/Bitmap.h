//A bitmap that stores image pixels as vectors
//with doubleing-point color values (0-infinity)
//and maps them to 0..255 when exporting.
#ifndef BITMAP_H
#define BITMAP_H

#include "Vector.h"
#include <cstdio>

struct BitmapPixel{
    Vector color;
    double depth;
};

class Bitmap {
private:
    int width_;
    int height_;
    BitmapPixel** bits_;
public:
    //Specifying width and height is necessary.
    Bitmap(int width, int height);
	Bitmap(char* filepath);
    ~Bitmap();

    //Getters, setters
    void setPixel(int x, int y, Vector pixel, double depth);
    Vector getPixel(int x, int y);
    double getDepth(int x, int y);
    int getWidth() {return width_;}
    int getHeight() {return height_;}

    void saveToFile(char* filename);

    void foreach(BitmapPixel(*callback)(BitmapPixel));
};

#endif
