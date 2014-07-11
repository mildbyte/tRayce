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
	void reinhardMap();

    void foreach(BitmapPixel(*callback)(BitmapPixel));

	//Performs bilinear texture sampling, adapted from the Wikipedia article
	Vector textureSample(double u, double v) {
		u = u * width_ - 0.5;
		v = v * height_ - 0.5;
		int x = floor(u);
		int y = floor(v);
		double u_ratio = u - x;
		double v_ratio = v - y;
		double u_opposite = 1 - u_ratio;
		double v_opposite = 1 - v_ratio;

		//Wrap pixels outside of the range around the texture
		Vector result = (getPixel(x % width_, y % height_) * u_opposite + getPixel((x + 1) % width_,  y % height_) * u_ratio) * v_opposite +
			(getPixel(x % width_, (y + 1) % height_) * u_opposite + getPixel((x + 1) % width_, (y + 1) % height_) * u_ratio) * v_ratio;
		return result;
	}
};

#endif
