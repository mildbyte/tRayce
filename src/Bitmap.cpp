#include "Bitmap.h"

unsigned char clamp (double value) {
    unsigned int result = (unsigned int)(255 * value);
    if (result > 255)
        result = 255;

    return (unsigned char)result;
}

Bitmap::Bitmap(int width, int height) {
    //Set the variables, allocate memory space.
    width_ = width;
    height_ = height;

    bits_ = new BitmapPixel* [width + 1];
    for (register int i = 0; i < width; i++) {bits_[i] = new BitmapPixel[height+1];}
}

//Loads a BMP image from file
Bitmap::Bitmap(char* filepath) {
	//Read the header

#pragma warning(disable: 4996)
	FILE* input = fopen(filepath, "rb");
#pragma warning(default: 4996)

	unsigned char bfheader[14];
	unsigned char biheader[40];

	fread(&bfheader, 1, sizeof(bfheader), input);
	fread(&biheader, 1, sizeof(biheader), input);

	width_ = 0;
	height_ = 0;
	
	for (int i = 7; i >= 4; i--) width_ = (width_ << 8) + biheader[i];
	for (int i = 11; i >= 8; i--) height_ = (height_ << 8) + biheader[i];

	bits_ = new BitmapPixel*[width_ + 1];
	for (register int i = 0; i < width_; i++) { bits_[i] = new BitmapPixel[height_ + 1]; }

	int rowSize = (width_ * 3);
	unsigned char* row = new unsigned char[rowSize];

	for (int i = height_-1; i >= 0; i--) {
		fread(row, 1, rowSize, input);
		for (int j = 0; j < width_; j++) {
			Vector pix;
			for (int k = 0; k < 3; k++) pix[2 - k] = (row[j * 3 + k] / 255.0);
			setPixel(j, i, pix, 0);
		}
	}

	delete(row);
}

Bitmap::~Bitmap() {
    for (register int i = width_ - 1; i >= 0; i--) {delete[] bits_[i];}
    delete[] bits_;
}

void Bitmap::setPixel(int x, int y, Vector pixel, double depth)
{
    //Check if the address is out of bounds
    if ((x >= width_) || (x < 0) || (y >= height_) || (y < 0)) return;

    bits_[x][y].color = pixel;
    bits_[x][y].depth = depth;
}

Vector Bitmap::getPixel(int x, int y) {
	if ((x >= width_) || (x < 0) || (y >= height_) || (y < 0)) {
		return Vector(0, 0, 0);
	}

	return bits_[x][y].color;
}



void Bitmap::foreach(BitmapPixel (*callback)(BitmapPixel)) {
    for (int i = 0; i < width_; i++) {
        for (int j = 0; j < height_; j++) {
            bits_[i][j] = callback(bits_[i][j]);
        }
    }
}

void Bitmap::reinhardMap() {
	for (int i = 0; i < width_; i++)
		for (int j = 0; j < height_; j++)
		bits_[i][j].color.reinhardMap();
}

void Bitmap::saveToFile(char* filename) {
    //Converts the color values to 0..255, clamps them
    //and saves the bitmap to file (overwrites)

    //Bitmap saving code inspired by
    //http://stackoverflow.com/questions/2654480/writing-bmp-image-in-pure-c-c-without-other-libraries

    int filesize = 3*width_*height_;
    unsigned char bfheader [14] = {'B','M',0,0,0,0,0,0,0,0,54,0,0,0};
    unsigned char biheader [40] = {40,0,0,0,0,0,0,0,0,0,0,0,1,0,24,0};
    
    bfheader[ 2] = (unsigned char)(filesize    );
	bfheader[ 3] = (unsigned char)(filesize>> 8);
	bfheader[ 4] = (unsigned char)(filesize>>16);
	bfheader[ 5] = (unsigned char)(filesize>>24);

	biheader[ 4] = (unsigned char)(width_    );
	biheader[ 5] = (unsigned char)(width_>> 8);
	biheader[ 6] = (unsigned char)(width_>>16);
	biheader[ 7] = (unsigned char)(width_>>24);
	biheader[ 8] = (unsigned char)(height_    );
	biheader[ 9] = (unsigned char)(height_>> 8);
	biheader[10] = (unsigned char)(height_>>16);
	biheader[11] = (unsigned char)(height_>>24);
    

    //Open the output file and write the header
    FILE* output = fopen(filename, "wb");

    fwrite(&bfheader, 1, sizeof(bfheader), output);
    fwrite(&biheader, 1, sizeof(biheader), output);

    //Output the bitmap

    //BMP requires every row to be padded to 4 bytes
    int padding = 4 - ((width_ * 3) % 4);
    if (padding == 4) padding = 0; //side effect if width mod 4 = 0 :)

    for (int i = height_ - 1; i >= 0; i--) {
        for (int j = 0; j < width_; j++) {
            //Export the clamped and 0-255 mapped colours.
			//Assume they are already normalized to 0..1.
            
            //Write the B, G, R to the output
            unsigned char clamped;
            clamped = clamp(bits_[j][i].color.getZ());
            fwrite(&clamped, 1, 1, output);

            clamped = clamp(bits_[j][i].color.getY());
            fwrite(&clamped, 1, 1, output);

            clamped = clamp(bits_[j][i].color.getX());
            fwrite(&clamped, 1, 1, output);
        }
        
        //Pad the row to 4 bytes
        for (int p = 0; p < padding; p++) fputc(0, output);
    }

    fclose(output);
}
