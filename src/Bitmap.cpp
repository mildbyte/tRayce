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

Bitmap::~Bitmap() {
    for (register int i = width_ - 1; i >= 0; i--) {delete[] bits_[i];}
    delete[] bits_;
}

void Bitmap::setPixel(int x, int y, Vector pixel, double depth)
{
    //Check if the address is out of bounds
    if ((x > width_) || (x < 0) || (y > height_) || (y < 0)) return;

    bits_[x][y].color = pixel;
    bits_[x][y].depth = depth;
}

void Bitmap::foreach(BitmapPixel (*callback)(BitmapPixel)) {
    for (int i = 0; i < width_; i++) {
        for (int j = 0; j < height_; j++) {
            bits_[i][j] = callback(bits_[i][j]);
        }
    }
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
    FILE* output;
	
	fopen_s(&output, filename, "wb");

    fwrite(&bfheader, 1, sizeof(bfheader), output);
    fwrite(&biheader, 1, sizeof(biheader), output);

    //Output the bitmap

    //BMP requires every row to be padded to 4 bytes
    int padding = 4 - ((width_ * 3) % 4);
    if (padding == 4) padding = 0; //side effect if width mod 4 = 0 :)

    for (int i = height_ - 1; i >= 0; i--) {
        for (int j = 0; j < width_; j++) {
            //Normalize and clamp the colors
            bits_[j][i].color.reinhardMap();
            
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
