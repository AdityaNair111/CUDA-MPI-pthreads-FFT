//
// Created by brian on 11/20/18.
//

#pragma once

class Complex;

class InputImage {
public:

    InputImage(const char* filename);
    int get_width() const;
    int get_height() const;

    //returns a pointer to the image data.  Note the return is a 1D
    //array which represents a 2D image.  The data for row 1 is
    //immediately following the data for row 0 in the 1D array
    Complex* get_image_data() const;

    //use this to save output from forward DFT
    void save_image_data(const char* filename, Complex* d, int w, int h);
    //use this to save output from reverse DFT
    void save_image_data_real(const char* filename, Complex* d, int w, int h);

private:
    int w;
    int h;
    Complex* data;
};
