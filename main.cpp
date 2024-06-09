/*===============================================================
             Filtering in the Frequency-Domain
              (Lowpass and Highpass Filtering)
              for instructional purposes only
                    by Prof. Lan
                    Apr. 15, 2010

In Dev-C++, first create a project file.  Then add Filter_Freq_students2.cpp,
FFT1.c, FFT2.c, and bmp.cpp to this project.  Remember to put all these files
in the same directory.
===============================================================*/
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include "bmp.h"
#define PI 3.1415926

#include <tr1/unordered_map>
using namespace std::tr1;

using namespace std;
using std::cout;

void fft1(float data[], int nn, int isign);
void fft2(float data[], int nn, int isign);
void spectrum_shift(int mm);

int R[MaxBMPSizeX][MaxBMPSizeY]; // MaxBMPSizeX and MaxBMPSizeY are defined in "bmp.h"
int G[MaxBMPSizeX][MaxBMPSizeY];
int B[MaxBMPSizeX][MaxBMPSizeY];
int r[MaxBMPSizeX][MaxBMPSizeY];
int g[MaxBMPSizeX][MaxBMPSizeY];
int b[MaxBMPSizeX][MaxBMPSizeY];


// You can try various combinations of D and filt_order values.
const float D = 20; // the D parameter of a Butterworth filter
const int filt_order = 2; // order of the Butterworth filter

float data[MaxBMPSizeX * MaxBMPSizeY * 2]; // a long 1D array to keep the image and its spectrum
float sp_re[MaxBMPSizeX][MaxBMPSizeY]; // real part of the spectrum
float sp_im[MaxBMPSizeX][MaxBMPSizeY]; // imaginary part of the spectrum
float tmp[MaxBMPSizeX][MaxBMPSizeY];

int sample() {
    int width, height;
    int i, j;

    open_bmp("images/cameraman.bmp", R, R, R, width, height); // for gray images
    //open_bmp("lena.bmp", R, R, R, width, height); // for gray images
    printf("==========================================\n");

    /*----------------------------------------------------------------------------
          範例：（請將主程式放在此處）
    ----------------------------------------------------------------------------*/
    // convert to the long 1D array
    int ptr = 0; // pointer for the data array
    for (i = 0; i < height; i++) {
        for (j = 0; j < width; j++) {
            data[ptr] = R[j][height - i];     // real part of the input data
            data[ptr + 1] = 0;         // imaginary part of the input data
            ptr += 2;
        }
    }

    fft2(data, width, 1); // perform the forward fft2

    // fetch the 2D spectrum from the long 1-D spectrum array
    ptr = 0;
    for (i = 0; i < height; i++) {
        for (j = 0; j < width; j++) {
            sp_re[i][j] = data[ptr];    // real part of the 2D spectrum
            sp_im[i][j] = data[ptr + 1];  // imaginary part of the 2D spectrum
            ptr += 2;
        }
    }

    // shift the 2-D spectrum
    int mm = width / 2; // mm is one half of the original image width
    spectrum_shift(mm);

    // perform Ideal Filtering in this sample program
    // Note that you should write your own code to perform Butterworth Lowpass Filtering and Highpass
    // Filtering, or Gaussian LPF/HPF.
    // Define the H(u,v) function and perform F(u,v)H(u,v).

    float Dsq = D * D;
    float sq_dist;
    for (i = 0; i < height; i++) {
        for (j = 0; j < width; j++) {
            sq_dist = (i - mm) * (i - mm) + (j - mm) * (j - mm);
            //if (sq_dist > Dsq) {  // ideal lowpass filtering
            //    sp_re[i][j] = 0;    // Pay attention to the ringing effect when ideal filtering is applied.
            //    sp_im[i][j] = 0;
            //}
            if (sq_dist < Dsq) {  // ideal highpass filtering
                sp_re[i][j] = 0;
                sp_im[i][j] = 0;
            }
        }
    }

    // shift the 2-D spectrum back
    spectrum_shift(mm);

    // convert to the long 1D array
    ptr = 0; // pointer for the data array
    for (i = 0; i < height; i++) {
        for (j = 0; j < width; j++) {
            data[ptr] = sp_re[i][j];
            data[ptr + 1] = sp_im[i][j];
            ptr += 2;
        }
    }

    fft2(data, width, -1); // perform the inverse fft2

    // convert back to the 2D image
    ptr = 0; // pointer for the data array
    for (i = 0; i < height; i++) {
        for (j = 0; j < width; j++) {
            sp_re[i][j] = data[ptr];
            ptr += 2;
        }
    }

    // overflow and underflow handling (by linear scaling)
    float min = 1E99;
    float max = -1E99;
    for (i = 0; i < height; i++) {
        for (j = 0; j < width; j++) {
            if (sp_re[i][j] > max) max = sp_re[i][j];
            if (sp_re[i][j] < min) min = sp_re[i][j];
        }
    }
    float sf;
    sf = 255 / (max - min);
    printf("%f \t %f \t %f \n", max, min, sf);
    for (i = 0; i < height; i++) {
        for (j = 0; j < width; j++) {
            R[j][height - i] = (int)((sp_re[i][j] - min) * sf);
        }
    }
    //printf("%d \n", R[12][16]); // check the data value

    // 儲存處理結果至新的圖檔中
    //save_bmp("lena_new1.bmp", r, g, b); // for true color images
    save_bmp("new.bmp", R, R, R); // for gray images
    printf("Job Finished!\n");

    // 關閉 bmp 影像圖檔
    close_bmp();

    system("PAUSE"); /* so that the command window holds a while */
    return 0;
}

void spectrum_shift(int mm) {
    for (int i = 0; i < mm; i++) {
        for (int j = 0; j < mm; j++) {
            //------------------------ the real part
            tmp[i][j] = sp_re[i][j]; // upper-left <-> lower-right
            sp_re[i][j] = sp_re[mm + i][mm + j];
            sp_re[mm + i][mm + j] = tmp[i][j];
            tmp[i][j] = sp_re[i][mm + j]; // upper-right <-> lower-left
            sp_re[i][mm + j] = sp_re[mm + i][j];
            sp_re[mm + i][j] = tmp[i][j];
            //--------------------------- the imag. part
            tmp[i][j] = sp_im[i][j]; // upper-left <-> lower-right
            sp_im[i][j] = sp_im[mm + i][mm + j];
            sp_im[mm + i][mm + j] = tmp[i][j];
            tmp[i][j] = sp_im[i][mm + j]; // upper-right <-> lower-left
            sp_im[i][mm + j] = sp_im[mm + i][j];
            sp_im[mm + i][j] = tmp[i][j];
        }
    }
}

float** convolve(int** arr, float** mask, int arr_width, int arr_height, int mask_width, int mask_height) {
    int convolved_width = arr_width - mask_width + 1;
    int convolved_height = arr_height - mask_height + 1;
    // 創建二維數組，size為捲積後的大小
    float** convolved_arr = new float* [convolved_height];
    for (int i = 0; i < convolved_height; ++i) {
        convolved_arr[i] = new float[convolved_width];
    }
    // 初始化elements為0
    for (int i = 0; i < convolved_height; ++i) {
        for (int j = 0; j < convolved_width; ++j) {
            convolved_arr[i][j] = 0.0f;
        }
    }
    int _w = 0;
    int _h = 0;
    for (int i = mask_height / 2; i < arr_height - mask_height / 2; ++i) {
        for (int j = mask_width / 2; j < arr_width - mask_width / 2; ++j) {
            float sum = 0.0f;
            for (int u = 0; u < mask_height; ++u) {
                for (int v = 0; v < mask_width; ++v) {
                    int arr_x = i + u - mask_height / 2;
                    int arr_y = j + v - mask_width / 2;
                    sum += arr[arr_x][arr_y] * mask[u][v];
                }
            }
            if (sum < 0.0f) {
                sum = 0.0f;
            }
            else if (sum > 255.0f) {
                sum = 255.0f;
            }
            else {
                sum = round(sum);
            }
            convolved_arr[_h][_w] = sum;
            _w += 1;
        }
        _h += 1;
        _w = 0;
    }
    return convolved_arr;
}

int** zero_padding(int** arr, int width, int height, int mask_size) {
    int padding_size = mask_size / 2;
    int padded_width = width + 2 * padding_size;
    int padded_height = height + 2 * padding_size;

    // 創建二維數組，size為填充後的大小
    int** padded_arr = new int* [padded_height];
    for (int i = 0; i < padded_height; ++i) {
        padded_arr[i] = new int[padded_width];
    }
    // 填充0
    for (int i = 0; i < padded_height; ++i) {
        for (int j = 0; j < padded_width; ++j) {
            padded_arr[i][j] = 0;
        }
    }
    // 填充原始數據進去padded_arr
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            padded_arr[i + padding_size][j + padding_size] = arr[i][j];
        }
    }
    return padded_arr;
}

struct Line {
    double theta;
    double s;
    int count;
    int start_x;
    int start_y;
    int end_x;
    int end_y;
};

void line_info(int lines, vector<Line>& S_map) {
    if (S_map.empty()) {
        cout << "S_map is EMPTY!" << endl;
    }

    // 打印前 lines 個直線信息
    for (int i = 0; i < lines; ++i) {
        cout << "Theta: " << S_map[i].theta * 180 / PI << " S: " << S_map[i].s << " Count: " << S_map[i].count << endl;
        cout << "Start : (" << S_map[i].start_x << ", " << S_map[i].start_y << ")" << endl << endl;
    }
}

void line_detection(int lines) {
    int width, height;
    open_bmp("images/house.bmp", R, R, R, width, height);
    float sobel_x[3][3] = { {1,0,-1},{2,0,-2},{1,0,-1} };
    float** sobel_x_ptr = new float* [3];
    for (int i = 0; i < 3; ++i) {
        sobel_x_ptr[i] = sobel_x[i];
    }
    float sobel_y[3][3] = { {1,2,1},{0,0,0},{-1,-2,-1} };
    float** sobel_y_ptr = new float* [3];
    for (int i = 0; i < 3; ++i) {
        sobel_y_ptr[i] = sobel_y[i];
    }
    int** arr_ptr = new int* [height];
    for (int i = 0; i < height; ++i) {
        arr_ptr[i] = R[i];
    }
    int** padded_arr = zero_padding(arr_ptr, width, height, 3);
    float** Gx = convolve(padded_arr, sobel_x_ptr, width + 2, height + 2, 3, 3);
    float** Gy = convolve(padded_arr, sobel_y_ptr, width + 2, height + 2, 3, 3);
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            float value = sqrt(pow(Gx[i][j], 2) + pow(Gy[i][j], 2));
            if (value < 0.0f) {
                value = 0.0f;
            }
            else if (value > 255.0f) {
                value = 255.0f;
            }
            else {
                value = (int)round(value);
            }
            if (value > 128) {
                r[i][j] = value;
            }
            else {
                r[i][j] = 0;
            }
        }
    }
    save_bmp("output_images/line_detection_edges.bmp", r, r, r);
    vector<Line> S_map;
    //unordered_map<double, Line*> hash_map;
    vector<double> theta;
    for (int t = -90; t < 90; t += 5) {
        double _theta = -t * PI / 180.0f;
        theta.push_back(_theta);
    }
    double tol = 1e-2; //容許誤差值
    for (int i = 1; i < height; ++i) {
        for (int j = 1; j < width; ++j) {
            if (r[i][j] != 0) {
                for (double t : theta) {
                    double s = j * cos(t) + i * sin(t);
                    bool found = false;
                    bool exist = false;
                    for (Line& line : S_map) { //在S_map中找相同s, 若找到則set exist為true
                        if (line.s > 0.0 && s < 0.0) {
                            exist = false;
                        }
                        else if (line.s < 0.0 && s>0.0) {
                            exist = false;
                        }
                        else {
                            if (fabs(line.s - s) < tol) exist = true;
                        }

                        if (line.theta == t && exist) { //篩選相同theta, 若對應s值exist, 則該line票數+1
                            line.end_x = j;
                            line.end_y = i;
                            line.count++;
                            found = true;
                            break;
                        }

                        if (!found) S_map.push_back({ t, s, 1, j, i, j, i });//建立一個new line, initial count = 1
                    }
                }
            }
        }
    }
    // 冒泡排序sort
    for (int i = 0; i < S_map.size() - 1; ++i) {
        for (int j = 0; j < S_map.size() - i - 1; ++j) {
            if (S_map[j].count < S_map[j + 1].count) {
                swap(S_map[j], S_map[j + 1]);
            }
        }
    }
    int max_lines;
    if (lines <= S_map.size()) {
        max_lines = lines;
    }
    else {
        max_lines = S_map.size();
    }

    line_info(max_lines, S_map);

    for (int k = 0; k <= max_lines; ++k) {
        Line _line = S_map[k];
        double _s = _line.s;
        double sin_theta = sin(_line.theta);
        double cos_theta = cos(_line.theta);
        for (int x = 0; x < width; ++x) {
            double y = (_s - x * cos_theta) / sin_theta;
            int y_int = round(y);
            if (y_int >= 0 && y_int < height) {
                G[y_int][x] = 255;
            }
        }
    }


    save_bmp("output_images/output_with_lines.bmp", G, G, G);

}

float** butter_worth(int** arr, int arr_size, float cutoff_radius, int order) {
    float** res = new float* [arr_size];
    for (int i = 0; i < arr_size; ++i) {
        res[i] = new float[arr_size];
    }
    for (int u = 0; u < arr_size; ++u) {
        for (int v = 0; v < arr_size; ++v) {
            float distance = sqrt(pow((u - arr_size / 2), 2) + pow((v - arr_size / 2), 2));
            res[u][v] = 1 / (1 + pow((distance / cutoff_radius), 2 * order));
        }
    }
    return res;
}

void butterworth_lowpass(float cr, int order) {
    int height, width;
    open_bmp("images/cameraman.bmp", R, R, R, width, height);
    int** arr_ptr = new int* [height];
    for (int i = 0; i < height; ++i) {
        arr_ptr[i] = R[i];
    }
    float** res = new float* [height];
    for (int i = 0; i < height; ++i) {
        res[i] = new float[width];
    }
    int ptr = 0; // pointer for the data array
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            data[ptr] = R[j][height - i];
            data[ptr + 1] = 0;
            ptr += 2;
        }
    }
    fft2(data, width, 1);
    ptr = 0;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            sp_re[i][j] = data[ptr];    // real part of the 2D spectrum
            sp_im[i][j] = data[ptr + 1];  // imaginary part of the 2D spectrum
            ptr += 2;
        }
    }

    // shift the 2-D spectrum
    int mm = width / 2; // mm is one half of the original image width
    spectrum_shift(mm);
    float** mask = butter_worth(arr_ptr, width, cr, order);
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            float mask_coef = mask[i][j];
            sp_re[i][j] *= mask_coef;
            sp_im[i][j] *= mask_coef;
        }
    }

    spectrum_shift(mm);
    ptr = 0; // pointer for the data array
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            data[ptr] = sp_re[i][j];
            data[ptr + 1] = sp_im[i][j];
            ptr += 2;
        }
    }

    fft2(data, width, -1); // perform the inverse fft2

    // convert back to the 2D image
    ptr = 0; // pointer for the data array
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            sp_re[i][j] = data[ptr];
            ptr += 2;
        }
    }

    // overflow and underflow handling (by linear scaling)
    float min = 1E99;
    float max = -1E99;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if (sp_re[i][j] > max) max = sp_re[i][j];
            if (sp_re[i][j] < min) min = sp_re[i][j];
        }
    }
    float sf;
    sf = 255 / (max - min);
    //printf("%f \t %f \t %f \n", max, min, sf);
    for (int i = 0; i <= height; i++) {
        for (int j = 0; j < width; j++) {
            int value = (int)(sp_re[i][j]);
            if (value > 255) {
                value = 255;
            }
            else if (value < 0) {
                value = 0;
            }
            R[j][height - i] = value;
            //R[j][height - i] = (int)((sp_re[i][j] - min) * sf);
        }
    }

    save_bmp("output_images/output_butterworth_lowpass.bmp", R, R, R);
}

void butterworth_highpass(float cr, int order) {
    int height, width;
    open_bmp("images/cameraman.bmp", R, R, R, width, height);
    int** arr_ptr = new int* [height];
    for (int i = 0; i < height; ++i) {
        arr_ptr[i] = R[i];
    }
    float** res = new float* [height];
    for (int i = 0; i < height; ++i) {
        res[i] = new float[width];
    }
    int ptr = 0; // pointer for the data array
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            data[ptr] = R[j][height - i];
            data[ptr + 1] = 0;
            ptr += 2;
        }
    }
    fft2(data, width, 1);
    ptr = 0;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            sp_re[i][j] = data[ptr];    // real part of the 2D spectrum
            sp_im[i][j] = data[ptr + 1];  // imaginary part of the 2D spectrum
            ptr += 2;
        }
    }

    // shift the 2-D spectrum
    int mm = width / 2; // mm is one half of the original image width
    spectrum_shift(mm);
    float** mask = butter_worth(arr_ptr, width, cr, order);
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            float mask_coef = (1 - mask[i][j]);
            sp_re[i][j] *= mask_coef;
            sp_im[i][j] *= mask_coef;
        }
    }

    spectrum_shift(mm);
    ptr = 0; // pointer for the data array
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            data[ptr] = sp_re[i][j];
            data[ptr + 1] = sp_im[i][j];
            ptr += 2;
        }
    }

    fft2(data, width, -1); // perform the inverse fft2

    // convert back to the 2D image
    ptr = 0; // pointer for the data array
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            sp_re[i][j] = data[ptr];
            ptr += 2;
        }
    }

    // overflow and underflow handling (by linear scaling)
    float min = 1E99;
    float max = -1E99;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if (sp_re[i][j] > max) max = sp_re[i][j];
            if (sp_re[i][j] < min) min = sp_re[i][j];
        }
    }
    float sf;
    sf = 255 / (max - min);
    //printf("%f \t %f \t %f \n", max, min, sf);
    for (int i = 0; i <= height; i++) {
        for (int j = 0; j < width; j++) {
            int value = (int)(sp_re[i][j]);
            if (value > 255) {
                value = 255;
            }
            else if (value < 0) {
                value = 0;
            }
            R[j][height - i] = value;
            //R[j][height - i] = (int)((sp_re[i][j] - min) * sf);
        }
    }
    save_bmp("output_images/output_butterworth_highpass.bmp", R, R, R);
}

void dilation(int** arr, int** structuring_element, int arr_width, int arr_height, int element_size) {
    // 對arr做zero_padding, 方便convolve
    int** padded_arr = zero_padding(arr, arr_width, arr_height, element_size);
    int convolved_width = arr_width;
    int convolved_height = arr_height;
    for (int i = element_size / 2; i < arr_height + element_size / 2; ++i) {
        for (int j = element_size / 2; j < arr_width + element_size / 2; ++j) {
            int sum = 0;
            for (int u = 0; u < element_size; ++u) {
                for (int v = 0; v < element_size; ++v) {
                    int arr_x = j + v - element_size / 2;
                    int arr_y = i + u - element_size / 2;
                    sum += padded_arr[arr_y][arr_x] * structuring_element[u][v];
                }
            }
            if (sum > 0) {
                arr[i - element_size / 2][j - element_size / 2] = 255;
            }
            else {
                arr[i - element_size / 2][j - element_size / 2] = 0;
            }
            sum = 0;
        }
    }
}

void erosion(int** arr, int** structuring_element, int arr_width, int arr_height, int element_size) {
    // 對arr做zero_padding, 方便convolve
    int** padded_arr = zero_padding(arr, arr_width, arr_height, element_size);
    int convolved_width = arr_width;
    int convolved_height = arr_height;
    int sum_of_structuring_element = 0;
    for (int i = 0; i < element_size; ++i) {
        for (int j = 0; j < element_size; ++j) {
            if (structuring_element[i][j] > 0) sum_of_structuring_element += 1;
        }
    }
    for (int i = element_size / 2; i < arr_height + element_size / 2; ++i) {
        for (int j = element_size / 2; j < arr_width + element_size / 2; ++j) {
            int sum = 0;
            for (int u = 0; u < element_size; ++u) {
                for (int v = 0; v < element_size; ++v) {
                    int arr_x = j + v - element_size / 2;
                    int arr_y = i + u - element_size / 2;
                    sum += padded_arr[arr_y][arr_x] * structuring_element[u][v];
                }
            }
            if (sum == sum_of_structuring_element * 255) {
                arr[i - element_size / 2][j - element_size / 2] = 255;
            }
            else {
                arr[i - element_size / 2][j - element_size / 2] = 0;
            }
            sum = 0;
        }
    }
}


void morphological_noise_removal() {
    int element_size = 23;
    int** se = new int* [element_size];
    for (int i = 0; i < element_size; ++i) {
        se[i] = new int[element_size];
    }
    for (int i = 0; i < element_size; ++i) {
        for (int j = 0; j < element_size; ++j) {
            se[i][j] = 1;
        }
    }
    int width, height;
    open_bmp("images/noisy_rectangle.bmp", R, R, R, width, height);
    int** arr_ptr = new int* [height];
    for (int i = 0; i < height; ++i) {
        arr_ptr[i] = R[i];
    }
    //opening
    erosion(arr_ptr, se, width, height, element_size);
    dilation(arr_ptr, se, width, height, element_size);
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            r[i][j] = arr_ptr[i][j];
        }
    }

    for (int i = 0; i < height; ++i) {
        arr_ptr[i] = R[i];
    }
    //closing
    dilation(arr_ptr, se, width, height, element_size);
    erosion(arr_ptr, se, width, height, element_size);
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            r[i][j] = arr_ptr[i][j];
        }
    }

    save_bmp("output_images/output_morphological_noise_removal.bmp", r, r, r);
}

void edge_detection_sobel() {
    int width, height;
    open_bmp("images/house.bmp", R, R, R, width, height);
    float sobel_x[3][3] = { {1,0,-1},{2,0,-2},{1,0,-1} };
    float** sobel_x_ptr = new float* [3];
    for (int i = 0; i < 3; ++i) {
        sobel_x_ptr[i] = sobel_x[i];
    }
    float sobel_y[3][3] = { {1,2,1},{0,0,0},{-1,-2,-1} };
    float** sobel_y_ptr = new float* [3];
    for (int i = 0; i < 3; ++i) {
        sobel_y_ptr[i] = sobel_y[i];
    }
    int** arr_ptr = new int* [height];
    for (int i = 0; i < height; ++i) {
        arr_ptr[i] = R[i];
    }
    int** padded_arr = zero_padding(arr_ptr, width, height, 3);
    float** Gx = convolve(padded_arr, sobel_x_ptr, width + 2, height + 2, 3, 3);
    float** Gy = convolve(padded_arr, sobel_y_ptr, width + 2, height + 2, 3, 3);
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            float value = sqrt(pow(Gx[i][j], 2) + pow(Gy[i][j], 2));
            if (value < 0.0f) {
                value = 0.0f;
            }
            else if (value > 255.0f) {
                value = 255.0f;
            }
            else {
                value = (int)round(value);
            }
            if (value > 128) {
                R[i][j] = value;
            }
            else {
                R[i][j] = 0;
            }

        }
    }
    save_bmp("output_images/output_edge_detection_sobel.bmp", R, R, R);
    close_bmp();
}

void test_dilation() {
    int element_size = 23;
    int** se = new int* [element_size];
    for (int i = 0; i < element_size; ++i) {
        se[i] = new int[element_size];
    }
    for (int i = 0; i < element_size; ++i) {
        for (int j = 0; j < element_size; ++j) {
            se[i][j] = 1;
        }
    }
    int width, height;
    open_bmp("output_images/test_erosion.bmp", R, R, R, width, height);
    int** arr_ptr = new int* [height];
    for (int i = 0; i < height; ++i) {
        arr_ptr[i] = R[i];
    }

    dilation(arr_ptr, se, width, height, element_size);
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            r[i][j] = arr_ptr[i][j];
        }
    }
    save_bmp("output_images/test_dilation.bmp", r, r, r);
}

void test_erosion() {
    int element_size = 23;
    int** se = new int* [element_size];
    for (int i = 0; i < element_size; ++i) {
        se[i] = new int[element_size];
    }
    for (int i = 0; i < element_size; ++i) {
        for (int j = 0; j < element_size; ++j) {
            se[i][j] = 1;
        }
    }
    int width, height;
    open_bmp("images/noisy_rectangle.bmp", R, R, R, width, height);
    int** arr_ptr = new int* [height];
    for (int i = 0; i < height; ++i) {
        arr_ptr[i] = R[i];
    }

    erosion(arr_ptr, se, width, height, element_size);
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            r[i][j] = arr_ptr[i][j];
        }
    }
    save_bmp("output_images/test_erosion.bmp", r, r, r);
}

int main(int argc, char* argv[])
{
    //sample();
    //test_erosion();
    //test_dilation();

    //line_detection(10);
    //butterworth_lowpass(20, 2);
    //butterworth_highpass(20, 2);
    morphological_noise_removal();



    cout << "Job Finished!" << endl;
    return 0;
}
