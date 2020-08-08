/*
* test of self compiled fftw3
* 2D version
* @author maeda
* @date 2008/09/09
*/
 
#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>
#include <iostream>
#include <vector>
#define _USE_MATH_DEFINES

using namespace std;

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>


#include "fftw3.h"
#pragma comment( lib, "fftw3.lib" )
 
#define SIZEX 256
#define SIZEY 256
#define SIZEZ 1
#define SIZE (SIZEX*SIZEY*SIZEZ)


cv::Mat_<uchar> image = cv::imread("image/a.jpg", 0);

int DCexchange2D( fftw_complex *data, int cols, int rows, int depth )
{
	int i,j,k;
	int p1,p2,p3;    // point position
	int c2,r2,d2;    // temporary for cols/2,rows/2
	double re,im; // temporary

	if( data==NULL )       return false;
	if( (rows<0 || cols<0) || depth<0) return false;
 
	c2 = cols/2;
	r2 = rows/2;
    d2 = depth/2;
 
	for( k=0; k<depth; k++ ){
		for( j=0; j<r2; j++ ){
            for ( i=0; i<cols; i++ ){
                // exchange p1( i, j ) <-> p2( (cols/2+i)%cols, rows/2+j )
                p1 = j*cols + i + k*rows*cols;
                p2 = (r2+j)*cols + (c2+i)%cols + k*rows*cols;
                re = data[p1][0];
                im = data[p1][1];
                data[p1][0] = data[p2][0];
                data[p1][1] = data[p2][1];
                data[p2][0] = re;
                data[p2][1] = im;
            }
		}
	}

	for( k=0; k<d2; k++ ){
		for( j=0; j<rows; j++ ){
            for ( i=0; i<cols; i++ ){
                // exchange p1( i, j ) <-> p2( (cols/2+i)%cols, rows/2+j )
                p1 = j*cols + i + k*rows*cols;
                p2 = j*cols + i + (d2+k)*rows*cols;
                re = data[p1][0];
                im = data[p1][1];
                data[p1][0] = data[p2][0];
                data[p1][1] = data[p2][1];
                data[p2][0] = re;
                data[p2][1] = im;
            }
		}
	}

	return true;
}

int fft3d(void){

	fftw_complex *in  = NULL;
	fftw_complex *out = NULL;
	fftw_complex *out2 = NULL;
	fftw_plan p, ip   = NULL;
	int i,j,k,idx;
 
	size_t mem_size = sizeof(fftw_complex) * SIZE;
	in  = (fftw_complex*)fftw_malloc( mem_size );
	out = (fftw_complex*)fftw_malloc( mem_size );
	out2 = (fftw_complex*)fftw_malloc( mem_size );

	if( !in || !out ){
		fprintf( stderr, "failed to allocate %d[byte] memory(-.-)\n", (int)mem_size );
		return false;
	}
 
	// !! row-major alignment is recommended, but here, column-major.
	p = fftw_plan_dft_3d( SIZEY, SIZEX, SIZEZ, in, out, FFTW_FORWARD, FFTW_ESTIMATE );
	ip = fftw_plan_dft_3d( SIZEY, SIZEX, SIZEZ, out, out2, FFTW_BACKWARD, FFTW_ESTIMATE );

	// input data creation
	//printf("----- INPUT -----\n");
	for( k=0; k<SIZEZ; k++ ){
        for( j=0; j<SIZEY; j++ ){
            for( i=0; i<SIZEX; i++ ){
                idx = SIZEZ*k+SIZEX*j+i; // column-major alignment
                in[idx][0] = image[i][j];  //1 + 2*sin(2*M_PI*i/SIZEX) + sin(4*M_PI*j/SIZEY);
                in[idx][1] = 0;
            }
        }
    }
 
	fftw_execute(p);
    DCexchange2D(out, SIZEX, SIZEY, SIZEZ);
    //cv::Mat resultimage[SIZEX][SIZEY];
 
	// output is DC exchanged and scaled.
	double scale = 1. / SIZE;
	//printf("\n----- RESULT -----\n");
	for( j=0; j<SIZEY; j++ ){
		for( i=0; i<SIZEX; i++ ){
			idx = SIZEX*j+i;
			//printf("%d %d %lf %lf\n", i, j, out[idx][0]*scale, out[idx][1]*scale );
            image[i][j] = int(abs(out[idx][0]));
		}
	}
    cv::imwrite("Result1.png", image);
 
    //DCexchange2D(out2, SIZEX, SIZEY, SIZEZ);
	fftw_execute(ip);
	//printf("\n----- RESULT -----\n");
	for( j=0; j<SIZEY; j++ ){
		for( i=0; i<SIZEX; i++ ){
			idx = SIZEX*j+i;
			//printf("%d %d %lf %lf\n", i, j, out2[idx][0]*scale, out2[idx][1]*scale );
            image[i][j] = int(abs(out2[idx][0]*scale));
		}
	}
    cv::imwrite("Result2.png", image);

	if( ip   ) fftw_destroy_plan(ip);
	if( p   ) fftw_destroy_plan(p);
	if( in  ) fftw_free(in);
	if( out ) fftw_free(out);
	if( out2 ) fftw_free(out2);

    return true;
}



int main( void ){

    fft3d();
    //vector<cv::Mat> planes;
    //cv::split(resultimage, planes);

    //cout << "dims:" << resultimage.dims << endl;
    //cv::imwrite("Result.png", image);
 
	return true;
}