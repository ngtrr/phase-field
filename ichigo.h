#pragma once

#include <stdio.h>
#include <stdlib.h>
#include<string>
#include <time.h> 
#include <math.h>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <opencv2/opencv.hpp>
#include <complex.h>
#include <fftw3.h>
#include <Eigen/Core>
#include <Eigen/LU>
#include <sys/stat.h>


#define DRND(x) ((double)(x)/RAND_MAX*rand())//乱数の関数設定

#define ND 256			//差分計算における計算領域一辺の分割数(高速フーリエ変換を用いるため２のべき乗)
#define IG 8				//2^IG=ND
#define SIZEX (ND)
#define SIZEY (ND)
#define SIZEZ 1
#define SIZE (SIZEX*SIZEY*SIZEZ)

	/*extern int nd=ND, ndm=ND-1; 	//計算領域の一辺の差分分割数(差分ブロック数)、ND-1を定義
	extern int nd2=ND/2;				 	//ND/2を定義：高速フ−リエ変換で使用
	extern int ig=IG;						//2^ig=ND
	extern double time1;					//計算カウント数(時間に比例)


	extern int name_id = 0;
	extern int num_id = 0;

	extern double filter[3][3][3];


	extern string m_path;
	extern string p_path;
	extern double time1max = 100000;
	extern int graph_step = 100;
	extern int name_max = 0;
	extern int num_max = 0;
	extern double Ms = 6.02E+5;
	extern double alpha=1.5;
  	extern double K1 = 2.7E+3, K2 = -6.1E+3;
	extern double Ks = 1.65E+8;
	extern double stretch = 0.034;
  	extern double ram100 = 0.0E-2, ram111 = 1.64E-3;
  	extern double c11 = 1.6E+11, c12 = 1.52E+11, c44 = 0.43E+11;
	extern double A = 4.0E-9;
	extern double delt = 0.1;
	extern double mu0 = 1.0;
	extern double ld = 18E-9;
	extern double G_para = 8.0E-11;
	extern double B = 2.00E+11;
	extern double smob = 4.0E-9;
	extern double ds_fac = 0.01;

  	extern double Astar;
	extern double G = G_para/(ld*ld);
	extern double xf[SIZEX];
	extern double yf[SIZEY];
	extern double zf[SIZEZ];
	extern double fai[ND][ND];
	extern double faifour[ND][ND];
	extern double faifour_i[ND][ND];
	extern double m_ave[3];
	extern double N[3];
	extern double m[ND][ND][3];
	extern double mstar1[ND][ND][3];
	extern double mstar2[ND][ND][3];
	extern double mfour[ND][ND][3];
	extern double mstarfour[ND][ND][3];
	extern double mstar2four[ND][ND][3];
	extern double mfour_i[ND][ND][3];
	extern double mstarfour_i[ND][ND][3];
	extern double mstar2four_i[ND][ND][3];

	extern double g[ND][ND][3];
	extern double gstar[ND][ND][3];
	extern double gfour[ND][ND][3];
	extern double gstarfour[ND][ND][3];
	extern double gfour_i[ND][ND][3];
	extern double gstarfour_i[ND][ND][3];
	extern double h[ND][ND][3];
	extern double hfour[ND][ND][3];
	extern double hfour_i[ND][ND][3];
	extern double Hanis[ND][ND][3];
	extern double Hms[ND][ND][3];
	extern double Hexternal[ND][ND][3];
	extern double Helastic[ND][ND][3];
	extern double Hme[ND][ND][3];
    extern double p[ND][ND][2];
    extern double pfour[ND][ND][2];
    extern double pfour_i[ND][ND][2];
    extern double p_prefour[ND][ND][2];
    extern double p_prefour_i[ND][ND][2];
    extern double p_postfour[ND][ND][2];
    extern double p_postfour_i[ND][ND][2];
	extern double Dr[ND][ND][3];
	extern double Drfour[ND][ND][3];
	extern double Drfour_i[ND][ND][3];
	extern double Dr_prefour[ND][ND][3];
	extern double Dr_prefour_i[ND][ND][3];

	extern double Dlandau[ND][ND][2];
	extern double Dgradient[ND][ND][2];
	extern double Dme[ND][ND][2];
	extern double Delastic[ND][ND][2];
	extern double Eanis[ND][ND];
	extern double Eexch[ND][ND];
	extern double Eexternal[ND][ND];
	extern double Ems[ND][ND];
	extern double Eme[ND][ND];
	extern double Eanis_ave;
	extern double Eexch_ave;
	extern double Eexternal_ave;
	extern double Ems_ave;
	extern double Eme_ave;
	extern double fourier_output[ND][ND];
	extern double fourier_output_i[ND][ND];
	extern double fourier_input[ND][ND];
	extern double filter_output[ND][ND];
	extern double epsilon_zero[ND][ND][3][3];
	extern double epsilon_zerofour[ND][ND][3][3];
	extern double epsilon_zerofour_i[ND][ND][3][3];
	extern double epsilon_zero_prefour[ND][ND][3][3];
	extern double epsilon_zero_prefour_i[ND][ND][3][3];
	extern double epsilon_zero_postfour[ND][ND][3][3];
	extern double epsilon_zero_postfour_i[ND][ND][3][3];
	extern double epsilon_phase[3][3];
	extern double epsilon_eigen[ND][ND][3][3];

	extern double epsilon_sum[3][3];
	extern double epsilon_homo[3][3];
	extern double eta[ND][ND][3][3];
	extern double epsilon_zero_grad[ND][ND][3][3][2];
	extern double epsilon_zero_grad_m[ND][ND][3][3][3];
	extern double eta_grad[ND][ND][3][3][3];
	extern double c[3][3][3][3];
	extern double s[3][3][3][3];
	extern MatrixXd c_matrix(6,6);
	extern MatrixXd s_matrix(6,6);
	extern double Dfour[ND][ND];
	extern double Dfour_i[ND][ND];
	extern double Nfour[ND][ND][3][3];
	extern double Nfour_i[ND][ND][3][3];
	extern double u[ND][ND][3];
	extern double ufour[ND][ND][3];
	extern double ufour_i[ND][ND][3];
	extern double u_grad[ND][ND][3][3];
	extern double sigma_a[3][3];*/


    int evolution(void);
	void load_txt(string file_path);
	void ini000();			//初期場の設定サブル−チン
	void apply_stress(double step_num, string file_path);			//計算途中で条件を変更(応力印加など)
	void graph_s1();		//組織描画サブル−チン

	int DCexchange2D();
	int fft3d();
	int ifft3d();
	int convolution3D(int switch_num);
	void grad_fai();
	void grad_u();
	int four_axis(int k,int i,int j);

