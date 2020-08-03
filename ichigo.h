#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string>
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


/*#define DRND(x) ((double)(x)/RAND_MAX*rand())//乱数の関数設定

#define ND 256			//差分計算における計算領域一辺の分割数(高速フーリエ変換を用いるため２のべき乗)
#define IG 8				//2^IG=ND
#define SIZEX (256)
#define SIZEY (256)
#define SIZEZ 1
#define SIZE (SIZEX*SIZEY*SIZEZ)*/

	extern int nd, ndm; 	//計算領域の一辺の差分分割数(差分ブロック数)、ND-1を定義
	extern int nd2;				 	//ND/2を定義：高速フ−リエ変換で使用
	extern int ig;						//2^ig=ND
	extern double time1;					//計算カウント数(時間に比例)


	extern int name_id;
	extern int num_id;
	extern double filter[3][3][3];

	extern std::string m_path;
	extern std::string p_path;
	extern double time1max;
	extern int graph_step;
	extern int name_max;
	extern int num_max;
	extern double Ms;
	extern double alpha;
  	extern double K1, K2;
	extern double Ks;
	extern double stretch;
  	extern double ram100, ram111;
  	extern double c11, c12, c44;
	extern double A;
	extern double delt;
	extern double mu0;
	extern double ld;
	extern double G_para;
	extern double B;
	extern double smob;
	extern double ds_fac;

  	extern double Astar;
	extern double G;
	extern double xf[256];
	extern double yf[256];
	extern double zf[1];
	extern double fai[256][256];
	extern double faifour[256][256];
	extern double faifour_i[256][256];
	extern double m_ave[3];
	extern double N[3];
	extern double m[256][256][3];
	extern double mstar1[256][256][3];
	extern double mstar2[256][256][3];
	extern double mfour[256][256][3];
	extern double mstarfour[256][256][3];
	extern double mstar2four[256][256][3];
	extern double mfour_i[256][256][3];
	extern double mstarfour_i[256][256][3];
	extern double mstar2four_i[256][256][3];

	extern double g[256][256][3];
	extern double gstar[256][256][3];
	extern double gfour[256][256][3];
	extern double gstarfour[256][256][3];
	extern double gfour_i[256][256][3];
	extern double gstarfour_i[256][256][3];
	extern double h[256][256][3];
	extern double hfour[256][256][3];
	extern double hfour_i[256][256][3];
	extern double Hanis[256][256][3];
	extern double Hms[256][256][3];
	extern double Hexternal[256][256][3];
	extern double Helastic[256][256][3];
	extern double Hme[256][256][3];
    extern double p[256][256][2];
    extern double pfour[256][256][2];
    extern double pfour_i[256][256][2];
    extern double p_prefour[256][256][2];
    extern double p_prefour_i[256][256][2];
    extern double p_postfour[256][256][2];
    extern double p_postfour_i[256][256][2];
	extern double Dr[256][256][3];
	extern double Drfour[256][256][3];
	extern double Drfour_i[256][256][3];
	extern double Dr_prefour[256][256][3];
	extern double Dr_prefour_i[256][256][3];

	extern double Dla256au[256][256][2];
	extern double Dgradient[256][256][2];
	extern double Dme[256][256][2];
	extern double Delastic[256][256][2];
	extern double Eanis[256][256];
	extern double Eexch[256][256];
	extern double Eexternal[256][256];
	extern double Ems[256][256];
	extern double Eme[256][256];
	extern double Eanis_ave;
	extern double Eexch_ave;
	extern double Eexternal_ave;
	extern double Ems_ave;
	extern double Eme_ave;
	extern double fourier_output[256][256];
	extern double fourier_output_i[256][256];
	extern double fourier_input[256][256];
	extern double filter_output[256][256];
	extern double epsilon_zero[256][256][3][3];
	extern double epsilon_zerofour[256][256][3][3];
	extern double epsilon_zerofour_i[256][256][3][3];
	extern double epsilon_zero_prefour[256][256][3][3];
	extern double epsilon_zero_prefour_i[256][256][3][3];
	extern double epsilon_zero_postfour[256][256][3][3];
	extern double epsilon_zero_postfour_i[256][256][3][3];
	extern double epsilon_phase[3][3];
	extern double epsilon_eigen[256][256][3][3];

	extern double epsilon_sum[3][3];
	extern double epsilon_homo[3][3];
	extern double eta[256][256][3][3];
	extern double epsilon_zero_grad[256][256][3][3][2];
	extern double epsilon_zero_grad_m[256][256][3][3][3];
	extern double eta_grad[256][256][3][3][3];
	extern double c[3][3][3][3];
	extern double s[3][3][3][3];
	//extern Eigen::MatrixXd c_matrix(6,6);
	//extern Eigen::MatrixXd s_matrix(6,6);
	extern double Dfour[256][256];
	extern double Dfour_i[256][256];
	extern double Nfour[256][256][3][3];
	extern double Nfour_i[256][256][3][3];
	extern double u[256][256][3];
	extern double ufour[256][256][3];
	extern double ufour_i[256][256][3];
	extern double u_grad[256][256][3][3];
	extern double sigma_a[3][3];


    int evolution(std::string file_path);
	void load_txt(std::string file_path);
	void ini000(std::string file_path);			//初期場の設定サブル−チン
	void apply_stress(double step_num, std::string file_path);			//計算途中で条件を変更(応力印加など)
	void graph_s1();		//組織描画サブル−チン

	int DCexchange2D();
	int fft3d();
	int ifft3d();
	int convolution3D(int switch_num);
	void grad_fai();
	void grad_u();
	int four_axis(int k,int i,int j);

