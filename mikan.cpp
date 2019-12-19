//g++ -o samplepleple example.cpp `pkg-config --cflags opencv4` `pkg-config --libs opencv4`
//for magenetostriction model by llg equation with phase-field simulation only 2d

/*
plot "LLG_Galfenol_300_m_2d.png" binary filetype=png with rgbimage
replot "check300_2d.txt" with vector lc rgb "#000000"
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <opencv2/opencv.hpp>
#include <complex.h>
#include <fftw3.h>
#include <Eigen/Core>
#include <Eigen/LU>

//#include "wingxa.h"


using namespace std;
using namespace Eigen;

#define DRND(x) ((double)(x)/RAND_MAX*rand())//乱数の関数設定

#define ND 256			//差分計算における計算領域一辺の分割数(高速フーリエ変換を用いるため２のべき乗)
#define IG 8				//2^IG=ND
#define SIZEX (ND)
#define SIZEY (ND)
#define SIZEZ 1
#define SIZE (SIZEX*SIZEY*SIZEZ)

	int nd=ND, ndm=ND-1; 	//計算領域の一辺の差分分割数(差分ブロック数)、ND-1を定義
	int nd2=ND/2;				 	//ND/2を定義：高速フ−リエ変換で使用
	int ig=IG;						//2^ig=ND
	double alpha=0.5;
	double time1;					//計算カウント数(時間に比例)
	double time1max = 10000;

	double filter[3][3][3];


    //*******************  FePd  ************************
	/*double Ms = 6.02E+5;
  	double K1 = 2.7E+3, K2 = -6.1E+3;
  	double ram100 = 0.0E+4, ram111 = 1.64E-3;
  	double c11 = 1.495E+11, c12 = 1.43E+12, c44 = 7.05E+10;
	double A = 2.0E-11;
  	double Astar;
	double delt = 0.01;
	double mu0 = 1.0;
	double ld = 1.8E-8;*/

    //******************* NiMnGa ************************
	double Ms = 6.02E+5;
  	double K1 = 2.7E+3, K2 = -6.1E+3;
  	double ram100 = 0.0E+4, ram111 = 1.64E-3;
  	double c11 = 1.6E+11, c12 = 1.52E+11, c44 = 0.43E+11;
	double A = 2.0E-11;
  	double Astar;
	double delt = 0.01;
	double mu0 = 1.0;
	double ld = 18E-9;
	double G = 1.0E-3/(250*ld*ld);

	double Q1 = 2.32E+10;
	double Q2 = -0.63E+8;
	double Q3 = 0.40E+10;
	double Q4 = 7.50E+10;
	//double Q1 = 2.4;
	//double Q2 = -1.0;
	//double Q3 = 0.83;
	//double Q4 = 0.04;

	double B = 4.00E+6;


	double smob = 1.0E-11;
	double ds_fac = 0.001;

	double root3 = 1.73205080757;
	double root2 = 1.41421356237;
	double root6 = 2.44948974278;

	double xf[SIZEX];
	double yf[SIZEY];
	double zf[SIZEZ];

	double fai[ND][ND];
	double faifour[ND][ND];
	double faifour_i[ND][ND];

	double m_ave[3];

	double N[3];

	double m[ND][ND][3];
	double mstar1[ND][ND][3];
	double mstar2[ND][ND][3];
	double mfour[ND][ND][3];
	double mstarfour[ND][ND][3];
	double mstar2four[ND][ND][3];
	double mfour_i[ND][ND][3];
	double mstarfour_i[ND][ND][3];
	double mstar2four_i[ND][ND][3];

	double g[ND][ND][3];
	double gstar[ND][ND][3];
	double gfour[ND][ND][3];
	double gstarfour[ND][ND][3];
	double gfour_i[ND][ND][3];
	double gstarfour_i[ND][ND][3];

	double h[ND][ND][3];
	double hfour[ND][ND][3];
	double hfour_i[ND][ND][3];

	double Hanis[ND][ND][3];
	double Hms[ND][ND][3];
	double Hexternal[ND][ND][3];
	double Helastic[ND][ND][3];
	double Hme[ND][ND][3];

	double P[ND][ND][3];
	double Pfour[ND][ND][3];
	double Pfour_i[ND][ND][3];

	double P_prefour[ND][ND][3];
	double P_prefour_i[ND][ND][3];

	double Plandau[ND][ND][3];
	double Pgradient[ND][ND][3];
	double Pme[ND][ND][3];
	double Pelastic[ND][ND][3];

	double fourier_output[ND][ND];
	double fourier_output_i[ND][ND];
	double fourier_input[ND][ND];

	double filter_output[ND][ND];

	double epsilon_zero[ND][ND][3][3];
	double epsilon_zerofour[ND][ND][3][3];
	double epsilon_zerofour_i[ND][ND][3][3];

	double epsilon_zero_prefour[ND][ND][3][3];
	double epsilon_zero_prefour_i[ND][ND][3][3];

	double epsilon_zero_postfour[ND][ND][3][3];
	double epsilon_zero_postfour_i[ND][ND][3][3];

	double epsilon_sum[3][3];
	double epsilon_homo[3][3];

	double eta[ND][ND][3][3];

	double epsilon_zero_grad[3][3][3];
	double eta_grad[ND][ND][3][3][3];

	double c[3][3][3][3];
	double s[3][3][3][3];
	MatrixXd c_matrix(6,6);
	MatrixXd s_matrix(6,6);

	double Dfour[ND][ND];
	double Dfour_i[ND][ND];
	double Nfour[ND][ND][3][3];
	double Nfour_i[ND][ND][3][3];

	double u[ND][ND][3];
	double ufour[ND][ND][3];
	double ufour_i[ND][ND][3];
	double u_grad[ND][ND][3][3];

	double sigma_a[3][3];

	double e1[ND][ND];
	double e2[ND][ND];
	double e3[ND][ND];
	double e4[ND][ND];
	double e5[ND][ND];
	double e6[ND][ND];

	double e1_grad[3];
	double e2_grad[3];
	double e3_grad[3];
	double e4_grad[3];
	double e5_grad[3];
	double e6_grad[3];


	void ini000();			//初期場の設定サブル−チン
	void graph_s1();		//組織描画サブル−チン
	void table();				//sinとcosのテーブルとビット反転テーブルの作成サブル−チン
	void fft();					//１次元高速フーリエ変換
	void rcfft();				//２次元高速フーリエ変換
	void datsave();			//デ−タ保存サブル−チン

	int DCexchange2D();
	int fft3d();
	int ifft3d();
	int convolution3D(int switch_num);
	void grad_fai();
	void grad_u();
	int four_axis(int k,int i,int j);


//******* メインプログラム ******************************************
int main(void)
{

	int   i, j, k, l, ii, jj, kk, ll, iii, jjj;		//整数
	int   ip, im, jp, jm, Nstep;							//整数
	double mlength;

	srand(time(NULL));

//****** 計算条件および物質定数の設定 ****************************************
	Astar = (2 * A)/(mu0 * Ms * Ms * ld * ld);
	Astar = 0.0625 ;
	cout << "Astar : " << Astar << endl;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<3;k++){
				Hms[i][j][k] = 0;//init
				Helastic[i][j][k] = 0;//ok
			}
			Hexternal[i][j][0] = 0.0E+6;//ok
			Hexternal[i][j][1] = 0.0E+6;//ok
			Hexternal[i][j][2] = 0.0E+6;//ok

			//Hanis[i][j][0] = 0;//(4 * K1)/(3 );// * 1.0E+7;//ok
			//Hanis[i][j][1] = 0;//(4 * K1)/(3 );// * 1.0E+7;//ok
			//Hanis[i][j][2] = 0;//(4 * K1)/(3 );// * 1.0E+7;//ok

		}
	}

	sigma_a[0][0] = 0;
	sigma_a[1][1] = 0;
	sigma_a[2][2] = 0;

	sigma_a[0][1] = 0;
	sigma_a[0][2] = 0;
	sigma_a[1][0] = 0;
	sigma_a[1][2] = 0;
	sigma_a[2][0] = 0;
	sigma_a[2][1] = 0;

	for(i=0;i<=ndm;i++){
		xf[i] = i - nd2;
		yf[i] = i - nd2;
	}

	N[0] = 0.333;
	N[1] = 0.333;
	N[2] = 0.333;

	c[0][0][0][0] = c[1][1][1][1] = c[2][2][2][2] = c11;
	c[1][2][1][2] = c[0][2][0][2] = c[0][1][0][1] = c44;
	c[0][0][1][1] = c[0][0][2][2] = c[1][1][2][2] = c[1][1][0][0] = c[2][2][0][0] = c[2][2][1][1] = c12;

	c_matrix(0,0) = c_matrix(1,1) = c_matrix(2,2) = c11;
	c_matrix(3,3) = c_matrix(4,4) = c_matrix(5,5) = c44;
	c_matrix(0,1) = c_matrix(0,2) = c_matrix(1,2) = c_matrix(1,0) = c_matrix(2,0) = c_matrix(2,1) = c12;

	s_matrix = c_matrix.inverse();

	s[0][0][0][0] = s[1][1][1][1] = s[2][2][2][2] = s_matrix(0,0);
	s[1][2][1][2] = s[0][2][0][2] = s[0][1][0][1] = s_matrix(3,3);
	s[0][0][1][1] = s[0][0][2][2] = s[1][1][2][2] = s[1][1][0][0] = s[2][2][0][0] = s[2][2][1][1] = s_matrix(0,1);


//*** sinおよびcosテ−ブル、ビット反転テーブル、および初期場の設定 ***************

	ini000();		//初期場の設定
 	//gwinsize(INXY,INXY); ginit(2); gsetorg(0,0);//描画Window表示

//**** シミュレーションスタート ******************************
start: ;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			Hanis[i][j][0] = -1/(mu0 * Ms) * (K1*(2 * m[i][j][0] * m[i][j][1] * m[i][j][1] + 2 * m[i][j][0] * m[i][j][2] * m[i][j][2]) + 2*K2 * m[i][j][0] * m[i][j][1] * m[i][j][1] * m[i][j][2] * m[i][j][2]);//ok
			Hanis[i][j][1] = -1/(mu0 * Ms) * (K1*(2 * m[i][j][1] * m[i][j][2] * m[i][j][2] + 2 * m[i][j][1] * m[i][j][0] * m[i][j][0]) + 2*K2 * m[i][j][1] * m[i][j][2] * m[i][j][2] * m[i][j][0] * m[i][j][0]);//ok
			Hanis[i][j][2] = -1/(mu0 * Ms) * (K1*(2 * m[i][j][2] * m[i][j][0] * m[i][j][0] + 2 * m[i][j][2] * m[i][j][1] * m[i][j][1]) + 2*K2 * m[i][j][2] * m[i][j][0] * m[i][j][0] * m[i][j][1] * m[i][j][1]);//ok
		}
	}

	//if(time1<=100.){Nstep=10;} else{Nstep=200;}		//データ保存する時間間隔の変更
	//if((((int)(time1) % Nstep)==0)) {datsave();} 	//一定繰返しカウント毎に組織データを保存
	if((((int)(time1) % 10)==0)) {graph_s1();} 		//一定繰返しカウント毎に組織を表示
	//if((((int)(time1) % 100)==0)) {datsave();} 		//一定繰返しカウント毎にデータを保存


	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			e1[i][j] = (epsilon_zero[i][j][0][0] + epsilon_zero[i][j][1][1] + epsilon_zero[i][j][2][2])/root3;
			e2[i][j] = (epsilon_zero[i][j][0][0] - epsilon_zero[i][j][1][1])/root2;
			e3[i][j] = (2*epsilon_zero[i][j][2][2] - epsilon_zero[i][j][1][1] - epsilon_zero[i][j][0][0])/root6;
			e4[i][j] = epsilon_zero[i][j][1][2];
			e5[i][j] = epsilon_zero[i][j][0][2];
			e6[i][j] = epsilon_zero[i][j][0][1];
		}
	}
	
	e1_grad[0] = 1/root3;
	e1_grad[1] = 1/root3;
	e1_grad[2] = 1/root3;
	e2_grad[0] = 1/root2;
	e2_grad[1] = -1/root2;
	e2_grad[2] = 0;
	e3_grad[0] = -1/root6;
	e3_grad[1] = -1/root6;
	e3_grad[2] = 2/root6;
	e4_grad[0] = 0;
	e4_grad[1] = 0;
	e4_grad[2] = 0;
	e5_grad[0] = 0;
	e5_grad[1] = 0;
	e5_grad[2] = 0;
	e6_grad[0] = 0;
	e6_grad[1] = 0;
	e6_grad[2] = 0;

//***** 化学ポテンシャル ************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<3;k++){
				Plandau[i][j][0] = Q1 * (2*e1[i][j]*e1_grad[0]) + Q2*(2*e2[i][j]*e2_grad[0] + 2*e3[i][j]*e3_grad[0]) + Q3*(3*e3[i][j]*e3[i][j]*e3_grad[0] - 3*(e3_grad[0]*e2[i][j]*e2[i][j] + e3[i][j]*2*e2[i][j]*e2_grad[0])) + Q4*(4*e2[i][j]*e2[i][j]*e2[i][j]*e2_grad[0] + 2*(2*e2[i][j]*e2_grad[0]*e3[i][j]*e3[i][j] + 2*e3[i][j]*e3_grad[0]*e2[i][j]*e2[i][j]) + 4*e3[i][j]*e3[i][j]*e3[i][j]*e3_grad[0]);
				Plandau[i][j][1] = Q1 * (2*e1[i][j]*e1_grad[1]) + Q2*(2*e2[i][j]*e2_grad[1] + 2*e3[i][j]*e3_grad[1]) + Q3*(3*e3[i][j]*e3[i][j]*e3_grad[1] - 3*(e3_grad[1]*e2[i][j]*e2[i][j] + e3[i][j]*2*e2[i][j]*e2_grad[1])) + Q4*(4*e2[i][j]*e2[i][j]*e2[i][j]*e2_grad[1] + 2*(2*e2[i][j]*e2_grad[1]*e3[i][j]*e3[i][j] + 2*e3[i][j]*e3_grad[1]*e2[i][j]*e2[i][j]) + 4*e3[i][j]*e3[i][j]*e3[i][j]*e3_grad[1]);
				Plandau[i][j][2] = Q1 * (2*e1[i][j]*e1_grad[2]) + Q2*(2*e2[i][j]*e2_grad[2] + 2*e3[i][j]*e3_grad[2]) + Q3*(3*e3[i][j]*e3[i][j]*e3_grad[2] - 3*(e3_grad[2]*e2[i][j]*e2[i][j] + e3[i][j]*2*e2[i][j]*e2_grad[2])) + Q4*(4*e2[i][j]*e2[i][j]*e2[i][j]*e2_grad[2] + 2*(2*e2[i][j]*e2_grad[2]*e3[i][j]*e3[i][j] + 2*e3[i][j]*e3_grad[2]*e2[i][j]*e2[i][j]) + 4*e3[i][j]*e3[i][j]*e3[i][j]*e3_grad[2]);
			}
		}
	}



//***** 勾配ポテンシャル ***********************
	/*for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1; jp=j+1; jm=j-1;
			if(i==ndm){ip=0;}  if(i==0){im=ndm;}	//周期的境界条件
			if(j==ndm){jp=0;}  if(j==0){jm=ndm;}
			for(k=0;k<3;k++){
				Pgradient[i][j][k] = -1 * G * (epsilon_zero[ip][j][k][k] + epsilon_zero[im][j][k][k] + epsilon_zero[i][jp][k][k] + epsilon_zero[i][jm][k][k] - 4.0*epsilon_zero[i][j][k][k]);
			}
		}
	}*/


//***** 磁気弾性ポテンシャル ***********************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<3;k++){
				Pme[i][j][k] = B * (m[i][j][k]*m[i][j][k] - 1/3);
			}
		}
	}


	for(k=0;k<3;k++){
		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
				if(i==j){
					if(i==k){
						epsilon_zero_grad[i][j][k] = 1;
					}else{
						epsilon_zero_grad[i][j][k] = 0;
					}
				}else{
					epsilon_zero_grad[i][j][k] = 0;
				}
			}
		}
	}


	for(ii=0;ii<3;ii++){
		for(jj=0;jj<3;jj++){
			epsilon_sum[ii][jj] = 0;
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					epsilon_sum[ii][jj] += epsilon_zero[i][j][ii][jj];
				}
			}
		}
	}


	epsilon_homo[0][0] = s[0][0][0][0] * sigma_a[0][0] + s[0][0][1][1] * (sigma_a[1][1] + sigma_a[2][2]) + epsilon_sum[0][0]/SIZE;
	epsilon_homo[1][1] = s[0][0][0][0] * sigma_a[1][1] + s[0][0][1][1] * (sigma_a[0][0] + sigma_a[2][2]) + epsilon_sum[1][1]/SIZE;
	epsilon_homo[2][2] = s[0][0][0][0] * sigma_a[2][2] + s[0][0][1][1] * (sigma_a[0][0] + sigma_a[1][1]) + epsilon_sum[2][2]/SIZE;
	//epsilon_homo[0][1] = epsilon_homo[1][0] = 1/2 * s[0][1][0][1] * sigma_a[0][1] + epsilon_sum[0][1]/SIZE;
	//epsilon_homo[1][2] = epsilon_homo[2][1] = 1/2 * s[0][1][0][1] * sigma_a[1][2] + epsilon_sum[1][2]/SIZE;
	//epsilon_homo[2][0] = epsilon_homo[0][2] = 1/2 * s[0][1][0][1] * sigma_a[2][0] + epsilon_sum[2][0]/SIZE;

	for(l=0;l<3;l++){
		for(k=0;k<3;k++){
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					fourier_input[i][j] = epsilon_zero[i][j][k][l];
				}
			}
			fft3d();
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					epsilon_zerofour[i][j][k][l] = fourier_output[i][j];
					epsilon_zerofour_i[i][j][k][l] = fourier_output_i[i][j];
				}
			}
		}
	}


	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			Dfour[i][j] = c44 * c44 * c11 * (xf[i] * xf[i] + yf[j] * yf[j] + zf[k] * zf[k]) * (xf[i] * xf[i] + yf[j] * yf[j] + zf[k] * zf[k]) * (xf[i] * xf[i] + yf[j] * yf[j] + zf[k] * zf[k]) + c44 * (c11-c12-2*c44) * (c11+c12) * (xf[i] * xf[i] + yf[j] * yf[j] + zf[k] * zf[k]) * (xf[i] * xf[i] * yf[j] * yf[j] + zf[k] * zf[k] * yf[j] * yf[j] + xf[i] * xf[i] * zf[k] * zf[k] ) + (c11 - c12 - 2*c44) * (c11 - c12 - 2*c44) * (c11+2*c12+c44) * xf[i] * xf[i] * yf[j] * yf[j] * zf[k] * zf[k];
			Nfour[i][j][0][0] = c44 * c44 + c44 * (c11 - c44) * (xf[i] * xf[i] + yf[j] * yf[j] + zf[k] * zf[k]) * (yf[j] * yf[j] + zf[k] * zf[k]) + (c11-c12-2*c44) * (c11+c12) * yf[j] * yf[j] * zf[k] * zf[k];
			Nfour[i][j][1][1] = c44 * c44 + c44 * (c11 - c44) * (xf[i] * xf[i] + yf[j] * yf[j] + zf[k] * zf[k]) * (xf[i] * xf[i] + zf[k] * zf[k]) + (c11-c12-2*c44) * (c11+c12) * xf[i] * xf[i] * zf[k] * zf[k];
			Nfour[i][j][2][2] = c44 * c44 + c44 * (c11 - c44) * (xf[i] * xf[i] + yf[j] * yf[j] + zf[k] * zf[k]) * (yf[j] * yf[j] + xf[i] * xf[i]) + (c11-c12-2*c44) * (c11+c12) * yf[j] * yf[j] * xf[i] * xf[i];
			Nfour[i][j][0][1] = Nfour[i][j][1][0] = -1 * (c12 + c44) * xf[i] * yf[j] * (c44 * (xf[i] * xf[i] + yf[j] * yf[j] + zf[k] * zf[k]) + (c11-c12-2*c44) * zf[k] * zf[k]);
			Nfour[i][j][0][2] = Nfour[i][j][2][0] = -1 * (c12 + c44) * zf[k] * xf[i] * (c44 * (xf[i] * xf[i] + yf[j] * yf[j] + zf[k] * zf[k]) + (c11-c12-2*c44) * yf[j] * yf[j]);
			Nfour[i][j][1][2] = Nfour[i][j][2][1] = -1 * (c12 + c44) * zf[k] * yf[j] * (c44 * (xf[i] * xf[i] + yf[j] * yf[j] + zf[k] * zf[k]) + (c11-c12-2*c44) * xf[i] * xf[i]);
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(kk=0;kk<3;kk++){
				ufour[i][j][kk] = 0;
				ufour_i[i][j][kk] = 0;
				for(jj=0;jj<3;jj++){
					for(ii=0;ii<3;ii++){
						for(ll=0;ll<3;ll++){
							if(c[ii][jj][kk][ll] * four_axis(jj, i, j) * four_axis(ll, i, j) == 0 ){
								ufour[i][j][kk] += 0;
								ufour_i[i][j][kk] += 0;
							}else{
								//ufour[i][j][kk] += 1 / (c[ii][jj][kk][ll] * four_axis(jj, i, j) * four_axis(ll, i, j)) * four_axis(jj, i, j) * c[ii][jj][kk][ll] * epsilon_zerofour_i[i][j][kk][ll];
								//ufour_i[i][j][kk] += -1 / (c[ii][jj][kk][ll] * four_axis(jj, i, j) * four_axis(ll, i, j)) * four_axis(jj, i, j) * c[ii][jj][kk][ll] * epsilon_zerofour[i][j][kk][ll];
								ufour[i][j][kk] += 1 * four_axis(jj, i, j) * c[ii][jj][kk][ll] * epsilon_zerofour_i[i][j][kk][ll] * Nfour[i][j][kk][jj] / Dfour[i][j];
								ufour_i[i][j][kk] += -1 * four_axis(jj, i, j) * c[ii][jj][kk][ll] * epsilon_zerofour[i][j][kk][ll] * Nfour[i][j][kk][jj] / Dfour[i][j];
								//cout << "ufour  : "<< ufour[i][j][kk] << ufour_i[i][j][kk] << endl;
							}
						}
					}
				}
			}
		}
	}

	for(k=0;k<3;k++){
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				fourier_output[i][j] = ufour[i][j][k];
				fourier_output_i[i][j] = ufour_i[i][j][k];
			}
		}
		ifft3d();
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				u[i][j][k] = fourier_input[i][j];
			}
		}
	}

	grad_u();

	for(l=0;l<3;l++){
		for(k=0;k<3;k++){
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					eta[i][j][k][l] = 0.5 * (u_grad[i][j][k][l] + u_grad[i][j][l][k]);
				}
			}
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<3;k++){
				Pelastic[i][j][k] = 0;

				for(ii=0;ii<3;ii++){
					for(jj=0;jj<3;jj++){
						for(kk=0;kk<3;kk++){
							for(ll=0;ll<3;ll++){
								Pelastic[i][j][k] += c[ii][jj][kk][ll] * (epsilon_homo[ii][jj] + eta[i][j][ii][jj] - epsilon_zero[i][j][ii][jj]) * epsilon_zero_grad[kk][ll][k];
							}
						}
					}
				}

			}
		}
	}



	cout << " ********************************************************** " << endl;
	cout << "Plandau - x   " << dec << -1 * Plandau[10][10][0] << endl;
	cout << "Pme - x   " << dec << -1 * Pme[10][10][0] << endl;
	cout << "Pelastic - x   " << dec << -1 * Pelastic[10][10][0] << endl;

	/*for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<3;k++){
				P[i][j][k]= -1 * smob * (Plandau[i][j][k] + Pgradient[i][j][k] + Pme[i][j][k] + Pelastic[i][j][k]);
				epsilon_zero[i][j][k][k] += ( P[i][j][k] + ds_fac*( 2.0*DRND(1.) - 1.0) )*delt;
				if(epsilon_zero[i][j][k][k]>=0.043){epsilon_zero[i][j][k][k]=0.043;}  if(epsilon_zero[i][j][k][k]<=-0.043){epsilon_zero[i][j][k][k]=-0.043;}
			}
		}
	}*/


	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<3;k++){
				P[i][j][k] = Plandau[i][j][k] + Pelastic[i][j][k];
			}
		}
	}

	if(time1 == 0){

		for(k=0;k<3;k++){
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					fourier_input[i][j] = P[i][j][k];
				}
			}
			fft3d();
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					Pfour[i][j][k] = fourier_output[i][j];
					Pfour_i[i][j][k] = fourier_output_i[i][j];
				}
			}
		}


		for(k=0;k<3;k++){
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					fourier_input[i][j] = epsilon_zero[i][j][k][k];
				}
			}
			fft3d();
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					epsilon_zerofour[i][j][k][k] = fourier_output[i][j];
					epsilon_zerofour_i[i][j][k][k] = fourier_output_i[i][j];
				}
			}
		}

		for(k=0;k<3;k++){
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					epsilon_zero_postfour[i][j][k][k] = (epsilon_zerofour[i][j][k][k] + -1*smob*delt*Pfour[i][j][k])/(1 + G*(four_axis(0,i,j)*four_axis(0,i,j) + four_axis(1,i,j)*four_axis(1,i,j) + four_axis(2,i,j)*four_axis(2,i,j))*smob*delt);
					epsilon_zero_postfour_i[i][j][k][k] = (epsilon_zerofour_i[i][j][k][k] + -1*smob*delt*Pfour_i[i][j][k])/(1 + G*(four_axis(0,i,j)*four_axis(0,i,j) + four_axis(1,i,j)*four_axis(1,i,j) + four_axis(2,i,j)*four_axis(2,i,j))*smob*delt);
				}
			}
		}

	}else{


		for(k=0;k<3;k++){
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					P_prefour[i][j][k] = Pfour[i][j][k];
					P_prefour_i[i][j][k] = Pfour_i[i][j][k];

					epsilon_zero_prefour[i][j][k][k] = epsilon_zerofour[i][j][k][k];
					epsilon_zero_prefour_i[i][j][k][k] = epsilon_zerofour_i[i][j][k][k];
				}
			}
		}



		for(k=0;k<3;k++){
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					fourier_input[i][j] = P[i][j][k];
				}
			}
			fft3d();
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					Pfour[i][j][k] = fourier_output[i][j];
					Pfour_i[i][j][k] = fourier_output_i[i][j];
				}
			}
		}


		for(k=0;k<3;k++){
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					fourier_input[i][j] = epsilon_zero[i][j][k][k];
				}
			}
			fft3d();
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					epsilon_zerofour[i][j][k][k] = fourier_output[i][j];
					epsilon_zerofour_i[i][j][k][k] = fourier_output_i[i][j];
				}
			}
		}

		for(k=0;k<3;k++){
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					epsilon_zero_postfour[i][j][k][k] = ((4*epsilon_zero_prefour[i][j][k][k] - epsilon_zerofour[i][j][k][k]) + 2*-1*smob*delt*(2*Pfour[i][j][k] - P_prefour[i][j][k]))/(3 + 2*G*(four_axis(0,i,j)*four_axis(0,i,j) + four_axis(1,i,j)*four_axis(1,i,j) + four_axis(2,i,j)*four_axis(2,i,j))*smob*delt);
					epsilon_zero_postfour_i[i][j][k][k] = ((4*epsilon_zero_prefour_i[i][j][k][k] - epsilon_zerofour_i[i][j][k][k]) + 2*-1*smob*delt*(2*Pfour_i[i][j][k] - P_prefour_i[i][j][k]))/(3 + 2*G*(four_axis(0,i,j)*four_axis(0,i,j) + four_axis(1,i,j)*four_axis(1,i,j) + four_axis(2,i,j)*four_axis(2,i,j))*smob*delt);
				}
			}
		}


	}


	for(k=0;k<3;k++){
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				fourier_output[i][j] = epsilon_zero_postfour[i][j][k][k];
				fourier_output_i[i][j] = epsilon_zero_postfour_i[i][j][k][k];
			}
		}
		ifft3d();
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				epsilon_zero[i][j][k][k] = fourier_input[i][j];
			}
		}
	}



	cout << "epsilon_zero - x   " << dec << epsilon_zero[10][10][0][0] << endl;
	cout << "epsilon_zero - y   " << dec << epsilon_zero[10][10][1][1] << endl;
	cout << "epsilon_zero - z   " << dec << epsilon_zero[10][10][2][2] << endl;
	/*cout << "epsilon_zero - xy  " << dec << epsilon_zero[10][10][0][1] << endl;
	cout << "epsilon_zero - yz  " << dec << epsilon_zero[10][10][1][2] << endl;
	cout << "epsilon_zero - zx  " << dec << epsilon_zero[10][10][2][0] << endl;*/
	cout << "m  -  x   " << dec << m[10][10][0] << endl;
	cout << "m  -  y   " << dec << m[10][10][1] << endl;
	cout << "m  -  z   " << dec << m[10][10][2] << endl;


	for(k=0;k<3;k++){
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				fourier_input[i][j] = m[i][j][k];
			}
		}
		fft3d();
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				mfour[i][j][k] = fourier_output[i][j];
				mfour_i[i][j][k] = fourier_output_i[i][j];
				if (isinf(mfour[i][j][k]) == 1){
					cout << "mfour        " << i << " : " << j << "   -    " << dec << mfour[i][j][k] << endl;
				}
			}
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//cout << "faifour        " << i << " : " << j << "   -    " << dec << faifour[i][j] << endl;
			if(xf[i]*xf[i] + yf[j]*yf[j] == 0){
				faifour[i][j] = 0;
				faifour_i[i][j] = 0;
				//faifour[i][j] = Ms*(mfour_i[i][j][0]*xf[i] + mfour_i[i][j][1]*yf[j] + mfour_i[i][j][2]*0 )/1;
				//faifour_i[i][j] = -1*Ms*(mfour[i][j][0]*xf[i] + mfour[i][j][1]*yf[j] + mfour[i][j][2]*0 )/1;
			}else{
				faifour[i][j] = Ms*(mfour_i[i][j][0]*xf[i] + mfour_i[i][j][1]*yf[j] + mfour_i[i][j][2]*0 )/(xf[i]*xf[i] + yf[j]*yf[j] + 0);
				faifour_i[i][j] = -1*Ms*(mfour[i][j][0]*xf[i] + mfour[i][j][1]*yf[j] + mfour[i][j][2]*0 )/(xf[i]*xf[i] + yf[j]*yf[j] + 0);
			}
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			fourier_output[i][j] = faifour[i][j];
			fourier_output_i[i][j] = faifour_i[i][j];
		}
	}
	ifft3d();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			fai[i][j] = fourier_input[i][j];
		}
	}

	grad_fai();

	for(k=0;k<3;k++){
		m_ave[k] = 0;
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				m_ave[k] += m[i][j][k] / SIZE;
			}
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<3;k++){
				//cout << "   " << endl;
				//cout << "Hms   " << Hms[i][j][k] << endl;
				Hms[i][j][k] = -1 * m_ave[k] * Ms * N[k];
				//cout << "Hms   " << Hms[i][j][k] << endl;
			}
		}
	}


	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			Hme[i][j][0] = -1/(mu0 * Ms) * 2 * B * epsilon_zero[i][j][0][0] * m[i][j][0];
			Hme[i][j][1] = -1/(mu0 * Ms) * 2 * B * epsilon_zero[i][j][1][1] * m[i][j][1];
			Hme[i][j][2] = -1/(mu0 * Ms) * 2 * B * epsilon_zero[i][j][2][2] * m[i][j][2];
		}
	}


	//*********************************  STEP 1  ******************************************************

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<3;k++){
				h[i][j][k] = (Hanis[i][j][k]  + Hms[i][j][k] + Hexternal[i][j][k] + Hme[i][j][k])/Ms;
				//cout << "h   " << h[i][j][k] * Ms << endl;
				//cout << "Hel   " << Helastic[i][j][k] * Ms << endl;
			}
		}
	}

	//cout << "m  :   " << m[100][100][0]  << " : " << m[100][100][1]  << " : " << m[100][100][2] << endl;
	//cout << "h  :   " << h[100][100][0]  << " : " << h[100][100][1]  << " : " << h[100][100][2] << endl;

	// hfour mfour　の計算 (fft)
	for(k=0;k<3;k++){
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				fourier_input[i][j] = h[i][j][k];
			}
		}
		fft3d();
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				hfour[i][j][k] = fourier_output[i][j];
				hfour_i[i][j][k] = fourier_output_i[i][j];
			}
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<3;k++){
				gfour[i][j][k] = (mfour[i][j][k] + delt*hfour[i][j][k])/(1+(xf[i]*xf[i] + yf[j]*yf[j] + 0)*Astar*delt);
				gfour_i[i][j][k] = (mfour_i[i][j][k] + delt*hfour_i[i][j][k])/(1+(xf[i]*xf[i] + yf[j]*yf[j] + 0)*Astar*delt);
			}
		}
	}
	

	for(k=0;k<3;k++){
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				fourier_output[i][j] = gfour[i][j][k];
				fourier_output_i[i][j] = gfour_i[i][j][k];
			}
		}
		ifft3d();
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				g[i][j][k] = fourier_input[i][j];
			}
		}
	}



	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			mstar1[i][j][0] = m[i][j][0] + (g[i][j][1] * m[i][j][2] - g[i][j][2] * m[i][j][1] );
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			fourier_input[i][j] = mstar1[i][j][0];
		}
	}
	fft3d();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			mstarfour[i][j][0] = fourier_output[i][j];
			mstarfour_i[i][j][0] = fourier_output_i[i][j];
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			gstarfour[i][j][0] = (mstarfour[i][j][0] + delt*hfour[i][j][0])/(1+(xf[i]*xf[i] + yf[j]*yf[j] + 0)*Astar*delt);
			gstarfour_i[i][j][0] = (mstarfour_i[i][j][0] + delt*hfour_i[i][j][0])/(1+(xf[i]*xf[i] + yf[j]*yf[j] + 0)*Astar*delt);
		}
	}


	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			fourier_output[i][j] = gstarfour[i][j][0];
			fourier_output_i[i][j] = gstarfour_i[i][j][0];
		}
	}
	ifft3d();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			gstar[i][j][0] = fourier_input[i][j];
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			mstar1[i][j][1] = m[i][j][1] + (g[i][j][2] * mstar1[i][j][0] - gstar[i][j][0] * m[i][j][2] );
		}
	}


	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			fourier_input[i][j] = mstar1[i][j][1];
		}
	}
	fft3d();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			mstarfour[i][j][1] = fourier_output[i][j];
			mstarfour_i[i][j][1] = fourier_output_i[i][j];
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			gstarfour[i][j][1] = (mstarfour[i][j][1] + delt*hfour[i][j][1])/(1+(xf[i]*xf[i] + yf[j]*yf[j] + 0)*Astar*delt);
			gstarfour_i[i][j][1] = (mstarfour_i[i][j][1] + delt*hfour_i[i][j][1])/(1+(xf[i]*xf[i] + yf[j]*yf[j] + 0)*Astar*delt);
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			fourier_output[i][j] = gstarfour[i][j][1];
			fourier_output_i[i][j] = gstarfour_i[i][j][1];
		}
	}
	ifft3d();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			gstar[i][j][1] = fourier_input[i][j];
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			mstar1[i][j][2] = m[i][j][2] + (gstar[i][j][0] * mstar1[i][j][1] - gstar[i][j][1] * mstar1[i][j][0] );
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			fourier_input[i][j] = mstar1[i][j][2];
		}
	}
	fft3d();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			mstarfour[i][j][2] = fourier_output[i][j];
			mstarfour_i[i][j][2] = fourier_output_i[i][j];
		}
	}


	//*********************************  STEP 2  ******************************************************

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<3;k++){
				mstar2four[i][j][k] = (mstarfour[i][j][k] + alpha*delt*hfour[i][j][k])/(1+(xf[i]*xf[i] + yf[j]*yf[j] + 0)*Astar*alpha*delt);
				mstar2four_i[i][j][k] = (mstarfour_i[i][j][k] + alpha*delt*hfour_i[i][j][k])/(1+(xf[i]*xf[i] + yf[j]*yf[j] + 0)*Astar*alpha*delt);
			}
		}
	}

	for(k=0;k<3;k++){
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				fourier_output[i][j] = mstar2four[i][j][k];
				fourier_output_i[i][j] = mstar2four_i[i][j][k];
			}
		}
		ifft3d();
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				mstar2[i][j][k] = fourier_input[i][j];
			}
		}
	}

	//*********************************  STEP 3  ******************************************************

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			mlength = sqrt( mstar2[i][j][0] * mstar2[i][j][0] + mstar2[i][j][1] * mstar2[i][j][1] + mstar2[i][j][2] * mstar2[i][j][2] );
			for(k=0;k<3;k++){
				m[i][j][k] = mstar2[i][j][k] / mlength;
			}
		}
	}


	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			mlength = sqrt( m[i][j][0] * m[i][j][0] + m[i][j][1] * m[i][j][1] + m[i][j][2] * m[i][j][2] );
			for(k=0;k<3;k++){
				m[i][j][k] = m[i][j][k] / mlength;
			}
		}
	}



	//if(keypress()){return 0;}//キー待ち状態

	time1=time1+1.0;								//計算カウント数の加算
	if(time1<time1max){goto start;}	//最大カウント数に到達したかどうかの判断

end:;
  return 0;
}

//************ 初期場の設定サブル−チン *************
void ini000()
{
	int i, j ,k;
	double mlength;
	int c;
	
	/*
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			s1h[i][j]=0.2*DRND(1.0); s2h[i][j]=0.2*DRND(1.0);//場を最大20%の乱数にて設定
			//if(abs(j-nd2)<(nd/40)){s1h[i][j]=DRND(1.0); s2h[i][j]=DRND(1.0);}
		}
	}*/

	cv::Mat_<uchar> image;
	if(SIZEX == 256){
		image = cv::imread("a.jpg" ,0);
	}else if(SIZEX == 512){
		image = cv::imread("b.jpg" ,0);
	}else{
		cout << "error : no image to load" << endl;
	}

	srand(time(NULL)); // 乱数初期化

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<3;k++){
				m[i][j][k] = rand() % 201 - 100;
				//epsilon_zero[i][j][k][k] = DRND(0.01);
			}

			c = rand() % 3;
			if(c==0){
				epsilon_zero[i][j][0][0] = -1 * DRND(0.01);
				epsilon_zero[i][j][1][1] = epsilon_zero[i][j][2][2] = -1/2 * epsilon_zero[i][j][0][0];
			}else if(c==1){
				epsilon_zero[i][j][1][1] = -1 * DRND(0.01);
				epsilon_zero[i][j][0][0] = epsilon_zero[i][j][2][2] = -1/2 * epsilon_zero[i][j][1][1];
			}else{
				epsilon_zero[i][j][2][2] = -1 * DRND(0.01);
				epsilon_zero[i][j][0][0] = epsilon_zero[i][j][1][1] = -1/2 * epsilon_zero[i][j][2][2];
			}


			//m[i][j][0] = int(image[i][j]);
			//m[i][j][1] = int(256-image[i][j]);
			//m[i][j][2] = int(100-image[i][j]/2);
			mlength = sqrt( m[i][j][0] * m[i][j][0] + m[i][j][1] * m[i][j][1] + m[i][j][2] * m[i][j][2] );
			for(k=0;k<3;k++){
				m[i][j][k] = m[i][j][k] / mlength;
			}
		}
	}
}

//******* 組織の描画サブルーチン ***************************************
void graph_s1()
{
	int i, j, ii, jj;													//整数
	double col, col_R, col_G, col_B, col_RG;	//色
	int ixmin=0, iymin=0, igx, igy, irad0;		//スクリーン座標系の設定
	double c, x, xmax, xmin, y, ymax, ymin, rad0, dia0;//規格化座標系の設定

	printf("time %f\n",time1);//計算カウント数の表示
    

	//差分ブロックの半分の長さ	//スクリーン座標系に変換（+1は整数化時の切捨て補正）
	cv::Mat chann(cv::Size(nd, nd), CV_8UC3, cv::Scalar(255, 255, 255));

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			col_R=epsilon_zero[i][j][0][0];//場の色をRGBにて設定
			col_G=epsilon_zero[i][j][1][1];
			col_B=epsilon_zero[i][j][2][2];
			col_R *= 100000;
			col_G *= 100000;
			col_B *= 100000;
			col_R += 128;
			col_G += 128;
			col_B += 128;

			chann.at<cv::Vec3b>(i,j) = cv::Vec3b(abs(int(col_B)), abs(int(col_G)), abs(int(col_R)));
		}
	}
	cv::imwrite("LLG_Terfenol_" + std::to_string(int(time1)) + "_strain_2d.png", chann);

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			col_R=eta[i][j][0][0];//場の色をRGBにて設定
			col_G=eta[i][j][1][1];
			col_B=eta[i][j][2][2];
			cout << u[i][j][1] << endl;
			col_R *= 10000000;
			col_G *= 10000000;
			col_B *= 10000000;
			col_R += 128;
			col_G += 128;
			col_B += 128;

			chann.at<cv::Vec3b>(i,j) = cv::Vec3b(abs(int(col_B)), abs(int(col_G)), abs(int(col_R)));
		}
	}
	cv::imwrite("LLG_Terfenol_" + std::to_string(int(time1)) + "_u_2d.png", chann);

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			col_R=-P[i][j][0];//場の色をRGBにて設定
			col_G=-P[i][j][1];
			col_B=-P[i][j][2];
			col_R *= 0.0000000001;
			col_G *= 0.0000000001;
			col_B *= 0.0000000001;
			col_R += 128;
			col_G += 128;
			col_B += 128;

			chann.at<cv::Vec3b>(i,j) = cv::Vec3b(abs(int(col_B)), abs(int(col_G)), abs(int(col_R)));
		}
	}
	cv::imwrite("LLG_Terfenol_" + std::to_string(int(time1)) + "_p_2d.png", chann);

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			col_R=-Pelastic[i][j][0];//場の色をRGBにて設定
			col_G=-Pelastic[i][j][1];
			col_B=-Pelastic[i][j][2];
			col_R *= 0.0000000001;
			col_G *= 0.0000000001;
			col_B *= 0.0000000001;
			col_R += 128;
			col_G += 128;
			col_B += 128;

			chann.at<cv::Vec3b>(i,j) = cv::Vec3b(abs(int(col_B)), abs(int(col_G)), abs(int(col_R)));
		}
	}
	cv::imwrite("LLG_Terfenol_" + std::to_string(int(time1)) + "_pelas_2d.png", chann);

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			col_R=m[i][j][0];//場の色をRGBにて設定
			col_G=m[i][j][1];
			col_B=m[i][j][2];
			col_R *= 100;
			col_G *= 100;
			col_B *= 100;
			col_R += 128;
			col_G += 128;
			col_B += 128;

			chann.at<cv::Vec3b>(i,j) = cv::Vec3b(abs(int(col_B)), abs(int(col_G)), abs(int(col_R)));
		}
	}
	cv::imwrite("LLG_Terfenol_" + std::to_string(int(time1)) + "_m_2d.png", chann);

	/*ofstream outputfile("check_Terfenol_" + std::to_string(int(time1)) + "_2d.txt");
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			if(i%(nd/16)==0 && j%(nd/16)==0){
				col_R=m[i][j][0];//場の色をRGBにて設定
				col_G=m[i][j][1];
				col_B=m[i][j][2];
				col_R *= 100;
				col_G *= 100;
				col_B *= 100;
				col_R += 128;
				col_G += 128;
				col_B += 128;
				outputfile <<  j << " " << ndm-i << " " <<  0.5*(nd/16)*m[i][j][1] << " " <<  0.5*(nd/16)*m[i][j][0] << endl;
			}
		}
	}
	outputfile.close();*/
}


int DCexchange2D( fftw_complex *data, int cols, int rows, int depth )
{
	int i,j,k;
	int p1,p2;    // point position
	int c2,r2,d2;    // temporary for cols/2,rows/2
	double re,im; // temporary

	if( data==NULL )       return false;
	if( rows<0 || cols<0) return false;
 
	c2 = cols/2;
	r2 = rows/2;
    //d2 = depth/2;
 

		for( j=0; j<r2; j++ ){
            for ( i=0; i<cols; i++ ){
                // exchange p1( i, j ) <-> p2( (cols/2+i)%cols, rows/2+j )
                p1 = j*cols + i;
                p2 = (r2+j)*cols + (c2+i)%cols;
                re = data[p1][0];
                im = data[p1][1];
                data[p1][0] = data[p2][0];
                data[p1][1] = data[p2][1];
                data[p2][0] = re;
                data[p2][1] = im;
            }
		}

	/*for( k=0; k<d2; k++ ){
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
	}*/

	return true;
}

int fft3d(void){

	fftw_complex *in  = NULL;
	fftw_complex *out = NULL;
	fftw_plan p       = NULL;
	int i,j,k,idx;
 
	size_t mem_size = sizeof(fftw_complex) * SIZE;
	in  = (fftw_complex*)fftw_malloc( mem_size );
	out = (fftw_complex*)fftw_malloc( mem_size );
 
	if( !in || !out ){
		fprintf( stderr, "failed to allocate %d[byte] memory(-.-)\n", (int)mem_size );
		return false;
	}
 
	// !! row-major alignment is recommended, but here, column-major.
	p = fftw_plan_dft_2d( SIZEY, SIZEX, in, out, FFTW_FORWARD, FFTW_ESTIMATE );
 
	// input data creation
	//printf("----- INPUT -----\n");
        for( j=0; j<SIZEY; j++ ){
            for( i=0; i<SIZEX; i++ ){
                idx = SIZEX*j+i; // column-major alignment
                in[idx][0] = fourier_input[i][j];
                in[idx][1] = 0;
            }
        }
 
	fftw_execute(p);
    DCexchange2D(out, SIZEX, SIZEY, SIZEZ);
 
	// output is DC exchanged and scaled.
	double scale = 1. / SIZE;
	//printf("\n----- RESULT -----\n");
	for( j=0; j<SIZEY; j++ ){
		for( i=0; i<SIZEX; i++ ){
			idx = SIZEX*j+i;
			//printf("fft :  %d %d %lf %lf\n", i, j, out[idx][0]*scale, out[idx][1]*scale );
			fourier_output[i][j] = out[idx][0];// * scale;
			fourier_output_i[i][j] = out[idx][1];// * scale;
		}
	}
 
	if( p   ) fftw_destroy_plan(p);
	if( in  ) fftw_free(in);
	if( out ) fftw_free(out);

    return true;
}


int ifft3d(void){

	fftw_complex *in2  = NULL;
	fftw_complex *out2 = NULL;
	fftw_plan ip       = NULL;
	int i,j,k,idx;
 
	size_t mem_size = sizeof(fftw_complex) * SIZE;
	in2  = (fftw_complex*)fftw_malloc( mem_size );
	out2 = (fftw_complex*)fftw_malloc( mem_size );
 
	if( !in2 || !out2 ){
		fprintf( stderr, "failed to allocate %d[byte] memory(-.-)\n", (int)mem_size );
		return false;
	}
 
	// !! row-major alignment is recommended, but here, column-major.
	ip = fftw_plan_dft_2d( SIZEY, SIZEX, in2, out2, FFTW_BACKWARD, FFTW_ESTIMATE );
 
	// input data creation
	//printf("----- INPUT -----\n");
        for( j=0; j<SIZEY; j++ ){
            for( i=0; i<SIZEX; i++ ){
                idx = SIZEX*j+i; // column-major alignment
                in2[idx][0] = fourier_output[i][j];
                in2[idx][1] = fourier_output_i[i][j];
            }
        }
 
	DCexchange2D(in2, SIZEX, SIZEY, SIZEZ);
	fftw_execute(ip);
 
	// output is DC exchanged and scaled.
	double scale = 1. / SIZE;
	//printf("\n----- RESULT -----\n");
	for( j=0; j<SIZEY; j++ ){
		for( i=0; i<SIZEX; i++ ){
			idx = SIZEX*j+i;

			//printf("ifft :  %d %d %lf %lf\n", i, j, out2[idx][0]*scale, out2[idx][1]*scale );
			fourier_input[i][j] = out2[idx][0] * scale;
		}
	}
	//cout << "out  :  " << out2[10000][0] << endl;
 
	if( ip   ) fftw_destroy_plan(ip);
	if( in2  ) fftw_free(in2);
	if( out2 ) fftw_free(out2);

    return true;
}

int convolution3D(int switch_num){
	int i, j, k, ii, jj, kk;
	double switch_cal[3][3][3];

	switch (switch_num){
	case 0:
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				for(k=0;k<SIZEZ;k++){

					for(ii=0;ii<3;ii++){
						for(jj=0;jj<3;jj++){
							for(kk=0;kk<3;kk++){
								if(i+ii-1 > ndm){
									if(jj == 1 && kk == 1){
										switch_cal[ii][jj][kk] = fai[0][j + jj-1] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(i+ii-1 < 0){
									if(jj == 1 && kk == 1){
										switch_cal[ii][jj][kk] = fai[ndm][j + jj-1] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(j+jj-1 > ndm){
									if(ii == 1 && kk == 1){
										switch_cal[ii][jj][kk] = fai[i + ii-1][0] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(j+jj-1 < 0){
									if(ii == 1 && kk == 1){
										switch_cal[ii][jj][kk] = fai[i + ii-1][ndm] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(k+kk-1 > SIZEZ-1){
									if(ii == 1 && jj == 1){
										switch_cal[ii][jj][kk] = fai[i + ii-1][j + jj-1] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(k+kk-1 < 0){
									if(ii == 1 && jj == 1){
										switch_cal[ii][jj][kk] = fai[i + ii-1][j + jj-1] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else {
									switch_cal[ii][jj][kk] = fai[i + ii-1][j + jj-1] * filter[ii][jj][kk];
								}
								//filter[ii][jj][kk] = 0;
							}
						}
					}

					filter_output[i][j] = 0;
					for(ii=0;ii<3;ii++){
						for(jj=0;jj<3;jj++){
							for(kk=0;kk<3;kk++){
								filter_output[i][j] += switch_cal[ii][jj][kk];
							}
						}
					}

				}
			}
		}
		break;
	
	case 1:
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				for(k=0;k<SIZEZ;k++){

					for(ii=0;ii<3;ii++){
						for(jj=0;jj<3;jj++){
							for(kk=0;kk<3;kk++){
								if(i+ii-1 > ndm){
									if(jj == 1 && kk == 1){
										switch_cal[ii][jj][kk] = u[0][j + jj-1][0] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(i+ii-1 < 0){
									if(jj == 1 && kk == 1){
										switch_cal[ii][jj][kk] = u[ndm][j + jj-1][0] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(j+jj-1 > ndm){
									if(ii == 1 && kk == 1){
										switch_cal[ii][jj][kk] = u[i + ii-1][0][0] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(j+jj-1 < 0){
									if(ii == 1 && kk == 1){
										switch_cal[ii][jj][kk] = u[i + ii-1][ndm][0] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(k+kk-1 > SIZEZ-1){
									if(ii == 1 && jj == 1){
										switch_cal[ii][jj][kk] = u[i + ii-1][j + jj-1][0] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(k+kk-1 < 0){
									if(ii == 1 && jj == 1){
										switch_cal[ii][jj][kk] = u[i + ii-1][j + jj-1][0] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else {
									switch_cal[ii][jj][kk] = u[i + ii-1][j + jj-1][0] * filter[ii][jj][kk];
								}
								//filter[ii][jj][kk] = 0;
							}
						}
					}

					filter_output[i][j] = 0;
					for(ii=0;ii<3;ii++){
						for(jj=0;jj<3;jj++){
							for(kk=0;kk<3;kk++){
								filter_output[i][j] += switch_cal[ii][jj][kk];
							}
						}
					}

				}
			}
		}
		break;

	case 2:
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				for(k=0;k<SIZEZ;k++){

					for(ii=0;ii<3;ii++){
						for(jj=0;jj<3;jj++){
							for(kk=0;kk<3;kk++){
								if(i+ii-1 > ndm){
									if(jj == 1 && kk == 1){
										switch_cal[ii][jj][kk] = u[0][j + jj-1][1] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(i+ii-1 < 0){
									if(jj == 1 && kk == 1){
										switch_cal[ii][jj][kk] = u[ndm][j + jj-1][1] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(j+jj-1 > ndm){
									if(ii == 1 && kk == 1){
										switch_cal[ii][jj][kk] = u[i + ii-1][0][1] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(j+jj-1 < 0){
									if(ii == 1 && kk == 1){
										switch_cal[ii][jj][kk] = u[i + ii-1][ndm][1] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(k+kk-1 > SIZEZ-1){
									if(ii == 1 && jj == 1){
										switch_cal[ii][jj][kk] = u[i + ii-1][j + jj-1][1] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(k+kk-1 < 0){
									if(ii == 1 && jj == 1){
										switch_cal[ii][jj][kk] = u[i + ii-1][j + jj-1][1] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else {
									switch_cal[ii][jj][kk] = u[i + ii-1][j + jj-1][1] * filter[ii][jj][kk];
								}
								//filter[ii][jj][kk] = 0;
							}
						}
					}

					filter_output[i][j] = 0;
					for(ii=0;ii<3;ii++){
						for(jj=0;jj<3;jj++){
							for(kk=0;kk<3;kk++){
								filter_output[i][j] += switch_cal[ii][jj][kk];
							}
						}
					}

				}
			}
		}
		break;
	case 3:
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				for(k=0;k<SIZEZ;k++){

					for(ii=0;ii<3;ii++){
						for(jj=0;jj<3;jj++){
							for(kk=0;kk<3;kk++){
								if(i+ii-1 > ndm){
									if(jj == 1 && kk == 1){
										switch_cal[ii][jj][kk] = u[0][j + jj-1][2] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(i+ii-1 < 0){
									if(jj == 1 && kk == 1){
										switch_cal[ii][jj][kk] = u[ndm][j + jj-1][2] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(j+jj-1 > ndm){
									if(ii == 1 && kk == 1){
										switch_cal[ii][jj][kk] = u[i + ii-1][0][2] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(j+jj-1 < 0){
									if(ii == 1 && kk == 1){
										switch_cal[ii][jj][kk] = u[i + ii-1][ndm][2] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(k+kk-1 > SIZEZ-1){
									if(ii == 1 && jj == 1){
										switch_cal[ii][jj][kk] = u[i + ii-1][j + jj-1][2] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(k+kk-1 < 0){
									if(ii == 1 && jj == 1){
										switch_cal[ii][jj][kk] = u[i + ii-1][j + jj-1][2] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else {
									switch_cal[ii][jj][kk] = u[i + ii-1][j + jj-1][2] * filter[ii][jj][kk];
								}
								//filter[ii][jj][kk] = 0;
							}
						}
					}

					filter_output[i][j] = 0;
					for(ii=0;ii<3;ii++){
						for(jj=0;jj<3;jj++){
							for(kk=0;kk<3;kk++){
								filter_output[i][j] += switch_cal[ii][jj][kk];
							}
						}
					}

				}
			}
		}
		break;
	default:
		break;
	}
	return 1;
}
/*
void laplacian(){
	int i, j, k;

	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			for(k=0;k<3;k++){
				filter[i][j][k] = 0;
			}
		}
	}

	filter[0][1][1] = 1;
	filter[1][0][1] = 1;
	filter[1][1][0] = 1;
	filter[1][1][1] = -6;
	filter[1][1][2] = 1;
	filter[1][2][1] = 1;
	filter[2][1][1] = 1;

	convolution3D();
}
*/
void grad_fai(){
	int i, j, k;

	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			for(k=0;k<3;k++){
				filter[i][j][k] = 0;
			}
		}
	}
	filter[0][1][1] = 0.5;
	filter[2][1][1] = -0.5;
	convolution3D(0);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			Hms[i][j][0] = filter_output[i][j];
		}
	}

	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			for(k=0;k<3;k++){
				filter[i][j][k] = 0;
			}
		}
	}
	filter[1][0][1] = 0.5;
	filter[1][2][1] = -0.5;
	convolution3D(0);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			Hms[i][j][1] = filter_output[i][j];
		}
	}

	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			for(k=0;k<3;k++){
				filter[i][j][k] = 0;
			}
		}
	}
	filter[1][1][0] = 0.5;
	filter[1][1][2] = -0.5;
	convolution3D(0);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			Hms[i][j][2] = filter_output[i][j];
		}
	}
}


void grad_u(){
	int i, j, k;

	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			for(k=0;k<3;k++){
				filter[i][j][k] = 0;
			}
		}
	}
	filter[0][1][1] = -0.5;
	filter[2][1][1] = 0.5;
	convolution3D(1);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			u_grad[i][j][0][k] = filter_output[i][j];
		}
	}

	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			for(k=0;k<3;k++){
				filter[i][j][k] = 0;
			}
		}
	}
	filter[1][0][1] = -0.5;
	filter[1][2][1] = 0.5;
	convolution3D(2);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			u_grad[i][j][1][k] = filter_output[i][j];
		}
	}

	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			for(k=0;k<3;k++){
				filter[i][j][k] = 0;
			}
		}
	}
	filter[1][1][0] = -0.5;
	filter[1][1][2] = 0.5;
	convolution3D(3);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			u_grad[i][j][2][k] = filter_output[i][j];
		}
	}
}

int four_axis(int k,int i,int j){
	if(k == 0){
		return xf[i];
	}else if(k == 1){
		return yf[j];
	}else if(k == 2){
		return 0;
	}else{
		return 0;
	}
	return 0;
}