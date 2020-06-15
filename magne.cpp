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

	//**************************	Terfenol-D	**************************************
	double Ms = 8.0E+5;							//飽和磁化
  	double K1 = -6.0E+4, K2 = 0.0E+4;					//異方性定数
  	double ram100 = 0.0E+4, ram111 = 1.64E-3;				//磁歪定数（ラムダ）
  	double c11 = 1.41E+11, c12 = 6.48E+10, c44 = 4.87E+10;			//立方晶の弾性係数テンソル成分
	double A = 9.0E-12;							//交換定数
  	double Astar;								//交換定数を計算ステップ用の変数として変換
	double delt = 0.1;							//時間きざみステップ
	double mu0 = 1.0;							//真空の透磁率
	double ld = 1.0E-9;							//セルサイズ

	//**************************	Galfenol	**************************************
	/*double Ms = 1.432E+6;
  	double K1 = 2.0E+4, K2 = -4.5E+4;
  	double ram100 = 1.32E-4, ram111 = 0;
  	double c11 = 1.96E+11, c12 = 1.56E+11, c44 = 1.23E+11;
	double A = 1.3E-11;
  	double Astar;
	double delt = 0.1;	
	double mu0 = 1.0;
	double ld = 1.0E-06;*/

	double root3 = 1.73205080757;
	double root2 = 1.41421356237;
	double root6 = 2.44948974278;

	double xf[SIZEX];		//フーリエ空間座標
	double yf[SIZEY];		//フーリエ空間座標
	double zf[SIZEZ];		//フーリエ空間座標

	double fai[ND][ND];		//フィールド関数
	double faifour[ND][ND];		//フィールド関数のフーリエ変換（実数部分）
	double faifour_i[ND][ND];	//フィールド関数のフーリエ変換（虚数部分）

	double m_ave[3];		//磁気モーメントの成分の平均
	double m2_ave[3];		//磁気モーメントの各成分の二乗の平均
	double mm_ave[3];		//磁気モーメントの各成分をx→y、y→z、z→xと変換したときの平均

	double N[3];			//反磁場係数（形状に依存）

	double m[ND][ND][3];		//各領域の磁気モーメント
	double mstar1[ND][ND][3];	//計算ステップ用の変数（ガウスザイデル）
	double mstar2[ND][ND][3];	//計算ステップ用の変数（ガウスザイデル）
	double mfour[ND][ND][3];	//フーリエ変換（実数部分）
	double mstarfour[ND][ND][3];	//フーリエ変換（虚数部分）
	double mstar2four[ND][ND][3];	//フーリエ変換（実数部分）
	double mfour_i[ND][ND][3];	//フーリエ変換（虚数部分）
	double mstarfour_i[ND][ND][3];	//フーリエ変換（実数部分）
	double mstar2four_i[ND][ND][3];	//フーリエ変換（虚数部分）

	double g[ND][ND][3];		//計算ステップ用の変数（ガウスザイデル）
	double gstar[ND][ND][3];	//計算ステップ用の変数（ガウスザイデル）
	double gfour[ND][ND][3];	//フーリエ変換（実数部分）
	double gstarfour[ND][ND][3];	//フーリエ変換（虚数部分）
	double gfour_i[ND][ND][3];	//フーリエ変換（実数部分）
	double gstarfour_i[ND][ND][3];	//フーリエ変換（虚数部分）

	double h[ND][ND][3];		//磁界の総和（磁化モーメント基準）
	double hfour[ND][ND][3];	//磁界の総和（磁化モーメント基準）のフーリエ変換（実数部分）
	double hfour_i[ND][ND][3];	//磁界の総和（磁化モーメント基準）のフーリエ変換（虚数部分）

	double Hanis[ND][ND][3];	//結晶磁気異方性による磁界
	double Hms[ND][ND][3];		//静磁エネルギーによる磁界
	double Hexternal[ND][ND][3];	//外部磁界
	double Helastic[ND][ND][3];	//弾性エネルギーによる磁界
	double Hlandau[ND][ND][3];	//自由エネルギーによる磁界
	double Hgradient[ND][ND][3];	//勾配エネルギーによる磁界
	double Hme[ND][ND][3];		//磁気弾性エネルギーによる磁界

	double fourier_output[ND][ND];		//フーリエ変換のアウトプットの置き場所（実数部分）
	double fourier_output_i[ND][ND];	//フーリエ変換のアウトプットの置き場所（虚数部分）
	double fourier_input[ND][ND];		//フーリエ変換のインプットの置き場所

	double filter_output[ND][ND];

	double epsilon_zero[ND][ND][3][3];		//固有ひずみ（アイゲンひずみ）
	double epsilon_zerofour[ND][ND][3][3];		//固有ひずみのフーリエ変換（実数部分）
	double epsilon_zerofour_i[ND][ND][3][3];	//固有ひずみのフーリエ変換（虚数部分）

	double epsilon_homo[3][3];	//均一ひずみ

	double eta[ND][ND][3][3];	//不均一ひずみ

	double epsilon_zero_grad[ND][ND][3][3][3];	//固有ひずみの偏微分
	double epsilon_homo_grad[3][3][3];		//均一ひずみの偏微分
	double eta_grad[ND][ND][3][3][3];		//不均一ひずみの偏微分

	double c[3][3][3][3];		//弾性係数テンソル
	double s[3][3][3][3];		//弾性コンプライアンス係数テンソル
	MatrixXd c_matrix(6,6);		//弾性係数テンソルのアイゲンへの受け渡し変数
	MatrixXd s_matrix(6,6);		//弾性コンプライアンス係数テンソルのアイゲンへの受け渡し変数

	double Dfour[ND][ND];		//K行列の補因子のフーリエ変換（実数部分）
	double Dfour_i[ND][ND];		//K行列の補因子のフーリエ変換（虚数部分）
	double Nfour[ND][ND][3][3];	//K行列の行列式のフーリエ変換（実数部分）
	double Nfour_i[ND][ND][3][3];	//K行列の行列式のフーリエ変換（虚数部分）

	double u[ND][ND][3];		//変位
	double ufour[ND][ND][3];	//変位のフーリエ変換（実数部分）
	double ufour_i[ND][ND][3];	//変位のフーリエ変換（虚数部分）
	double u_grad[ND][ND][3][3];	//変位の偏微分

	double sigma_a[3][3];　　//応力

	double e1[ND][ND];	//弾性ひずみ
	double e2[ND][ND];	//弾性ひずみ
	double e3[ND][ND];	//弾性ひずみ
	double e4[ND][ND];	//弾性ひずみ
	double e5[ND][ND];	//弾性ひずみ
	double e6[ND][ND];	//弾性ひずみ

	double e1_grad[ND][ND][3];	//弾性ひずみの偏微分
	double e2_grad[ND][ND][3];	//弾性ひずみの偏微分
	double e3_grad[ND][ND][3];	//弾性ひずみの偏微分
	double e4_grad[ND][ND][3];	//弾性ひずみの偏微分
	double e5_grad[ND][ND][3];	//弾性ひずみの偏微分
	double e6_grad[ND][ND][3];	//弾性ひずみの偏微分

	double Q1 ,Q2 ,Q3 ,Q4 ,Q5;	
	double B;

	void ini000();			//初期場の設定サブル−チン
	void graph_s1();		//組織描画サブル−チン
	void graph_h();		//組織描画サブル−チン
	void graph_mstar1();		//組織描画サブル−チン
	int DCexchange2D();
	int fft3d();
	int ifft3d();
	int convolution3D(int switch_num);
	void grad_fai();
	void grad_u();
	int four_axis(int k,int i,int j);

int main(void){

	int i;		//x方向領域
	int j;		//y方向領域
	int k;		//方向の指定（0:x,1:y,2:z）
	int l;
	int ii;
	int jj;
	int kk;
	int ll;

	double mlength;

	srand(time(NULL));	//乱数初期化

	Astar = (2 * A)/(mu0 * Ms * Ms * ld * ld);
	//Astar = 0.0625 ;
	cout << "Astar : " << Astar << endl;

	for(i=0;i<=ndm;i++){			//初期設定
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

			sigma_a[0][0] = 0;
			sigma_a[1][1] = 0;
			sigma_a[2][2] = 0;

			sigma_a[0][1] = 0;
			sigma_a[0][2] = 0;
			sigma_a[1][0] = 0;
			sigma_a[1][2] = 0;
			sigma_a[2][0] = 0;
			sigma_a[2][1] = 0;

		}
	}

	for(i=0;i<=ndm;i++){	//フーリエ変換の周期的境界条件のために設定
		xf[i] = i - nd2;
		yf[i] = i - nd2;
	}

	N[0] = 0.333;		//球状の形状を仮定
	N[1] = 0.333;
	N[2] = 0.333;

	c[0][0][0][0] = c[1][1][1][1] = c[2][2][2][2] = c11;		//立方体の対称性を仮定
	c[1][2][1][2] = c[0][2][0][2] = c[0][1][0][1] = c44;
	c[0][0][1][1] = c[0][0][2][2] = c[1][1][2][2] = c[1][1][0][0] = c[2][2][0][0] = c[2][2][1][1] = c12;

	c_matrix(0,0) = c_matrix(1,1) = c_matrix(2,2) = c11;		//アイゲンに渡す
	c_matrix(3,3) = c_matrix(4,4) = c_matrix(5,5) = c44;
	c_matrix(0,1) = c_matrix(0,2) = c_matrix(1,2) = c_matrix(1,0) = c_matrix(2,0) = c_matrix(2,1) = c12;

	s_matrix = c_matrix.inverse();		//逆行列として代入

	s[0][0][0][0] = s[1][1][1][1] = s[2][2][2][2] = s_matrix(0,0);	//アイゲンから逆行列情報をもらう
	s[1][2][1][2] = s[0][2][0][2] = s[0][1][0][1] = s_matrix(3,3);
	s[0][0][1][1] = s[0][0][2][2] = s[1][1][2][2] = s[1][1][0][0] = s[2][2][0][0] = s[2][2][1][1] = s_matrix(0,1);

	//*** 初期場の設定 ***************
	ini000();		//初期場の設定


	//**** シミュレーションスタート ******************************
	start: ;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){		//各領域について磁界の成分を計算
			Hanis[i][j][0] = -1/(mu0 * Ms) * (K1*(2 * m[i][j][0] * m[i][j][1] * m[i][j][1] + 2 * m[i][j][0] * m[i][j][2] * m[i][j][2]) + 2*K2 * m[i][j][0] * m[i][j][1] * m[i][j][1] * m[i][j][2] * m[i][j][2]);//ok
			Hanis[i][j][1] = -1/(mu0 * Ms) * (K1*(2 * m[i][j][1] * m[i][j][2] * m[i][j][2] + 2 * m[i][j][1] * m[i][j][0] * m[i][j][0]) + 2*K2 * m[i][j][1] * m[i][j][2] * m[i][j][2] * m[i][j][0] * m[i][j][0]);//ok
			Hanis[i][j][2] = -1/(mu0 * Ms) * (K1*(2 * m[i][j][2] * m[i][j][0] * m[i][j][0] + 2 * m[i][j][2] * m[i][j][1] * m[i][j][1]) + 2*K2 * m[i][j][2] * m[i][j][0] * m[i][j][0] * m[i][j][1] * m[i][j][1]);//ok
		}
	}

	//if(time1<=100.){Nstep=10;} else{Nstep=200;}		//データ保存する時間間隔の変更
	//if((((int)(time1) % Nstep)==0)) {datsave();} 	//一定繰返しカウント毎に組織データを保存
	if((((int)(time1) % 100)==0)) {graph_s1();}//graph_fai();graph_h();graph_mstar1();} 		//一定繰返しカウント毎に組織を表示
	//if((((int)(time1) % 100)==0)) {datsave();} 		//一定繰返しカウント毎にデータを保存


	for(k=0;k<3;k++){
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){	//各領域の磁気モーメントの成分をフーリエ変換用の変数に渡す
				fourier_input[i][j] = m[i][j][k];
			}
		}
		fft3d();　//フーリエ変換（outputに算出）
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				mfour[i][j][k] = fourier_output[i][j];		//フーリエ変換の計算結果を保存
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
			if(xf[i]*xf[i] + yf[j]*yf[j] == 0){		//原点の場合
				faifour[i][j] = 0;
				faifour_i[i][j] = 0;
				//faifour[i][j] = Ms*(mfour_i[i][j][0]*xf[i] + mfour_i[i][j][1]*yf[j] + mfour_i[i][j][2]*0 )/1;
				//faifour_i[i][j] = -1*Ms*(mfour[i][j][0]*xf[i] + mfour[i][j][1]*yf[j] + mfour[i][j][2]*0 )/1;
			}else{						//フィールド関数のフーリエ変換結果の計算（磁気モーメントのみフーリエ変換されている）
				faifour[i][j] = Ms*(mfour_i[i][j][0]*xf[i] + mfour_i[i][j][1]*yf[j] + mfour_i[i][j][2]*0 )/(xf[i]*xf[i] + yf[j]*yf[j] + 0);
				faifour_i[i][j] = -1*Ms*(mfour[i][j][0]*xf[i] + mfour[i][j][1]*yf[j] + mfour[i][j][2]*0 )/(xf[i]*xf[i] + yf[j]*yf[j] + 0);
			}
		}
	}

	for(i=0;i<=ndm;i++){　　//フィールド関数のフーリエ変換結果をフーリエ変換用変数に渡す
		for(j=0;j<=ndm;j++){
			fourier_output[i][j] = faifour[i][j];
			fourier_output_i[i][j] = faifour_i[i][j];
		}
	}
	ifft3d();		//逆フーリエ変換
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){	//逆フーリエ変換結果をフィールド関数に代入（フィールド関数の計算終了）
			fai[i][j] = fourier_input[i][j];
		}
	}

	grad_fai();

	for(k=0;k<3;k++){
		m_ave[k] = 0;		//平均値の保存場用のクリア
		m2_ave[k] = 0;
		mm_ave[k] = 0;
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){		//各平均の算出
				m_ave[k] += m[i][j][k] / SIZE;
				m2_ave[k] += m[i][j][k] * m[i][j][k] / SIZE;
				mm_ave[k] += m[i][j][(k + 1) % 3] * m[i][j][(k + 2) % 3] / SIZE;
			}
		}
	}

	//cout << "m_ave  :   " << m_ave[0] << " : " << m_ave[1] << " : "  << m_ave[2] << endl;
	for(i=0;i<=ndm;i++){		//反磁界の追加項の計算
		for(j=0;j<=ndm;j++){
			for(k=0;k<3;k++){
				//cout << "   " << endl;
				//cout << "Hms   " << Hms[i][j][k] << endl;
				Hms[i][j][k] = -1 * m_ave[k] * Ms * N[k];
				//cout << "Hms   " << Hms[i][j][k] << endl;
			}
		}
	}
	


	for(i=0;i<=ndm;i++){		//各領域の固有ひずみの計算
		for(j=0;j<=ndm;j++){
			epsilon_zero[i][j][0][0] = 3/2 * ram100 * (m[i][j][0] * m[i][j][0] - 1/3);
			epsilon_zero[i][j][1][1] = 3/2 * ram100 * (m[i][j][1] * m[i][j][1] - 1/3);
			epsilon_zero[i][j][2][2] = 3/2 * ram100 * (m[i][j][2] * m[i][j][2] - 1/3);
			epsilon_zero[i][j][0][1] = epsilon_zero[i][j][1][0] = 3/2 * ram111 * m[i][j][0] * m[i][j][1];
			epsilon_zero[i][j][0][2] = epsilon_zero[i][j][2][0] = 3/2 * ram111 * m[i][j][0] * m[i][j][2];
			epsilon_zero[i][j][1][2] = epsilon_zero[i][j][2][1] = 3/2 * ram111 * m[i][j][1] * m[i][j][2];

					//各領域の固有ひずみの偏微分を計算
			epsilon_zero_grad[i][j][0][0][0] = 3/2 * ram100 * (2 * m[i][j][0]);
			epsilon_zero_grad[i][j][0][1][0] = epsilon_zero_grad[i][j][1][0][0] = 3/2 * ram111 * m[i][j][1];
			epsilon_zero_grad[i][j][0][2][0] = epsilon_zero_grad[i][j][2][0][0] = 3/2 * ram111 * m[i][j][2];

			epsilon_zero_grad[i][j][1][1][1] = 3/2 * ram100 * (2 * m[i][j][1]);
			epsilon_zero_grad[i][j][0][1][1] = epsilon_zero_grad[i][j][1][0][1] = 3/2 * ram111 * m[i][j][0];
			epsilon_zero_grad[i][j][1][2][1] = epsilon_zero_grad[i][j][2][1][1] = 3/2 * ram111 * m[i][j][2];

			epsilon_zero_grad[i][j][2][2][2] = 3/2 * ram100 * (2 * m[i][j][2]);
			epsilon_zero_grad[i][j][0][2][2] = epsilon_zero_grad[i][j][2][0][2] = 3/2 * ram111 * m[i][j][0];
			epsilon_zero_grad[i][j][1][2][2] = epsilon_zero_grad[i][j][2][1][2] = 3/2 * ram111 * m[i][j][1];
		}
	}
	//均一ひずみを計算
	epsilon_homo[0][0] = s[0][0][0][0] * sigma_a[0][0] + s[0][0][1][1] * (sigma_a[1][1] + sigma_a[2][2]) + 3/2 * ram100 * (m2_ave[0] - 1/3);
	epsilon_homo[1][1] = s[0][0][0][0] * sigma_a[1][1] + s[0][0][1][1] * (sigma_a[0][0] + sigma_a[2][2]) + 3/2 * ram100 * (m2_ave[1] - 1/3);
	epsilon_homo[2][2] = s[0][0][0][0] * sigma_a[2][2] + s[0][0][1][1] * (sigma_a[0][0] + sigma_a[1][1]) + 3/2 * ram100 * (m2_ave[2] - 1/3);
	epsilon_homo[0][1] = epsilon_homo[1][0] = 1/2 * s[0][1][0][1] * sigma_a[0][1] + 3/2 * ram111 * mm_ave[2];
	epsilon_homo[1][2] = epsilon_homo[2][1] = 1/2 * s[0][1][0][1] * sigma_a[1][2] + 3/2 * ram111 * mm_ave[0];
	epsilon_homo[2][0] = epsilon_homo[0][2] = 1/2 * s[0][1][0][1] * sigma_a[2][0] + 3/2 * ram111 * mm_ave[1];

	for(l=0;l<3;l++){		//固有ひずみをフーリエ変換
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


	for(i=0;i<=ndm;i++){		//行列Kに対する補因子Nと行列式Dのフーリエ変換した結果を代入
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
						for(ll=0;ll<3;ll++){		//変位のフーリエ変換した結果を代入
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

	for(k=0;k<3;k++){		//変位の逆フーリエ変換
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

	grad_u();		//変位を偏微分

	for(l=0;l<3;l++){		//不均一ひずみを計算
		for(k=0;k<3;k++){
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					eta[i][j][k][l] = 0.5 * (u_grad[i][j][k][l] + u_grad[i][j][l][k]);
				}
			}
		}
	}

	for(i=0;i<=ndm;i++){		//弾性エネルギーによる磁界を計算
		for(j=0;j<=ndm;j++){
			for(k=0;k<3;k++){
				Helastic[i][j][k] = 0;

				for(ii=0;ii<3;ii++){
					for(jj=0;jj<3;jj++){
						for(kk=0;kk<3;kk++){
							for(ll=0;ll<3;ll++){
								Helastic[i][j][k] += -1/ (mu0 * Ms) * (c[ii][jj][kk][ll] * (eta[i][j][ii][jj] - epsilon_zero[i][j][ii][jj]) * epsilon_zero_grad[i][j][ii][jj][k]);
							}
						}
					}
				}


			}
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){		//弾性ひずみを計算
			e1[i][j] = (epsilon_zero[i][j][0][0] + epsilon_zero[i][j][1][1] + epsilon_zero[i][j][2][2]) / root3;
			e2[i][j] = (epsilon_zero[i][j][0][0] - epsilon_zero[i][j][1][1]) / root2;
			e3[i][j] = (2 * epsilon_zero[i][j][2][2] - epsilon_zero[i][j][1][1] - epsilon_zero[i][j][0][0]) / root6;
			e4[i][j] = epsilon_zero[i][j][1][2];
			e5[i][j] = epsilon_zero[i][j][0][2];
			e6[i][j] = epsilon_zero[i][j][0][1];

			//弾性ひずみの偏微分を計算
			e1_grad[i][j][0] = root3 * ram100 * m[i][j][0];
			e1_grad[i][j][1] = root3 * ram100 * m[i][j][1];
			e1_grad[i][j][2] = root3 * ram100 * m[i][j][2];
			e2_grad[i][j][0] = 3 / root2 * ram100 * m[i][j][0];
			e2_grad[i][j][1] = -3 / root2 * ram100 * m[i][j][1];
			e2_grad[i][j][2] = 0;
			e3_grad[i][j][0] = -1 * root3 / root2 * ram100 * m[i][j][0];
			e3_grad[i][j][1] = -1 * root3 / root2 * ram100 * m[i][j][1];
			e3_grad[i][j][2] = root6  * ram100 * m[i][j][2];
			e4_grad[i][j][0] = epsilon_zero_grad[i][j][1][2][0];
			e4_grad[i][j][1] = epsilon_zero_grad[i][j][1][2][1];
			e4_grad[i][j][2] = epsilon_zero_grad[i][j][1][2][2];
			e5_grad[i][j][0] = epsilon_zero_grad[i][j][0][2][0];
			e5_grad[i][j][1] = epsilon_zero_grad[i][j][0][2][1];
			e5_grad[i][j][2] = epsilon_zero_grad[i][j][0][2][2];
			e6_grad[i][j][0] = epsilon_zero_grad[i][j][0][1][0];
			e6_grad[i][j][1] = epsilon_zero_grad[i][j][0][1][1];
			e6_grad[i][j][2] = epsilon_zero_grad[i][j][0][1][2];
		}
	}


	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			Hlandau[i][j][0] = 2 * Q1 * e1[i][j] * e1_grad[i][j][0] + Q5 * (2 * e5[i][j] + 2 * e6[i][j]) * e5_grad[i][j][0] * e6_grad[i][j][0];
			Hlandau[i][j][1] = 2 * Q2 * e2[i][j] * e2_grad[i][j][1] - 6 * Q3 * e3[i][j] * e2[i][j] * e2_grad[i][j][1] + 2 * Q4 * (e2[i][j] * e2[i][j] + e3[i][j] * e3[i][j]) * 2 * e2[i][j] * e2_grad[i][j][1] + Q5 * (2 * e4[i][j] + 2 * e6[i][j]) * e4_grad[i][j][1] * e6_grad[i][j][1];
			Hlandau[i][j][2] = 2 * Q2 * e3[i][j] * e3_grad[i][j][2] + 3 * Q3 * e3[i][j] * e3[i][j] * e3_grad[i][j][2] + 2 * Q4 * (e2[i][j] * e2[i][j] + e3[i][j] * e3[i][j]) * 2 * e3[i][j] * e3_grad[i][j][2] + Q5 * (2 * e4[i][j] + 2 * e5[i][j]) * e4_grad[i][j][2] * e5_grad[i][j][2];

			Hme[i][j][0] = 6 * B * ram100 * (m[i][j][0] * m[i][j][0] * m[i][j][0] - 1/3 * m[i][j][0]);
			Hme[i][j][1] = 6 * B * ram100 * (m[i][j][1] * m[i][j][1] * m[i][j][1] - 1/3 * m[i][j][1]);
			Hme[i][j][2] = 6 * B * ram100 * (m[i][j][2] * m[i][j][2] * m[i][j][2] - 1/3 * m[i][j][2]);
		}
	}


	//*********************************  STEP 1  ******************************************************

	for(i=0;i<=ndm;i++){		//磁界の総和（磁気モーメント換算）を算出
		for(j=0;j<=ndm;j++){
			for(k=0;k<3;k++){
				h[i][j][k] = (Hanis[i][j][k]  + Hms[i][j][k] + Hexternal[i][j][k] + Helastic[i][j][k] + Hlandau[i][j][k] + Hme[i][j][k])/Ms;
				//cout << "h   " << h[i][j][k] * Ms << endl;
				//cout << "Hel   " << Helastic[i][j][k] * Ms << endl;
			}
		}
	}

	//cout << "m  :   " << m[100][100][0]  << " : " << m[100][100][1]  << " : " << m[100][100][2] << endl;
	//cout << "h  :   " << h[100][100][0]  << " : " << h[100][100][1]  << " : " << h[100][100][2] << endl;

	// hfour mfour　の計算 (fft)
	for(k=0;k<3;k++){		//磁化の総和をフーリエ変換
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

	for(i=0;i<=ndm;i++){		//計算ステップ用の変数gのフーリエ変換の答えを計算
		for(j=0;j<=ndm;j++){
			for(k=0;k<3;k++){
				gfour[i][j][k] = (mfour[i][j][k] + delt*hfour[i][j][k])/(1+(xf[i]*xf[i] + yf[j]*yf[j] + 0)*Astar*delt);
				gfour_i[i][j][k] = (mfour_i[i][j][k] + delt*hfour_i[i][j][k])/(1+(xf[i]*xf[i] + yf[j]*yf[j] + 0)*Astar*delt);
			}
		}
	}
	

	for(k=0;k<3;k++){		//gの逆フーリエ変換
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



	for(i=0;i<=ndm;i++){		//m*[x]の計算（ガウスザイデル）
		for(j=0;j<=ndm;j++){
			mstar1[i][j][0] = m[i][j][0] + (g[i][j][1] * m[i][j][2] - g[i][j][2] * m[i][j][1] );
		}
	}

	for(i=0;i<=ndm;i++){		//フーリエ変換
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

	for(i=0;i<=ndm;i++){		//g*[x]のフーリエ変換の結果を代入
		for(j=0;j<=ndm;j++){
			gstarfour[i][j][0] = (mstarfour[i][j][0] + delt*hfour[i][j][0])/(1+(xf[i]*xf[i] + yf[j]*yf[j] + 0)*Astar*delt);
			gstarfour_i[i][j][0] = (mstarfour_i[i][j][0] + delt*hfour_i[i][j][0])/(1+(xf[i]*xf[i] + yf[j]*yf[j] + 0)*Astar*delt);
		}
	}


	for(i=0;i<=ndm;i++){		//g*[x]の逆フーリエ変換
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

	for(i=0;i<=ndm;i++){		//m*[y]の計算（ガウスザイデル）↓xと同様
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

	for(i=0;i<=ndm;i++){		//m*[z]の計算（ガウスザイデル）
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

	for(i=0;i<=ndm;i++){		//m**の計算
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

	for(i=0;i<=ndm;i++){		//m**の二乗平均を用いて磁気モーメントを算出
		for(j=0;j<=ndm;j++){
			mlength = sqrt( mstar2[i][j][0] * mstar2[i][j][0] + mstar2[i][j][1] * mstar2[i][j][1] + mstar2[i][j][2] * mstar2[i][j][2] );
			for(k=0;k<3;k++){
				m[i][j][k] = mstar2[i][j][k] / mlength;
			}
		}
	}


	for(i=0;i<=ndm;i++){		//磁気モーメントの成分ををm[x]あたりに変換
		for(j=0;j<=ndm;j++){
			mlength = sqrt( m[i][j][0] * m[i][j][0] + m[i][j][1] * m[i][j][1] + m[i][j][2] * m[i][j][2] );
			for(k=0;k<3;k++){
				m[i][j][k] = m[i][j][k] / mlength;
			}
		}
	}

	//cout << "m  :   " << m[100][100][0]  << m[100][100][1]  << m[100][100][2] << endl;

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
	int i, j;													//整数
	double col_R, col_G, col_B;	//色

	printf("time %d\n",int(time1));//計算カウント数の表示
	//差分ブロックの半分の長さ	//スクリーン座標系に変換（+1は整数化時の切捨て補正）
	cv::Mat chann(cv::Size(nd, nd), CV_8UC3, cv::Scalar(255, 255, 255));

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

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			col_R=u[i][j][0];//場の色をRGBにて設定
			col_G=u[i][j][1];
			col_B=u[i][j][2];
			col_R *= 1000000;
			col_G *= 1000000;
			col_B *= 1000000;
			col_R += 128;
			col_G += 128;
			col_B += 128;

			chann.at<cv::Vec3b>(i,j) = cv::Vec3b(abs(int(col_B)), abs(int(col_G)), abs(int(col_R)));
		}
	}
	cv::imwrite("LLG_Terfenol_" + std::to_string(int(time1)) + "_u_2d.png", chann);

	ofstream outputfile("check_Terfenol_" + std::to_string(int(time1)) + "_2d.txt");
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
	outputfile.close();
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