//for magenetostriction model by llg equation with phase-field simulation

#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>
#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>
#include <complex.h>
#include <fftw3.h>
#include <Eigen/Core>
#include <Eigen/LU>

using namespace std;
using namespace Eigen;

#define DRND(x) ((double)(x)/RAND_MAX*rand())//乱数の関数設定

#define ND 64			//差分計算における計算領域一辺の分割数(高速フーリエ変換を用いるため２のべき乗)
#define IG 8				//2^IG=ND
#define SIZEX (ND)
#define SIZEY (ND)
#define SIZEZ (ND)
#define SIZE (SIZEX*SIZEY*SIZEZ)

	int nd=ND, ndm=ND-1; 	//計算領域の一辺の差分分割数(差分ブロック数)、ND-1を定義
	int nd2=ND/2;				 	//ND/2を定義：高速フ−リエ変換で使用
	int ig=IG;						//2^ig=ND
	double alpha=0.5;
	double time1;					//計算カウント数(時間に比例)
	double time1max = 10000;

	double filter[3][3][3];

	double Ms = 1.432E+6;
  	double K1 = 2.0E+4, K2 = -4.5E+4;
  	double ram100 = 1.32E-4, ram111 = 0;
  	double c11 = 1.96E+11, c12 = 1.56E+11, c44 = 1.23E+11;
	double A = 1.3E-11;
  	double Astar;
	double delt = 0.1;
	double mu0 = 1.0;
	double ld = 1.0E-06;

	double xf[SIZEX];
	double yf[SIZEY];
	double zf[SIZEZ];

	double fai[ND][ND][SIZEZ];
	double faifour[ND][ND][SIZEZ];
	double faifour_i[ND][ND][SIZEZ];

	double m_ave[3];
	double m2_ave[3];
	double mm_ave[3];

	double N[3];

	double m[ND][ND][SIZEZ][3];
	double mstar1[ND][ND][SIZEZ][3];
	double mstar2[ND][ND][SIZEZ][3];
	double mfour[ND][ND][SIZEZ][3];
	double mstarfour[ND][ND][SIZEZ][3];
	double mstar2four[ND][ND][SIZEZ][3];
	double mfour_i[ND][ND][SIZEZ][3];
	double mstarfour_i[ND][ND][SIZEZ][3];
	double mstar2four_i[ND][ND][SIZEZ][3];

	double g[ND][ND][SIZEZ][3];
	double gstar[ND][ND][SIZEZ][3];
	double gfour[ND][ND][SIZEZ][3];
	double gstarfour[ND][ND][SIZEZ][3];
	double gfour_i[ND][ND][SIZEZ][3];
	double gstarfour_i[ND][ND][SIZEZ][3];

	double h[ND][ND][SIZEZ][3];
	double hfour[ND][ND][SIZEZ][3];
	double hfour_i[ND][ND][SIZEZ][3];

	double Hanis[ND][ND][SIZEZ][3];
	double Hms[ND][ND][SIZEZ][3];
	double Hexternal[ND][ND][SIZEZ][3];
	double Helastic[ND][ND][SIZEZ][3];

	double fourier_output[ND][ND][SIZEZ];
	double fourier_output_i[ND][ND][SIZEZ];
	double fourier_input[ND][ND][SIZEZ];

	double filter_output[ND][ND][SIZEZ];

	double epsilon_zero[ND][ND][SIZEZ][3][3];
	double epsilon_zerofour[ND][ND][SIZEZ][3][3];
	double epsilon_zerofour_i[ND][ND][SIZEZ][3][3];

	double epsilon_homo[3][3];

	double eta[ND][ND][SIZEZ][3][3];

	double epsilon_zero_grad[ND][ND][SIZEZ][3][3][3];
	double epsilon_homo_grad[3][3][3];
	double eta_grad[ND][ND][SIZEZ][3][3][3];

	double c[3][3][3][3];
	double s[3][3][3][3];
	MatrixXd c_matrix(6,6);
	MatrixXd s_matrix(6,6);

	double Dfour[ND][ND][SIZEZ];
	double Dfour_i[ND][ND][SIZEZ];
	double Nfour[ND][ND][SIZEZ][3][3];
	double Nfour_i[ND][ND][SIZEZ][3][3];

	double u[ND][ND][SIZEZ][3];
	double ufour[ND][ND][SIZEZ][3];
	double ufour_i[ND][ND][SIZEZ][3];
	double u_grad[ND][ND][SIZEZ][3][3];

	double sigma_a[3][3];

	void ini000();			//初期場の設定サブル−チン
	void graph_s1();		//組織描画サブル−チン
	void graph_h();		//組織描画サブル−チン
	void graph_mstar1();		//組織描画サブル−チン
	int DCexchange3D();
	int fft3d();
	int ifft3d();
	int convolution3D(int switch_num);
	void grad_fai();
	void grad_u();
	int four_axis(int v,int i,int j, int k);

int main(void){

	int i;
	int j;
	int k;
	int l;
	int v;
	int ii;
	int jj;
	int kk;
	int ll;

	double mlength;

	srand(time(NULL));

	//Astar = (2 * A)/(mu0 * Ms * Ms * ld * ld);
	Astar = 0.0625/16;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<SIZEZ;k++){
				Hexternal[i][j][k][0] = 0.0E+3;//ok
				Hexternal[i][j][k][1] = 0.0E+6;//ok
				Hexternal[i][j][k][2] = 0.0E+6;//ok


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
	}

	for(i=0;i<=ndm;i++){
		xf[i] = i - nd2;
		yf[i] = i - nd2;
	}
	for(k=0;k<SIZEZ;k++){
		zf[k] = k - SIZEZ/2;
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

	//*** 初期場の設定 ***************
	ini000();		//初期場の設定


	//**** シミュレーションスタート ******************************
	start: ;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<SIZEZ;k++){
				Hanis[i][j][k][0] = -1/(mu0 * Ms) * (K1*(2 * m[i][j][k][0] * m[i][j][k][1] * m[i][j][k][1] + 2 * m[i][j][k][0] * m[i][j][k][2] * m[i][j][k][2]) + 2*K2 * m[i][j][k][0] * m[i][j][k][1] * m[i][j][k][1] * m[i][j][k][2] * m[i][j][k][2]);//ok
				Hanis[i][j][k][1] = -1/(mu0 * Ms) * (K1*(2 * m[i][j][k][1] * m[i][j][k][2] * m[i][j][k][2] + 2 * m[i][j][k][1] * m[i][j][k][0] * m[i][j][k][0]) + 2*K2 * m[i][j][k][1] * m[i][j][k][2] * m[i][j][k][2] * m[i][j][k][0] * m[i][j][k][0]);//ok
				Hanis[i][j][k][2] = -1/(mu0 * Ms) * (K1*(2 * m[i][j][k][2] * m[i][j][k][0] * m[i][j][k][0] + 2 * m[i][j][k][2] * m[i][j][k][1] * m[i][j][k][1]) + 2*K2 * m[i][j][k][2] * m[i][j][k][0] * m[i][j][k][0] * m[i][j][k][1] * m[i][j][k][1]);//ok
			}
		}
	}

	//if(time1<=100.){Nstep=10;} else{Nstep=200;}		//データ保存する時間間隔の変更
	//if((((int)(time1) % Nstep)==0)) {datsave();} 	//一定繰返しカウント毎に組織データを保存
	if((((int)(time1) % 100 )==0)) {graph_s1();}//graph_fai();graph_h();graph_mstar1();} 		//一定繰返しカウント毎に組織を表示
	//if((((int)(time1) % 100)==0)) {datsave();} 		//一定繰返しカウント毎にデータを保存


	for(v=0;v<3;v++){
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				for(k=0;k<SIZEZ;k++){
					fourier_input[i][j][k] = m[i][j][k][v];
				}
			}
		}
		fft3d();
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				for(k=0;k<SIZEZ;k++){
					mfour[i][j][k][v] = fourier_output[i][j][k];
					mfour_i[i][j][k][v] = fourier_output_i[i][j][k];
					if (isinf(mfour[i][j][k][v]) == 1){
						cout << "mfour        " << i << " : " << j << "   -    " << dec << mfour[i][j][k][v] << endl;
					}
				}
			}
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<SIZEZ;k++){
				//cout << "faifour        " << i << " : " << j << "   -    " << dec << faifour[i][j] << endl;
				if(xf[i]*xf[i] + yf[j]*yf[j] + zf[k]*zf[k] == 0){
					faifour[i][j][k] = 0;
					faifour_i[i][j][k] = 0;
					//faifour[i][j] = Ms*(mfour_i[i][j][0]*xf[i] + mfour_i[i][j][1]*yf[j] + mfour_i[i][j][2]*0 )/1;
					//faifour_i[i][j] = -1*Ms*(mfour[i][j][0]*xf[i] + mfour[i][j][1]*yf[j] + mfour[i][j][2]*0 )/1;
				}else{
					faifour[i][j][k] = Ms*(mfour_i[i][j][k][0]*xf[i] + mfour_i[i][j][k][1]*yf[j] + mfour_i[i][j][k][2]*0 )/(xf[i]*xf[i] + yf[j]*yf[j] + zf[k]*zf[k]);
					faifour_i[i][j][k] = -1*Ms*(mfour[i][j][k][0]*xf[i] + mfour[i][j][k][1]*yf[j] + mfour[i][j][k][2]*0 )/(xf[i]*xf[i] + yf[j]*yf[j] + zf[k]*zf[k]);
				}
			}
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<SIZEZ;k++){
				fourier_output[i][j][k] = faifour[i][j][k];
				fourier_output_i[i][j][k] = faifour_i[i][j][k];
			}
		}
	}
	ifft3d();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<SIZEZ;k++){
				fai[i][j][k] = fourier_input[i][j][k];
			}
		}
	}

	grad_fai();

	for(v=0;v<3;v++){
		m_ave[v] = 0;
		m2_ave[v] = 0;
		mm_ave[v] = 0;
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				for(k=0;k<SIZEZ;k++){
					m_ave[v] += m[i][j][k][v] / SIZE;
					m2_ave[v] += m[i][j][k][v] * m[i][j][k][v] / SIZE;
					mm_ave[v] += m[i][j][k][(v + 1) % 3] * m[i][j][k][(v + 2) % 3] / SIZE;
				}
			}
		}
	}

	//cout << "m_ave  :   " << m_ave[0] << " : " << m_ave[1] << " : "  << m_ave[2] << endl;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<SIZEZ;k++){
				for(v=0;v<3;v++){
					Hms[i][j][k][v] = -1 * m_ave[v] * Ms * N[v];
				}
			}
		}
	}
	


	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<SIZEZ;k++){
				epsilon_zero[i][j][k][0][0] = 3/2 * ram100 * (m[i][j][k][0] * m[i][j][k][0] - 1/3);
				epsilon_zero[i][j][k][1][1] = 3/2 * ram100 * (m[i][j][k][1] * m[i][j][k][1] - 1/3);
				epsilon_zero[i][j][k][2][2] = 3/2 * ram100 * (m[i][j][k][2] * m[i][j][k][2] - 1/3);
				epsilon_zero[i][j][k][0][1] = epsilon_zero[i][j][k][1][0] = 3/2 * ram111 * m[i][j][k][0] * m[i][j][k][1];
				epsilon_zero[i][j][k][0][2] = epsilon_zero[i][j][k][2][0] = 3/2 * ram111 * m[i][j][k][0] * m[i][j][k][2];
				epsilon_zero[i][j][k][1][2] = epsilon_zero[i][j][k][2][1] = 3/2 * ram111 * m[i][j][k][1] * m[i][j][k][2];


				epsilon_zero_grad[i][j][k][0][0][0] = 3/2 * ram100 * (2 * m[i][j][k][0]);
				epsilon_zero_grad[i][j][k][0][1][0] = epsilon_zero_grad[i][j][k][1][0][0] = 3/2 * ram111 * m[i][j][k][1];
				epsilon_zero_grad[i][j][k][0][2][0] = epsilon_zero_grad[i][j][k][2][0][0] = 3/2 * ram111 * m[i][j][k][2];

				epsilon_zero_grad[i][j][k][1][1][1] = 3/2 * ram100 * (2 * m[i][j][k][1]);
				epsilon_zero_grad[i][j][k][0][1][1] = epsilon_zero_grad[i][j][k][1][0][1] = 3/2 * ram111 * m[i][j][k][0];
				epsilon_zero_grad[i][j][k][1][2][1] = epsilon_zero_grad[i][j][k][2][1][1] = 3/2 * ram111 * m[i][j][k][2];

				epsilon_zero_grad[i][j][k][2][2][2] = 3/2 * ram100 * (2 * m[i][j][k][2]);
				epsilon_zero_grad[i][j][k][0][2][2] = epsilon_zero_grad[i][j][k][2][0][2] = 3/2 * ram111 * m[i][j][k][0];
				epsilon_zero_grad[i][j][k][1][2][2] = epsilon_zero_grad[i][j][k][2][1][2] = 3/2 * ram111 * m[i][j][k][1];
			}
		}
	}


	epsilon_homo[0][0] = s[0][0][0][0] * sigma_a[0][0] + s[0][0][1][1] * (sigma_a[1][1] + sigma_a[2][2]) + 3/2 * ram100 * (m2_ave[0] - 1/3);
	epsilon_homo[1][1] = s[0][0][0][0] * sigma_a[1][1] + s[0][0][1][1] * (sigma_a[0][0] + sigma_a[2][2]) + 3/2 * ram100 * (m2_ave[1] - 1/3);
	epsilon_homo[2][2] = s[0][0][0][0] * sigma_a[2][2] + s[0][0][1][1] * (sigma_a[0][0] + sigma_a[1][1]) + 3/2 * ram100 * (m2_ave[2] - 1/3);
	epsilon_homo[0][1] = epsilon_homo[1][0] = 1/2 * s[0][1][0][1] * sigma_a[0][1] + 3/2 * ram111 * mm_ave[2];
	epsilon_homo[1][2] = epsilon_homo[2][1] = 1/2 * s[0][1][0][1] * sigma_a[1][2] + 3/2 * ram111 * mm_ave[0];
	epsilon_homo[2][0] = epsilon_homo[0][2] = 1/2 * s[0][1][0][1] * sigma_a[2][0] + 3/2 * ram111 * mm_ave[1];

	for(l=0;l<3;l++){
		for(v=0;v<3;v++){
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					for(k=0;k<SIZEZ;k++){
						fourier_input[i][j][k] = epsilon_zero[i][j][k][v][l];
					}
				}
			}
			fft3d();
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					for(k=0;k<SIZEZ;k++){
						epsilon_zerofour[i][j][k][v][l] = fourier_output[i][j][k];
						epsilon_zerofour_i[i][j][k][v][l] = fourier_output_i[i][j][k];
					}
				}
			}
		}
	}



	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<SIZEZ;k++){
				for(kk=0;kk<3;kk++){
					ufour[i][j][k][kk] = 0;
					ufour_i[i][j][k][kk] = 0;
					for(jj=0;jj<3;jj++){
						for(ii=0;ii<3;ii++){
							for(ll=0;ll<3;ll++){
								if(c[ii][jj][kk][ll] * four_axis(jj, i, j, k) * four_axis(ll, i, j, k) == 0 ){
									ufour[i][j][k][kk] += 0;
									ufour_i[i][j][k][kk] += 0;
								}else{
									ufour[i][j][k][kk] += 1 / (c[ii][jj][kk][ll] * four_axis(jj, i, j, k) * four_axis(ll, i, j, k)) * four_axis(jj, i, j, k) * c[ii][jj][kk][ll] * epsilon_zerofour_i[i][j][k][kk][ll];
									ufour_i[i][j][k][kk] += -1 / (c[ii][jj][kk][ll] * four_axis(jj, i, j, k) * four_axis(ll, i, j, k)) * four_axis(jj, i, j, k) * c[ii][jj][kk][ll] * epsilon_zerofour[i][j][k][kk][ll];
									//cout << "ufour  : (" << i << ", " << j << ", " << kk << ") = " << ufour[i][j][k][kk] << " : " << ufour_i[i][j][k][kk] << endl;
								}
							}
						}
					}
				}
			}
		}
	}

	for(v=0;v<3;v++){
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				for(k=0;k<SIZEZ;k++){
					fourier_output[i][j][k] = ufour[i][j][k][v];
					fourier_output_i[i][j][k] = ufour_i[i][j][k][v];
				}
			}
		}
		ifft3d();
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				for(k=0;k<SIZEZ;k++){
					u[i][j][k][v] = fourier_input[i][j][k];
					//cout << "u  : (" << i << ", " << j << ", " << k << ") =  " << u[i][j][k][v] << endl;
				}
			}
		}
	}

	grad_u();

	for(l=0;l<3;l++){
		for(v=0;v<3;v++){
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					for(k=0;k<SIZEZ;k++){
						eta[i][j][k][v][l] = 0.5 * (u_grad[i][j][k][v][l] + u_grad[i][j][k][l][v]);
					}
				}
			}
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<SIZEZ;k++){
				for(v=0;v<3;v++){
					Helastic[i][j][k][v] = 0;

					for(ii=0;ii<3;ii++){
						for(jj=0;jj<3;jj++){
							for(kk=0;kk<3;kk++){
								for(ll=0;ll<3;ll++){
									Helastic[i][j][k][v] += -1/ (mu0 * Ms) * (c[ii][jj][kk][ll] * (eta[i][j][k][ii][jj] + epsilon_homo[ii][jj] - epsilon_zero[i][j][k][ii][jj]) * epsilon_zero_grad[i][j][k][ii][jj][k]);
								}
							}
						}
					}


				}
			}
		}
	}

	//*********************************  STEP 1  ******************************************************

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<SIZEZ;k++){
				for(v=0;v<3;v++){
					h[i][j][k][v] = (Hanis[i][j][k][v]  + Hms[i][j][k][v] + Hexternal[i][j][k][v] + Helastic[i][j][k][v])/Ms;
					//cout << "h   " << h[i][j][k] * Ms << endl;
					//cout << "Hel   " << Helastic[i][j][k] * Ms << endl;
				}
			}
		}
	}


	// hfour mfour　の計算 (fft)
	for(v=0;v<3;v++){
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				for(k=0;k<SIZEZ;k++){
					fourier_input[i][j][k] = h[i][j][k][v];
				}
			}
		}
		fft3d();
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				for(k=0;k<SIZEZ;k++){
					hfour[i][j][k][v] = fourier_output[i][j][k];
					hfour_i[i][j][k][v] = fourier_output_i[i][j][k];
				}
			}
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<SIZEZ;k++){
				for(v=0;v<3;v++){
					gfour[i][j][k][v] = (mfour[i][j][k][v] + delt*hfour[i][j][k][v])/(1+(xf[i]*xf[i] + yf[j]*yf[j] + zf[k]*zf[k])*Astar*delt);
					gfour_i[i][j][k][v] = (mfour_i[i][j][k][v] + delt*hfour_i[i][j][k][v])/(1+(xf[i]*xf[i] + yf[j]*yf[j] + zf[k]*zf[k])*Astar*delt);
				}
			}
		}
	}
	

	for(v=0;v<3;v++){
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				for(k=0;k<SIZEZ;k++){
					fourier_output[i][j][k] = gfour[i][j][k][v];
					fourier_output_i[i][j][k] = gfour_i[i][j][k][v];
				}
			}
		}
		ifft3d();
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				for(k=0;k<SIZEZ;k++){
					g[i][j][k][v] = fourier_input[i][j][k];
				}
			}
		}
	}


	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<SIZEZ;k++){
				mstar1[i][j][k][0] = m[i][j][k][0] + (g[i][j][k][1] * m[i][j][k][2] - g[i][j][k][2] * m[i][j][k][1] );
			}
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<SIZEZ;k++){
				fourier_input[i][j][k] = mstar1[i][j][k][0];
			}
		}
	}
	fft3d();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<SIZEZ;k++){
				mstarfour[i][j][k][0] = fourier_output[i][j][k];
				mstarfour_i[i][j][k][0] = fourier_output_i[i][j][k];
			}
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<SIZEZ;k++){
				gstarfour[i][j][k][0] = (mstarfour[i][j][k][0] + delt*hfour[i][j][k][0])/(1+(xf[i]*xf[i] + yf[j]*yf[j] + zf[k]*zf[k])*Astar*delt);
				gstarfour_i[i][j][k][0] = (mstarfour_i[i][j][k][0] + delt*hfour_i[i][j][k][0])/(1+(xf[i]*xf[i] + yf[j]*yf[j] + zf[k]*zf[k])*Astar*delt);
			}
		}
	}


	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<SIZEZ;k++){
				fourier_output[i][j][k] = gstarfour[i][j][k][0];
				fourier_output_i[i][j][k] = gstarfour_i[i][j][k][0];
			}
		}
	}
	ifft3d();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<SIZEZ;k++){
				gstar[i][j][k][0] = fourier_input[i][j][k];
			}
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<SIZEZ;k++){
				mstar1[i][j][k][1] = m[i][j][k][1] + (g[i][j][k][2] * mstar1[i][j][k][0] - gstar[i][j][k][0] * m[i][j][k][2] );
			}
		}
	}


	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<SIZEZ;k++){
				fourier_input[i][j][k] = mstar1[i][j][k][1];
			}
		}
	}
	fft3d();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<SIZEZ;k++){
				mstarfour[i][j][k][1] = fourier_output[i][j][k];
				mstarfour_i[i][j][k][1] = fourier_output_i[i][j][k];
			}
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<SIZEZ;k++){
				gstarfour[i][j][k][1] = (mstarfour[i][j][k][1] + delt*hfour[i][j][k][1])/(1+(xf[i]*xf[i] + yf[j]*yf[j] + zf[k]*zf[k])*Astar*delt);
				gstarfour_i[i][j][k][1] = (mstarfour_i[i][j][k][1] + delt*hfour_i[i][j][k][1])/(1+(xf[i]*xf[i] + yf[j]*yf[j] + zf[k]*zf[k])*Astar*delt);
			}
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<SIZEZ;k++){
				fourier_output[i][j][k] = gstarfour[i][j][k][1];
				fourier_output_i[i][j][k] = gstarfour_i[i][j][k][1];
			}
		}
	}
	ifft3d();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<SIZEZ;k++){
				gstar[i][j][k][1] = fourier_input[i][j][k];
			}
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<SIZEZ;k++){
				mstar1[i][j][k][2] = m[i][j][k][2] + (gstar[i][j][k][0] * mstar1[i][j][k][1] - gstar[i][j][k][1] * mstar1[i][j][k][0] );
			}
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<SIZEZ;k++){
				fourier_input[i][j][k] = mstar1[i][j][k][2];
			}
		}
	}
	fft3d();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<SIZEZ;k++){
				mstarfour[i][j][k][2] = fourier_output[i][j][k];
				mstarfour_i[i][j][k][2] = fourier_output_i[i][j][k];
			}
		}
	}


	//*********************************  STEP 2  ******************************************************

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<SIZEZ;k++){
				for(v=0;v<3;v++){
					mstar2four[i][j][k][v] = (mstarfour[i][j][k][v] + alpha*delt*hfour[i][j][k][v])/(1+(xf[i]*xf[i] + yf[j]*yf[j] + zf[k]*zf[k])*Astar*alpha*delt);
					mstar2four_i[i][j][k][v] = (mstarfour_i[i][j][k][v] + alpha*delt*hfour_i[i][j][k][v])/(1+(xf[i]*xf[i] + yf[j]*yf[j] + zf[k]*zf[k])*Astar*alpha*delt);
				}
			}
		}
	}

	for(v=0;v<3;v++){
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				for(k=0;k<SIZEZ;k++){
					fourier_output[i][j][k] = mstar2four[i][j][k][v];
					fourier_output_i[i][j][k] = mstar2four_i[i][j][k][v];
				}
			}
		}
		ifft3d();
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				for(k=0;k<SIZEZ;k++){
					mstar2[i][j][k][v] = fourier_input[i][j][k];
				}
			}
		}
	}

	//*********************************  STEP 3  ******************************************************

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<SIZEZ;k++){
				mlength = sqrt( mstar2[i][j][k][0] * mstar2[i][j][k][0] + mstar2[i][j][k][1] * mstar2[i][j][k][1] + mstar2[i][j][k][2] * mstar2[i][j][k][2] );
				for(v=0;v<3;v++){
					m[i][j][k][v] = mstar2[i][j][k][v] / mlength;
				}
			}
		}
	}


	/*for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<SIZEZ;k++){
				mlength = sqrt( m[i][j][k][0] * m[i][j][k][0] + m[i][j][k][1] * m[i][j][k][1] + m[i][j][k][2] * m[i][j][k][2] );
				for(v=0;v<3;v++){
					m[i][j][k][v] = m[i][j][k][v] / mlength;
				}
			}
		}
	}*/

	time1=time1+1.0;								//計算カウント数の加算
	if(time1<time1max){goto start;}	//最大カウント数に到達したかどうかの判断

	end:;
  return 0;

}


//************ 初期場の設定サブル−チン *************
void ini000()
{
	int i, j, k, v;
	double mlength;
	
	cv::Mat_<uchar> image;
	if(SIZEX == 256){
		image = cv::imread("a.jpg" ,0);
	}else if(SIZEX == 512){
		image = cv::imread("c.jpg" ,0);
	}else{
		cout << "error : no image to load" << endl;
	}

	srand(time(NULL)); // 乱数初期化

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<SIZEZ;k++){
				for(v=0;v<3;v++){
					m[i][j][k][v] = rand();
				}
				//m[i][j][0] = int(image[i][j]);
				//m[i][j][1] = int(256-image[i][j]);
				//m[i][j][2] = int(100-image[i][j]/2);
				mlength = sqrt( m[i][j][k][0] * m[i][j][k][0] + m[i][j][k][1] * m[i][j][k][1] + m[i][j][k][2] * m[i][j][k][2] );
				for(v=0;v<3;v++){
					m[i][j][k][v] = m[i][j][k][v] / mlength;
				}
			}
		}
	}
}

//******* 組織の描画サブルーチン ***************************************
void graph_s1()
{
	int i, j, k;													//整数
	double col_R, col_G, col_B;	//色

	printf("time %d\n",int(time1));//計算カウント数の表示
	//差分ブロックの半分の長さ	//スクリーン座標系に変換（+1は整数化時の切捨て補正）
	cv::Mat chann(cv::Size(nd, nd), CV_8UC3, cv::Scalar(255, 255, 255));

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			col_R=m[i][j][0][0];//場の色をRGBにて設定
			col_G=m[i][j][0][1];
			col_B=m[i][j][0][2];
			col_R *= 100;
			col_G *= 100;
			col_B *= 100;
			col_R += 128;
			col_G += 128;
			col_B += 128;

			chann.at<cv::Vec3b>(i,j) = cv::Vec3b(abs(int(col_B)), abs(int(col_G)), abs(int(col_R)));
		}
	}
	cv::imwrite("LLG_permalloy_" + std::to_string(int(time1)) + "_m_xy.png", chann);

	/*for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			col_R=m[0][j][i][0];//場の色をRGBにて設定
			col_G=m[0][j][i][1];
			col_B=m[0][j][i][2];
			col_R *= 255;
			col_G *= 255;
			col_B *= 255;
			//col_R += 128;
			//col_G += 128;
			//col_B += 128;

			chann.at<cv::Vec3b>(i,j) = cv::Vec3b(abs(int(col_B)), abs(int(col_G)), abs(int(col_R)));
		}
	}
	cv::imwrite("LLG_permalloy_" + std::to_string(int(time1)) + "_m_zy.png", chann);

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			col_R=m[i][0][j][0];//場の色をRGBにて設定
			col_G=m[i][0][j][1];
			col_B=m[i][0][j][2];
			col_R *= 255;
			col_G *= 255;
			col_B *= 255;
			//col_R += 128;
			//col_G += 128;
			//col_B += 128;

			chann.at<cv::Vec3b>(i,j) = cv::Vec3b(abs(int(col_B)), abs(int(col_G)), abs(int(col_R)));
		}
	}
	cv::imwrite("LLG_permalloy_" + std::to_string(int(time1)) + "_m_xz.png", chann);*/

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			col_R=u[i][j][0][0];//場の色をRGBにて設定
			col_G=u[i][j][0][1];
			col_B=u[i][j][0][2];
			col_R *= 100000;
			col_G *= 100000;
			col_B *= 100000;
			col_R += 128;
			col_G += 128;
			col_B += 128;

			chann.at<cv::Vec3b>(i,j) = cv::Vec3b(abs(int(col_B)), abs(int(col_G)), abs(int(col_R)));
		}
	}
	cv::imwrite("LLG_permalloy_" + std::to_string(int(time1)) + "_u.png", chann);
}

int DCexchange3D( fftw_complex *data, int cols, int rows, int depth )
{
	int i,j,k;
	int p1,p2;    // point position
	int c2,r2,d2;    // temporary for cols/2,rows/2
	double re,im; // temporary

	if( data==NULL )       return false;
	if( rows<0 || cols<0) return false;
	if( depth<0) return false;
 
	c2 = cols/2;
	r2 = rows/2;
    d2 = depth/2;
 

	for( k=0; k<depth; k++ ){
		for( j=0; j<r2; j++ ){
            for ( i=0; i<cols; i++ ){
                // exchange p1( i, j ) <-> p2( (cols/2+i)%cols, rows/2+j )
                p1 = j*cols + i + k*rows*cols;;
                p2 = (r2+j)*cols + (c2+i)%cols + k*rows*cols;;
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
	p = fftw_plan_dft_3d( SIZEY, SIZEX, SIZEZ, in, out, FFTW_FORWARD, FFTW_ESTIMATE );
 
	// input data creation
    for( k=0; k<SIZEZ; k++ ){
        for( j=0; j<SIZEY; j++ ){
            for( i=0; i<SIZEX; i++ ){
                idx = SIZEX*SIZEY*k + SIZEX*j+i; // column-major alignment
                in[idx][0] = fourier_input[i][j][k];
                in[idx][1] = 0;
            }
		}
	}
 
	fftw_execute(p);
    DCexchange3D(out, SIZEX, SIZEY, SIZEZ);
 
	// output is DC exchanged and scaled.
	double scale = 1. / SIZE;
    for( k=0; k<SIZEZ; k++ ){
		for( j=0; j<SIZEY; j++ ){
			for( i=0; i<SIZEX; i++ ){
				idx = SIZEX*SIZEY*k + SIZEX*j+i;
				//printf("fft :  %d %d %lf %lf\n", i, j, out[idx][0]*scale, out[idx][1]*scale );
				fourier_output[i][j][k] = out[idx][0];// * scale;
				fourier_output_i[i][j][k] = out[idx][1];// * scale;
			}
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
	ip = fftw_plan_dft_3d( SIZEY, SIZEX, SIZEZ, in2, out2, FFTW_BACKWARD, FFTW_ESTIMATE );
 
	// input data creation
    for( k=0; k<SIZEZ; k++ ){
        for( j=0; j<SIZEY; j++ ){
            for( i=0; i<SIZEX; i++ ){
                idx = SIZEX*SIZEY*k + SIZEX*j+i; // column-major alignment
                in2[idx][0] = fourier_output[i][j][k];
                in2[idx][1] = fourier_output_i[i][j][k];
            }
        }
	}

	DCexchange3D(in2, SIZEX, SIZEY, SIZEZ);
	fftw_execute(ip);
 
	// output is DC exchanged and scaled.
	double scale = 1. / SIZE;
    for( k=0; k<SIZEZ; k++ ){
		for( j=0; j<SIZEY; j++ ){
			for( i=0; i<SIZEX; i++ ){
				idx = SIZEX*SIZEY*k + SIZEX*j+i;

				//printf("ifft :  %d %d %lf %lf\n", i, j, out2[idx][0]*scale, out2[idx][1]*scale );
				fourier_input[i][j][k] = out2[idx][0] * scale;
			}
		}
	}
 
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
										switch_cal[ii][jj][kk] = fai[0][j + jj-1][k + kk-1] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(i+ii-1 < 0){
									if(jj == 1 && kk == 1){
										switch_cal[ii][jj][kk] = fai[ndm][j + jj-1][k + kk-1] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(j+jj-1 > ndm){
									if(ii == 1 && kk == 1){
										switch_cal[ii][jj][kk] = fai[i + ii-1][0][k + kk-1] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(j+jj-1 < 0){
									if(ii == 1 && kk == 1){
										switch_cal[ii][jj][kk] = fai[i + ii-1][ndm][k + kk-1] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(k+kk-1 > SIZEZ-1){
									if(ii == 1 && jj == 1){
										switch_cal[ii][jj][kk] = fai[i + ii-1][j + jj-1][0] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(k+kk-1 < 0){
									if(ii == 1 && jj == 1){
										switch_cal[ii][jj][kk] = fai[i + ii-1][j + jj-1][ndm] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else {
									switch_cal[ii][jj][kk] = fai[i + ii-1][j + jj-1][k + kk-1] * filter[ii][jj][kk];
								}
								//filter[ii][jj][kk] = 0;
							}
						}
					}

					filter_output[i][j][k] = 0;
					for(ii=0;ii<3;ii++){
						for(jj=0;jj<3;jj++){
							for(kk=0;kk<3;kk++){
								filter_output[i][j][k] += switch_cal[ii][jj][kk];
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
										switch_cal[ii][jj][kk] = u[0][j + jj-1][k + kk-1][0] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(i+ii-1 < 0){
									if(jj == 1 && kk == 1){
										switch_cal[ii][jj][kk] = u[ndm][j + jj-1][k + kk-1][0] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(j+jj-1 > ndm){
									if(ii == 1 && kk == 1){
										switch_cal[ii][jj][kk] = u[i + ii-1][0][k + kk-1][0] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(j+jj-1 < 0){
									if(ii == 1 && kk == 1){
										switch_cal[ii][jj][kk] = u[i + ii-1][ndm][k + kk-1][0] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(k+kk-1 > SIZEZ-1){
									if(ii == 1 && jj == 1){
										switch_cal[ii][jj][kk] = u[i + ii-1][j + jj-1][0][0] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(k+kk-1 < 0){
									if(ii == 1 && jj == 1){
										switch_cal[ii][jj][kk] = u[i + ii-1][j + jj-1][ndm][0] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else {
									switch_cal[ii][jj][kk] = u[i + ii-1][j + jj-1][k + kk-1][0] * filter[ii][jj][kk];
								}
								//filter[ii][jj][kk] = 0;
							}
						}
					}

					filter_output[i][j][k] = 0;
					for(ii=0;ii<3;ii++){
						for(jj=0;jj<3;jj++){
							for(kk=0;kk<3;kk++){
								filter_output[i][j][k] += switch_cal[ii][jj][kk];
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
										switch_cal[ii][jj][kk] = u[0][j + jj-1][k + kk-1][1] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(i+ii-1 < 0){
									if(jj == 1 && kk == 1){
										switch_cal[ii][jj][kk] = u[ndm][j + jj-1][k + kk-1][1] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(j+jj-1 > ndm){
									if(ii == 1 && kk == 1){
										switch_cal[ii][jj][kk] = u[i + ii-1][0][k + kk-1][1] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(j+jj-1 < 0){
									if(ii == 1 && kk == 1){
										switch_cal[ii][jj][kk] = u[i + ii-1][ndm][k + kk-1][1] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(k+kk-1 > SIZEZ-1){
									if(ii == 1 && jj == 1){
										switch_cal[ii][jj][kk] = u[i + ii-1][j + jj-1][0][1] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(k+kk-1 < 0){
									if(ii == 1 && jj == 1){
										switch_cal[ii][jj][kk] = u[i + ii-1][j + jj-1][ndm][1] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else {
									switch_cal[ii][jj][kk] = u[i + ii-1][j + jj-1][k + kk-1][1] * filter[ii][jj][kk];
								}
								//filter[ii][jj][kk] = 0;
							}
						}
					}

					filter_output[i][j][k] = 0;
					for(ii=0;ii<3;ii++){
						for(jj=0;jj<3;jj++){
							for(kk=0;kk<3;kk++){
								filter_output[i][j][k] += switch_cal[ii][jj][kk];
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
										switch_cal[ii][jj][kk] = u[0][j + jj-1][k + kk-1][2] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(i+ii-1 < 0){
									if(jj == 1 && kk == 1){
										switch_cal[ii][jj][kk] = u[ndm][j + jj-1][k + kk-1][2] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(j+jj-1 > ndm){
									if(ii == 1 && kk == 1){
										switch_cal[ii][jj][kk] = u[i + ii-1][0][k + kk-1][2] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(j+jj-1 < 0){
									if(ii == 1 && kk == 1){
										switch_cal[ii][jj][kk] = u[i + ii-1][ndm][k + kk-1][2] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(k+kk-1 > SIZEZ-1){
									if(ii == 1 && jj == 1){
										switch_cal[ii][jj][kk] = u[i + ii-1][j + jj-1][0][2] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else if(k+kk-1 < 0){
									if(ii == 1 && jj == 1){
										switch_cal[ii][jj][kk] = u[i + ii-1][j + jj-1][ndm][2] * filter[ii][jj][kk];
									}else{
										switch_cal[ii][jj][kk] = 0;
									}
								}else {
									switch_cal[ii][jj][kk] = u[i + ii-1][j + jj-1][k + kk-1][2] * filter[ii][jj][kk];
								}
								//filter[ii][jj][kk] = 0;
							}
						}
					}

					filter_output[i][j][k] = 0;
					for(ii=0;ii<3;ii++){
						for(jj=0;jj<3;jj++){
							for(kk=0;kk<3;kk++){
								filter_output[i][j][k] += switch_cal[ii][jj][kk];
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
			for(k=0;k<SIZEZ;k++){
				Hms[i][j][k][0] = filter_output[i][j][k];
			}
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
			for(k=0;k<SIZEZ;k++){
				Hms[i][j][k][1] = filter_output[i][j][k];
			}
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
			for(k=0;k<SIZEZ;k++){
				Hms[i][j][k][2] = filter_output[i][j][k];
			}
		}
	}
}


void grad_u(){
	int i, j, k, v;

	for(v=0;v<3;v++){
		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
				for(k=0;k<3;k++){
					filter[i][j][k] = 0;
				}
			}
		}
		filter[0][1][1] = -0.5;
		filter[2][1][1] = 0.5;
		convolution3D(v+1);
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				for(k=0;k<SIZEZ;k++){
					u_grad[i][j][k][0][v] = filter_output[i][j][k];
				}
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
		convolution3D(v+1);
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				for(k=0;k<SIZEZ;k++){
					u_grad[i][j][k][1][v] = filter_output[i][j][k];
				}
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
				for(k=0;k<SIZEZ;k++){
					u_grad[i][j][k][2][k] = filter_output[i][j][k];
				}
			}
		}
	}
}

int four_axis(int v,int i,int j, int k){
	if(v == 0){
		return xf[i];
	}else if(v == 1){
		return yf[j];
	}else if(v == 2){
		return zf[k];
	}else{
		return 0;
	}
	return 0;
}