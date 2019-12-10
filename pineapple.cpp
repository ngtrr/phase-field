//g++ -o samplepleple example.cpp `pkg-config --cflags opencv4` `pkg-config --libs opencv4`

#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <complex.h>
#include <fftw3.h>

#include <opencv2/opencv.hpp>

using namespace std;

#define DRND(x) ((double)(x)/RAND_MAX*rand())//乱数の関数設定//same

#define ND 256			//差分計算における計算領域一辺の分割数(高速フーリエ変換を用いるため２のべき乗)//same
#define IG 8				//2^IG=ND//same
#define INXY 400		//描画window１辺のピクセルサイズ(正方形の描画領域)//erace
#define SIZEX (ND)
#define SIZEY (ND)
#define SIZEZ 1
#define SIZE (SIZEX*SIZEY*SIZEZ)

	int nd=ND, ndm=ND-1; 	//計算領域の一辺の差分分割数(差分ブロック数)、ND-1を定義//same
	int nd2=ND/2;				 	//ND/2を定義：高速フ−リエ変換で使用//adopt
	int ig=IG;						//2^ig=ND//same
	double alpha=0.5;
	double rr=8.3145;			//ガス定数//adopt
	double time1;					//計算カウント数(時間に比例)//same


	double filter[3][3][3];



	//**************************	FePd	**************************************
	double Ms = 8.0E+5;
  	double K1 = -6.0E+4, K2 = 0.0E+4;
	double A = 1.0E-6;
  	double Astar;
	double myu0 = 1.0;
	double ld = 18.0E-9;
	double B = 4.0E+8;




	double s1h[ND][ND], s2h[ND][ND];		//マルテンサイトのフェーズフィールド//adopt

	double qs;					//フ−リエ変換(qs:-1)とフ−リエ逆変換(qs:1)の区別//adopt
	double xi[ND][ND], xr[ND][ND], xif[ND], xrf[ND];//フ−リエ変換の実部・虚部配列//adopt
	double s[ND],c[ND];			//sinとcosのテーブル//adopt
	int ik[ND];					//ビット反転テーブル//adopt
	double s1, s2;										//マルテンサイトのフェーズフィールド//adopt
	double ep11h0[ND][ND], ep22h0[ND][ND];				//組織内の変態歪//adopt
	double ep11qrh0[ND][ND],	ep11qih0[ND][ND];		//拘束歪変動量のフーリエ変換//adopt
	double ep22qrh0[ND][ND],	ep22qih0[ND][ND];		//拘束歪変動量のフーリエ変換//adopt
	double s1k_chem, s1k_str, s1k_su[ND][ND];			//ポテンシャル//adopt
	double s2k_chem, s2k_str, s2k_su[ND][ND];			//ポテンシャル//adopt
	double c11, c12, c44, lam0, mu0, nu0; 				//弾性定数//unify//comp
	double eta_s1[4][4], eta_s2[4][4];						//アイゲン歪成分//adopt
	double ec11[ND][ND], ec22[ND][ND];					//拘束歪変動量（実空間）//adopt
	double ep11T, ep22T;								//?
	double ep11_0, ep22_0;									//組織内の変態歪の平均値//adopt
	double ep11_a, ep22_a, ep12_a, ep21_a;				//外力に起因する歪//adopt
	double sig11_a, sig22_a;								//外力//adopt
	double Z11ep, Z12ep, Z21ep, Z22ep;					//フーリエ逆変換時の係数//?
	double sum11, sum22;										//s1とs2の空間積分
	double s1ddtt, s2ddtt;									//s1とs2の時間変化量（発展方程式の左辺）//adopt

	double temp, delt;										//計算領域、温度、時間きざみ//unify adopt same
	double time1max;													//計算カウント数の最大値（計算終了カウント）//same
	double smob;															//モビリティー（結晶変態の緩和係数）//adopt
	double nxx, nyy, nxy, alnn;								//フーリエ空間の基本ベクトルの積、ノルム//erace

	double AA0, AA1, AA2, AA3;								//ギズブエネルギー内の係数//adopt
	double kappa_s1, kappa_s2;								//勾配エネルギ−係数//adopt
	double ds_fac;										//結晶変態の揺らぎ係数//adopt






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
	double Hme[ND][ND][3];

	double fourier_output[ND][ND];
	double fourier_output_i[ND][ND];
	double fourier_input[ND][ND];

	double filter_output[ND][ND];



	void ini000();			//初期場の設定サブル−チン//unify
	void graph_s1();		//組織描画サブル−チン//unify
	void table();				//sinとcosのテーブルとビット反転テーブルの作成サブル−チン//adopt
	void fft();					//１次元高速フーリエ変換//adopt
	void rcfft();				//２次元高速フーリエ変換//adopt

	int DCexchange2D();
	int fft3d();
	int ifft3d();
	int convolution3D(int switch_num);
	void grad_fai();
	int four_axis(int k,int i,int j);

//******* メインプログラム ******************************************
int main(void)
{

	int   i, j, k, l, ii, jj, kk, ll;					//整数//unify//comp
	int   ip, im, jp, jm;							//整数//unify
	double mlength;


	srand(time(NULL));

	Astar = (2 * A)/(myu0 * Ms * Ms * ld * ld);
	//Astar = 0.0625 ;
	cout << "Astar : " << Astar << endl;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<3;k++){
				Hms[i][j][k] = 0;//init
			}
			Hexternal[i][j][0] = 0.0E+6;//ok
			Hexternal[i][j][1] = 0.0E+6;//ok
			Hexternal[i][j][2] = 0.0E+6;//ok

		}
	}

	for(i=0;i<=ndm;i++){
		xf[i] = i - nd2;
		yf[i] = i - nd2;
	}

	N[0] = 0.333;
	N[1] = 0.333;
	N[2] = 0.333;


//****** 計算条件および物質定数の設定 ****************************************

	delt=0.025;					//時間きざみ入力

	temp=300.0;						//温度(K)
	

	time1=0.0;						//初期計算カウント数の設定
	time1max=1.0+1.0e+07;	//最大計算カウント数の設定

	smob=1.0;							//モビリティー（結晶変態の緩和係数）
	ds_fac=0.01;					//結晶変態の揺らぎ係数

	AA0=1000.0/rr/temp;		//マルテンサイト変態の化学的駆動力
	AA1=1.0;  AA2=3.0*AA1+12.0;  AA3=2.0*AA1+12.0;	//ギズブエネルギー内の係数

	kappa_s1=kappa_s2=5.0e-12/rr/temp/ld/ld;				//勾配エネルギ−係数

//*** s1場のアイゲン歪の設定 ***************
	eta_s1[1][1]=0.083; eta_s1[2][2]=-0.083;
	eta_s1[3][3]=0.;
	eta_s1[1][2]=eta_s1[2][1]=eta_s1[1][3]=eta_s1[3][1]=eta_s1[2][3]=eta_s1[3][2]=0.;

//*** s2場のアイゲン歪の設定 ***************
	eta_s2[1][1]=eta_s1[2][2];
	eta_s2[2][2]=eta_s1[1][1];
	eta_s2[3][3]=0.;
	eta_s2[1][2]=eta_s2[2][1]=eta_s2[1][3]=eta_s2[3][1]=eta_s2[2][3]=eta_s2[3][2]=0.;

//***** FePdの弾性定数 ****************************
  c11=149.5;
  c44=70.5;
  c12=143.6;
  printf("c11= %f  \n", c11);
  printf("c12= %f  \n", c12);
  printf("c44= %f  \n", c44);
  //c12=c11-2.0*c44;
	lam0=c12;		mu0=c44;//cijのデーからラーメ定数を設定
	nu0=lam0/2.0/(lam0+mu0);//ポアソン比
	printf("nu0= %f  \n", nu0);//ポアソン比の値を表示（ラーメ定数の妥当性の確認のため）

//*** 外力の設定 *******************************
 	sig22_a=0.;//本計算では、外力を考慮していないので０を設定
	ep11_a=-0.5*lam0/mu0/(3.0*lam0+2.0*mu0)*sig22_a;
	ep22_a=(lam0+mu0)/mu0/(3.0*lam0+2.0*mu0)*sig22_a;
	ep12_a=ep21_a=0.0;

//*** sinおよびcosテ−ブル、ビット反転テーブル、および初期場の設定 ***************

	table();		//sinおよびcosテ−ブルとビット反転テーブルの設定

	ini000();		//初期場の設定
 	//gwinsize(INXY,INXY); ginit(2); gsetorg(0,0);//描画Window表示

//**** シミュレーションスタート ******************************
start: ;

	//if(time1<=100.){Nstep=10;} else{Nstep=200;}		//データ保存する時間間隔の変更
	//if((((int)(time1) % Nstep)==0)) {datsave();} 	//一定繰返しカウント毎に組織データを保存
	if((((int)(time1) % 100)==0)) {graph_s1();} 		//一定繰返しカウント毎に組織を表示
	//if((((int)(time1) % 100)==0)) {datsave();} 		//一定繰返しカウント毎にデータを保存

//***** 勾配ポテンシャル ***********************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1; jp=j+1; jm=j-1;
			if(i==ndm){ip=0;}  if(i==0){im=ndm;}	//周期的境界条件
			if(j==ndm){jp=0;}  if(j==0){jm=ndm;}
			s1k_su[i][j]=-kappa_s1*(s1h[ip][j]+s1h[im][j]+s1h[i][jp]+s1h[i][jm]-4.0*s1h[i][j]);			//式(4.6)
			s2k_su[i][j]=-kappa_s2*(s2h[ip][j]+s2h[im][j]+s2h[i][jp]+s2h[i][jm]-4.0*s2h[i][j]);
		}
	}

//**** アイゲン歪場[式(4.7)]のフ−リエ変換 ep11 ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			xr[i][j]=ep11h0[i][j]=eta_s1[1][1]*s1h[i][j]+eta_s2[1][1]*s2h[i][j];
			xi[i][j]=0.;
		}
	}
	qs=-1.; rcfft();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){ 
			ep11qrh0[i][j]=xr[i][j];
			ep11qih0[i][j]=xi[i][j];
		}
	}
	ep11qrh0[0][0]=ep11qih0[0][0]=0.;


//**** アイゲン歪場[式(4.7)]のフ−リエ変換 ep22 ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			xr[i][j]=ep22h0[i][j]=eta_s1[2][2]*s1h[i][j]+eta_s2[2][2]*s2h[i][j];
			xi[i][j]=0.;
		}
	}
	qs=-1.; rcfft();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){ 
			ep22qrh0[i][j]=xr[i][j];
			ep22qih0[i][j]=xi[i][j];
		}
	}
	ep22qrh0[0][0]=ep22qih0[0][0]=0.;

//*** アイゲン歪場の平均値の算出 ***
	sum11=sum22=0.;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){ sum11+=ep11h0[i][j];  sum22+=ep22h0[i][j]; }
	}
  ep11_0=sum11/nd/nd;  ep22_0=sum22/nd/nd;

//***** 拘束歪変動量ec11の計算[式(4.9)] *************************************
	for(i=0;i<=ndm;i++){
		if(i<=nd2-1){ii=i;}  if(i>=nd2){ii=i-nd;}
		for(j=0;j<=ndm;j++){
			if(j<=nd2-1){jj=j;}  if(j>=nd2){jj=j-nd;}
			alnn=sqrt((double)ii*(double)ii+(double)jj*(double)jj);
			if(alnn==0.){alnn=1.;}
			nxx=(double)ii/alnn*(double)ii/alnn;
			nyy=(double)jj/alnn*(double)jj/alnn;
			Z11ep=nxx*(2.0*(1.0-nu0)-nxx-nu0/(1.0-nu0)*nyy)/(1.0-2.0*nu0);
			Z12ep=nxx*(2.0*nu0      -nyy-nu0/(1.0-nu0)*nxx)/(1.0-2.0*nu0);
			xr[i][j]=Z11ep*ep11qrh0[i][j]+Z12ep*ep22qrh0[i][j];
			xi[i][j]=Z11ep*ep11qih0[i][j]+Z12ep*ep22qih0[i][j];
	 }
	}
	qs=1.; rcfft();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ec11[i][j]=xr[i][j];
		}
	}

//***** 拘束歪変動量ec22の計算[式(4.9)] *****************************
	for(i=0;i<=ndm;i++){
		if(i<=nd2-1){ii=i;}  if(i>=nd2){ii=i-nd;}
		for(j=0;j<=ndm;j++){
			if(j<=nd2-1){jj=j;}  if(j>=nd2){jj=j-nd;}
			alnn=sqrt((double)ii*(double)ii+(double)jj*(double)jj);
			if(alnn==0.){alnn=1.;}
			nxx=(double)ii/alnn*(double)ii/alnn;
			nyy=(double)jj/alnn*(double)jj/alnn;
			Z21ep=nyy*(2.0*nu0      -nxx-nu0/(1.0-nu0)*nyy)/(1.0-2.0*nu0);
			Z22ep=nyy*(2.0*(1.0-nu0)-nyy-nu0/(1.0-nu0)*nxx)/(1.0-2.0*nu0);
			xr[i][j]=Z21ep*ep11qrh0[i][j]+Z22ep*ep22qrh0[i][j];
			xi[i][j]=Z21ep*ep11qih0[i][j]+Z22ep*ep22qih0[i][j];
	 }
	}
	qs=1.; rcfft();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ec22[i][j]=xr[i][j];
		}
	}

//******  ポテンシャルの計算 ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){

			s1=s1h[i][j];  	s2=s2h[i][j];

//******  化学ポテンシャルの計算[式(4.4)] ********************************
			s1k_chem=AA0*s1*(AA1-AA2*s1+AA3*(s1*s1+s2*s2));
			s2k_chem=AA0*s2*(AA1-AA2*s2+AA3*(s1*s1+s2*s2));

//******  弾性ポテンシャルの計算[式(4.8)] ********************************

			ep11T=ep11h0[i][j]-ep11_0-ec11[i][j]-ep11_a;
			ep22T=ep22h0[i][j]-ep22_0-ec22[i][j]-ep22_a;

			s1k_str=ep11T*((lam0+2.0*mu0)*eta_s1[1][1]+lam0*eta_s1[2][2])
						 +ep22T*((lam0+2.0*mu0)*eta_s1[2][2]+lam0*eta_s1[1][1]);
			s2k_str=ep11T*((lam0+2.0*mu0)*eta_s2[1][1]+lam0*eta_s2[2][2])
						 +ep22T*((lam0+2.0*mu0)*eta_s2[2][2]+lam0*eta_s2[1][1]);

//****** フェーズフィールドの時間発展の計算[式(4.10)] ********************************
			s1ddtt=-smob*(s1k_chem+s1k_su[i][j]+s1k_str);
			s2ddtt=-smob*(s2k_chem+s2k_su[i][j]+s2k_str);
			s1h[i][j]=s1h[i][j]+( s1ddtt+ds_fac*(2.0*DRND(1.)-1.0) )*delt;//陽解法
			s2h[i][j]=s2h[i][j]+( s2ddtt+ds_fac*(2.0*DRND(1.)-1.0) )*delt;

//*** sの変域(0<=s<=1)の補正 ***
			if(s1h[i][j]>=1.0){s1h[i][j]=1.0;}  if(s1h[i][j]<=0.0){s1h[i][j]=0.0;}
			if(s2h[i][j]>=1.0){s2h[i][j]=1.0;}  if(s2h[i][j]<=0.0){s2h[i][j]=0.0;}
		}
	}

















	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			Hanis[i][j][0] = -1/(myu0 * Ms) * (K1*(2 * m[i][j][0] * m[i][j][1] * m[i][j][1] + 2 * m[i][j][0] * m[i][j][2] * m[i][j][2]) + 2*K2 * m[i][j][0] * m[i][j][1] * m[i][j][1] * m[i][j][2] * m[i][j][2]);//ok
			Hanis[i][j][1] = -1/(myu0 * Ms) * (K1*(2 * m[i][j][1] * m[i][j][2] * m[i][j][2] + 2 * m[i][j][1] * m[i][j][0] * m[i][j][0]) + 2*K2 * m[i][j][1] * m[i][j][2] * m[i][j][2] * m[i][j][0] * m[i][j][0]);//ok
			Hanis[i][j][2] = -1/(myu0 * Ms) * (K1*(2 * m[i][j][2] * m[i][j][0] * m[i][j][0] + 2 * m[i][j][2] * m[i][j][1] * m[i][j][1]) + 2*K2 * m[i][j][2] * m[i][j][0] * m[i][j][0] * m[i][j][1] * m[i][j][1]);//ok
		}
	}

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

	//cout << "m_ave  :   " << m_ave[0] << " : " << m_ave[1] << " : "  << m_ave[2] << endl;
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
			Hme[i][j][0] = -1/(myu0 * Ms) * 2*B*(ep11h0[i][j]*m[i][j][0]);
			Hme[i][j][1] = -1/(myu0 * Ms) * 2*B*(ep22h0[i][j]*m[i][j][1]);
			Hme[i][j][2] = -1/(myu0 * Ms) * 0;
		}
	}

	//*********************************  STEP 1  ******************************************************

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<3;k++){
				h[i][j][k] = (Hanis[i][j][k]  + Hms[i][j][k] + Hexternal[i][j][k] + Hme[i][j][k])/Ms;
				//cout << "h   " << h[i][j][k] * Ms << endl;
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

	//cout << "m  :   " << m[100][100][0]  << m[100][100][1]  << m[100][100][2] << endl;




















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

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			s1h[i][j]=0.2*DRND(1.0); s2h[i][j]=0.2*DRND(1.0);//場を最大20%の乱数にて設定
			//if(abs(j-nd2)<(nd/40)){s1h[i][j]=DRND(1.0); s2h[i][j]=DRND(1.0);}
		}
	}
	
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
	int i, j, ii, jj;													//整数
	double col, col_R, col_G, col_B, col_RG;	//色
	int ixmin=0, iymin=0, igx, igy, irad0;		//スクリーン座標系の設定
	double c, x, xmax, xmin, y, ymax, ymin, rad0, dia0;//規格化座標系の設定
	int ixmax=INXY, iymax=INXY;								//描画Window範囲

	//差分ブロックの半分の長さ	//スクリーン座標系に変換（+1は整数化時の切捨て補正）
	cv::Mat chann(cv::Size(256, 256), CV_8UC3, cv::Scalar(255, 255, 255));
	cout << "time " << (int)time1 << endl;

	for(i=0;i<=nd;i++){
		for(j=0;j<=nd;j++){
			x=rad0+dia0*i;  igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
			y=rad0+dia0*j;  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
			//座標計算			//スクリーン座標系に変換
			ii=i; jj=j; if(i==nd){ii=0;} if(j==nd){jj=0;}//周期的境界条件

			col_R=s1h[ii][jj];//場の色をRGBにて設定
			col_G=s2h[ii][jj];
			col_RG=col_R+col_G;  if(col_RG>1.){col_RG=1.;}  col_B=1.-col_RG;
			if(col_R>=0.999){col_R=1.;} if(col_R<=0.001){col_R=0.;}//RGBの変域補正
			if(col_G>=0.999){col_G=1.;} if(col_G<=0.001){col_G=0.;}
			if(col_B>=0.999){col_B=1.;} if(col_B<=0.001){col_B=0.;}
			col_R *= 255;
			col_G *= 255;
			col_B *= 255;

			chann.at<cv::Vec3b>(ii,jj) = cv::Vec3b(int(col_B), int(col_G), int(col_R));

		}
	}
	cv::imwrite("test" + std::to_string(time1) + ".png", chann);


	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			col_R = m[i][j][0];//場の色をRGBにて設定
			col_G = m[i][j][1];
			col_B = m[i][j][2];
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
			col_R = ep11h0[i][j];//場の色をRGBにて設定
			col_G = ep22h0[i][j];
			col_B = 0;
			col_R *= 100;
			col_G *= 100;
			col_B *= 100;
			col_R += 128;
			col_G += 128;
			col_B += 128;

			chann.at<cv::Vec3b>(i,j) = cv::Vec3b(abs(int(col_B)), abs(int(col_G)), abs(int(col_R)));
		}
	}
	cv::imwrite("A_Terfenol_" + std::to_string(int(time1)) + "_epsilon.png", chann);
}

//******* Sin, Cos のテーブルおよびビット反転テーブルの設定 ***************
void table()
{
	int it, it1, it2, mc, mn;
	double q;

	q=2.0*M_PI/nd;
	for(it=0;it<=nd2-1;it++){ c[it]=cos(q*it); s[it]=sin(q*it); }//Sin, Cos のテーブル

	ik[0]=0; mn=nd2; mc=1;
	for(it1=1;it1<=ig;it1++){
		for(it2=0;it2<=mc-1;it2++){
			ik[it2+mc]=ik[it2]+mn;				//ビット反転テーブル
		}
		mn=mn/2; mc=2*mc;
	}
}

//********** １次元高速フーリエ変換 **************************************
void fft()
{
	int ix, ka, kb, l2, lf, mf, n2, nf;
	double tj, tr;

	l2=1;
	for(lf=1;lf<=ig;lf++){
		n2=nd2/l2;
		for(mf=1;mf<=l2;mf++){
			for(nf=0;nf<=n2-1;nf++){
				ix=nf*l2;
				ka=nf+2*n2*(mf-1);
				kb=ka+n2;
				tr=xrf[ka]-xrf[kb];  					tj=xif[ka]-xif[kb];
				xrf[ka]=xrf[ka]+xrf[kb]; 			xif[ka]=xif[ka]+xif[kb];
				xrf[kb]=tr*c[ix]-tj*qs*s[ix];	xif[kb]=tj*c[ix]+tr*qs*s[ix];
			}
		}
		l2=l2*2;
	}

}

//************ ２次元高速フーリエ変換 ***********************************
void rcfft()
{
	int i, ic, ir, j;

	for(ir=0;ir<=ndm;ir++){
		for(ic=0;ic<=ndm;ic++){
			xrf[ic]=xr[ir][ic];	xif[ic]=xi[ir][ic];
		}
		fft();
		for(ic=0;ic<=ndm;ic++){
			xr[ir][ic]=xrf[ik[ic]];	xi[ir][ic]=xif[ik[ic]];
		}
	}
	for(ic=0;ic<=ndm;ic++){
		for(ir=0;ir<=ndm;ir++){
			xrf[ir]=xr[ir][ic];	xif[ir]=xi[ir][ic];
		}
		fft();
		for(ir=0;ir<=ndm;ir++){
			xr[ir][ic]=xrf[ik[ir]];	xi[ir][ic]=xif[ik[ir]];
		}
	}
	if(qs>0.){return;}
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			xr[i][j]=xr[i][j]/nd/nd;	xi[i][j]=xi[i][j]/nd/nd;
		}
	}

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