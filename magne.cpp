//for magenetostriction model by llg equation with phase-field simulation


#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>
#include <iostream>
#include <opencv2/opencv.hpp>


#define DRND(x) ((double)(x)/RAND_MAX*rand())//乱数の関数設定

#define ND 256			//差分計算における計算領域一辺の分割数(高速フーリエ変換を用いるため２のべき乗)
#define IG 8				//2^IG=ND
#define INXY 400		//描画window１辺のピクセルサイズ(正方形の描画領域)

	int nd=ND, ndm=ND-1; 	//計算領域の一辺の差分分割数(差分ブロック数)、ND-1を定義
	int nd2=ND/2;				 	//ND/2を定義：高速フ−リエ変換で使用
	int ig=IG;						//2^ig=ND
	double PI=3.14159;		//円周率
	double rr=8.3145;			//ガス定数
	double time1;					//計算カウント数(時間に比例)


    double Ms = 1.432E+6;
    double K1 = 2E+4, K2 = -4.5E+4;
    double ram100 = 2.64E-4, ram111 = 0;
    double c11 = 1.96E+11, c12 = 1.56E+11, c44 = 1.23E+11;
    double Astar = 0.0625;
		double delt = 0.1;


	void ini000();			//初期場の設定サブル−チン
	void graph_s1();		//組織描画サブル−チン
	void table();				//sinとcosのテーブルとビット反転テーブルの作成サブル−チン
	void fft();					//１次元高速フーリエ変換
	void rcfft();				//２次元高速フーリエ変換

int main(void){
  double E, Eanis, Eexch, Ems, Eexternal, Eelastic;
  double m[ND][ND][3];
	double 

	//*** sinおよびcosテ−ブル、ビット反転テーブル、および初期場の設定 ***************
	table();		//sinおよびcosテ−ブルとビット反転テーブルの設定
	ini000();		//初期場の設定

	//**** シミュレーションスタート ******************************
	start: ;

	//if(time1<=100.){Nstep=10;} else{Nstep=200;}		//データ保存する時間間隔の変更
	//if((((int)(time1) % Nstep)==0)) {datsave();} 	//一定繰返しカウント毎に組織データを保存
	if((((int)(time1) % 10)==0)) {graph_s1();} 		//一定繰返しカウント毎に組織を表示
	//if((((int)(time1) % 100)==0)) {datsave();} 		//一定繰返しカウント毎にデータを保存





























	time1=time1+1.0;								//計算カウント数の加算
	if(time1<time1max){goto start;}	//最大カウント数に到達したかどうかの判断

	end:;
  return 0;

}



//************ 初期場の設定サブル−チン *************
void ini000()
{
	int i, j;
  //srand(time(NULL)); // 乱数初期化

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			s1h[i][j]=0.2*DRND(1.0); s2h[i][j]=0.2*DRND(1.0);//場を最大20%の乱数にて設定
			//if(abs(j-nd2)<(nd/40)){s1h[i][j]=DRND(1.0); s2h[i][j]=DRND(1.0);}
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

	//gcls(); //画面クリア
	xmin=0.; xmax=1.; ymin=0.; ymax=1.;//描画領域（規格化されている）

	printf("time %f\n",time1);//計算カウント数の表示
	dia0=1.0/nd;
	rad0=dia0/2.0;   						irad0=(ixmax-ixmin)/(xmax-xmin)*rad0+1;
	//差分ブロックの半分の長さ	//スクリーン座標系に変換（+1は整数化時の切捨て補正）
	cv::Mat chann(cv::Size(256, 256), CV_8UC3, cv::Scalar(255, 255, 255));

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
}


//******* Sin, Cos のテーブルおよびビット反転テーブルの設定 ***************
void table()
{
	int it, it1, it2, mc, mn;
	double q;

	q=2.0*PI/nd;
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

void laplacian(int i, int j){

}



