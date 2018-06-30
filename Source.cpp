//------------------------------------------------------------
// 標準的な差分進化のアルゴリズムを実装したプログラム
//------------------------------------------------------------
//------------------------------------------------------------
//「error C2065: 'M_PI' : 定義されていない識別子です」への対応
//　VisualStudioの場合は以下のコメントを外してコンパイル
//------------------------------------------------------------
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
//#include <fstream>
//#include <random>
//#include "random.h"

using namespace std;

// 以下はメルセンヌツイスターで使用する
#include "MT.h"
//------------------------------------------------------------
// メルセンヌツイスターで使用する関数の概説
// genrand_int32() //符号なし32ビット長整数
// genrand_int31() //符号なし31ビット長整数
// genrand_real1() //一様実乱数[0,1] (32ビット精度)
// genrand_real2() //一様実乱数[0,1) (32ビット精度)
// genrand_real3() //一様実乱数(0,1) (32ビット精度)
// genrand_res53() //一様実乱数[0,1) (53ビット精度)
//------------------------------------------------------------
//------------------------------------------------------------
// VisualStudioの場合は以下のコメントを外してコンパイル
//------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE
#pragma warning(disable:4996)
//------------------------------------------------------------
// 設定パラメータ 実験時にこの部分を書き換える
//------------------------------------------------------------
#define MAXPSIZE		30			//最大個体数
#define MAXGSIZE		100			//最大次元数（問題の次元数）
#define FUNCNO		1					//最適化問題の種類
#define RANGE		100				//最適化問題の定義域
#define MRATE			0.9				//突然変異率
#define CRATE			0.6				//交叉率
#define MAXGENERATIO	1500	//最大繰り返し回数
#define DE_ALGORITHM_NO	5			//DEのアルゴリズム
#define EXTIME			1				//試行回数
#define TERMINATE		1.0e-70	//終了条件
#define M_PI	3.14159265359
//------------------------------------------------------------
double cVect[MAXPSIZE][MAXGSIZE];				//時刻Tの個体（ベクトル）
double cFitness[MAXPSIZE];							//時刻Tの個体の評価値
double pBestVector[MAXPSIZE][MAXGSIZE];	//時刻Tにおける各個体の最良解の履歴
double pBestFitness[MAXPSIZE];			//時刻Tにおける各個体の最良解評価値の履歴
double gBestVector[MAXGSIZE];			//時刻Tにおける集団全体での最良解の履歴
double gBestFitness;							//時刻Tにおける集団全体での最良解評価値の履歴
double pVect1[MAXGSIZE];					//親ベクトル１
double pVect2[MAXGSIZE];					//親ベクトル２
double pVect3[MAXGSIZE];					//親ベクトル３
double nVect[MAXPSIZE][MAXGSIZE];		//時刻T+1の個体（ベクトル）
double nFitness[MAXPSIZE];					//時刻T+1の個体の評価値
double gBestHistory[EXTIME][MAXGENERATIO];	//各試行におけるgBestの履歴
double pDiversity[EXTIME][MAXGENERATIO];		//各試行における集団の多様性
double Mrate[MAXPSIZE], cMrate = 0.0;
double Crate[MAXPSIZE], cCrate = 0.0;
double C = 0.1;	//JADE のパラメータ

int gTable[EXTIME];					//最適解発見時の世代数
int sRate[EXTIME];					//最適解の発見率
double gBestTable[EXTIME];		//終了時点でのgBestの値
//------------------------------------------------------------
//------------------------------------------------------------
// [min, max]の一様乱数
//------------------------------------------------------------
double uniform(double min, double max) {
	return min + (max - min)*genrand_real1();
}
//------------------------------------------------------------
// テスト関数
// F1 Sphere関数
// [-5.12, +5.12] x=(0,..,0) F1=0 
//------------------------------------------------------------
double Sphere(double *x) {
	int i;
	double sum = 0.0;
	for (i = 0; i<MAXGSIZE; i++) {
		sum += x[i] * x[i];
	}
	return sum;
}
//------------------------------------------------------------
// F2 Rosenbrock関数 Chain型
// [-2.048, +2.048] x=(1,..,1) F2=0 
//------------------------------------------------------------
double Rosenbrock(double *x) {
	int i;
	double sum = 0.0;
	/*
	for (i = 0; i<MAXGSIZE - 1; i++){
	sum += 100 * (x[i + 1] - x[i] * x[i])*(x[i + 1] - x[i] * x[i]) + (x[i] - 1.0)*(x[i] - 1.0);
	}
	*/
	// UNDXの文献の式に統一
	// Rosenbrock関数 Star型
	for (i = 2; i<MAXGSIZE; i++) {
		sum += (100 * (x[0] - x[i] * x[i])*(x[0] - x[i] * x[i]) + (x[1] - 1.0)*(x[1] - 1.0));
	}
	return sum;
}
//------------------------------------------------------------
// F3 Rastrigin関数
//------------------------------------------------------------
double Rastrigin(double *x) {
	int i;
	double sum1 = 0, sum2 = 0;
	for (i = 0; i<MAXGSIZE; i++) {
		sum1 += x[i] * x[i] - 10 * cos(2 * M_PI*x[i]);
	}
	return (10 * MAXGSIZE + sum1);
}
//------------------------------------------------------------
// F4 Griewank関数
//------------------------------------------------------------
double Griewank(double *x) {
	int i;
	double sum1 = 0.0, sum2 = 1.0;
	for (i = 0; i<MAXGSIZE; i++) {
		sum1 += x[i] * x[i];
		sum2 *= cos(x[i] / sqrt((i + 1)*1.0));
	}
	return ((sum1 / 4000) - sum2 + 1);
}
//------------------------------------------------------------
// F5 Ackley関数
//------------------------------------------------------------
double Ackley(double *x) {
	int i;
	double sum1, sum2;
	for (i = 0, sum1 = 0.0; i<MAXGSIZE; i++)sum1 += x[i] * x[i];
	for (i = 0, sum2 = 0.0; i<MAXGSIZE; i++)sum2 += cos(2.0*M_PI*MAXGSIZE);
	return (20.0 + M_E - 20.0 * exp(-0.2*sqrt(sum1 / MAXGSIZE)) - exp(sum2 / MAXGSIZE));
}
//------------------------------------------------------------
// F6 Schwefel関数
//------------------------------------------------------------
double Schwefel(double *x) {
	int i;
	double sum;
	for (i = 0, sum = 0.0; i<MAXGSIZE; i++) {
		sum += x[i] * sin(sqrt(fabs(x[i])));
	}
	return (sum + 418.9828872724338 * MAXGSIZE);
}
//------------------------------------------------------------
// F7 Ridge関数
//------------------------------------------------------------
double Ridge(double *x) {
	int i, j;
	double sum1, sum2;
	for (i = 0, sum1 = 0; i<MAXGSIZE; i++) {
		for (j = 0, sum2 = 0; j <= i; j++) {
			sum2 += x[j];
		}
		sum1 += sum2 * sum2;
	}
	return sum1;
}
//------------------------------------------------------------
// F8 Bohachevsky関数
//------------------------------------------------------------
double Bohachevsky(double *x) {
	int i;
	double sum;
	for (i = 0, sum = 0; i<MAXGSIZE - 1; i++) {
		sum += x[i] * x[i] + 2 * x[i + 1] * x[i + 1] - 0.3*cos(3 * M_PI*x[i]) - 0.4*cos(4 * M_PI*x[i + 1]) + 0.7;
	}
	return sum;
}
//------------------------------------------------------------
// F9 Schaffer関数
//------------------------------------------------------------
double Schaffer(double *x) {
	int i;
	double tp1 = 0, tp2 = 0, tp3 = 0, tp4 = 0, tp5 = 0;
	for (i = 0, tp5 = 0; i<MAXGSIZE - 1; i++) {
		tp1 = x[i] * x[i] + x[i + 1] * x[i + 1];
		tp2 = pow(tp1, 0.25);
		tp3 = pow(tp1, 0.1);
		tp4 = sin(50 * tp3) * sin(50 * tp3);
		tp5 += tp2 * (tp4 + 1.0);
	}
	return tp5;
}
//------------------------------------------------------------
// F10 Ellipsoid関数
//------------------------------------------------------------
double Ellipsoid(double *x) {
	int i;
	double sum;
	for (i = 0, sum = 0; i<MAXGSIZE; i++) {
		sum += (pow(1000, (double)i / (double)(MAXGSIZE - 1)) * x[i]) * (pow(1000, (double)i / (double)(MAXGSIZE - 1)) * x[i]);
	}
	return sum;
}
//------------------------------------------------------------
// F11 k-tablet関数
//------------------------------------------------------------
double K_Tablet(double *x) {
	int i;
	double sum1, sum2;
	for (i = 0, sum1 = 0; i<(int)MAXGSIZE / 2; i++) {
		sum1 += x[i] * x[i];
	}
	for (i = (int)MAXGSIZE / 2 + 1, sum2 = 0; i<MAXGSIZE; i++) {
		sum2 += (100 * x[i]) * (100 * x[i]);
	}
	return (sum1 + sum2);
}
//------------------------------------------------------------
// F12 Shifted-Rastrigin関数
//------------------------------------------------------------
double Shifted_Rastrigin(double *x) {
	int i;
	double sum1 = 0.0, sum2 = 0.0;
	for (i = 0; i<MAXGSIZE; i++) sum1 += (x[i] - 1) * (x[i] - 1) - 10 * cos(2 * M_PI*(x[i] - 1));
	sum2 = 10 * MAXGSIZE + sum1;
	return sum2;
}
//------------------------------------------------------------
// F13 Cigar関数
//------------------------------------------------------------
double Cigar(double *x) {
	int i;
	double sum;
	for (i = 1, sum = 0; i<MAXGSIZE; i++) {
		sum += (1000 * x[i]) * (1000 * x[i]);
	}
	return (x[0] * x[0] + sum);
}
//------------------------------------------------------------
// F13 2n-minima関数
//------------------------------------------------------------
double N2_Minima(double *x) {
	int i;
	double sum;
	for (i = 0, sum = 0; i<MAXGSIZE; i++) {
		sum += (pow(x[i], 4) - 16 * pow(x[i], 2) + 5 * x[i]);
	}
	return (sum / 2);
}
//------------------------------------------------------------
// 目的関数値を計算
//------------------------------------------------------------
double Calc_Objective_Function(double *x)
{
	if (FUNCNO == 1)return  Sphere(x);
	else if (FUNCNO == 2)return  Rosenbrock(x);
	else if (FUNCNO == 3)return  Rastrigin(x);
	else if (FUNCNO == 4)return  Griewank(x);
	else if (FUNCNO == 5)return  Ackley(x);
	else if (FUNCNO == 6)return  Schwefel(x);
	else if (FUNCNO == 7)return  Ridge(x);
	else if (FUNCNO == 8)return  Bohachevsky(x);
	else if (FUNCNO == 9)return  Schaffer(x);
	else if (FUNCNO == 10)return Ellipsoid(x);
	else if (FUNCNO == 11)return K_Tablet(x);
	else if (FUNCNO == 12)return Shifted_Rastrigin(x);
	else if (FUNCNO == 13)return Cigar(x);
	else if (FUNCNO == 14)return N2_Minima(x);
	else exit(0);
}
//------------------------------------------------------------
//初期個体群の生成
//------------------------------------------------------------
void Init_Vector(void)
{
	int i, j;		//繰り返し用変数
	double r;		//定義域用変数
	r = RANGE;
	for (i = 0; i<MAXPSIZE; i++) {
		for (j = 0; j<MAXGSIZE; j++) {
			cVect[i][j] = RANGE*(genrand_real1() * 2 - 1);
		}
	}
	//初期ベクトルをpBestVectorに保存する
	for (i = 0; i<MAXPSIZE; i++) {
		for (j = 0; j<MAXGSIZE; j++) {
			pBestVector[i][j] = cVect[i][j];
		}
	}
}
//------------------------------------------------------------
// 初期集団の評価（目的関数値を適応度として利用）
//------------------------------------------------------------
void Evaluate_Init_Vector(void) {
	int i;
	for (i = 0; i<MAXPSIZE; i++) {
		cFitness[i] = Calc_Objective_Function(cVect[i]);
		//初期値を初期pBestとして保存
		pBestFitness[i] = cFitness[i];
	}
}
//------------------------------------------------------------
//新規ベクトルの生成のための親ベクトルの選択
//------------------------------------------------------------
void Select_pVector(int pop1)
{
	register int i;
	int pop2, pop3, pop4;
	do {
		pop2 = (int)MAXPSIZE*genrand_real1();
		pop3 = (int)MAXPSIZE*genrand_real1();
		pop4 = (int)MAXPSIZE*genrand_real1();
	} while (pop1 == pop2 || pop1 == pop3 || pop1 == pop4 || pop2 == pop3 || pop2 == pop4 || pop3 == pop4);
	for (i = 0; i<MAXGSIZE; i++)pVect1[i] = cVect[pop2][i];
	for (i = 0; i<MAXGSIZE; i++)pVect2[i] = cVect[pop3][i];
	for (i = 0; i<MAXGSIZE; i++)pVect3[i] = cVect[pop4][i];
}
//------------------------------------------------------------
//DEの操作
//------------------------------------------------------------
void DE_Operation(int pop1)
{
	register int i;
	int N = 0, L = 0;
	//DE/rand/1/exp
	if (DE_ALGORITHM_NO == 1) {
		for (i = 0; i<MAXGSIZE; i++) nVect[pop1][i] = cVect[pop1][i];
		N = (int)(genrand_real1()*MAXGSIZE);
		L = 0;
		do {
			nVect[pop1][N] = pVect1[N] + MRATE * (pVect2[N] - pVect3[N]);
			if (nVect[pop1][N] < -RANGE) nVect[pop1][N] = pVect1[N] + genrand_real1() * (-RANGE - pVect1[N]);
			if (nVect[pop1][N] >  RANGE) nVect[pop1][N] = pVect1[N] + genrand_real1() * (RANGE - pVect1[N]);
			N = (N + 1) % MAXGSIZE;
			L++;
		} while (genrand_real1() < CRATE && L < MAXGSIZE);
	}
	//DE/best/1/exp
	else if (DE_ALGORITHM_NO == 2) {
		for (i = 0; i<MAXGSIZE; i++) nVect[pop1][i] = cVect[pop1][i];
		N = (int)(genrand_real1()*MAXGSIZE);
		L = 0;
		do {
			nVect[pop1][N] = gBestVector[N] + MRATE * (pVect2[N] - pVect3[N]);
			if (nVect[pop1][N] < -RANGE) nVect[pop1][N] = pVect1[N] + genrand_real1() * (-RANGE - pVect1[N]);
			if (nVect[pop1][N] >  RANGE) nVect[pop1][N] = pVect1[N] + genrand_real1() * (RANGE - pVect1[N]);
			N = (N + 1) % MAXGSIZE;
			L++;
		} while (genrand_real1() < CRATE && L < MAXGSIZE);
	}
	//DE/rand/1/bin
	else if (DE_ALGORITHM_NO == 3) {
		for (i = 0; i<MAXGSIZE; i++) nVect[pop1][i] = cVect[pop1][i];
		N = (int)(genrand_real1()*MAXGSIZE);
		for (L = 0; L<MAXGSIZE; L++) {
			if (L == 0 || genrand_real1() < CRATE) {
				nVect[pop1][N] = pVect1[N] + MRATE * (pVect2[N] - pVect3[N]);
				if (nVect[pop1][N] < -RANGE) nVect[pop1][N] = pVect1[N] + genrand_real1() * (-RANGE - pVect1[N]);
				if (nVect[pop1][N] >  RANGE) nVect[pop1][N] = pVect1[N] + genrand_real1() * (RANGE - pVect1[N]);
			}
			else {
				nVect[pop1][N] = cVect[pop1][N];
			}
			N = (N + 1) % MAXGSIZE;
		}
	}
	//DE/best/1/bin
	else if (DE_ALGORITHM_NO == 4) {
		for (i = 0; i<MAXGSIZE; i++) nVect[pop1][i] = cVect[pop1][i];
		N = (int)(genrand_real1()*MAXGSIZE);
		for (L = 0; L<MAXGSIZE; L++) {
			if (L == 0 || genrand_real1() < CRATE) {
				nVect[pop1][N] = gBestVector[N] + MRATE * (pVect2[N] - pVect3[N]);
				if (nVect[pop1][N] < -RANGE) nVect[pop1][N] = pVect1[N] + genrand_real1() * (-RANGE - pVect1[N]);
				if (nVect[pop1][N] >  RANGE) nVect[pop1][N] = pVect1[N] + genrand_real1() * (RANGE - pVect1[N]);
			}
			else {
				nVect[pop1][N] = cVect[pop1][N];
			}
			N = (N + 1) % MAXGSIZE;
		}
	}
	//JADE	DE/curreny-to-pbest/1
	else if (DE_ALGORITHM_NO == 5) {
		for (i = 0; i<MAXGSIZE; i++) nVect[pop1][i] = cVect[pop1][i];
		N = (int)(genrand_real1()*MAXGSIZE);
		for (L = 0; L<MAXGSIZE; L++) {
			if (L == 0 || genrand_real1() < CRATE) {
				nVect[pop1][N] = pVect1[pop1] + Mrate[pop1] * (gBestVector[N] - pVect1[N]) + Mrate[pop1] * (pVect1[pop1] - pVect2[pop1]);	//変更中
				if (nVect[pop1][N] < -RANGE) nVect[pop1][N] = pVect1[N] + genrand_real1() * (-RANGE - pVect1[N]);
				if (nVect[pop1][N] >  RANGE) nVect[pop1][N] = pVect1[N] + genrand_real1() * (RANGE - pVect1[N]);
			}
			else {
				nVect[pop1][N] = cVect[pop1][N];
			}
			N = (N + 1) % MAXGSIZE;
		}
	}
	else exit(0);
}

//------------------------------------------------------------
//新しいベクトルの評価
//------------------------------------------------------------
void Evaluate_New_Vector(int pop1)
{
	nFitness[pop1] = Calc_Objective_Function(nVect[pop1]);

}
//------------------------------------------------------------
//ベクトルの比較
//------------------------------------------------------------
void Compare_Vector(void)
{
	int i, j;				//繰返し用変数.
	for (i = 0; i<MAXPSIZE; i++) {
		//新しいベクトルが良ければ置き換え操作を行う
		if (nFitness[i] < cFitness[i]) {
			cFitness[i] = nFitness[i];
			for (j = 0; j<MAXGSIZE; j++)cVect[i][j] = nVect[i][j];
		}
		else continue;
	}
}
//------------------------------------------------------------
//エリート選択
//------------------------------------------------------------
void Select_Elite_Vector(int itime, int gtime)
{
	int i;						//繰返し用変数.
	int num;				//添字
	double best;			//一時保存用
	for (i = 0, num = 0, best = cFitness[0]; i<MAXPSIZE; i++) {
		if (cFitness[i]<best) {
			best = cFitness[i];
			num = i;
		}
	}
	for (i = 0; i<MAXGSIZE; i++)gBestVector[i] = cVect[num][i];
	gBestFitness = cFitness[num];
	gBestHistory[itime][gtime] = gBestFitness;
	printf("%.70lf\n", gBestFitness);
}


//分布の生成
double rand_cauchy(double mu, double gamma) {
	return mu + gamma * tan(M_PI*(genrand_real1() - 0.5));
}
double rand_normal(double mu, double sigma) {
	double z = sqrt(-2.0*log(genrand_real1())) * sin(2.0*M_PI*genrand_real3());
	return mu + sigma * z;
}


//JADE 作成中
/*
void Pameter_Filter(double Sf, double Scr) {
	int i;
	//	std::random_device seed_gen;//乱数の生成
	//	std::default_random_engine engine(seed_gen());//擬似乱数の生成
	// 位置母数cMrate、尺度母数Sfで分布させる
	//	std::cauchy_distribution<> cauchy_data(cMrate, Sf); //Sfが0超過でないとエラーが発生する
	//	std::cauchy_distribution<> cauchy_data(0.0, 1.0);
	// 平均cCrate、標準偏差Scr^2で分布させる
	//	std::normal_distribution<> nomal_data(cCrate, Scr*Scr);//Scrが0より大きくないとエラーが発生する
	//	std::normal_distribution<> nomal_data(0, 1.0);

	for (i = 0; i < MAXPSIZE; i++) {
	reMrate:
		//Mrate[i] = cauchy_data((unsigned int)genrand_int32);if (Mrate[i]>1 || Mrate[i]<0) goto	reMrate;//生成方法が不明
		//			Mrate[i] = cauchy_data(engine);if (Mrate[i]>1 || Mrate[i]<0) goto	reMrate;//生成方法が不明
		Mrate[i] = rand_normal(cMrate, Sf); if (Mrate[i]>1 || Mrate[i]<0) goto	reMrate;//生成方法が不明
	reCrate:
		//		Crate[i] = nomal_data(engine); if (Crate[i]>1 || Crate[i]<0) goto	reCrate;//生成方法が不明
		Crate[i] = rand_normal(cCrate, Scr*Scr); if (Crate[i]>1 || Crate[i]<0) goto	reCrate;//生成方法が不明
	}
}
*/
void Pameter_Filter(double Sf, double Scr) {
	int i;
	for (i = 0; i < MAXPSIZE; i++) {
	reMrate:
		Mrate[i] = rand_cauchy(cMrate, Sf);; if (Mrate[i]>1 || Mrate[i]<0) goto	reMrate;
	reCrate:
		Crate[i] = rand_normal(cCrate, Scr*Scr); if (Mrate[i]>1 || Mrate[i]<0) goto	reMrate;
	}
}


//パラメータの更新	JADE 作成中
void New_parameter() {
	int i;
	double Sn = 0, Sf = 1.0, Sf2 = 0, Scr = 1.0;
	for (i = 0; i < MAXPSIZE; i++) {
		if (fabs(nVect[i][MAXPSIZE])	>	fabs(cVect[i][MAXPSIZE])) {
			Sn += 1;
			Sf2 += Mrate[i] * Mrate[i];
			Scr += Crate[i];
		}
	}
	if (Sn != 0) {
		cMrate = (1 - C)*cMrate + (C*Sf2 / Sf);
		cCrate = (1 - C)*cCrate + ((C*Scr) / Sn);
	}
	Pameter_Filter(Sf, Scr);
}




//------------------------------------------------------------
//ファイル出力1　最良値の推移
//------------------------------------------------------------
void Output_To_File1(void)
{
	int i, j;							//繰返し用変数.
	FILE *fp;						//ファイルポインタ
	char filename[50];		//ファイル名
	time_t timer;				//時間計測用
	struct tm *t_st;			//時間計測用
	time(&timer);				//時間の取得
	t_st = localtime(&timer);	//時間の変換
	sprintf_s(filename, "DE_gBestHistory%04d%02d%02d%02d%02d.txt",
		t_st->tm_year + 1900, t_st->tm_mon + 1,
		t_st->tm_mday, t_st->tm_hour, t_st->tm_min);
	fp = fopen(filename, "a");
	for (i = 0; i<MAXGENERATIO; i++) {
		for (j = 0; j<EXTIME; j++) {
			fprintf(fp, "%15.70lf", gBestHistory[j][i]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
}
//------------------------------------------------------------
//集団の多様性の評価
//------------------------------------------------------------
void Calc_Diversity(int itime, int gtime)
{
	int i, j;					//繰返し用変数.
	double sum[MAXGSIZE];		//合計
	double ave[MAXGSIZE];		//平均
	double div;					//分散
	//初期化
	for (i = 0; i<MAXGSIZE; i++) {
		sum[i] = 0.0;
		ave[i] = 0.0;
	}
	//合計
	for (i = 0; i<MAXPSIZE; i++) {
		for (j = 0; j<MAXGSIZE; j++) {
			sum[j] += cVect[i][j];
		}
	}
	//平均
	for (i = 0; i<MAXGSIZE; i++) {
		ave[i] = sum[i] / MAXPSIZE;
	}
	//分散
	for (i = 0, div = 0.0; i<MAXPSIZE; i++) {
		for (j = 0; j<MAXGSIZE; j++) {
			div += (ave[j] - cVect[i][j])*(ave[j] - cVect[i][j]);
		}
	}
	pDiversity[itime][gtime] = div / MAXPSIZE;
}
//------------------------------------------------------------
//ファイル出力2　分散値の推移
//------------------------------------------------------------
void Output_To_File2(void)
{
	int i, j;					//繰返し用変数.
	FILE *fp;					//ファイルポインタ
	char filename[50];			//ファイル名
	time_t timer;				//時間計測用
	struct tm *t_st;			//時間計測用
	time(&timer);				//時間の取得
	t_st = localtime(&timer);	//時間の変換
	sprintf_s(filename, "DE_pDiversity%04d%02d%02d%02d%02d.txt",
		t_st->tm_year + 1900, t_st->tm_mon + 1,
		t_st->tm_mday, t_st->tm_hour, t_st->tm_min);
	fp = fopen(filename, "a");
	for (i = 0; i<MAXGENERATIO; i++) {
		for (j = 0; j<EXTIME; j++) {
			fprintf(fp, "%15.70lf", pDiversity[j][i]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
}
//------------------------------------------------------------
//ファイル出力3　最適解発見世代，収束率
//------------------------------------------------------------
void Output_To_File3(void)
{
	register int i;
	FILE *fp;
	char filename[50];
	time_t timer;
	struct tm *t_st;
	time(&timer);
	t_st = localtime(&timer);
	sprintf_s(filename, "DE_sRate_gTable%04d%02d%02d%02d%02d.txt",
		t_st->tm_year + 1900, t_st->tm_mon + 1,
		t_st->tm_mday, t_st->tm_hour, t_st->tm_min);
	fp = fopen(filename, "a");
	for (i = 0; i<EXTIME; i++) {
		fprintf(fp, "%6d %d\n", gTable[i], sRate[i]);
	}
}
//------------------------------------------------------------
//メイン関数
//------------------------------------------------------------
int main(void)
{
	int episode;
	int iteration;
	int pop;
	init_genrand((unsigned)time(NULL));	//MTの初期化
	for (iteration = 0; iteration<EXTIME; iteration++) {
		episode = 0;
		Init_Vector();
		Evaluate_Init_Vector();
		while (episode<MAXGENERATIO) {
			Select_Elite_Vector(iteration, episode);
			Calc_Diversity(iteration, episode);
			for (pop = 0; pop<MAXPSIZE; pop++) {
				Select_pVector(pop);
				DE_Operation(pop);
				New_parameter(); //JADE関数の中に入れるべき
				Evaluate_New_Vector(pop);
			}
			Compare_Vector();
			episode++;
		}
		gTable[iteration] = episode;
		if (gBestFitness<TERMINATE)sRate[iteration] = 1;
		else sRate[iteration] = 0;
		gBestTable[iteration] = gBestFitness;
	}
	Output_To_File1();
	Output_To_File2();
	Output_To_File3();
	getchar();
	return 0;
}