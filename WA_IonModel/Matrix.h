/***************************Matrix.h**********************************************
*矩阵类
*
*Contact:
*		wuguanbin@shao.ac.cn
*声明：
*		Matrix A(2,3);				生成一个2x3的0零矩阵
*		Matrix B(A);				生成一个与A相同的矩阵
*		Matrix C;					声明C是一个矩阵，未初始化
*
*功能：
*		A[i]=k;						对矩阵A的第i(int)个元素赋值k(int or double)
*		A(i,j)=k;					对矩阵A的第i(int，从0开始计数)行第j(int，从0开始计数)列元素赋值k(int or double)
*		B=A;						用A矩阵向B矩阵赋值
*		B==A;						判断B与A是否相等，相等返回true，否则返回false
*		A+B;						A、B矩阵相加，可赋值给C矩阵，如：C=A+B;
*		A-B;						A、B矩阵相减，可赋值给C矩阵，如：C=A-B;
*		A*B;						A、B矩阵相乘，可赋值给C矩阵，如：C=A*B;
*		A*a;						a(int or double)与矩阵A数乘，可赋值给C矩阵，如：C=A*a;
*		A.trans();					矩阵A的转置，可赋值给矩阵C，如：C=A.trans();
*		A.inv();					矩阵A的逆矩阵，可赋值给矩阵C，如：C=A.inv();
*		A.getRow();					获取矩阵A的行数
*		A.getCol();					获取矩阵A的列数
*		A.print();					将矩阵A输出到屏幕
*		erase=EraseZerosCol(A,B);	清除矩阵A中所有元素为0的列，并将清除后的结果赋给矩阵B，返回矩阵A中被清除的列号给erase(set<int>)
*		erase=EraseZerosRow(A,B);	清除矩阵A中所有元素为0的行，并将清除后的结果赋给矩阵B，返回矩阵A中被清除的行号给erase(set<int>)
*		CreateDiag(A,sig);			给定一个 vector<double> 类型的sig，生成以sig为对角元素的对角阵A
*		fileprint(fp);				给定一个FILE*fp,将矩阵输入该文件
*		EraseRowAccordToRowNum(Matrix&a,set<int>& b,Matrix& output)	根据行号删除矩阵对应的行
*
*Last Update:
*		2019 05 24:finished the whole class
*		2019 05	26:新增 EraseZerosRow 、CreateDiag函数
*		2019 05	27:新增 fileprint 函数
*		2019 05	29:新增 EraseRowAccordToRowNum 函数
**********************************************************************************/
#pragma once
#include<cmath>
#include<cstring>
#include<stdio.h>
#include<new>
#include<iostream>
#include<set>
#include<vector>
#include<stdio.h>
#define POS(i,j,k) ((i)*j+k)
using std::vector;
using std::set;
class Matrix {
public:
	Matrix(int n, int m);
	Matrix();
	Matrix(Matrix &a);
	double operator=(double a);
	friend Matrix operator*(Matrix &a, Matrix &b);
	Matrix operator*(double a);
	friend Matrix operator*(double a, Matrix&b);
	Matrix operator+(Matrix &a);
	Matrix operator-(Matrix &a);
	Matrix trans();
	bool operator==(Matrix &a);
	double& operator[](int a);
	Matrix &operator=(Matrix& other);
	double& operator()(int i, int j);
	double getRow();
	double getCol();
	Matrix inv();
	friend std::set<int> EraseZerosCol(Matrix &input, Matrix &output);
	friend std::set<int> EraseZerosRow(Matrix &input, Matrix &output);
	friend void CreateDiag(Matrix &a, vector<double> sig);
	friend void EraseRowAccordToRowNum(Matrix&a, set<int>& b, Matrix& output);
	void print();
	void fileprint(FILE*fp);
	~Matrix();


private:
	int row, col;
	double *elements;
	void initMatrix();
	int matinv(double *A, int n);
	int *imat(int n, int m);
	double *mat(int n, int m);
	void matcpy(double *A, const double *B, int n, int m);
	int ludcmp(double *A, int n, int *indx, double *d);
	void lubksb(const double *A, int n, const int *indx, double *b);
};