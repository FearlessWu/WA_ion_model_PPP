/***************************Matrix.h**********************************************
*������
*
*Contact:
*		wuguanbin@shao.ac.cn
*������
*		Matrix A(2,3);				����һ��2x3��0�����
*		Matrix B(A);				����һ����A��ͬ�ľ���
*		Matrix C;					����C��һ������δ��ʼ��
*
*���ܣ�
*		A[i]=k;						�Ծ���A�ĵ�i(int)��Ԫ�ظ�ֵk(int or double)
*		A(i,j)=k;					�Ծ���A�ĵ�i(int����0��ʼ����)�е�j(int����0��ʼ����)��Ԫ�ظ�ֵk(int or double)
*		B=A;						��A������B����ֵ
*		B==A;						�ж�B��A�Ƿ���ȣ���ȷ���true�����򷵻�false
*		A+B;						A��B������ӣ��ɸ�ֵ��C�����磺C=A+B;
*		A-B;						A��B����������ɸ�ֵ��C�����磺C=A-B;
*		A*B;						A��B������ˣ��ɸ�ֵ��C�����磺C=A*B;
*		A*a;						a(int or double)�����A���ˣ��ɸ�ֵ��C�����磺C=A*a;
*		A.trans();					����A��ת�ã��ɸ�ֵ������C���磺C=A.trans();
*		A.inv();					����A������󣬿ɸ�ֵ������C���磺C=A.inv();
*		A.getRow();					��ȡ����A������
*		A.getCol();					��ȡ����A������
*		A.print();					������A�������Ļ
*		erase=EraseZerosCol(A,B);	�������A������Ԫ��Ϊ0���У����������Ľ����������B�����ؾ���A�б�������кŸ�erase(set<int>)
*		erase=EraseZerosRow(A,B);	�������A������Ԫ��Ϊ0���У����������Ľ����������B�����ؾ���A�б�������кŸ�erase(set<int>)
*		CreateDiag(A,sig);			����һ�� vector<double> ���͵�sig��������sigΪ�Խ�Ԫ�صĶԽ���A
*		fileprint(fp);				����һ��FILE*fp,������������ļ�
*		EraseRowAccordToRowNum(Matrix&a,set<int>& b,Matrix& output)	�����к�ɾ�������Ӧ����
*
*Last Update:
*		2019 05 24:finished the whole class
*		2019 05	26:���� EraseZerosRow ��CreateDiag����
*		2019 05	27:���� fileprint ����
*		2019 05	29:���� EraseRowAccordToRowNum ����
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