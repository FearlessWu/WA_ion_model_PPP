//***************************Matrix.cpp***********************************
#include"Matrix.h"
//���캯�������������������� n x m 0����
Matrix::Matrix(int n, int m) {
	
	if (n <= 0 || m <= 0)
		throw("Matrix error: ��������������");
	elements = new double[n*m];
	this->row = n;
	this->col = m;
	initMatrix();
}
//���캯����0x0 ����
Matrix::Matrix() :row(0), col(0), elements(nullptr) {
	
}
//���캯���������ʼ��
Matrix::Matrix(Matrix&a)
{	
	if (col + row > 0)
		delete[] elements;
	row = a.row;
	col = a.col;
	elements = new double[row*col];
	for (int i = 0;i < row*col;++i)
		elements[i] = a.elements[i];

}
//��������
Matrix::~Matrix() {
	
	if (row + col > 0)
	{
		delete[]elements;
		row = 0;
		col = 0;
	}
	else
	{
		delete elements;
		row = 0;
		col = 0;
	}
}
//���������Ļ
void Matrix::print() {
	
	printf("\n");
	for (int i = 0;i<row;++i)
		for (int j = 0;j<col;++j)
		{
			if (j == col - 1)
				printf(" %10.8f\n", elements[POS(i, col, j)]);
			else
				printf(" %10.8f ", elements[POS(i, col, j)]);
		}
	printf("\n");
}
//����ת�ú���
Matrix Matrix::trans() {
	
	Matrix C(col, row);
	for (int i = 0;i<row;++i)
		for (int j = 0;j<col;++j)
		{
			C.elements[POS(j, row, i)] = elements[POS(i, col, j)];
		}
	return C;
}
//��þ�������
double Matrix::getRow() {
	
	return row;
}
//��þ�������
double Matrix::getCol() {
	
	return col;
}
//����=����������ڶԾ���Ԫ�ظ�ֵ
double  Matrix::operator=(double a) {
	
	return a;
}
//����=����������ھ���Ծ���ֵ
Matrix &Matrix::operator=(Matrix &other) {
	
	if (row != 0 || col != 0)
	{
		if (row != other.row || col != other.col)
			throw("Matrix error: ����ά�Ȳ�һ�²��ܸ�ֵ��");
	}
	if (row == 0 && col == 0) {
		row = other.row;
		col = other.col;
		elements = new double[row*col];
	}
	//elements = new double[row*col];
	for (int i = 0;i < row*col;++i)
		elements[i] = other.elements[i];
	return *this;
}
//����*��������������˾���
Matrix Matrix::operator*(double a) {
	
	Matrix C(*this);
	for (int i = 0;i<row*col;++i)
		C.elements[i] = elements[i] * a;
	return C;
}
//����*��������������˾���
Matrix operator*(double a, Matrix &b) {
	
	Matrix C(b);
	for (int i = 0;i<C.row*C.col;++i)
		C.elements[i] = b.elements[i] * a;
	return C;
}
//����+����������ھ������
Matrix Matrix::operator+(Matrix &a) {
	
	if (row != a.row || col != a.col)
		throw("Matrix error: ��ͬά�Ⱦ��������");
	Matrix C(a);
	for (int i = 0;i<row*col;++i)
		C.elements[i] = elements[i] + a[i];
	return C;
}
//����+����������ھ������
Matrix Matrix::operator-(Matrix &a) {
	
	if (row != a.row || col != a.col)
		throw("Matrix error: ��ͬά�Ⱦ��������");
	Matrix C(a);
	for (int i = 0;i<row*col;++i)
		C.elements[i] = elements[i] - a[i];
	return C;
}
//����==��������жϾ����Ƿ����
bool Matrix::operator==(Matrix &a) {
	
	if (row != a.row || col != a.col)
		return false;

	for (int i = 0;i<row*col;i++)
	{
		if (fabs(elements[i] - a.elements[i])>1e-14)
			return false;
	}
	return true;
}
//����()����������ڻ�ȡ��i�е�j��Ԫ��,i j����0��ʼ����
double& Matrix::operator()(int i, int j) {
	
	if (i > row || j > col || i < 0 || j < 0)
		throw("Matrix error: ��������к��±겻�Ϸ���");
	return elements[POS(i, col, j)];
}
//����*����������ھ������
Matrix operator*(Matrix &a, Matrix &b) {
	
	if (a.col != b.row)
		throw("Matrix error: �������ά�ȱ���һ��");
	else if (a.row <= 0 || a.col <= 0 || b.row <= 0 || b.col <= 0)
	{
		throw("Matrix error: ����ά�ȱ������0");
	}
	Matrix C(a.row, b.col);

	double tmp = 0;

	int k = b.col;
	for (int i = 0;i<a.row;i++)
	{
		for (int j = 0;j<k;++j)
		{
			for (int l = 0;l<a.col;++l)
			{
				tmp += a.elements[POS(i, a.col, l)] * b.elements[POS(l, k, j)];
			}
			if (fabs(tmp)<1e-12)
				C.elements[POS(i, k, j)] = 0;
			else
				C.elements[POS(i, k, j)] = tmp;
			tmp = 0;
		}
	}
	return C;
}
//����[]����������ڻ�ȡ�����i��Ԫ��
double& Matrix::operator[](int i) {
	
	if (i<0 || i >= row*col)
		throw("Matrix error: []�±���������");
	return elements[i];
}
//��������
Matrix Matrix::inv() {
	
	if (row != col || row <= 0)
		throw("Matrix error: ����������Ϊ����");
	Matrix tmp1(row, row), tmp2(row, row), tmp3(row, row);
	int flag;
	tmp1 = trans();
	if (!matinv(tmp1.elements, row))
		throw("Matrix error: ���������������������");
	tmp2 = tmp1.trans();
	for (int i = 0;i<row*row;++i)
		tmp3.elements[i] = tmp2.elements[i];
	return tmp3;
}
//���Ԫ��ȫΪ����У����ر���� set<int> �к�
std::set<int> EraseZerosRow(Matrix &input, Matrix &output)
{
	
	if (input.row <= 0 || input.col <= 0)
		throw("Matrix error: �������ά�Ȳ������ڵ���0��");
	if (output.row + output.col>0)
	{
		delete[]output.elements;
		output.row = 0;
		output.col = 0;
	}
	std::set<int> erase;

	for (int i = 0;i < input.row;++i) {
		int flag = 0;
		for (int j = 0;j < input.col;++j) {
			if (input(i, j) != 0.0)
			{
				flag = 1;
				break;
			}
		}
		if (flag == 0)
			erase.insert(i);
	}
	int row = input.row - (int)erase.size();
	int col = input.col;
	if (row == 0) {
		return erase;
	}
	output.row = row;
	output.col = col;
	output.elements = new double[row*col];
	std::set<int>::iterator iter;
	int k = 0;
	int r, c;
	for (int i = 0;i < input.row;++i)
	{
		iter = erase.find(i);
		if (iter != erase.end())
			continue;
		for (int j = 0;j < input.col;++j) {

			output(k, j) = input(i, j);
		}
		k++;
	}
	return erase;
}
//���Ԫ��ȫΪ����У����ر���� set<int> �к�
std::set<int> EraseZerosCol(Matrix& input, Matrix& output) {
	
	if (input.row <= 0 || input.col <= 0)
		throw("Matrix error: �������ά�Ȳ������ڵ���0��");
	if (output.row + output.col>0)
	{
		delete[]output.elements;
		output.row = 0;
		output.col = 0;
	}
	std::set<int> erase;

	for (int i = 0;i < input.col;++i) {
		int flag = 0;
		for (int j = 0;j < input.row;++j) {
			if (input(j, i) != 0.0)
			{
				flag = 1;
				break;
			}
		}
		if (flag == 0)
			erase.insert(i);
	}
	int row = input.row;
	int col = input.col - (int)erase.size();
	if (col == 0) {
		return erase;
	}
	output.row = row;
	output.col = col;
	output.elements = new double[row*col];
	std::set<int>::iterator iter;
	int k = 0;
	int r, c;
	for (int j = 0;j < input.col;++j)
	{
		iter = erase.find(j);
		if (iter != erase.end())
			continue;
		for (int i = 0;i < input.row;++i) {

			output(i, k) = input(i, j);
		}
		k++;
	}
	return erase;

}
//��sigΪ�Խ�Ԫ�����ɶԽ���
void CreateDiag(Matrix &a, vector<double> sig) {
	if (a.row + a.col>0)
		delete[]a.elements;
	if (a.row < 0 || a.col < 0)
		throw("Matrix error: �����ά�Ȳ�������0");
	a.row = sig.size();
	a.col = a.row;
	a.elements = new double[a.row*a.col];
	memset(a.elements, 0, sizeof(double)*a.row*a.col);
	vector<double>::iterator iter = sig.begin();
	for (int i = 0;i < sig.size();++i)
	{
		a(i, i) = *iter;
		iter++;
	}
}
//����������ļ�
void Matrix::fileprint(FILE*fp)
{
	Matrix t(*this);
	if (col + row <= 0)
		throw("Matrix error: �ļ����ά�Ȳ������ڵ���0��");
	for (int i = 0;i<row;++i)
		for (int j = 0;j<col;++j)
		{
			if (j == col - 1)
				fprintf(fp, "	%10.5f\n", t(i, j));
			else
				fprintf(fp, "	%10.5f ", t(i, j));
		}
	fprintf(fp, "\n");
}
//�����к���������Ӧ����
void EraseRowAccordToRowNum(Matrix&input, set<int>& b, Matrix& output)
{
	if (output.row + output.col>0)
		delete[] output.elements;
	if (input.row - b.size()<0)
		throw("Matrix error: set<>������кŹ���");
	if ((input.row - b.size())<1e-12)
	{
		output.col = 0;
		output.row = 0;
		return;
	}
	output.row = input.row - b.size();
	output.col = input.col;
	output.elements = new double[output.row*output.col];
	int i = 0;
	for (int j = 0;j<input.row;++j)
	{
		if (b.count(j))
			continue;
		for (int k = 0;k<input.col;++k)
		{
			output(i, k) = input(j, k);
		}
		i++;
	}

}
//****************************************  Private Function  ********************************************
int Matrix::matinv(double *A, int n)
{

	double d, *B;
	int i, j, *indx;

	indx = imat(n, 1); B = mat(n, n); matcpy(B, A, n, n);
	if (ludcmp(B, n, indx, &d)) { free(indx); free(B); return 0; }
	for (j = 0;j<n;j++) {
		for (i = 0;i<n;i++) A[i + j*n] = 0.0; A[j + j*n] = 1.0;
		lubksb(B, n, indx, A + j*n);
	}
	free(indx); free(B);
	return 1;
}
int* Matrix::imat(int n, int m)
{
	int *p;

	if (n <= 0 || m <= 0) return NULL;
	if (!(p = (int *)malloc(sizeof(int)*n*m))) {
		throw("integer matrix memory allocation error: n=%d,m=%d\n", n, m);
	}
	return p;
}
double* Matrix::mat(int n, int m)
{
	double *p;

	if (n <= 0 || m <= 0) return NULL;
	if (!(p = (double *)malloc(sizeof(double)*n*m))) {
		throw("matrix memory allocation error: n=%d,m=%d\n", n, m);
	}
	return p;
}
void Matrix::matcpy(double *A, const double *B, int n, int m)
{
	memcpy(A, B, sizeof(double)*n*m);
}
int Matrix::ludcmp(double *A, int n, int *indx, double *d)
{
	double big, s, tmp, *vv = mat(n, 1);
	int i, imax = 0, j, k;

	*d = 1.0;
	for (i = 0;i<n;i++) {
		big = 0.0; for (j = 0;j<n;j++) if ((tmp = fabs(A[i + j*n]))>big) big = tmp;
		if (big>0.0) vv[i] = 1.0 / big; else { free(vv); return -1; }
	}
	for (j = 0;j<n;j++) {
		for (i = 0;i<j;i++) {
			s = A[i + j*n]; for (k = 0;k<i;k++) s -= A[i + k*n] * A[k + j*n]; A[i + j*n] = s;
		}
		big = 0.0;
		for (i = j;i<n;i++) {
			s = A[i + j*n]; for (k = 0;k<j;k++) s -= A[i + k*n] * A[k + j*n]; A[i + j*n] = s;
			if ((tmp = vv[i] * fabs(s)) >= big) { big = tmp; imax = i; }
		}
		if (j != imax) {
			for (k = 0;k<n;k++) {
				tmp = A[imax + k*n]; A[imax + k*n] = A[j + k*n]; A[j + k*n] = tmp;
			}
			*d = -(*d); vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (A[j + j*n] == 0.0) { free(vv); return -1; }
		if (j != n - 1) {
			tmp = 1.0 / A[j + j*n]; for (i = j + 1;i<n;i++) A[i + j*n] *= tmp;
		}
	}
	free(vv);
	return 0;
}
void Matrix::lubksb(const double *A, int n, const int *indx, double *b)
{
	double s;
	int i, ii = -1, ip, j;

	for (i = 0;i<n;i++) {
		ip = indx[i]; s = b[ip]; b[ip] = b[i];
		if (ii >= 0) for (j = ii;j<i;j++) s -= A[i + j*n] * b[j]; else if (s) ii = i;
		b[i] = s;
	}
	for (i = n - 1;i >= 0;i--) {
		s = b[i]; for (j = i + 1;j<n;j++) s -= A[i + j*n] * b[j]; b[i] = s / A[i + i*n];
	}
}
void Matrix::initMatrix() {

	memset(elements, 0, sizeof(double)*row*col);
}