#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>
#include <vector>
#include <conio.h>

// ���������� ���������� �� ����������
int	intervalsX, intervalsY, intervalsZ;

// ���� �����
double	hx = 1,
hy = 1,
hz = 1;

// ���������� ������ � ����� �������
double	begX = 0, endX = 6,
begY = 0, endY = 6,
begZ = 0, endZ = 6;

// ������ �����
int dim;

struct Point {
	double x, y, z;
	int globalNumber;
	void set(double x1, double y1, double z1, int glN)
	{
		x = x1;
		y = y1;
		z = z1;
		globalNumber = glN;
	};
};

class SLAE
{
	
	// ������������ ���-�� ��������
	int maxiter = 10000;
	double eps = 1e-6;

public:
	SLAE();
	~SLAE();
	// ���-�� ������� ����������
	int m, m2;
	// ����������� �������
	int n;
	// ���-�� ��-��� � ���������� al2, au2 - n1 � al3, au3 - n2
	int n1, n2;
	std::vector <double> di; // ������� ���������
							 // ���������, �������� � �������
	std::vector <double> al1, au1;
	// ���������, ������������� �� m ��-��� ����/���� ������� ��������� � �������� � ���
	std::vector <double> al2, au2;
	// ���������, ������������� �� m2-��� ����/���� ������� ��������� � �������� � ���
	std::vector <double> al3, au3;
	// ������ ������ �����
	std::vector <double> f;
	// ������� �������
	std::vector <double> x;
	// ���������� ������
	std::vector <double> x_pred;

	// ����� �����
	void Jakobi();
	// ��������� ������ �� ������
	void Mult_Ax(std::vector <double> &vo);
	// ������������� �������
	double Residual();
	double summ(int i, std::vector <double> &x_next, std::vector <double> x_pred);

	double CalcF1BC(double x, double y, double z);
	double CalcF(double x, double y, double z);
};
