#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <conio.h>

// ���������� ���������� �� ����������
const int	intervalsX = 5,
	intervalsY = 5,
	intervalsZ = 5;

// ���� �����
double hx, hy, hz;

// ���������� ������ � ����� �������
double	begX = 0, endX = 6,
	begY = 0, endY = 6,
	begZ = 0, endZ = 6;

// ������ �����
const int dim = (intervalsX + 1) * (intervalsY + 1) * (intervalsZ + 1);

struct Point {
	double x, y, z;
	void set(double x1, double y1, double z1)
	{
		x = x1;
		y = y1;
		z = z1;
	};
};

class Matrix
{
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
	// ������������ ���-�� ��������
	int maxiter = 10000; 
	double eps = 1e-6;

public:
	Matrix(); 
	~Matrix();
	// ����� �����
	void Jakobi();
	// ��������� ������ �� ������
	void Mult_Ax(std::vector <double> &vo);
	// ������������� �������
	double Residual(); 
	double summ(int i, std::vector <double> &x_next, std::vector <double> x_pred);
	
};
