#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>
#include <vector>
#include <conio.h>

// Количество интервалов по координате
int	intervalsX, intervalsY, intervalsZ;

// Шаги сетки
double	hx = 1,
hy = 1,
hz = 1;

// Координаты начала и конца области
double	begX = 0, endX = 6,
begY = 0, endY = 6,
begZ = 0, endZ = 6;

// Размер сетки
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
	
	// Максимальное кол-во итераций
	int maxiter = 10000;
	double eps = 1e-6;

public:
	SLAE();
	~SLAE();
	// Кол-во нулевых диагоналей
	int m, m2;
	// Размерность матрицы
	int n;
	// Кол-во эл-тов в диагоналях al2, au2 - n1 и al3, au3 - n2
	int n1, n2;
	std::vector <double> di; // главная диагональ
							 // Диагонали, прижатые к главной
	std::vector <double> al1, au1;
	// Диагонали, расположенные на m эл-тов выше/ниже главной диагонали и прижатых к ней
	std::vector <double> al2, au2;
	// Диагонали, расположенные на m2-тов выше/ниже главной диагонали и прижатых к ней
	std::vector <double> al3, au3;
	// Вектор правой части
	std::vector <double> f;
	// Решение системы
	std::vector <double> x;
	// Предыдущий вектор
	std::vector <double> x_pred;

	// Метод Якоби
	void Jakobi();
	// Умножение мтрицы на вектор
	void Mult_Ax(std::vector <double> &vo);
	// Относительная невязка
	double Residual();
	double summ(int i, std::vector <double> &x_next, std::vector <double> x_pred);

	double CalcF1BC(double x, double y, double z);
	double CalcF(double x, double y, double z);
};
