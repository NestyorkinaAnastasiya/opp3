#include "slae.h"

SLAE::SLAE()
{
	FILE *fn;
	fn = fopen("features.txt", "r");
	fscanf(fn, "%d %d %d %d %lf", &n, &m, &m2, &maxiter, &eps);
	fclose(fn);
	n = dim;
	m = intervalsX - 1;
	m2 = (intervalsX + 1)*(intervalsY + 1) - m - 3;
	n1 = n - m - 2;
	n2 = n - m - m2 - 3;
	x.resize(n);

	// Выделение памяти для диагоналей матрицы
	au3.resize(n2);
	au2.resize(n1);
	au1.resize(n - 1);
	di.resize(n);
	al1.resize(n - 1);
	al2.resize(n1);
	al3.resize(n2);
	x_pred.resize(n);
	f.resize(n);

	for (int i = 0; i < n; i++)
		di[i] = 1;
}

SLAE::~SLAE() {}

// Относительная невязка
double SLAE::Residual()
{
	double res = 0, sum = 0.0;

	Mult_Ax(x_pred);
	for (int i = 0; i < n; i++)
		res += (f[i] - x_pred[i])*(f[i] - x_pred[i]);

	for (int i = 0; i < n; i++)
		sum += f[i] * f[i];

	res = sqrt(res / sum);
	return res;
}

double SLAE::summ(int i, std::vector <double> &x_next, std::vector <double> x_pred)
{
	double sum;
	sum = di[i] * x_pred[i];
	// Нижний треугольник
	// Если в столбце есть эл-ты al1
	if (i >= 1)
	{	// Смещаемся по строке
		int j = i - 1;
		sum += al1[j] * x_next[j];
		// Если в столбце есть эл-ты al2
		if (i >= n - n1)
		{	// Смещаемся по строке
			j -= m + 1;
			sum += al2[j] * x_next[j];
			// Если в столбце есть эл-ты al3
			if (i >= n - n2)
			{	// Смещаемся по строке
				j -= m2 + 1;
				sum += al3[j] * x_next[j];
			}
		}
	}
	// Верхний треугольник
	// Если в строке есть эл-ты au1
	if (i < n - 1)
	{	// Смещаемся по столбцу
		int j = i + 1;
		sum += au1[i] * x_pred[j];
		// Если в строке есть эл-ты au2
		if (i < n1)
		{	// Смещаемся по столбцу
			j += m + 1;
			sum += au2[i] * x_pred[j];
			// Если в строке есть эл-ты au3
			if (i < n2)
			{	// Смещаемся по столбцу
				j += m2 + 1;
				sum += au3[i] * x_pred[j];
			}
		}
	}
	return sum;
}

// Метод Якоби
void SLAE::Jakobi()
{
	double res, var, sum;
	int k;
	for (int i = 0; i < n; i++)
		x[i] = 0;
	res = Residual();
	for (k = 0; res >= eps && k < maxiter; k++)
	{
		for (int i = 0; i < n; i++)
			x_pred[i] = x[i];

		for (int i = 0; i < n; i++)
		{
			var = 1 / di[i];
			sum = summ(i, x_pred, x_pred);
			x[i] += var * (f[i] - sum);
		}

		res = Residual();
	}
}

// Умножение матрицы на вектор
void SLAE::Mult_Ax(std::vector <double> &vo)
{
	// Умножение главной диагонали
	for (int i = 0; i < n; i++)
		vo[i] = di[i] * x[i];
	// Умножение диагоналей, прилегающих к главной
	for (int i = 0, j = 1; i < n - 1; i++, j++)
	{
		vo[i] += au1[i] * x[j];
		vo[j] += al1[i] * x[i];
	}
	// Умножение диагоналей, смещённых на m от главной и прилежащих
	for (int i = 0, j = m + 2; i < n1; i++, j++)
	{
		vo[i] += au2[i] * x[j];
		vo[j] += al2[i] * x[i];
	}
	// Умножение диагоналей, смещённых на m+m2+2 от главной и прилежащих
	for (int i = 0, j = m + m2 + 3; i < n2; i++, j++)
	{
		vo[i] += au3[i] * x[j];
		vo[j] += al3[i] * x[i];
	}
}

double SLAE::CalcF1BC(double x, double y, double z)
{
	return 4;
}

double SLAE::CalcF(double x, double y, double z)
{
	return 4;
}