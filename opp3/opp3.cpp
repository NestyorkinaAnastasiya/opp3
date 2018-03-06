#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>
#include <math.h>
#include <ctime>
#include <vector>

int size = 1, rank = 0;
// Шаги сетки
double	hx = 0.5,
hy = 0.5,
hz = 0.5;

// Координаты начала и конца области
double	begX = 0, endX = 2,
begY = 0, endY = 2,
begZ = 0, endZ = 2;

struct Point 
{
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

struct Grid 
{
	// Максимальное кол-во итераций
	int maxiter = 1000;
	double eps = 1e-14;

	// Количество интервалов по координате
	int	intervalsX, intervalsY, intervalsZ;
	// Количество интервалов по х на процесс
	int tasksPerProcess;

	// Координаты начала и конца области для процесса
	double	firstX, lastX;

	// Размер сетки
	int dim;

	Grid();
	std::vector<Point> points;
	std::vector<double> oldU, newU, f;
	double *oldData, *newData;
	std::vector<int> numbersOfLeftShadowBorders;
	std::vector<int> numbersOfRightShadowBorders;
	
	void Calculate1Node(int i);
	void CalculateU();
	bool BelongToShadowBorders(int node);
	double CalcF(double x, double y, double z);
	double CalcF1BC(double x, double y, double z);

} grid;


// Построение сетки в зависимости от координат начала и конца
// и количества интервалов по соответствующим координатам
Grid::Grid()
{
	double	l;
	double x, y, z;

	grid.points.reserve(dim);

	// Бьём задачи по координате х
	l = abs(endX - begX);
	intervalsX = l / hx;
	tasksPerProcess = intervalsX / size;

	// Если количество интервалов не кратно числу процессов
	// то распределяем оставшиеся задачи по первым процессам
	int residue = intervalsX % size;
	if (residue)
		for (int i = 0; i < residue && rank == i; i++)
			tasksPerProcess++;

	l = abs(endY - begY);
	intervalsY = l / hy;
	l = abs(endZ - begZ);
	intervalsZ = l / hz;
	dim = (intervalsX + 1) * (intervalsY + 1) * (intervalsZ + 1);
	// Глобальный номер первого элемента сетки
	int number_beg = 0;
	firstX = begX;
	if (size != 1)
	{
		// Если процесс не входит в область задач с лишними подзадачами
		if (residue <= rank)
		{
			// Смещаем по области задач с лишними подзадачами
			number_beg = residue * (tasksPerProcess + 1);
			// -//- по остальным
			number_beg += (rank - residue)*tasksPerProcess;

			firstX = number_beg * hx;
		}
		else // Если входит, то количество задач на процесс совпадает с предыдущими
		{
			number_beg = rank * tasksPerProcess;
			firstX = number_beg * hx;
		}

		//	Нам нужны теневые границы
		// Для первого и последнего процесса +1 теневая грань
		if ((rank == 0 || rank == size - 1))
		{
			if (rank == size - 1) 
			{ 
				firstX -= hx;
				number_beg--;
			}
			tasksPerProcess++;
		}
		else // для остальных +2
		{
			tasksPerProcess += 2;
			firstX -= hx;
			number_beg --;
		}
	}
	int number;
	// Построение сетки для процесса
	z = begZ;
	// Начало и конец области с теневыми
	lastX = firstX + hx * tasksPerProcess;
	for (int i = 0; i < intervalsZ+1; i++)
	{
		y = begY;
		for (int j = 0; j < intervalsY+1; j++)
		{
			x = firstX;
			for (int k = 0; k < tasksPerProcess+1; k++)
			{
				Point p;
				number = number_beg + k + (intervalsX + 1)*j + (intervalsX + 1)*(intervalsY + 1)*i;
				p.set(x, y, z, number);
				points.push_back(p);
				oldU.push_back(1);

				// Расчёт правых частей
				if ((abs(x - begX) < eps) || (abs(x - endX) < eps) ||
					(abs(y - begY) < eps) || (abs(y - endY) < eps) ||
					(abs(z - begZ) < eps) || (abs(z - endZ) < eps))
					f.push_back(CalcF1BC(x, y, z));
				else
				{
					f.push_back(CalcF(x, y, z));

					// Формирование теневых границ
					if (k == 0)
						numbersOfLeftShadowBorders.push_back(oldU.size() - 1);
					if (k == tasksPerProcess)
						numbersOfRightShadowBorders.push_back(oldU.size() - 1);
				}

				x += hx;
			}
			y += hy;
		}
		z += hz;
	}
	newU.resize(grid.oldU.size());
}

// i - номер узла, для которого происходит расчёт
void Grid::Calculate1Node(int i)
{
	double result;
	// Смещения, для расчёта соседей по y и по z
	int offsetY = tasksPerProcess + 1, offsetZ = (tasksPerProcess + 1)*(intervalsY + 1);

	// Если узел принадлежит к границе, то накладываются первые краевые условия
	if ((abs(grid.points[i].x - begX) < eps) || (abs(grid.points[i].x - endX) < eps) ||
		(abs(grid.points[i].y - begY) < eps) || (abs(grid.points[i].y - endY) < eps) ||
		(abs(grid.points[i].z - begZ) < eps) || (abs(grid.points[i].z - endZ) < eps))
		newData[i] = grid.f[i];
	else // иначе если узел не принадлежит тененевой границе, 
		 // то производим расчёт значения в данном узле
		if (!BelongToShadowBorders(i))
		{
			result = (oldData[i - 1] + oldData[i + 1]) / pow(hx, 2) +
				(oldData[i - offsetY] + oldData[i + offsetY]) / pow(hy, 2) +
				(oldData[i - offsetZ] + oldData[i + offsetZ]) / pow(hz, 2) - f[i];
			result *= pow(hx, 2)*pow(hy, 2)*pow(hz, 2) / 2 / (pow(hy, 2)*pow(hz, 2) +
				pow(hz, 2)*pow(hx, 2) + pow(hx, 2)*pow(hy, 2));
			newData[i] = result;
		}
}

void Grid::CalculateU()
{
	double residual = 1;
	std::vector<double> leftBorder(numbersOfLeftShadowBorders.size());
	std::vector<double> rightBorder(numbersOfLeftShadowBorders.size());
	MPI_Request sendReqLeft, recvReqLeft, sendReqRight, recvReqRight;
	int iteration;
	for (iteration = 0; iteration < maxiter && residual > eps; iteration++)
	{	// Если итерация чётная, то в  oldU старые значения, а в newU - новые
		// иначе наоборот

		if (iteration % 2 == 0)
		{
			oldData = grid.oldU.data();
			newData = grid.newU.data();
		}
		else
		{
			oldData = grid.newU.data();
			newData = grid.oldU.data();
		}
		
		if (rank != 0)
		{	
			// Получение левой теневой границы (без неё не возможен счёт)
			MPI_Recv(leftBorder.data(), leftBorder.size(), MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			for (int i = 0; i < numbersOfLeftShadowBorders.size(); i++)
			{
				// Записали полученную левую теневую границу
				oldData[numbersOfLeftShadowBorders[i]] = leftBorder[i];
				Calculate1Node(numbersOfLeftShadowBorders[i] + 1);
				// Сформировали левую границу области самого процесса
				leftBorder[i] = newData[numbersOfLeftShadowBorders[i] + 1];
			}
			// Передача левой границы области процесса как теневой для левого соседнего процесса
			MPI_Isend(leftBorder.data(), leftBorder.size(), MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &sendReqLeft);
		}

		if (rank != size - 1)
		{
			MPI_Recv(rightBorder.data(), rightBorder.size(), MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			for (int i = 0; i < numbersOfRightShadowBorders.size(); i++)
			{
				Calculate1Node(numbersOfRightShadowBorders[i] - 1);
				rightBorder[i] = newData[numbersOfRightShadowBorders[i] - 1];
			}
			MPI_Isend(rightBorder.data(), rightBorder.size(), MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &sendReqRight);
		}

		// Проходимся по сетке (OpenMP)
		for (int i = 0; i < grid.oldU.size() && abs(grid.points[i].x - firstX) < eps 
			&& abs(grid.points[i].x - lastX) < eps; i++)
			Calculate1Node(i);

		// Расчёт абсолютной погрешности между текущим и предыдущим решениями
		residual = 0;
		for (int i = 0; i < oldU.size(); i++)
			residual += pow(newData[i] - oldData[i], 2);
		residual = sqrt(residual);
		if(rank != 0 && rank != size - 1)
		{
			MPI_Wait(&sendReqLeft, MPI_STATUS_IGNORE);
			MPI_Wait(&sendReqRight, MPI_STATUS_IGNORE);
		}
		else if (rank == 0)
			MPI_Wait(&sendReqRight, MPI_STATUS_IGNORE);
		else
			MPI_Wait(&sendReqLeft, MPI_STATUS_IGNORE);
		// ОБМЕНЯТЬСЯ СООБЩЕНИЯМИ про невязку (флаг)18+
		
	}
	// Если новые значения оказались в старом векторе, копируем результат новый вектор
	if ((iteration - 1) % 2 != 0)
		newU = oldU;
}

// Функция принадлежности узла node теневой границе
bool Grid::BelongToShadowBorders(int node)
{
	for (auto i : numbersOfLeftShadowBorders)
		if (i == node) return true;
	for (auto i : numbersOfRightShadowBorders)
		if (i == node) return true;
	return false;
}

// Учёт первых краевых условий
double Grid::CalcF1BC(double x, double y, double z)
{
	return 2;
}

// Правая часть уравнения Пуассона
double Grid::CalcF(double x, double y, double z)
{
	return 0;
}

/*
// Функция вычислительного потока
void* worker(void* me)
{	//Текущая задача
Task *currTask;
int id = *((int*)me);
unsigned int beginTime = clock();
// Пока есть свои задачи - выполняем свои
while (GetTask(&currTask))
{
currTask->Run();
delete[] currTask;
countOfTasks[id] = countOfTasks[id] + 1;
}
unsigned int endTime = clock();
timeOfWork[id] = (endTime - beginTime) * 1000 / CLOCKS_PER_SEC;
return 0;
}*/

int main(int argc, char **argv)
{
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//	MPI_Comm_size(MPI_COMM_WORLD, &size);
//	MPI_Status st;


	unsigned int beginTime = clock();
	grid.CalculateU();
	unsigned int endTime = clock();

	/*FILE *fo;
	fopen_s(&fo, "result56.txt", "a+b");
	fprintf(fo, "------------------------------------------------------------\n\n n: %d\tdimBlock: %d\tcountOfThreads: %d\n\n", n, m1_blocki, countOfWorkers);
	unsigned int sumOfTime = 0;
	for (int i = 0; i < countOfWorkers; i++)
	sumOfTime += timeOfWork[i];
	fprintf(fo, "Time: %d\n\n", (endTime - beginTime) * 1000 / CLOCKS_PER_SEC);
	for (int i = 0; i < countOfWorkers; i++)
	fprintf(fo, "process: %d\ttime: %d\tcountOfTasks: %d\n", ids[i], timeOfWork[i], countOfTasks[i]);
	fprintf(fo, "\n------------------------------------------------------------\n");
	fclose(fo);*/

	//MPI_Reduce(&globalRes, &resultOfGlobalTask, 1, MPI_INT, MPI_SUM, 0, currentComm);

	/*	if (rank == 0)
	std::cerr << "Global result = " << resultOfGlobalTask << "\n";*/
//	MPI_Finalize();
	return 0;
}