#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>
#include <ctime>
#include <vector>
#include "slae.h"


std::vector <Point> grid;
SLAE poissonSLAE;
int tasksPerProcess;
std::vector <int> globalNumber;


// Построение сетки в зависимости от координат начала и конца
// и количества интервалов по соответствующим координатам
void BuildGrid()
{
	double	l;
	int x, y, z;

	grid.reserve(dim);

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

	// Глобальный номер первого элемента сетки
	int number_beg = 0;

	if (size != 1)
	{
		// Если процесс не входит в область задач с лишними подзадачами
		if (residue <= rank)
		{
			// Смещаем по области задач с лишними подзадачами
			number_beg = residue * (tasksPerProcess + 1);
			// -//- по остальным
			number_beg += (rank - residue)*tasksPerProcess;

			l += number_beg * hx;
		}
		else // Если входит, то количество задач на процесс совпадает с предыдущими
		{
			number_beg = rank * tasksPerProcess;
			l += number_beg * hx;
		}

		//	Нам нужны теневые границы
		// Для первого и последнего процесса +1 теневая грань
		if ((rank == 0 || rank == size - 1))
		{
			if (rank == size - 1) 
			{ 
				l -= hx; 
				number_beg--;
			}
			tasksPerProcess++;
		}
		else // для остальных +2
		{
			tasksPerProcess += 2;
			l -= hx;
			number_beg --;
		}
	}
	int number;
	// Построение сетки для процесса
	z = begZ;
	for (int i = 0; i < intervalsZ+1; i++)
	{
		y = begY;
		for (int j = 0; j < intervalsY+1; j++)
		{
			x = l;
			for (int k = 0; k < tasksPerProcess+1; k++)
			{
				Point p;
				number = number_beg + k + (intervalsX + 1)*j + (intervalsX + 1)*(intervalsY + 1)*i;
				p.set(x, y, z, number);
				grid.push_back(p);
				x += hx;
			}
			y += hy;
		}
		z += hz;
	}
}

void CreateMatrix()
{
	double res;
	for (int i = 0; i < grid.size(); i++)
	{
		if (grid[i].x == begX || grid[i].x == endX ||
			grid[i].y == begY || grid[i].y == endY ||
			grid[i].z == begZ || grid[i].z == endZ)
			poissonSLAE.CalcF1BC(grid[i].x, grid[i].y, grid[i].z);
		else 
		{
			poissonSLAE.di[grid[i].globalNumber] = -2 / pow(hx,2) - 2 / pow(hy, 2) - 2 / pow(hz, 2);
			res = 1 / pow(hx, 2) + 1 / pow(hy, 2) + 1 / pow(hz, 2);

			poissonSLAE.al1[grid[i].globalNumber - 1] = res;
			poissonSLAE.al2[grid[i].globalNumber - poissonSLAE.m - 2] = res;
			poissonSLAE.al3[grid[i].globalNumber - poissonSLAE.m - poissonSLAE.m2 - 3] = res;

			poissonSLAE.au1[grid[i].globalNumber] = poissonSLAE.au2[grid[i].globalNumber] =
				poissonSLAE.au3[grid[i].globalNumber] = res;
			poissonSLAE.f[grid[i].globalNumber] = poissonSLAE.CalcF(grid[i].x, grid[i].y, grid[i].z);
		}
	}
}

int globalRes = 0, finalRes = 0;
int rank, size;

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

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Status st;

	unsigned int beginTime = clock();

	BuildGrid();

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
	MPI_Finalize();
	return 0;
}