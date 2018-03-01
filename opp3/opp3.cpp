#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>
#include <ctime>
#include <vector>
#include "matrix.h"


std::vector <Point> grid;

// Построение сетки в зависимости от координат начала и конца
// и количества интервалов по соответствующим координатам

void init()
{
	double	l;
	int x, y, z;
	grid.reserve(dim);

	l = abs(endX - begX);
	hx = l / intervalsX;
	l = abs(endY - begY);
	hy = l / intervalsY;
	l = abs(endZ - begZ);
	hz = l / intervalsZ;

	z = begZ;
	for (int i = 0; i < intervalsZ; i++)
	{
		y = begY;
		for (int j = 0; j < intervalsY; j++)
		{
			x = begX;
			for (int k = 0; k < intervalsX; k++)
			{
				Point p;
				p.set(x, y, z);
				grid.push_back(p);
				x += hx;
			}
			y += hy;
		}
		z += hz;
	}

}



int globalRes = 0, finalRes = 0;

int rank, size;


class Task
{
	// индексы начальных и конечных строк и столбцов в левой матрице
	int m1_begi, m1_endi, m1_begj, m1_endj;
	// -//- в правой матрице
	int m2_begi, m2_endi, m2_begj, m2_endj;
public:
	void SetM1(int bi, int ei, int bj, int ej)
	{
		m1_begi = bi;
		m1_endi = ei;
		m1_begj = bj;
		m1_endj = ej;
	};
	void SetM2(int bi, int ei, int bj, int ej)
	{
		m2_begi = bi;
		m2_endi = ei;
		m2_begj = bj;
		m2_endj = ej;
	};
	void Run();
};

void Task::Run()
{
	
}

// Очередь задач
/*std::queue <Task*> allTasks;
bool GetTask(Task **currTask)
{
	// Блокируем доступ других потоков для избежания ошибок
	// вследствие некорректной работы с очередью
	//pthread_mutex_lock(&mutex);

	// Если очередь задач пуста
	if (allTasks.empty())
	{
		// Снимаем замок
		//pthread_mutex_unlock(&mutex);
		return false;
	}
	else
	{
		// Достаём задачу из очереди
		*currTask = allTasks.front();
		allTasks.pop();
	}

	//pthread_mutex_unlock(&mutex);
	return true;
}
*/
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
	// Формирование очереди
	for (int m1_i = 0; m1_i < n / m1_blocki; m1_i++)
		for (int m1_j = 0; m1_j < m / m1_blockj; m1_j++)
			for (int m2_j = 0; m2_j < r / m2_blockj; m2_j++)
			{
				Task *t = new Task();
				t->SetM1(m1_i*m1_blocki, m1_i*m1_blocki + m1_blocki - 1,
					m1_j*m1_blockj, m1_j*m1_blockj + m1_blockj - 1);
				t->SetM2(m1_j*m1_blockj, m1_j*m1_blockj + m1_blockj - 1,
					m2_j*m2_blockj, m2_j*m2_blockj + m2_blockj - 1); // строки должны совпадать со столбцами в матрице m1
				allTasks.push(t);
			}
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Status st;

	unsigned int beginTime = clock();

	init();
	
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