#include <stdio.h> 
#include <omp.h>
#include <ctime>
#include <cmath>
#include <corecrt_math_defines.h>
#include <iostream>
#include <fstream>

int WorkSum(int currI, int n) {
	int sum = 0;
	for (int i = currI; i < n; ++i)
		sum += 1;
	return sum;
}

void Task1(int N) {
	int sum = 0;
#pragma omp parallel num_threads(2) 
	{
		int currThread = omp_get_thread_num();
		//sum+=Work(omp_get_thread_num()*N/omp_get_num_threads();(omp_get_thread_num()+1)*N/omp_get_num_threads()); для k нитей
		if (currThread == 0) {
			sum += WorkSum(0, N / 2);
			printf("[%d]: Sum=%d\n", currThread, sum);
		}
		if (currThread == 1) {
			sum += WorkSum(N / 2, N);
			printf("[%d]: Sum=%d\n", currThread, sum);
		}
	}
	printf("Sum = %d\n", sum);
}

void Task2(int N, int k) {
	int sum = 0;
#pragma omp parallel num_threads(k)
	{
		int localSum = 0;
		//sum += WorkSum(omp_get_thread_num()*N/omp_get_num_threads(),(omp_get_thread_num()+1)*N/omp_get_num_threads());
		int currThread = omp_get_thread_num();
		localSum += WorkSum(currThread * N / k, (currThread + 1) * N / k);
		printf("[%d]: Sum=%d\n", currThread, localSum);
		sum += localSum;
	}
	printf("Sum = %d\n", sum);
}

void Task3(int N, int k) {
	int sum = 0;
#pragma omp parallel num_threads(k) reduction(+:sum)
	{
#pragma omp for 
		for (int i = 0; i < N; ++i)
			sum += 1;
		printf("[%d]: Sum=%d\n", omp_get_thread_num(), sum);
	}
	printf("Sum=%d\n", sum);
}

void Task4() {
	int k = 4, N = 10;
	int sum = 0;
	printf("schedule static\n");
#pragma omp parallel num_threads(k) reduction(+:sum)
	{
#pragma omp for schedule(static)
		for (int i = 1; i <= N; ++i) {
			sum += 1;
			printf("Iteraion [%d] thread [%d]\n", i, omp_get_thread_num());
		}
		printf("[%d]: Sum=%d\n", omp_get_thread_num(), sum);
	}
	printf("Sum=%d\n", sum);

	sum = 0;
	printf("schedule static 1\n");
#pragma omp parallel num_threads(k) reduction(+:sum)
	{
#pragma omp for schedule(static,1)
		for (int i = 1; i <= N; ++i) {
			sum += 1;
			printf("Iteraion [%d] thread [%d]\n", i, omp_get_thread_num());
		}
		printf("[%d]: Sum=%d\n", omp_get_thread_num(), sum);
	}
	printf("Sum=%d\n", sum);

	sum = 0;
	printf("schedule static 2\n");
#pragma omp parallel num_threads(k) reduction(+:sum)
	{
#pragma omp for schedule(static,2)
		for (int i = 1; i <= N; ++i) {
			sum += 1;
			printf("Iteraion [%d] thread [%d]\n", i, omp_get_thread_num());
		}
		printf("[%d]: Sum=%d\n", omp_get_thread_num(), sum);
	}
	printf("Sum=%d\n", sum);

	sum = 0;
	printf("schedule dynamic\n");
#pragma omp parallel num_threads(k) reduction(+:sum)
	{
#pragma omp for schedule(dynamic)
		for (int i = 1; i <= N; ++i) {
			sum += 1;
			printf("Iteraion [%d] thread [%d]\n", i, omp_get_thread_num());
		}
		printf("[%d]: Sum=%d\n", omp_get_thread_num(), sum);
	}
	printf("Sum=%d\n", sum);

	sum = 0;
	printf("schedule dynamic 2\n");
#pragma omp parallel num_threads(k) reduction(+:sum)
	{
#pragma omp for schedule(dynamic,2)
		for (int i = 1; i <= N; ++i) {
			sum += 1;
			printf("Iteraion [%d] thread [%d]\n", i, omp_get_thread_num());
		}
		printf("[%d]: Sum=%d\n", omp_get_thread_num(), sum);
	}
	printf("Sum=%d\n", sum);

	sum = 0;
	printf("schedule guided\n");
#pragma omp parallel num_threads(k) reduction(+:sum)
	{
#pragma omp for schedule(guided)
		for (int i = 1; i <= N; ++i) {
			sum += 1;
			printf("Iteraion [%d] thread [%d]\n", i, omp_get_thread_num());
		}
		printf("[%d]: Sum=%d\n", omp_get_thread_num(), sum);
	}
	printf("Sum=%d\n", sum);

	sum = 0;
	printf("schedule guided 2\n");
#pragma omp parallel num_threads(k) reduction(+:sum)
	{
#pragma omp for schedule(guided,2)
		for (int i = 1; i <= N; ++i) {
			sum += 1;
			printf("Iteraion [%d] thread [%d]\n", i, omp_get_thread_num());
		}
		printf("[%d]: Sum=%d\n", omp_get_thread_num(), sum);
	}
	printf("Sum=%d\n", sum);
}

double Task5(long long n, int k) {
	double sum = 0, x, step = 1. / n;//a=0,b=1 границы интеграла;(b-a)/(double)n в общем виде
#pragma omp parallel for num_threads(k) reduction (+:sum) private(x)
	for (long long i = 0; i < n - 1; ++i) {//n-1 для левых
		x = i * step;//i+0.5 если средних; a+i*step но здесь частный случай для pi 0, сумма не нужна 
		sum += 4 / (1 + x * x);
	}
	return sum / n;
}

//double Task5Serial(long long n) {
//	double sum = 0, x, step = 1. / n;//a=0,b=1 границы интеграла;(b-a)/(double)n в общем виде
//	for (long long i = 0; i < n - 1; ++i) {//n-1 для левых
//		x = i * step;//i+0.5 если средних; a+i*step но здесь частный случай для pi 0, сумма не нужна 
//		sum += 4 / (1 + x * x);
//	}
//	return sum / n;
//}

int main() {

	int N = 100000000, sum = 0;
	const int arrNumbersSize = 7;
	int arrNumbers[arrNumbersSize] = { 1,2,5,10,100,1000,N };

	clock_t start, end;

	printf("Task1\n");
	for (int i = 0; i < arrNumbersSize; ++i) {
		start = clock();
		Task1(arrNumbers[i]);
		end = clock();
		printf("time %f\n", (double)(end - start) / CLOCKS_PER_SEC);
	}

	const int  arrThreadsSize = 4;
	int arrThreads[arrThreadsSize] = { 2,4,8,16 };

	printf("Task2\n");
	for (int k = 0; k < arrThreadsSize; ++k) {
		for (int i = 0; i < arrNumbersSize; ++i) {
			start = clock();
			Task2(arrNumbers[i], arrThreads[k]);
			end = clock();
			printf("time %f\n", (double)(end - start) / CLOCKS_PER_SEC);
		}
	}

	printf("Task3\n");
	for (int k = 0; k < arrThreadsSize; ++k) {
		for (int i = 0; i < arrNumbersSize; ++i) {
			start = clock();
			Task3(arrNumbers[i], arrThreads[k]);
			end = clock();
		}
		printf("time %f\n", (double)(end - start) / CLOCKS_PER_SEC);
	}

	printf("Task4\n");
	Task4();

	printf("Task5\n");
	std::ofstream OutputData;
	OutputData.open("OutputData.csv");
	long long n = 2;
	int i = 0;
	double eps = 0.001, time;
	printf("\t");
	OutputData << ';';
	for (int k = 1; k <= 20; ++k) {
		printf("\t%d\t", k);
		OutputData << ';' << k;
	}
	do {
		printf("\n%.10f\t", eps);
		OutputData << '\n' << eps << ';';
		for (int k = 1; k <= 20;) {
			start = clock();
			double myPi = Task5(n, k);
			end = clock();
			time = (double)(end - start) / CLOCKS_PER_SEC;
			if (fabs(M_PI - myPi) > eps)n *= 2;
			else {
				printf("%.10f\t", time);
				OutputData << time << ';';
				if (time > 120)break;
				++k;
			}
		}
		eps *= 0.1;
	} while (time <= 120);

	sum = 0;
	start = clock();
	for (int i = 0; i < N; ++i)
		sum += 1;
	end = clock();
	printf("%d\ntime SERIAL %f", sum, (double)(end - start) / CLOCKS_PER_SEC);
	return 0;
}