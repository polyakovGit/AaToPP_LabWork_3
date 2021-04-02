#include <stdio.h> 
#include <omp.h>
#include <ctime>
#include <cmath>
#include <corecrt_math_defines.h>

int WorkSum(int currI, int n) {
	int sum = 0;
	for (int i = currI; i < n; ++i)
		sum += 1;
	return sum;
}

void Task1(int N) {
	int sum = 0;
#pragma omp parallel num_threads(2) reduction(+:sum)
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
#pragma omp parallel num_threads(k) reduction(+:sum)
	{
		//sum += WorkSum(omp_get_thread_num()*N/omp_get_num_threads(),(omp_get_thread_num()+1)*N/omp_get_num_threads());
		int currThread = omp_get_thread_num();
		sum += WorkSum(currThread * N / k, (currThread + 1) * N / k);
		printf("[%d]: Sum=%d\n", currThread, sum);
	}
	printf("Sum = %d\n", sum);
}

void Task3(int N, int k) {
	int sum = 0;
#pragma omp parallel num_threads(k)
	{
#pragma omp for reduction(+:sum)
		for (int i = 0; i < N; ++i)
			sum += 1;
		printf("[%d]: Sum=%d\n", omp_get_thread_num(), sum);
	}
	printf("Sum=%d\n", sum);
}

void Task4() {
	int k = 4, N = 10;
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

double Task5Serial(long long n) {
	double sum = 0, x, step = 1. / n;//a=0,b=1 границы интеграла;(b-a)/(double)n в общем виде
	for (long long i = 0; i < n - 1; ++i) {//n-1 для левых
		x = i * step;//i+0.5 если средних; a+i*step но здесь частный случай для pi 0, сумма не нужна 
		sum += 4 / (1 + x * x);
	}
	return sum / n;
}

int main() {

	int N = 100000000, sum = 0;
	const int arrNumbersSize = 7;
	int arrNumbers[arrNumbersSize] = { 1,2,5,10,100,1000,N };

	clock_t start, end;

	//printf("Task1\n");
	//for (int i = 0; i < arrNumbersSize; ++i) {
	//	start = clock();
	//	Task1(arrNumbers[i]);
	//	end = clock();
	//	printf("time %f\n", (double)(end - start) / CLOCKS_PER_SEC);
	//}

	//start = clock();
	//Task2(N, 4);
	//end = clock();
	//printf("time %f\n", (double)(end - start) / CLOCKS_PER_SEC);
	//start = clock();
	//Task3(N, 4);
	//end = clock();
	//printf("time %f\n", (double)(end - start) / CLOCKS_PER_SEC);

	//const int  arrThreadsSize = 4;
	//int arrThreads[arrThreadsSize] = { 2,4,8,16 };
	//printf("Task2\n");
	//for (int k = 0; k < arrThreadsSize; ++k) {
	//	for (int i = 0; i < arrNumbersSize; ++i) {
	//		start = clock();
	//		Task2(arrNumbers[i], arrThreads[k]);
	//		end = clock();
	//		printf("time %f\n", (double)(end - start) / CLOCKS_PER_SEC);
	//	}
	//}

	//printf("Task3\n");
	//for (int k = 0; k < arrThreadsSize; ++k) {
	//	for (int i = 0; i < arrNumbersSize; ++i) {
	//		start = clock();
	//		Task3(arrNumbers[i], arrThreads[k]);
	//		end = clock();
	//	}
	//	printf("time %f\n", (double)(end - start) / CLOCKS_PER_SEC);
	//}

	//printf("Task4\n");
	//Task4();

	printf("Task5\n");
	long long n = 2;
	int i = 0;
	double eps = 0.001;
	for (int k = 1; k < 20; ++k) {
		do {
			start = clock();
			double myPi = Task5(n, k);
			end = clock();
			if (fabs(M_PI - myPi) > eps)n *= 2;
			else {
				//printf("ITERATION %d\nN=%lld\neps   %.10f\n FABS %.10f\n   PI %.10f\nvalue %.10f\ntime %f\n\n",i++, n, eps,
				//	fabs(M_PI - myPi), M_PI, myPi, (double)(end - start) / CLOCKS_PER_SEC);
				printf("%.10f\nk= %d\n%.10f\n\n", myPi,k, Task5Serial(n));//проверить с разным количеством нитей значение pi
				eps *= 0.1;
			}
		} while ((double)(end - start) / CLOCKS_PER_SEC < 120);
	}

	//добавить массив нитей, если нынешний тред==arrThreads[k]
//	for (int i = 0; i < arrSize; ++i) {
//		int sum = 0;
//		start = clock();
//#pragma omp parallel num_threads(2) reduction(+:sum)
//		{
//			int currThread = omp_get_thread_num();
//			if (currThread == 0) {
//				sum += WorkSum(0, arrNumbers[i] / 2);
//				printf("[%d]: Sum=%d\n", currThread, sum);
//			}
//			if (currThread == 1) {
//				sum += WorkSum(arrNumbers[i] / 2, arrNumbers[i]);
//				printf("[%d]: Sum=%d\n", currThread, sum);
//			}
//		}
//		printf("Sum = %d\n", sum);
//		end = clock();
//		printf("time %f\n", (double)(end - start) / CLOCKS_PER_SEC);
//	}


	sum = 0;
	start = clock();
	for (int i = 0; i < N; ++i)
		sum += 1;
	end = clock();
	printf("%d\ntime SERIAL %f", sum, (double)(end - start) / CLOCKS_PER_SEC);
	return 0;
}