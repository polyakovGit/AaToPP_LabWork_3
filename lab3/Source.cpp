#include <stdio.h> 
#include <omp.h>
#include <ctime>

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
		//sum+=Work(omp_get_thread_num()*N/omp_get_num_threads();++omp_get_thread_num()*N/omp_get_num_threads()); для k нитей
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