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
		//sum+=Work(omp_get_num_thread()*N/omp_get_threads_num();++omp_get_num_thread()*N/omp_get_threads_num()); для k нитей
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

int main() {
	int N = 1000000000, sum = 0, i;
	const int arrSize = 7;
	int arrNumbers[arrSize] = { 1,2,5,10,100,1000,N };

	clock_t start, end;

	printf("Task1");
	for (int i = 0; i < arrSize; ++i) {
		start = clock();
		Task1(arrNumbers[i]);
		end = clock();
		printf("time %f\n", (double)(end - start) / CLOCKS_PER_SEC);
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
	for (i = 0; i < N; ++i)
		sum += 1;
	end = clock();
	printf("%d\ntime SERIAL %f", sum, (double)(end - start) / CLOCKS_PER_SEC);
	return 0;
}