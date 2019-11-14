#include<stdio.h> 
#include<stdlib.h> 

#define max(x, y) (((x) > (y)) ? (x) : (y))
#define min(x, y) (((x) < (y)) ? (x) : (y))

void S1(int t, int i, int j, int k)
{
	static int c0Last = 0;
	static int superCounter = 0;
	static int counter = 0;
	if (c0Last!=t)
	{
		printf("^^^%d^^^   %d    ", counter, superCounter);
		printf(" t: %d\n", t);
		superCounter += counter;
		counter = 0;
	}
	//printf("c0: %d  :  c1: %d  :  c2: %d  :   x:%d\n", t, i, j, k);
	printf("(I%d%d) edge %s node {} (I%d%d)\n", i, k, (k==j) ? "[loop above]" : "[bend right, above]", i, j);
	printf("(I%d%d) edge %s node {} (I%d%d)\n", k, j, (k==i) ? "[loop above]" : "[bend right, above]", i, j);
	counter ++;
	c0Last = t;
}
int main(void)
{
	int n = 4;
	for (int c0 = 1; c0 < 3 * n - 1; c0 += 1) 
	{
		for (int c1 = 1; c1 <= -2 * n + c0 + 1; c1 += 1)
			for (int c2 = 1; c2 <= n; c2 += 1) 
			{
				if (2 * n + c1 + c2 >= c0 + 2) {
					S1(c0, c1, c2, -n + c0 - c1 - c2 + 2);
				} else
					S1(c0, c1, c2, -2 * n + c0 - c1 - c2 + 2);
			}
		for (int c1 = max(1, -2 * n + c0 + 2); c1 <= min(n, c0); c1 += 1)
			for (int c2 = 1; c2 <= min(n, c0 - c1 + 1); c2 += 1) 
			{
				if (n + c1 + c2 >= c0 + 2) {
					S1(c0, c1, c2, c0 - c1 - c2 + 2);
				} else
					S1(c0, c1, c2, -n + c0 - c1 - c2 + 2);
			}
	}
	printf("FOR1END\n");
	for (int c0 = 3 * n - 1; c0 < 4 * n - 1; c0 += 1) 
	{
		for (int c1 = 1; c1 <= -3 * n + c0 + 1; c1 += 1)
			for (int c2 = -3 * n + c0 - c1 + 2; c2 <= n; c2 += 1)
				S1(c0, c1, c2, -2 * n + c0 - c1 - c2 + 2);
			for (int c1 = -3 * n + c0 + 2; c1 <= n; c1 += 1)
				for (int c2 = 1; c2 <= n; c2 += 1) 
				{
					if (2 * n + c1 + c2 >= c0 + 2) {
						S1(c0, c1, c2, -n + c0 - c1 - c2 + 2);
					} else
						S1(c0, c1, c2, -2 * n + c0 - c1 - c2 + 2);
				}
	}
	printf("FOR2END\n");
	for (int c0 = 4 * n - 1; c0 < 5 * n - 1; c0 += 1)
		for (int c1 = -4 * n + c0 + 2; c1 <= n; c1 += 1)
			for (int c2 = -3 * n + c0 - c1 + 2; c2 <= n; c2 += 1)
				S1(c0, c1, c2, -2 * n + c0 - c1 - c2 + 2);
	printf("FOR3END\n");

	return 0;
}
