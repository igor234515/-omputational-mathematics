#include <iostream>
const int DIM1 = 3;
const int DIM2 = 5;

int ary[DIM1][DIM2];

int main()
{
	for (int i = 0; i < DIM1; ++i)
	{
		for (int j = 0; j < DIM2; j++)
		{
			std::cout << ary[i][j] << "\t";
		}
		std::cout << std::endl;
	}
	return 0;
}