#include<Eigen/Dense>
#include <iostream>


int main()
{
	int aray[2][2];
	aray[0][0] = 1;
	aray[0][1] = 2;
	for (int i = 0; i < 2; i++)
	{
		std::cout << aray[0][i] << "\t";
	}
	return 0;
}