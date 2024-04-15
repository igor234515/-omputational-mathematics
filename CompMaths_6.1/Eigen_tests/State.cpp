#include <Eigen/Dense>
#include <iostream>
int main()
{
	Eigen::Matrix3d f;
	f << 1, 2, 3,
		4, 5, 6,
		7, 8, 9;
	Eigen::Vector2f v;
	v << 1, 2;
	Eigen::Matrix<float, 2, 2> a;
	a << 3, -4,
		-5, -6;
	Eigen::Matrix2f r;
	r << 5, 6, 7, 8;
	Eigen::Matrix2f g;
	g << 10, 11, 12, 13;
	float cons = 10;
	r = a;

	Eigen::Matrix<double, 2, 3> m1;
	Eigen::Matrix<double, 2, 3> m2;
	m2 << 1, 2, 5,
		  6, 3, 4;

	m1.col(0) = m2.col(0) * 5;
	std::cout << m1;
	return 0;
}