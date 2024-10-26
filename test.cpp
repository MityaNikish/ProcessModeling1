#include "test.h"
#include "vector3D.h"
#include <iostream>

void TEST_vector(void)
{
	Vector3D V1{ 6, 2, 4 };
	Vector3D V2{ 3, 7, 0 };

	double multiply = V1 * (V2 / 0);
	
	std::cout << multiply << std::endl;

}



void TEST_all(void)
{
	TEST_vector();
}