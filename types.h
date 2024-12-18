#pragma once
#include "matrix.h"
#include "matrix3D.h"
#include "vector3D.h"

struct ExpanseGrid
{
	//	��� �� ������������
	double h = 0.01;
	//	���-�� ����� �����
	size_t nodes = 1001;
	//	��������� �����
	double starting_point = -4.0;
};


struct TimeGrid
{
	//	�������� ���
	double tau = 0.0001;
	//	���-�� ����� �����
	size_t nodes = 25001;
	//	��������� �����
	double starting_point = 0.0;
};


struct StartCondition
{
	//	��������� ������� ��� ���������
	std::function<double(double)> start_ro;
	//	��������� ������� ��� ��������
	std::function<double(double)> start_u;
	//	��������� ������� ��� ��������
	std::function<double(double)> start_p;
};


struct BorderlineCondition
{
	//	����� �� ��� ���������
	double left_borderline_ro;
	//	����� �� ��� ��������
	double left_borderline_u;
	//	����� �� ��� ��������
	double left_borderline_p;

	//	������ �� ��� ���������
	double right_borderline_ro;
	//	������ �� ��� ��������
	double right_borderline_u;
	//	������ �� ��� ��������
	double right_borderline_p;
};