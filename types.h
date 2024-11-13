#pragma once
#include "matrix.h"
#include "matrix3D.h"
#include "vector3D.h"

struct ExpanseGrid
{
	//	Шаг по пространству
	double h = 0.01;
	//	Кол-во узлов сетки
	size_t nodes = 1001;
	//	Начальная точка
	double starting_point = -4.0;
	//	Конечная точка
	double ending_point = 6.0;
};


struct TimeGrid
{
	//	Времяной шаг
	double tau = 0.0001;
	//	Кол-во узлов сетки
	size_t nodes = 25001;
	//	Стартовое время
	double starting_point = 0.0;
	//	Конечное время
	double ending_point = 2.5;
};


struct StartCondition
{
	//	Начальное условие для плотности
	std::function<double(double)> start_ro;
	//	Начальное условие для скорости
	std::function<double(double)> start_u;
	//	Начальное условие для давления
	std::function<double(double)> start_p;
};


struct BorderlineCondition
{
	//	Левое ГУ для плотности
	double left_borderline_ro;
	//	Левое ГУ для скорости
	double left_borderline_u;
	//	Левое ГУ для давления
	double left_borderline_p;

	//	Правое ГУ для плотности
	double right_borderline_ro;
	//	Правое ГУ для скорости
	double right_borderline_u;
	//	Правое ГУ для давления
	double right_borderline_p;
};