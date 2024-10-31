﻿#pragma once
#include <functional>
#include <filesystem>
#include "matrix.h"
#include "matrix3D.h"
#include "vector3D.h"


struct ExpanseGrid
{
	//	Начальная точка
	double starting_point = -4.0;
	//	Конечная точка
	double ending_point = 6.0;
	//	Шаг по пространству
	double h = 0.01;
	//	Кол-во узлов сетки
	size_t nodes = 1001;
};


struct TimeGrid
{
	//	Стартовое время
	double starting_point = 0.0;
	//	Конечное время
	double ending_point = 2.5;
	//	Времяной шаг
	double tau = 0.0001;
	//	Кол-во узлов сетки
	size_t nodes = 25001;
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
	std::function<double(double)> left_borderline_ro;
	//	Левое ГУ для скорости
	std::function<double(double)> left_borderline_u;
	//	Левое ГУ для давления
	std::function<double(double)> left_borderline_p;

	//	Правое ГУ для плотности
	std::function<double(double)> right_borderline_ro;
	//	Правое ГУ для скорости
	std::function<double(double)> right_borderline_u;
	//	Правое ГУ для давления
	std::function<double(double)> right_borderline_p;
};

//	Одношаговая схема Лакс-Вендроффа для моделирования квазиодномерного течения в канале


class GasDynamicsEquation final
{
private:
	//	Данные пространственной сетки
	const ExpanseGrid _expanse_grid;
	//	Данные времянной сетки
	const TimeGrid _time_grid;

	//	Начальные условия
	const StartCondition _start_condition;
	//	Граничные условия
	const BorderlineCondition _borderline_condition;

	//	Пространственно-временая плоскость плотности
	Matrix _ro;
	//	Пространственно-временая плоскость скорости
	Matrix _u;
	//	Пространственно-временая плоскость давления
	Matrix _p;

	//	Площадь поперечного сечения трубы
	const std::vector<double> _S;

	//	Отношение шага по времяни к шагу по пространству
	const double _alpha;

public:
	//	Показатель адиабаты
	static double gamma;
	//	Искусственная вязкость Лапидуса
	static double artificial_viscosity;

public:
	GasDynamicsEquation(const ExpanseGrid& expanse_grid, const TimeGrid& time_grid, StartCondition& start_condition, BorderlineCondition& borderline_condition, const std::vector<double>& S);
	GasDynamicsEquation(const ExpanseGrid& expanse_grid, const TimeGrid& time_grid, StartCondition& start_condition, BorderlineCondition& borderline_condition, const std::function<double(double)>& S_func);
	~GasDynamicsEquation() = default;

	GasDynamicsEquation(const GasDynamicsEquation&) = delete;
	GasDynamicsEquation(GasDynamicsEquation&&) = delete;

	GasDynamicsEquation& operator=(const GasDynamicsEquation&) = delete;
	GasDynamicsEquation& operator=(GasDynamicsEquation&&) = delete;

	//	Моделирует УГД
	void solving();

	//	Записывает массив данных плотности в файл
	void writeRO(const std::filesystem::path& file_path) const;
	//	Записывает массив данных скорости в файл
	void writeU(const std::filesystem::path& file_path) const;
	//	Записывает массив данных давления в файл
	void writeP(const std::filesystem::path& file_path) const;

private:
	void initConditions();
	void postProcessing();

	bool chekStopConditions(size_t n, size_t j, double eps) const;

	double E(size_t n, size_t j) const;
	double H(size_t n, size_t j) const;
	
	Vector3D U(size_t n, size_t j) const;
	Vector3D dU(size_t n, size_t j) const;
	Vector3D U_future(size_t n, size_t j) const;

	Vector3D F(size_t n, size_t j) const;
	Vector3D Q(size_t n, size_t j) const;

	Matrix3D A(size_t n, size_t j) const;
};