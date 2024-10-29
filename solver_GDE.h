#pragma once
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
	double h = 0.1;
	//	Кол-во узлов сетки
	size_t nodes = 101;
};


struct TimeGrid
{
	//	Стартовое время
	double starting_point = 0.0;
	//	Конечное время
	double ending_point = 2.5;
	//	Времяной шаг
	double tau = 0.01;
	//	Кол-во узлов сетки
	size_t nodes = 2501;
};


struct StartCondition
{
	//	Для плотности
	std::function<double(double)> start_ro;
	//	Для скорости
	std::function<double(double)> start_u;
	//	Для давления
	std::function<double(double)> start_p;
};


struct BorderlineCondition
{
	//	Начальное ГУ для плотности
	std::function<double(double)> start_borderline_ro;
	//	Начальное ГУ для скорости
	std::function<double(double)> start_borderline_u;
	//	Начальное ГУ для давления
	std::function<double(double)> start_borderline_p;

	//	Конечное ГУ для плотности
	std::function<double(double)> end_borderline_ro;
	//	Конечное ГУ для скорости
	std::function<double(double)> end_borderline_u;
	//	Конечное ГУ для давления
	std::function<double(double)> end_borderline_p;
};

//	Одношаговая схема Лакс-Вендроффа для моделирования квазиодномерного течения в канале
class GasDynamicsEquation
{
protected:
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

	//	Отношение шага по времяни к шагу по пространству
	const double _alpha;
	//	Показатель адиабаты
	const double _gamma = 1.4;
	//	Искусственная вязкость Лапидуса
	const double _artificial_viscosity = 2.0;

public:
	GasDynamicsEquation(const ExpanseGrid& expanse_grid, const TimeGrid& time_grid, StartCondition& start_condition, BorderlineCondition& borderline_condition);
	~GasDynamicsEquation() = default;

	GasDynamicsEquation(const GasDynamicsEquation&) = delete;
	GasDynamicsEquation(GasDynamicsEquation&&) = delete;

	GasDynamicsEquation& operator=(const GasDynamicsEquation&) = delete;
	GasDynamicsEquation& operator=(GasDynamicsEquation&&) = delete;

	virtual void solving();

	void writeRO(const std::filesystem::path& file_path) const;
	void writeU(const std::filesystem::path& file_path) const;
	void writeP(const std::filesystem::path& file_path) const;

protected:
	virtual void initConditions();
	virtual void postProcessing();

	virtual	bool chekStopConditions(size_t n, size_t j, double eps);

	double E(size_t n, size_t j) const;
	double H(size_t n, size_t j) const;
	
	Vector3D U(size_t n, size_t j) const;
	Vector3D dU(size_t n, size_t j) const;
	Vector3D U_future(size_t n, size_t j) const;

	Vector3D F(size_t n, size_t j) const;
	virtual Vector3D Q(size_t n, size_t j) const;

	Matrix3D A(size_t n, size_t j) const;
};