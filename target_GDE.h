#pragma once
#include <vector>
#include <filesystem>
#include "types.h"

struct Values
{
	//Давление в каналe
	double p;
	//Плотность в каналe
	double ro;
	//Скорость в каналe
	double u;
};

class TargetGDE
{
	//	Данные пространственной сетки
	const ExpanseGrid _expanse_grid;

	//	Граничные условия
	const BorderlineCondition _borderline_condition;

	//Площадь сечения канала
	const std::vector<double> _S;

	//	Результат вычислений
	std::vector<Values> _data;

	//Критическое сечение
	double _x_crit;
	//Площадь критического сечения
	double _S_crit;

	//Координаты скачка уплотнения
	double _x_s;

	//	Показатель адиабаты
	static double _gamma;

public:
	TargetGDE(const ExpanseGrid& expanse_grid, const BorderlineCondition& borderline_condition, const std::vector<double>& S, const double x_s = 1);
	TargetGDE(const ExpanseGrid& expanse_grid, const BorderlineCondition& borderline_condition, const std::function<double(double)>& S_func, const double x_s = 1);

	//	Моделирует УГД
	void solving();

	//	Записывает массив данных плотности в файл
	void writeRO(const std::filesystem::path& file_path) const;
	//	Записывает массив данных скорости в файл
	void writeU(const std::filesystem::path& file_path) const;
	//	Записывает массив данных давления в файл
	void writeP(const std::filesystem::path& file_path) const;

	//	Вернуть массив данных плотности
	std::vector<double> getRO() const;
	//	Вернуть массив данных скорости
	std::vector<double> getU() const;
	//	Вернуть массив данных давления
	std::vector<double> getP() const;

private:
	void initBottleneck();
};