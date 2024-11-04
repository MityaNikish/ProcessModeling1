#include "target_GDE.h"

#include <utility>
#include <functional>
#include <cfloat>
#include <fstream>

namespace
{
	static const size_t cycle_lock = 100000;
	static const double epsilon = 10e-6;

	double halfCutMethod(const std::function<double(double)>& func, const double left_boundary, const double right_boundary, const double epsilon, const size_t quantity_iteration)
	{
		double a = left_boundary;
		double b = right_boundary;
		double c;

		for (size_t i = 0; i < quantity_iteration; ++i)
		{
			double f_a = func(a);
			double f_b = func(b);

			c = (a + b) / 2;
			double f_c = func(c);

			a = (f_a * f_c) > 0 ? c : a;
			b = (f_b * f_c) > 0 ? c : b;

			if (abs(f_c) < epsilon)
			{
				break;
			}
		}
		return c;
	}

	void write(std::filesystem::path file_path, double* arr, size_t size)
	{
		std::ofstream fout(file_path, std::ios::trunc | std::ios::binary);
		fout.write((char*)arr, sizeof(double) * size);
	}
}

double TargetGDE::_gamma = 1.4;

TargetGDE::TargetGDE(const ExpanseGrid& expanse_grid, const BorderlineCondition& borderline_condition, const std::vector<double>& S, const double x_s) :
	_expanse_grid(expanse_grid),
	_borderline_condition(borderline_condition),
	_S(S),
	_data(expanse_grid.nodes),
	_x_s(x_s)
{
	initBottleneck();
}

TargetGDE::TargetGDE(const ExpanseGrid& expanse_grid, const BorderlineCondition& borderline_condition, const std::function<double(double)>& S_func, const double x_s) :
	_expanse_grid(expanse_grid),
	_borderline_condition(borderline_condition),
	_S(expanse_grid.nodes),
	_data(expanse_grid.nodes),
	_x_s(x_s)
{
	for (size_t i = 0; i < _expanse_grid.nodes; i++)
	{
		(*const_cast<std::vector<double>*>(&_S))[i] = S_func(i * _expanse_grid.h);
	}
	initBottleneck();
}


//	Моделирует УГД
void TargetGDE::solving()
{
	const double p0 = _borderline_condition.left_borderline_p(_expanse_grid.starting_point);
	const double ro0 = _borderline_condition.left_borderline_ro(_expanse_grid.starting_point);
	const double u0 = pow(_gamma * p0 / ro0, 0.5);

	const double S_crit = _S_crit;
	const double S_0 = _S[0];

	std::function<double(double)> func = [S_0, S_crit, M_crit = 1.0](double M)
		{
			const double Q = (1 + (_gamma - 1) / 2 * M_crit * M_crit) / (1 + (_gamma - 1) / 2 * M * M);
			return M_crit / M * pow(Q, (_gamma + 1) / (1 - _gamma) / 2) - S_0 / S_crit;
		};

	const double M_0 = halfCutMethod(func, 0.0001, 1, epsilon, cycle_lock);

	size_t index_crit = (_expanse_grid.starting_point + _x_crit) / _expanse_grid.h;
	size_t index_s = (_expanse_grid.starting_point + _x_s) / _expanse_grid.h + 1;
	index_s = index_s > _expanse_grid.nodes ? _expanse_grid.nodes : index_s;

	for (size_t i = 0; i < index_s; ++i)
	{
		const double S_i = _S[i];
		func = [S_i, S_0, M_0](double M)
			{
				const double Q = (1 + (_gamma - 1) / 2 * M_0 * M_0) / (1 + (_gamma - 1) / 2 * M * M);
				return M_0 / M * pow(Q, (_gamma + 1) / (1 - _gamma) / 2) - S_i / S_0;
			};

		double M = 1;
		if (i < index_crit)
		{
			M = halfCutMethod(func, 0.0001, 1.0, epsilon, cycle_lock);
		}
		else if (i > index_crit)
		{
			M = halfCutMethod(func, 1.0, 20, epsilon, cycle_lock);
		}

		const double Q = (1 + (_gamma - 1) / 2 * M_0 * M_0) / (1 + (_gamma - 1) / 2 * M * M);

		_data[i].ro = ro0 * pow(Q, 1 / (_gamma - 1));
		_data[i].p = p0 * pow(Q, _gamma / (_gamma - 1));
		_data[i].u = u0 * pow(Q, 0.5) * M / M_0;
	}

	if (index_s >= _expanse_grid.nodes) return;


	Values left = _data[index_s - 1];
	const double M_L = left.u / pow(_gamma * left.p / left.ro, 0.5);

	_data[index_s] = {
		left.ro * (_gamma + 1) / 2 * M_L * M_L / (1 + (_gamma - 1) / 2 * M_L * M_L),
		left.p * (2 * M_L * M_L * _gamma / (_gamma + 1) - (_gamma - 1) / (_gamma + 1)),
		left.u * (_gamma - 1) / (_gamma + 1) * (1 + 2 / (_gamma - 1) / (M_L * M_L))
	};
	Values right = _data[index_s];
	const double M_R = pow((1 + M_L * M_L * (_gamma - 1) / 2) / (_gamma * M_L * M_L - (_gamma - 1) / 2), 0.5);
	
	const double S_s = _S[index_s];


	for (size_t i = index_s + 1; i < _S.size(); ++i)
	{
		const double S_i = _S[i];
		func = [S_i, S_s, M_R](double M)
			{
				const double Q = (1 + (_gamma - 1) / 2 * M_R * M_R) / (1 + (_gamma - 1) / 2 * M * M);
				return M_R / M * pow(Q, (_gamma + 1) / (1 - _gamma) / 2) - S_i / S_s;
			};

		const double M = halfCutMethod(func, 0.0001, 1, epsilon, cycle_lock);
		const double Q = (1 + (_gamma - 1) / 2 * M_R * M_R) / (1 + (_gamma - 1) / 2 * M * M);

		_data[i].ro = right.ro * pow(Q, 1 / (_gamma - 1));
		_data[i].p = right.p * pow(Q, _gamma / (_gamma - 1));
		_data[i].u = right.u * pow(Q, 0.5) * M / M_R;
	}
	return;
}


//	Записывает массив данных плотности в файл
void TargetGDE::writeRO(const std::filesystem::path& file_path)  const
{
	write(file_path, getRO().data(), _expanse_grid.nodes);
}

//	Записывает массив данных скорости в файл
void TargetGDE::writeU(const std::filesystem::path& file_path)  const
{
	write(file_path, getU().data(), _expanse_grid.nodes);
}

//	Записывает массив данных давления в файл
void TargetGDE::writeP(const std::filesystem::path& file_path)  const
{
	write(file_path, getP().data(), _expanse_grid.nodes);
}


//	Вернуть массив данных плотности
std::vector<double> TargetGDE::getRO() const
{
	std::vector<double> res(_expanse_grid.nodes);
	for (size_t i = 0; i < _expanse_grid.nodes; i++)
	{
		res[i] = _data[i].ro;
	}
	return res;
}

//	Вернуть массив данных скорости
std::vector<double> TargetGDE::getU() const
{
	std::vector<double> res(_expanse_grid.nodes);
	for (size_t i = 0; i < _expanse_grid.nodes; i++)
	{
		res[i] = _data[i].u;
	}
	return res;
}

//	Вернуть массив данных давления
std::vector<double> TargetGDE::getP() const
{
	std::vector<double> res(_expanse_grid.nodes);
	for (size_t i = 0; i < _expanse_grid.nodes; i++)
	{
		res[i] = _data[i].p;
	}
	return res;
}


//	Вычисляет позицию и площадь наименьшего сечения трубы
void TargetGDE::initBottleneck()
{
	double A_min = DBL_MAX;
	size_t index_min = 0;

	for (size_t i = 0; i < _expanse_grid.nodes; ++i)
	{
		if (_S[i] < A_min)
		{
			index_min = i;
			A_min = _S[i];
		}
	}
	_x_crit = _expanse_grid.starting_point + static_cast<double>(index_min) * _expanse_grid.h;
	_S_crit = A_min;
}