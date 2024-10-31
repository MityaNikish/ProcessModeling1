#include "solver_VCSP_GDE.h"

#define NDEBUG
#include <cassert>
//#include <iostream>


namespace
{
	using enum Num;
}

VariableCrossSectionPipeGDE::VariableCrossSectionPipeGDE(const ExpanseGrid& expanse_grid, const TimeGrid& time_grid, StartCondition& start_condition, BorderlineCondition& borderline_condition, const std::vector<double>& S) :
	GasDynamicsEquation(expanse_grid, time_grid, start_condition, borderline_condition),
	_S(S)
{
	assert(_S.size() == expanse_grid.nodes && "The cross-section array is not of the right lenght!");
}

VariableCrossSectionPipeGDE::VariableCrossSectionPipeGDE(const ExpanseGrid& expanse_grid, const TimeGrid& time_grid, StartCondition& start_condition, BorderlineCondition& borderline_condition, const std::function<double(double)>& S_func) : 
	GasDynamicsEquation(expanse_grid, time_grid, start_condition, borderline_condition),
	_S(_expanse_grid.nodes)
{
	for (size_t i = 0; i < _expanse_grid.nodes; i++)
	{
		(*const_cast<std::vector<double>*>(&_S))[i] = S_func(i * _expanse_grid.h);
		//std::cout << _S[i] << "\n";
	}
}


//	Инициализирует сетки начальными и граничными условиями
void VariableCrossSectionPipeGDE::initConditions()
{
	for (size_t j = 0; j < _expanse_grid.nodes; ++j)
	{
		_ro.getElement(0, j) = _S[j] * _start_condition.start_ro(_expanse_grid.starting_point + _expanse_grid.h * static_cast<double>(j));
		_u.getElement(0, j) = _start_condition.start_u(_expanse_grid.starting_point + _expanse_grid.h * static_cast<double>(j));
		_p.getElement(0, j) = _S[j] * _start_condition.start_p(_expanse_grid.starting_point + _expanse_grid.h * static_cast<double>(j));
	}

	for (size_t n = 0; n < 2; ++n)
	{
		_ro.getElement(n, 0) = _S[0] * _borderline_condition.left_borderline_ro(_time_grid.starting_point + _time_grid.tau * static_cast<double>(n));
		_u.getElement(n, 0) = _borderline_condition.left_borderline_u(_time_grid.starting_point + _time_grid.tau * static_cast<double>(n));
		_p.getElement(n, 0) = _S[0] * _borderline_condition.left_borderline_p(_time_grid.starting_point + _time_grid.tau * static_cast<double>(n));

		_ro.getElement(n, _expanse_grid.nodes - 1) = _S[_expanse_grid.nodes - 1] * _borderline_condition.right_borderline_ro(_time_grid.starting_point + _time_grid.tau * static_cast<double>(n));
		_u.getElement(n, _expanse_grid.nodes - 1) = _borderline_condition.right_borderline_u(_time_grid.starting_point + _time_grid.tau * static_cast<double>(n));
		_p.getElement(n, _expanse_grid.nodes - 1) = _S[_expanse_grid.nodes - 1] * _borderline_condition.right_borderline_p(_time_grid.starting_point + _time_grid.tau * static_cast<double>(n));
	}
}

//	Проводит послевычеслительные операции
void VariableCrossSectionPipeGDE::postProcessing()
{
	for (size_t n = 0; n < 2; ++n)
	{
		for (size_t j = 0; j < _expanse_grid.nodes; ++j)
		{
			_ro.getElement(n, j) /= _S[j];
			_p.getElement(n, j) /= _S[j];
		}
	}
}

//	Проверка на удовлетворение условий остановки вычислений
bool VariableCrossSectionPipeGDE::chekStopConditions(size_t n, size_t j, double eps) const
{
	const double ro_ = _ro.getElement(n, j);
	const double ro_next = _ro.getElement(n, j + 1);

	const double u_ = _ro.getElement(n, j);
	const double u_next = _ro.getElement(n, j + 1);

	const double div = (ro_next * u_next - ro_ * u_) / _expanse_grid.h;

	return abs(div) * _S[j] / ro_ < eps;
}


//	Источник
Vector3D VariableCrossSectionPipeGDE::Q(size_t n, size_t j) const
{
	const double p = _p.getElement(n, j);
	double dS;
	if (j == 0)
	{
		dS = (_S[2] - _S[0]) / 2 / _expanse_grid.h;
	}
	else if (j == _expanse_grid.nodes - 1)
	{
		dS = (_S[_expanse_grid.nodes - 1] - _S[_expanse_grid.nodes - 3]) / 2 / _expanse_grid.h;
	}
	else
	{
		dS = (_S[j + 1] - _S[j - 1]) / 2 / _expanse_grid.h;
	}

	return Vector3D(0, p / _S[j] * dS, 0);
}