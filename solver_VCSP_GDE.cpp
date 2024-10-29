#include "solver_VCSP_GDE.h"

//#define NDEBUG
#include <cassert>
#include <iostream>


namespace
{
	using enum Num;
}

VariableCrossSectionPipeGDE::VariableCrossSectionPipeGDE(const ExpanseGrid& expanse_grid, const TimeGrid& time_grid, StartCondition& start_condition, BorderlineCondition& borderline_condition, const std::vector<double>& S)
	: GasDynamicsEquation(expanse_grid, time_grid, start_condition, borderline_condition), _S(S)
{
	assert(_S.size() == expanse_grid.nodes && "The cross-section array is not of the right lenght!");
}


//void VariableCrossSectionPipeGDE::solving()
//{
//	for (size_t n = 0; n < _time_grid.nodes - 1; ++n)
//	{
//		for (size_t j = 1; j < _expanse_grid.nodes - 1; ++j)
//		{
//			const double ro = _ro.getElement(n, j);
//			const double u = _u.getElement(n, j);
//			const double p = _p.getElement(n, j);
//
//			const Vector3D U_pref = U(n, j - 1);
//			const Vector3D U_ = U(n, j);
//			const Vector3D U_next = U(n, j + 1);
//
//
//			//const Matrix3D A_pref = _A_const;
//			//const Matrix3D A_ = _A_const;
//			//const Matrix3D A_next = _A_const;
//
//			//const Matrix3D A_half_pref = _A_const;
//			//const Matrix3D A_half_next = _A_const;
//
//			const Matrix3D A_pref = A(n, j - 1);
//			const Matrix3D A_ = A(n, j);
//			const Matrix3D A_next = A(n, j + 1);
//
//			const Matrix3D A_half_pref = (A_ + A_next) / 2;
//			const Matrix3D A_half_next = (A_pref + A_) / 2;
//
//
//			const Vector3D F_pref = F(n, j - 1);
//			const Vector3D F_ = F(n, j);
//			const Vector3D F_next = F(n, j + 1);
//
//			const double alpha = _time_grid.tau / _expanse_grid.h;
//
//			const Vector3D F_half_pref = (F_pref + F_ - A_half_pref * A_half_pref * (U_ - U_pref) * alpha) / 2;
//			const Vector3D F_half_next = (F_ + F_next - A_half_next * A_half_pref * (U_next - U_) * alpha) / 2;
//
//			//const Vector3D F_half_pref = (F_pref + F_ - A_half_pref * (F_ - F_pref) * alpha) / 2;
//			//const Vector3D F_half_next = (F_ + F_next - A_half_next * (F_next - F_) * alpha) / 2;
//
//			const Vector3D Q_ = Q(n, j);
//
//			//const Vector3D U_star = U_ - (F_half_next - F_half_pref) * _time_grid.tau / _expanse_grid.h + (Q_future - Q_ - (Q_next - Q_) * _time_grid.tau / _expanse_grid.h) / 2;
//			//const Vector3D U_star = U_ - (F_half_next - F_half_pref) * _time_grid.tau / _expanse_grid.h + (Q_next - Q_) * _time_grid.tau / _expanse_grid.h / 2;
//			const Vector3D U_star = U_ - (F_half_next - F_half_pref) * alpha - Q_ * _time_grid.tau;
//
//			const double ro_new = U_star[one];
//			const double u_new = U_star[two] / ro_new;
//			const double E_new = U_star[three] / ro_new;
//			const double p_new = (E_new - u_new * u_new / 2) * ro_new * (_gamma - 1);
//
//			_ro.getElement(n + 1, j) = ro_new;
//			_u.getElement(n + 1, j) = u_new;
//			_p.getElement(n + 1, j) = p_new;
//		}
//	}
//	postProcessing();
//}


//	Инициализирует сетки начальными и граничными условиями
void VariableCrossSectionPipeGDE::initConditions()
{
	for (size_t j = 0; j < _expanse_grid.nodes; ++j)
	{
		_ro.getElement(0, j) = _S[j] * _start_condition.start_ro(_expanse_grid.starting_point + _expanse_grid.h * static_cast<double>(j));
		_u.getElement(0, j) = _start_condition.start_u(_expanse_grid.starting_point + _expanse_grid.h * static_cast<double>(j));
		_p.getElement(0, j) = _S[j] * _start_condition.start_p(_expanse_grid.starting_point + _expanse_grid.h * static_cast<double>(j));
	}

	for (size_t n = 0; n < _time_grid.nodes; ++n)
	{
		_ro.getElement(n, 0) = _S[0] * _borderline_condition.start_borderline_ro(_time_grid.starting_point + _time_grid.tau * static_cast<double>(n));
		_u.getElement(n, 0) = _borderline_condition.start_borderline_u(_time_grid.starting_point + _time_grid.tau * static_cast<double>(n));
		_p.getElement(n, 0) = _S[0] * _borderline_condition.start_borderline_p(_time_grid.starting_point + _time_grid.tau * static_cast<double>(n));

		_ro.getElement(n, _expanse_grid.nodes - 1) = _S[_expanse_grid.nodes - 1] * _borderline_condition.end_borderline_ro(_time_grid.starting_point + _time_grid.tau * static_cast<double>(n));
		_u.getElement(n, _expanse_grid.nodes - 1) = _borderline_condition.end_borderline_u(_time_grid.starting_point + _time_grid.tau * static_cast<double>(n));
		_p.getElement(n, _expanse_grid.nodes - 1) = _S[_expanse_grid.nodes - 1] * _borderline_condition.end_borderline_p(_time_grid.starting_point + _time_grid.tau * static_cast<double>(n));
	}
}

//	Проводит послевычеслительные операции
void VariableCrossSectionPipeGDE::postProcessing()
{
	for (size_t n = 0; n < _time_grid.nodes; ++n)
	{
		for (size_t j = 0; j < _expanse_grid.nodes; ++j)
		{
			_ro.getElement(0, j) /= _S[j];
			_p.getElement(0, j) /= _S[j];
		}
	}
}


//	Проверка на удовлетворение условий остановки вычислений
bool VariableCrossSectionPipeGDE::chekStopConditions(size_t n, size_t j, double eps)
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
