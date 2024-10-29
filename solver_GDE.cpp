#include "solver_GDE.h"
#include <fstream>
#include <iostream>

namespace 
{
	using enum Num;

	void write(std::filesystem::path file_path, double* arr, size_t size)
	{
		std::ofstream fout(file_path, std::ios::trunc | std::ios::binary);
		fout.write((char*)arr, sizeof(double) * size);
	}
}


GasDynamicsEquation::GasDynamicsEquation(const ExpanseGrid& expanse_grid, const TimeGrid& time_grid, StartCondition& start_condition, BorderlineCondition& borderline_condition) :
	_expanse_grid(expanse_grid),
	_time_grid(time_grid),
	_start_condition(std::move(start_condition)),
	_borderline_condition(std::move(borderline_condition)),
	_ro(2, _expanse_grid.nodes),
	_u(2, _expanse_grid.nodes),
	_p(2, _expanse_grid.nodes),
	_alpha(_time_grid.tau / _expanse_grid.h)
{ }


//	Моделирует УГД
void GasDynamicsEquation::solving()
{
	initConditions();

	for (size_t n = 0; n < _time_grid.nodes - 1; ++n)
	{
		bool not_satisfy_condition = false;
		for (size_t j = 1; j < _expanse_grid.nodes - 1; ++j)
		{
			Vector3D U_star = U_future(n, j);
			//Vector3D U_star = dU(n, j);

			const double ro_new = U_star[one];
			const double u_new = U_star[two] / ro_new;
			const double E_new = U_star[three] / ro_new;
			const double p_new = (E_new - u_new * u_new / 2) * ro_new * (_gamma - 1);

			_ro.getElement(n + 1, j) = ro_new;
			_u.getElement(n + 1, j) = u_new;
			_p.getElement(n + 1, j) = p_new;

			if (!not_satisfy_condition && !chekStopConditions(n, j, 1.0e-6))
			{
				not_satisfy_condition = !not_satisfy_condition;
			}
		}
		if (!not_satisfy_condition)
		{
			std::cout << n << std::endl;
			break;
		}
	}

	postProcessing();
}


//	Запиывает данных плотности в файл
void GasDynamicsEquation::writeRO(const std::filesystem::path& file_path)  const
{
	write(file_path, _ro[0].getData(), _expanse_grid.nodes);
}

//	Запиывает данных скорости в файл
void GasDynamicsEquation::writeU(const std::filesystem::path& file_path)  const
{
	write(file_path, _u[0].getData(), _expanse_grid.nodes);
}

//	Запиывает данных давления в файл
void GasDynamicsEquation::writeP(const std::filesystem::path& file_path)  const
{
	write(file_path, _p[0].getData(), _expanse_grid.nodes);
}


//	Инициализирует сетки начальными и граничными условиями
void GasDynamicsEquation::initConditions()
{
	for (size_t j = 0; j < _expanse_grid.nodes; ++j)
	{
		_ro.getElement(0, j) = _start_condition.start_ro(_expanse_grid.starting_point + _expanse_grid.h * static_cast<double>(j));
		_u.getElement(0, j) = _start_condition.start_u(_expanse_grid.starting_point + _expanse_grid.h * static_cast<double>(j));
		_p.getElement(0, j) = _start_condition.start_p(_expanse_grid.starting_point + _expanse_grid.h * static_cast<double>(j));
	}

	for (size_t n = 0; n < 2; ++n)
	{
		_ro.getElement(n, 0) = _borderline_condition.start_borderline_ro(_time_grid.starting_point + _time_grid.tau * static_cast<double>(n));
		_u.getElement(n, 0) = _borderline_condition.start_borderline_u(_time_grid.starting_point + _time_grid.tau * static_cast<double>(n));
		_p.getElement(n, 0) = _borderline_condition.start_borderline_p(_time_grid.starting_point + _time_grid.tau * static_cast<double>(n));

		_ro.getElement(n, _expanse_grid.nodes - 1) = _borderline_condition.end_borderline_ro(_time_grid.starting_point + _time_grid.tau * static_cast<double>(n));
		_u.getElement(n, _expanse_grid.nodes - 1) = _borderline_condition.end_borderline_u(_time_grid.starting_point + _time_grid.tau * static_cast<double>(n));
		_p.getElement(n, _expanse_grid.nodes - 1) = _borderline_condition.end_borderline_p(_time_grid.starting_point + _time_grid.tau * static_cast<double>(n));
	}
}

//	Проводит послевычеслительные операции
void GasDynamicsEquation::postProcessing()
{
	return;
}


//	Проверка на удовлетворение условий остановки вычислений
bool GasDynamicsEquation::chekStopConditions(size_t n, size_t j, double eps)
{
	const double ro_ = _ro.getElement(n, j);
	const double ro_next = _ro.getElement(n, j + 1);

	const double u_ = _ro.getElement(n, j);
	const double u_next = _ro.getElement(n, j + 1);

	const double div = (ro_next * u_next - ro_ * u_) / _expanse_grid.h;

	return abs(div) / ro_ < eps;
}


//	Удельная полная энергия
double GasDynamicsEquation::E(size_t n, size_t j) const
{
	const double ro = _ro.getElement(n, j);
	const double u = _u.getElement(n, j);
	const double p = _p.getElement(n, j);

	return p / ro / (_gamma - 1) + u * u / 2;
}

//	Удельная польная энтальпия
double GasDynamicsEquation::H(size_t n, size_t j) const
{
	const double ro = _ro.getElement(n, j);
	const double p = _p.getElement(n, j);

	return E(n, j) + p / ro;
}


//	Консервативная переменная
Vector3D GasDynamicsEquation::U(size_t n, size_t j) const
{
	const double ro = _ro.getElement(n, j);
	const double u = _u.getElement(n, j);

	return Vector3D (ro, ro * u, ro * E(n, j));
}

//	beta: Борьба с осциляциями при помощи искусственной вязкости
Vector3D GasDynamicsEquation::dU(size_t n, size_t j) const
{
	Vector3D U_future_pref = U_future(n, j - 1);
	Vector3D U_future_ = U_future(n, j);
	Vector3D U_future_next = U_future(n, j + 1);
	
	Vector3D dU_pref = U_future_ - U_future_pref;
	Vector3D dU_next = U_future_next - U_future_;

	Vector3D dU_half_pref = dU_pref * pow(dU_pref * dU_pref, 0.5);
	Vector3D dU_half_next = dU_next * pow(dU_next * dU_next, 0.5);

	return U_future_ + (dU_half_next - dU_half_pref) * (_artificial_viscosity * _alpha);
}

//	Проводит вычисление консервативной переменной, следующей по временной сетке
Vector3D GasDynamicsEquation::U_future(size_t n, size_t j) const
{
	const double ro = _ro.getElement(n, j);
	const double u = _u.getElement(n, j);
	const double p = _p.getElement(n, j);

	const Vector3D U_pref = U(n, j - 1);
	const Vector3D U_ = U(n, j);
	const Vector3D U_next = U(n, j + 1);

	const Matrix3D A_pref = A(n, j - 1);
	const Matrix3D A_ = A(n, j);
	const Matrix3D A_next = A(n, j + 1);

	const Matrix3D A_half_pref = (A_ + A_next) / 2;
	const Matrix3D A_half_next = (A_pref + A_) / 2;

	const Vector3D F_pref = F(n, j - 1);
	const Vector3D F_ = F(n, j);
	const Vector3D F_next = F(n, j + 1);

	//const Vector3D F_half_pref = (F_pref + F_ - A_half_pref * (F_ - F_pref) * _alpha) / 2;
	//const Vector3D F_half_next = (F_ + F_next - A_half_next * (F_next - F_) * _alpha) / 2;

	const Vector3D F_half_pref = (F_pref + F_ - A_half_pref * A_half_pref * (U_ - U_pref) * _alpha) / 2;
	const Vector3D F_half_next = (F_ + F_next - A_half_next * A_half_next * (U_next - U_) * _alpha) / 2;

	return U_ - (F_half_next - F_half_pref) * _alpha - Q(n, j) * _time_grid.tau;
}

//	Значение потока
Vector3D GasDynamicsEquation::F(size_t n, size_t j) const
{
	const double ro = _ro.getElement(n, j);
	const double u = _u.getElement(n, j);
	const double p = _p.getElement(n, j);

	double ro_u = ro * u;

	return Vector3D(ro_u, ro_u * u + p, ro_u * E(n, j) + u * p);
}

//	Источник
Vector3D GasDynamicsEquation::Q(size_t n, size_t j) const
{
	return Vector3D();
}


//	Матрица Якоби (A = dF / dU)
Matrix3D GasDynamicsEquation::A(size_t n, size_t j) const
{
	const double u = _u.getElement(n, j);

	double H_ = H(n, j);

	double A_1_0 = (_gamma - 3) * u * u / 2;
	double A_1_1 = (3 - _gamma) * u;
	double A_1_2 = _gamma - 1;

	double A_2_0 = u * ((_gamma - 1) * u * u / 2 - H_);
	double A_2_1 = H_ - (_gamma - 1) * u * u;
	double A_2_2 = _gamma * u;

	return Matrix3D( Vector3D{ 0, 1, 0 }, Vector3D{ A_1_0, A_1_1, A_1_2 }, Vector3D{ A_2_0, A_2_1, A_2_2 } );
}