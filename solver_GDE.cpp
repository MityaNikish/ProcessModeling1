#include "solver_GDE.h"
#include <fstream>
#include <iostream>

#define NDEBUG
#include <cassert>


namespace 
{
	using enum Num;

	void write(std::filesystem::path file_path, double* arr, size_t size)
	{
		std::ofstream fout(file_path, std::ios::trunc | std::ios::binary);
		fout.write((char*)arr, sizeof(double) * size);
	}
}


//	Показатель адиабаты
double GasDynamicsEquation::gamma = 1.4;	//	1.4 - Воздух

//	Искусственная вязкость Лапидуса
double GasDynamicsEquation::artificial_viscosity = 2.0;	//	2.0 - по умолчанию


GasDynamicsEquation::GasDynamicsEquation(const ExpanseGrid& expanse_grid, const TimeGrid& time_grid, StartCondition& start_condition, BorderlineCondition& borderline_condition, const std::vector<double>& S) :
	_expanse_grid(expanse_grid),
	_time_grid(time_grid),
	_start_condition(std::move(start_condition)),
	_borderline_condition(std::move(borderline_condition)),
	_ro(2, _expanse_grid.nodes),
	_u(2, _expanse_grid.nodes),
	_p(2, _expanse_grid.nodes),
	_S(S),
	_alpha(_time_grid.tau / _expanse_grid.h)
{ 
	assert(_S.size() == expanse_grid.nodes && "The cross-section array is not of the right lenght!");
}


GasDynamicsEquation::GasDynamicsEquation(const ExpanseGrid& expanse_grid, const TimeGrid& time_grid, StartCondition& start_condition, BorderlineCondition& borderline_condition, const std::function<double(double)>& S_func) :
	_expanse_grid(expanse_grid),
	_time_grid(time_grid),
	_start_condition(std::move(start_condition)),
	_borderline_condition(std::move(borderline_condition)),
	_ro(2, _expanse_grid.nodes),
	_u(2, _expanse_grid.nodes),
	_p(2, _expanse_grid.nodes),
	_S(_expanse_grid.nodes),
	_alpha(_time_grid.tau / _expanse_grid.h)
{
	for (size_t i = 0; i < _expanse_grid.nodes; i++)
	{
		(*const_cast<std::vector<double>*>(&_S))[i] = S_func(i * _expanse_grid.h);
		//std::cout << _S[i] << "\n";
	}
}


//	Моделирует УГД
void GasDynamicsEquation::solving()
{
	initConditions();

	for (size_t n = 0; n < _time_grid.nodes - 1; ++n)
	{
		bool not_satisfy_condition = false;
		for (size_t j = 1; j < _expanse_grid.nodes - 1; ++j)
		{
			Vector3D U_star = dU(n, j);

			const double ro_new = U_star[one];
			const double u_new = U_star[two] / ro_new;
			const double E_new = U_star[three] / ro_new;
			const double p_new = (E_new - u_new * u_new / 2) * ro_new * (gamma - 1);

			_ro.getElement(n + 1, j) = ro_new;
			_u.getElement(n + 1, j) = u_new;
			_p.getElement(n + 1, j) = p_new;

			if (!not_satisfy_condition && !chekStopConditions(n, j, 1.0e-1))
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


//	Записывает массив данных плотности в файл
void GasDynamicsEquation::writeRO(const std::filesystem::path& file_path)  const
{
	write(file_path, _ro[0].getData(), _expanse_grid.nodes);
}

//	Записывает массив данных скорости в файл
void GasDynamicsEquation::writeU(const std::filesystem::path& file_path)  const
{
	write(file_path, _u[0].getData(), _expanse_grid.nodes);
}

//	Записывает массив данных давления в файл
void GasDynamicsEquation::writeP(const std::filesystem::path& file_path)  const
{
	write(file_path, _p[0].getData(), _expanse_grid.nodes);
}


//	Вернуть массив данных плотности
std::vector<double> GasDynamicsEquation::getRO() const
{
	Matrix slice = _ro[0];
	return std::vector<double>(slice.getData(), slice.getData() + _expanse_grid.nodes);
}

//	Вернуть массив данных скорости
std::vector<double> GasDynamicsEquation::getU() const
{
	Matrix slice = _u[0];
	return std::vector<double>(slice.getData(), slice.getData() + _expanse_grid.nodes);
}

//	Вернуть массив данных давления
std::vector<double> GasDynamicsEquation::getP() const
{
	Matrix slice = _p[0];
	return std::vector<double>(slice.getData(), slice.getData() + _expanse_grid.nodes);
}


//	Инициализирует сетки начальными и граничными условиями
void GasDynamicsEquation::initConditions()
{
	for (size_t j = 0; j < _expanse_grid.nodes; ++j)
	{
		_ro.getElement(static_cast<size_t>(0), j) = _S[j] * _start_condition.start_ro(_expanse_grid.starting_point + _expanse_grid.h * static_cast<double>(j));
		_u.getElement(static_cast<size_t>(0), j) = _start_condition.start_u(_expanse_grid.starting_point + _expanse_grid.h * static_cast<double>(j));
		_p.getElement(static_cast<size_t>(0), j) = _S[j] * _start_condition.start_p(_expanse_grid.starting_point + _expanse_grid.h * static_cast<double>(j));
	}

	for (size_t n = 0; n < 2; ++n)
	{
		_ro.getElement(n, static_cast<size_t>(0)) = _S[0] * _borderline_condition.left_borderline_ro;
		_u.getElement(n, static_cast<size_t>(0)) = _borderline_condition.left_borderline_u;
		_p.getElement(n, static_cast<size_t>(0)) = _S[0] * _borderline_condition.left_borderline_p;

		_ro.getElement(n, _expanse_grid.nodes - 1) = _S[_expanse_grid.nodes - 1] * _borderline_condition.right_borderline_ro;
		_u.getElement(n, _expanse_grid.nodes - 1) = _borderline_condition.right_borderline_u;
		_p.getElement(n, _expanse_grid.nodes - 1) = _S[_expanse_grid.nodes - 1] * _borderline_condition.right_borderline_p;
	}
}

//	Проводит послевычеслительные операции
void GasDynamicsEquation::postProcessing()
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
bool GasDynamicsEquation::chekStopConditions(size_t n, size_t j, double eps) const
{
	/*
	//	Условие остановки по t
	const double ro_ = _ro.getElement(n, j);
	const double ro_future = _ro.getElement(n + 1, j);

	const double div = (ro_future - ro_) / _time_grid.tau;

	return abs(div) * _S[j] / ro_ < eps;
	*/

	//	Условие остановки по x
	const double ro_ = _ro.getElement(n, j);
	const double ro_next = _ro.getElement(n, j + 1);

	const double u_ = _ro.getElement(n, j);
	const double u_next = _ro.getElement(n, j + 1);

	const double div = (ro_next * u_next - ro_ * u_) / _expanse_grid.h;

	return abs(div) * _S[j] / ro_ < eps;
	
	//return false;
}


//	Удельная полная энергия
double GasDynamicsEquation::E(size_t n, size_t j) const
{
	const double ro = _ro.getElement(n, j);
	const double u = _u.getElement(n, j);
	const double p = _p.getElement(n, j);

	return p / ro / (gamma - 1) + u * u / 2;
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

//	Борьба с осциляциями при помощи искусственной вязкости
Vector3D GasDynamicsEquation::dU(size_t n, size_t j) const
{
	//	Костыль
	Vector3D U_future_pref = j == 1 ? U_future(n, 1) : U_future(n, j - 1);
	Vector3D U_future_ = U_future(n, j);
	Vector3D U_future_next = j == _expanse_grid.nodes - 2 ? U_future(n, _expanse_grid.nodes - 2) : U_future(n, j + 1);

	/*
	Vector3D U_future_pref = U_future(n, j - 1);
	Vector3D U_future_ = U_future(n, j);
	Vector3D U_future_next = j == U_future(n, j + 1);
	*/
	
	Vector3D dU_pref = U_future_ - U_future_pref;
	Vector3D dU_next = U_future_next - U_future_;

	Vector3D dU_half_pref = dU_pref * dU_pref.normL2();
	Vector3D dU_half_next = dU_next * dU_next.normL2();

	return U_future_ + (dU_half_next - dU_half_pref) * (artificial_viscosity * _alpha);
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
	const Matrix3D _A = A(n, j);
	const Matrix3D A_next = A(n, j + 1);

	const Matrix3D A_half_pref = (_A + A_next) / 2;
	const Matrix3D A_half_next = (A_pref + _A) / 2;

	const Vector3D F_pref = F(n, j - 1);
	const Vector3D F_ = F(n, j);
	const Vector3D F_next = F(n, j + 1);

	const Vector3D F_half_pref = (F_pref + F_ - A_half_pref * (F_ - F_pref) * _alpha) / 2;
	const Vector3D F_half_next = (F_ + F_next - A_half_next * (F_next - F_) * _alpha) / 2;

	/*
	//	Альтернативная формула
	const Vector3D F_half_pref = (F_pref + F_ - A_half_pref * A_half_pref * (U_ - U_pref) * _alpha) / 2;
	const Vector3D F_half_next = (F_ + F_next - A_half_next * A_half_next * (U_next - U_) * _alpha) / 2;
	*/

	return U_ - (F_half_next - F_half_pref) * _alpha + Q(n, j) * _time_grid.tau;
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
	const double p = _p.getElement(n, j);
	double dS;
	if (j == 0)
	{
		dS = (_S[1] - _S[0]) / _expanse_grid.h;
	}
	else if (j == _expanse_grid.nodes - 1)
	{
		dS = (_S[_expanse_grid.nodes - 1] - _S[_expanse_grid.nodes - 2]) / _expanse_grid.h;
	}
	else
	{
		dS = (_S[j + 1] - _S[j - 1]) / 2 / _expanse_grid.h;
	}

	return Vector3D(0, p / _S[j] * dS, 0);
}


//	Матрица Якоби (A = dF / dU)
Matrix3D GasDynamicsEquation::A(size_t n, size_t j) const
{
	
	////	Константная матрица А
	//static Matrix3D A;
	//static bool create = false;

	//if (!create)
	//{
	//	const double u = _u.getElement(0, 0);

	//	double H_ = H(0, 0);

	//	double A_1_0 = (gamma - 3) * u * u / 2;
	//	double A_1_1 = (3 - gamma) * u;
	//	double A_1_2 = gamma - 1;

	//	double A_2_0 = u * ((gamma - 1) * u * u / 2 - H_);
	//	double A_2_1 = H_ - (gamma - 1) * u * u;
	//	double A_2_2 = gamma * u;

	//	create = !create;
	//	A = Matrix3D(Vector3D{ 0, 1, 0 }, Vector3D{ A_1_0, A_1_1, A_1_2 }, Vector3D{ A_2_0, A_2_1, A_2_2 });
	//	//std::cout << A.getMaxElem() << std::endl;
	//}
	//return A;
	//

	const double u = _u.getElement(n, j);

	double H_ = H(n, j);

	double A_1_0 = (gamma - 3) * u * u / 2;
	double A_1_1 = (3 - gamma) * u;
	double A_1_2 = gamma - 1;

	double A_2_0 = u * ((gamma - 1) * u * u / 2 - H_);
	double A_2_1 = H_ - (gamma - 1) * u * u;
	double A_2_2 = gamma * u;

	return Matrix3D( Vector3D{ 0, 1, 0 }, Vector3D{ A_1_0, A_1_1, A_1_2 }, Vector3D{ A_2_0, A_2_1, A_2_2 } );
}