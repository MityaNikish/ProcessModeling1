#include "test.h"
#include "solver_GDE.h"
#include "target_GDE.h"

#include <iostream>
#include <fstream>


namespace
{
	//	Инициализация путей к рабочим файлам
	static const std::filesystem::path file_path_ro = std::filesystem::current_path() / "Treatment";
	static const std::filesystem::path file_path_u = std::filesystem::current_path() / "Treatment";
	static const std::filesystem::path file_path_p = std::filesystem::current_path() / "Treatment";
	static const std::filesystem::path file_path_grid_info = std::filesystem::current_path() / "Treatment\\grid_info.txt";

	void write_grid_info(const std::filesystem::path& file_path, const ExpanseGrid& expanse_grid, const TimeGrid& time_grid)
	{
		std::ofstream fout(file_path);
		fout << expanse_grid.h << "\n";
		fout << expanse_grid.nodes << "\n";
		fout << time_grid.tau << "\n";
		fout << time_grid.nodes << "\n";
	}

	void read_grid_info(const std::filesystem::path& file_path, ExpanseGrid& expanse_grid, TimeGrid& time_grid)
	{
		std::ifstream fin(file_path);
		fin >> expanse_grid.h;
		fin >> expanse_grid.nodes;
		fin >> time_grid.tau;
		fin >> time_grid.nodes;
	}

	std::vector<double> operator-(const std::vector<double>& left, const std::vector<double>& right)
	{
		if (left.size() != right.size())
		{
			exit(1);
		}


		size_t size = left.size();
		std::vector<double> res(size);

		for (size_t i = 0; i < size; ++i)
		{
			res[i] = left[i] - right[i];
		}
		return res;
	}

	std::vector<double> operator-(const std::vector<double>&& left, const std::vector<double>&& right)
	{
		if (left.size() != right.size())
		{
			exit(1);
		}


		size_t size = left.size();
		std::vector<double> res(size);

		for (size_t i = 0; i < size; ++i)
		{
			res[i] = left[i] - right[i];
		}
		return res;
	}

	std::vector<double>& abs_vec(std::vector<double>& this_)
	{
		size_t size = this_.size();
		for (size_t i = 0; i < size; ++i)
		{
			this_[i] = abs(this_[i]);
		}
		return this_;
	}

	std::vector<double>& abs_vec(std::vector<double>&& this_)
	{
		size_t size = this_.size();
		for (size_t i = 0; i < size; ++i)
		{
			this_[i] = abs(this_[i]);
		}
		return this_;
	}


	void write_arr(std::filesystem::path file_path, double* arr, size_t size)
	{
		std::ofstream fout(file_path, std::ios::trunc | std::ios::binary);
		fout.write((char*)arr, sizeof(double) * size);
	}

	void read_arr(std::filesystem::path file_path, double* arr, size_t size)
	{
		std::ifstream fin(file_path, std::ios::binary);
		fin.read((char*)arr, sizeof(double) * size);
	}
}


//	Моделирование целевого решения (для разрывного квазиодномерного течения)
void modeling_target()
{
	//	Инициализация сетки
	ExpanseGrid expanse_grid{ 0.01, 101, 0.0, 1.01 };
	TimeGrid time_grid{ 0.0001, 10000, 0.0, 10.0 };
	////	Частные параметры сетки
	//read_grid_info(file_path_grid_info, expanse_grid, time_grid);
	//expanse_grid.starting_point = 0.0;
	//expanse_grid.ending_point = expanse_grid.starting_point + expanse_grid.h * (expanse_grid.nodes - 1);
	////

	//	Функция попреречного сечения трубы
	const double pos_bottleneck = 0.5;
	const double value_bottleneck = 0.5;
	std::function<double(double)> S_func = [pos_bottleneck, value_bottleneck](double arg) {
		const double part = (1 - arg / pos_bottleneck);
		return value_bottleneck + (1 - value_bottleneck) * part * part;
		};

	//	Граничные условия
	BorderlineCondition borderline_condition;
	borderline_condition.left_borderline_ro = 1.0;
	borderline_condition.left_borderline_u = 1.0237498;
	borderline_condition.left_borderline_p = 8.0;

	//	Моделирование целевого решения
	TargetGDE solver(expanse_grid, borderline_condition, S_func, 0.75);
	solver.solving();

	//	Запись результатов в файл 
	solver.writeRO(file_path_ro / "output_ro.raw");
	solver.writeU(file_path_u / "output_u.raw");
	solver.writeP(file_path_p / "output_p.raw");
}

//	Моделирование ударной трубы
void modeling1()
{
	//	Инициализация сетки
	//		h = 0.01; n_x = 1001; x = [-4; 6]; tau = 0.0001; n_t = 25000; t = 2.5; gamma = 1.4; artificial_viscosity = 2.5
	GasDynamicsEquation::gamma = 1.4;
	GasDynamicsEquation::artificial_viscosity = 2.5;
	ExpanseGrid expanse_grid{ 0.01, 1001, -4.0, 6.01 };
	TimeGrid time_grid{ 0.0001, 25000, 0.0, 2.5 };

	read_grid_info(file_path_grid_info, expanse_grid, time_grid);
	expanse_grid.starting_point = -4.0;
	expanse_grid.ending_point = expanse_grid.starting_point + expanse_grid.h * (expanse_grid.nodes - 1);
	time_grid.starting_point = 0.0;
	time_grid.ending_point = time_grid.starting_point + time_grid.tau * (time_grid.nodes - 1);

	//	Функция попреречного сечения трубы
	std::function<double(double)> S_func = [](double arg) { return 1; };

	//	Начальный условия
	StartCondition start_condition;
	start_condition.start_ro = std::function<double(double)>([](double arg) { return arg < 0 ? 1.0 : 0.125; });
	start_condition.start_u = std::function<double(double)>([](double arg) { return 0.0; });
	start_condition.start_p = std::function<double(double)>([](double arg) { return arg < 0 ? 1.0 : 0.1; });

	//	Граничные условия
	BorderlineCondition borderline_condition;
	borderline_condition.left_borderline_ro = 1.0;
	borderline_condition.left_borderline_u = 0.0;
	borderline_condition.left_borderline_p = 1.0;
	borderline_condition.right_borderline_ro = 0.125;
	borderline_condition.right_borderline_u = 0.0;
	borderline_condition.right_borderline_p = 0.1;

	//	Моделирование ударной трубы
	GasDynamicsEquation solver(expanse_grid, time_grid, start_condition, borderline_condition, S_func);
	solver.solving();

	//	Запись результатов в файл 
	solver.writeRO(file_path_ro / "output_ro.raw");
	solver.writeU(file_path_u / "output_u.raw");
	solver.writeP(file_path_p / "output_p.raw");
}

//	Моделирование квазиодномерного течения в канале
void modeling2()
{
	//	Инициализация сетки
	//		h = 0.01; n_x = 101; x = [0; 1]; tau = 0.0001; n_t = 100000; t = 10.0; gamma = 1.4; artificial_viscosity = 2.5
	GasDynamicsEquation::gamma = 1.4;
	GasDynamicsEquation::artificial_viscosity = 2.5;
	ExpanseGrid expanse_grid{ 0.01, 101, 0.0, 1.01 };
	TimeGrid time_grid{ 0.0001, 100000, 0.0, 10.0 };

	read_grid_info(file_path_grid_info, expanse_grid, time_grid);
	expanse_grid.starting_point = 0.0;
	expanse_grid.ending_point = expanse_grid.starting_point + expanse_grid.h * (expanse_grid.nodes - 1);
	time_grid.starting_point = 0.0;
	time_grid.ending_point = time_grid.starting_point + time_grid.tau * (time_grid.nodes - 1);


	//	Функция попреречного сечения трубы
	const double pos_bottleneck = 0.5;
	const double value_bottleneck = 0.5;
	std::function<double(double)> S_func = [pos_bottleneck, value_bottleneck](double arg) {
		const double part = (1 - arg / pos_bottleneck);
		return value_bottleneck + (1 - value_bottleneck) * part * part;
		};

	//	Начальный условия
	StartCondition start_condition;
	start_condition.start_ro = std::function<double(double)>([expanse_grid](double arg) { return (0.1933880 - 1.0) / (expanse_grid.ending_point - expanse_grid.starting_point) * arg + 1.0; });
	start_condition.start_u = std::function<double(double)>([expanse_grid](double arg) { return (5.2937598 - 1.0237498) / (expanse_grid.ending_point - expanse_grid.starting_point) * arg + 1.0237498; });
	start_condition.start_p = std::function<double(double)>([expanse_grid](double arg) { return (0.8018469 - 8.0) / (expanse_grid.ending_point - expanse_grid.starting_point) * arg + 8.0; });

	//	Граничные условия
	BorderlineCondition borderline_condition;
	borderline_condition.left_borderline_ro = 1.0;
	borderline_condition.left_borderline_u = 1.0237498;
	borderline_condition.left_borderline_p = 8.0;
	borderline_condition.right_borderline_ro = 0.1933880;
	borderline_condition.right_borderline_u = 5.2937598;
	borderline_condition.right_borderline_p = 0.8018469;

	//	Моделирование квазиодномерного течения в канале
	GasDynamicsEquation solver(expanse_grid, time_grid, start_condition, borderline_condition, S_func);
	solver.solving();

	//	Запись результатов в файл 
	solver.writeRO(file_path_ro / "output_ro.raw");
	solver.writeU(file_path_u / "output_u.raw");
	solver.writeP(file_path_p / "output_p.raw");
}

//	Моделирование разрывного квазиодномерного течения в канале
void modeling3()
{
	//	Инициализация сетки
	//		h = 0.01; n_x = 101; x = [0; 1]; tau = 0.0001; n_t = 100000; t = 10.0; gamma = 1.4; artificial_viscosity = 2.5
	GasDynamicsEquation::gamma = 1.4;
	GasDynamicsEquation::artificial_viscosity = 1.0;
	ExpanseGrid expanse_grid{ 0.01, 101, 0.0, 1.01 };
	TimeGrid time_grid{ 0.0001, 100000, 0.0, 10.0 };

	read_grid_info(file_path_grid_info, expanse_grid, time_grid);
	expanse_grid.starting_point = 0.0;
	expanse_grid.ending_point = expanse_grid.starting_point + expanse_grid.h * (expanse_grid.nodes - 1);
	time_grid.starting_point = 0.0;
	time_grid.ending_point = time_grid.starting_point + time_grid.tau * (time_grid.nodes - 1);

	//	Функция попреречного сечения трубы
	const double pos_bottleneck = 0.5;
	const double value_bottleneck = 0.5;
	std::function<double(double)> S_func = [pos_bottleneck, value_bottleneck](double arg) {
		const double part = (1 - arg / pos_bottleneck);
		return value_bottleneck + (1 - value_bottleneck) * part * part;
		};

	//	Начальный условия
	StartCondition start_condition;
	start_condition.start_ro = std::function<double(double)>([expanse_grid](double arg) { return (0.8835893 - 1.0) / (expanse_grid.ending_point - expanse_grid.starting_point) * arg + 1.0; });
	start_condition.start_u = std::function<double(double)>([expanse_grid](double arg) { return  (1.1586268 - 1.0237498) / (expanse_grid.ending_point - expanse_grid.starting_point) * arg + 1.0237498; });
	start_condition.start_p = std::function<double(double)>([expanse_grid](double arg) { return  (7.0315580 - 8.0) / (expanse_grid.ending_point - expanse_grid.starting_point) * arg + 8.0; });

	//	Граничные условия
	BorderlineCondition borderline_condition;
	borderline_condition.left_borderline_ro = 1.0;
	borderline_condition.left_borderline_u = 1.0237498;
	borderline_condition.left_borderline_p = 8.0;
	borderline_condition.right_borderline_ro = 0.8835893;
	borderline_condition.right_borderline_u = 1.1586268;
	borderline_condition.right_borderline_p = 7.0315580;

	//	Моделирование квазиодномерного течения в канале
	GasDynamicsEquation solver(expanse_grid, time_grid, start_condition, borderline_condition, S_func);
	solver.solving();

	//	Запись результатов в файл 
	solver.writeRO(file_path_ro / "output_ro.raw");
	solver.writeU(file_path_u / "output_u.raw");
	solver.writeP(file_path_p / "output_p.raw");
}

//	Вычисление невязки для решения не разрывного квазиодномерного течения в канале
void discrepancy_continuous_quasi1D_flow()
{
	//	Инициализация сетки
	//	h = 0.01; n_x = 101; x = [0; 1]; tau = 0.0001; n_t = 100000; t = 10.0; gamma = 1.4; artificial_viscosity = 2.5
	GasDynamicsEquation::gamma = 1.4;
	GasDynamicsEquation::artificial_viscosity = 2.5;
	ExpanseGrid expanse_grid{ 0.001, 1001, 0.0, 0.0 };
	TimeGrid time_grid{ 0.00001, 100000, 0.0, 0.0 };

	//	Функция попреречного сечения трубы
	const double pos_bottleneck = 0.5;
	const double value_bottleneck = 0.5;
	std::function<double(double)> S_func = [pos_bottleneck, value_bottleneck](double arg) {
		const double part = (1 - arg / pos_bottleneck);
		return value_bottleneck + (1 - value_bottleneck) * part * part;
		};

	//	Граничные условия
	BorderlineCondition borderline_condition;
	borderline_condition.left_borderline_ro = 1.0;
	borderline_condition.left_borderline_u = 1.0237498;
	borderline_condition.left_borderline_p = 8.0;
	borderline_condition.right_borderline_ro = 0.1933880;
	borderline_condition.right_borderline_u = 5.2937598;
	borderline_condition.right_borderline_p = 0.8018469;

	//	Начальный условия
	StartCondition start_condition;
	start_condition.start_ro = std::function<double(double)>([expanse_grid, borderline_condition](double arg) {
		return (borderline_condition.right_borderline_ro - borderline_condition.left_borderline_ro) / (expanse_grid.h - expanse_grid.nodes) * arg + borderline_condition.left_borderline_ro;
		});
	start_condition.start_u = std::function<double(double)>([expanse_grid, borderline_condition](double arg) {
		return (borderline_condition.right_borderline_u - borderline_condition.left_borderline_u) / (expanse_grid.h - expanse_grid.nodes) * arg + borderline_condition.left_borderline_u;
		});
	start_condition.start_p = std::function<double(double)>([expanse_grid, borderline_condition](double arg) {
		return (borderline_condition.right_borderline_p - borderline_condition.left_borderline_p) / (expanse_grid.h - expanse_grid.nodes) * arg + borderline_condition.left_borderline_p;
		});


	//	Моделирование не разрывного квазиодномерного течения в канале
	GasDynamicsEquation GDE(expanse_grid, time_grid, start_condition, borderline_condition, S_func);
	GDE.solving();

	GDE.writeRO(file_path_ro / "ro.raw");
	GDE.writeU(file_path_u / "u.raw");
	GDE.writeP(file_path_p / "p.raw");

	//	Моделирование целевого решения не разрывного квазиодномерного течения в канале
	TargetGDE GDE_target(expanse_grid, borderline_condition, S_func);
	GDE_target.solving();

	GDE_target.writeRO(file_path_ro / "t_ro.raw");
	GDE_target.writeU(file_path_u / "t_u.raw");
	GDE_target.writeP(file_path_p / "t_p.raw");

	//	Вычисление невязки
	std::vector<double> discrepancy_ro = abs_vec(GDE.getRO() - GDE_target.getRO());
	std::vector<double> discrepancy_u = abs_vec(GDE.getU() - GDE_target.getU());
	std::vector<double> discrepancy_p = abs_vec(GDE.getP() - GDE_target.getP());

	write_arr(file_path_ro / "discrep_ro.raw", discrepancy_ro.data(), discrepancy_ro.size());
	write_arr(file_path_u / "discrep_u.raw", discrepancy_u.data(), discrepancy_u.size());
	write_arr(file_path_p / "discrep_p.raw", discrepancy_p.data(), discrepancy_p.size());
}

//	Вычисление невязки для решения разрывного квазиодномерного течения в канале
void discrepancy_discontinuous_quasi1D_flow()
{
	//	Инициализация сетки
	//	h = 0.01; n_x = 101; x = [0; 1]; tau = 0.0001; n_t = 100000; t = 10.0; gamma = 1.4; artificial_viscosity = 2.5
	GasDynamicsEquation::gamma = 1.4;
	GasDynamicsEquation::artificial_viscosity = 2.5;
	ExpanseGrid expanse_grid{ 0.001, 1001, 0.0, 0.0 };
	TimeGrid time_grid{ 0.00001, 100000, 0.0, 0.0 };

	//	Функция попреречного сечения трубы
	const double pos_bottleneck = 0.5;
	const double value_bottleneck = 0.5;
	std::function<double(double)> S_func = [pos_bottleneck, value_bottleneck](double arg) {
		const double part = (1 - arg / pos_bottleneck);
		return value_bottleneck + (1 - value_bottleneck) * part * part;
		};

	//	Граничные условия
	BorderlineCondition borderline_condition;
	borderline_condition.left_borderline_ro = 1.0;
	borderline_condition.left_borderline_u = 1.0237498;
	borderline_condition.left_borderline_p = 8.0;
	borderline_condition.right_borderline_ro = 0.8835893;
	borderline_condition.right_borderline_u = 1.1586268;
	borderline_condition.right_borderline_p = 7.0315580;

	//	Начальный условия
	StartCondition start_condition;
	start_condition.start_ro = std::function<double(double)>([expanse_grid, borderline_condition](double arg) {
		return (borderline_condition.right_borderline_ro - borderline_condition.left_borderline_ro) / (expanse_grid.h - expanse_grid.nodes) * arg + borderline_condition.left_borderline_ro;
		});
	start_condition.start_u = std::function<double(double)>([expanse_grid, borderline_condition](double arg) {
		return (borderline_condition.right_borderline_u - borderline_condition.left_borderline_u) / (expanse_grid.h - expanse_grid.nodes) * arg + borderline_condition.left_borderline_u;
		});
	start_condition.start_p = std::function<double(double)>([expanse_grid, borderline_condition](double arg) {
		return (borderline_condition.right_borderline_p - borderline_condition.left_borderline_p) / (expanse_grid.h - expanse_grid.nodes) * arg + borderline_condition.left_borderline_p;
		});


	//	Моделирование разрывного квазиодномерного течения в канале
	GasDynamicsEquation GDE(expanse_grid, time_grid, start_condition, borderline_condition, S_func);
	GDE.solving();

	GDE.writeRO(file_path_ro / "ro.raw");
	GDE.writeU(file_path_u / "u.raw");
	GDE.writeP(file_path_p / "p.raw");

	//	Моделирование целевого решения разрывного квазиодномерного течения в канале
	TargetGDE GDE_target(expanse_grid, borderline_condition, S_func, 0.75);
	GDE_target.solving();

	GDE_target.writeRO(file_path_ro / "t_ro.raw");
	GDE_target.writeU(file_path_u / "t_u.raw");
	GDE_target.writeP(file_path_p / "t_p.raw");

	//	Вычисление невязки
	std::vector<double> discrepancy_ro = abs_vec(GDE.getRO() - GDE_target.getRO());
	std::vector<double> discrepancy_u = abs_vec(GDE.getU() - GDE_target.getU());
	std::vector<double> discrepancy_p = abs_vec(GDE.getP() - GDE_target.getP());

	write_arr(file_path_ro / "discrep_ro.raw", discrepancy_ro.data(), discrepancy_ro.size());
	write_arr(file_path_u / "discrep_u.raw", discrepancy_u.data(), discrepancy_u.size());
	write_arr(file_path_p / "discrep_p.raw", discrepancy_p.data(), discrepancy_p.size());
}

//	Вычисление разницы течений
void diff_flow()
{
	ExpanseGrid expanse_grid{ 0.001, 1001, 0.0, 0.0 };
	TimeGrid time_grid{ 0.00001, 100000, 0.0, 0.0 };

	std::vector<double> GDE_RO(expanse_grid.nodes);	read_arr(".\\Comparison\\ro.raw", GDE_RO.data(), GDE_RO.size());
	std::vector<double> GDE_U(expanse_grid.nodes);	read_arr(".\\Comparison\\u.raw", GDE_U.data(), GDE_U.size());
	std::vector<double> GDE_P(expanse_grid.nodes);	read_arr(".\\Comparison\\p.raw", GDE_P.data(), GDE_P.size());

	std::vector<double> GDE_TVD_RO(expanse_grid.nodes);	read_arr(".\\Comparison\\ro_tvd.raw", GDE_TVD_RO.data(), GDE_TVD_RO.size());
	std::vector<double> GDE_TVD_U(expanse_grid.nodes);	read_arr(".\\Comparison\\u_tvd.raw", GDE_TVD_U.data(), GDE_TVD_U.size());
	std::vector<double> GDE_TVD_P(expanse_grid.nodes);	read_arr(".\\Comparison\\p_tvd.raw", GDE_TVD_P.data(), GDE_TVD_P.size());

	//	Вычисление невязки
	std::vector<double> difference_ro = abs_vec(GDE_RO - GDE_TVD_RO);
	std::vector<double> difference_u = abs_vec(GDE_U - GDE_TVD_U);
	std::vector<double> difference_p = abs_vec(GDE_P - GDE_TVD_P);

	write_arr(".\\Comparison\\diff_ro.raw", difference_ro.data(), difference_ro.size());
	write_arr(".\\Comparison\\diff_u.raw", difference_u.data(), difference_u.size());
	write_arr(".\\Comparison\\diff_p.raw", difference_p.data(), difference_p.size());
}