#include "solver_VCSP_GDE.h"
#include <fstream>
#include "test.h"

#include <iostream>

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

int main()
{
	//	Инициализация путей к рабочим файлам
	std::filesystem::path file_path_ro = std::filesystem::current_path() / "Treatment\\output_ro.raw";
	std::filesystem::path file_path_u = std::filesystem::current_path() / "Treatment\\output_u.raw";
	std::filesystem::path file_path_p = std::filesystem::current_path() / "Treatment\\output_p.raw";
	std::filesystem::path file_path_grid_info = std::filesystem::current_path() / "grid_info.txt";

	//	Инициализация сетки
	ExpanseGrid expanse_grid;
	TimeGrid time_grid;

	////	Частные параметры сетки
	read_grid_info(file_path_grid_info, expanse_grid, time_grid);
	expanse_grid.starting_point = 0.0;
	expanse_grid.ending_point = expanse_grid.starting_point + expanse_grid.h * (expanse_grid.nodes - 1);
	time_grid.starting_point = 0.0;
	time_grid.ending_point = time_grid.starting_point + time_grid.tau * (time_grid.nodes - 1);
	////

	//expanse_grid.starting_point = 0;
	//expanse_grid.ending_point = 1;
	//expanse_grid.nodes = 101;
	//expanse_grid.h = (expanse_grid.ending_point - expanse_grid.starting_point) / (expanse_grid.nodes - 1);

	//time_grid.starting_point = 0;
	//time_grid.ending_point = 2.5;
	//time_grid.nodes = 101;
	//time_grid.tau = (time_grid.ending_point - time_grid.starting_point) / (time_grid.nodes - 1);
	

	//	Функция попреречного сечения трубы
	std::vector<double> S(expanse_grid.nodes);
	const double pos_bottleneck = 0.5;
	const double value_bottleneck = 0.5;

	for (size_t i = 0; i < expanse_grid.nodes; ++i)
	{
		const double part = (1 - i * expanse_grid.h / pos_bottleneck);
		S[i] = value_bottleneck + (1 - value_bottleneck) * part * part;
	}

	//	Начальный условия
	StartCondition start_condition;
	start_condition.start_ro = std::function<double(double)>([pos_bottleneck](double arg) { return arg < pos_bottleneck ? 1.0 : 0.1933880; });
	start_condition.start_u = std::function<double(double)>([pos_bottleneck](double arg) { return arg < pos_bottleneck ? 1.0237498 : 5.2937598; });
	start_condition.start_p = std::function<double(double)>([pos_bottleneck](double arg) { return arg < pos_bottleneck ? 8.0 : 0.8018469; });

	//	Граничные условия
	BorderlineCondition borderline_condition;
	borderline_condition.start_borderline_ro = std::function<double(double)>([](double arg) { return 1.0; });
	borderline_condition.start_borderline_u = std::function<double(double)>([](double arg) { return 1.0237498; });
	borderline_condition.start_borderline_p = std::function<double(double)>([](double arg) { return 8.0; });
	borderline_condition.end_borderline_ro = std::function<double(double)>([](double arg) { return 0.1933880; });
	borderline_condition.end_borderline_u = std::function<double(double)>([](double arg) { return 5.2937598; });
	borderline_condition.end_borderline_p = std::function<double(double)>([](double arg) { return 0.8018469; });

	//	Моделирование ударной трубы
	VariableCrossSectionPipeGDE solver(expanse_grid, time_grid, start_condition, borderline_condition, S);
	solver.solving();

	//	Запись результатов в файл 
	//solver.writeRO(file_path_ro);
	//solver.writeU(file_path_u);
	//solver.writeP(file_path_p);
	
	////	Тестирование некоторых функций
	//TEST_all();

	return 0;
}