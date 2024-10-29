#pragma once
#include "solver_GDE.h"

//	VariableCrossSectionPipeGasDynamicsEquation
class VariableCrossSectionPipeGDE final : public GasDynamicsEquation
{
	//	Площадь поперечного сечения трубы
	const std::vector<double> _S;
	//Matrix3D _A_const;

public:
	VariableCrossSectionPipeGDE(const ExpanseGrid& expanse_grid, const TimeGrid& time_grid, StartCondition& start_condition, BorderlineCondition& borderline_condition, const std::vector<double>& S);

	//void solving() override;

private:
	void initConditions() override;
	void postProcessing() override;

	bool chekStopConditions(size_t n, size_t j, double eps) override;

	Vector3D Q(size_t n, size_t j) const override;
};

