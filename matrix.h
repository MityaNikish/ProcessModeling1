#pragma once
#include <functional>

class Matrix
{
	double* ptr_;
	size_t row_;
	size_t col_;

public:
	Matrix();
	Matrix(size_t row, size_t col);
	~Matrix();

	Matrix(const Matrix& other);
	Matrix(Matrix&& other) noexcept;

	Matrix& operator=(const Matrix& other);
	Matrix& operator=(Matrix&& other) noexcept;

	Matrix operator[](const size_t index) const;
	Matrix& operator()(const std::function<double(double)>& func);

	bool operator==(const Matrix& other) const noexcept;

	size_t getQuantityRow(void) const noexcept;
	size_t getQuantityCal(void) const noexcept;
	double& getElement(int index_row, int index_col) const;
	Matrix getCut(size_t begin_index_row, size_t end_index_row, size_t begin_index_col, size_t end_index_col) const;
	double* getData(void) const;
};