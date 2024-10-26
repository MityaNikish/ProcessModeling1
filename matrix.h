#pragma once

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

	const Matrix operator+(const Matrix& other) const;
	void operator+=(const Matrix& other);
	const Matrix operator-(const Matrix& other) const;
	void operator-=(const Matrix& other);

	const Matrix operator-() const noexcept;

	const Matrix operator*(const Matrix& other) const;
	const Matrix operator*(const double arg) const noexcept;
	void operator*=(const double arg);

	Matrix operator[](const size_t index) const;
	Matrix operator()(const Matrix& other) const;

	bool operator==(const Matrix& other) const noexcept;

	size_t getQuantityRow(void) const noexcept;
	size_t getQuantityCal(void) const noexcept;
	double& getElement(size_t index_row, size_t index_col) const;
	Matrix getCut(size_t begin_index_row, size_t end_index_row, size_t begin_index_col, size_t end_index_col) const;
	double* getData(void) const;
};