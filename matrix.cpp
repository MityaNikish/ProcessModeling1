#include "matrix.h"

#include <memory>

Matrix::Matrix() : ptr_(nullptr), row_(0), col_(0) { }

Matrix::Matrix(size_t row, size_t col) : row_(row), col_(col)
{
	ptr_ = new double[row_ * col_];

	for (size_t i = 0; i < row_ * col_; ++i)
	{
		ptr_[i] = 0;
	}
}

Matrix::~Matrix()
{
	delete[] ptr_;
}


Matrix::Matrix(const Matrix& other)
{
	*this = other;
}

Matrix::Matrix(Matrix&& other) noexcept
{
	*this = other;
}


Matrix& Matrix::operator=(const Matrix& other)
{
	if (this == &other)
	{
		return *this;
	}
	                        
	delete[] ptr_;
	ptr_ = new double[other.col_ * other.row_];
	row_ = other.row_;
	col_ = other.col_;
	std::copy_n(other.ptr_, other.col_ * other.row_, ptr_);
	return *this;
}

Matrix& Matrix::operator=(Matrix&& other) noexcept
{
	delete[] ptr_;
	ptr_ = other.ptr_;
	other.ptr_ = nullptr;

	row_ = other.row_;
	other.row_ = 0;

	col_ = other.col_;
	other.col_ = 0;

	return *this;
}


const Matrix Matrix::operator+(const Matrix& other) const
{
	Matrix result(*this);

	if (this->row_ != other.row_ || this->col_ != other.col_)
	{
		throw "Исключение: Размерности матриц не совпадают.";
	}

	for (size_t i = 0; i < row_; ++i)
	{
		for (size_t j = 0; j < col_; ++j)
		{
			result.getElement(i, j) += other.getElement(i, j);
		}
	}

	return result;
}

void Matrix::operator+=(const Matrix& other)
{
	*this = *this + other;
}


const Matrix Matrix::operator-(const Matrix& other) const
{
	return Matrix(*this + -other);
}

void Matrix::operator-=(const Matrix& other)
{
	*this = *this - other;
}


const Matrix Matrix::operator-() const noexcept
{
	Matrix result(*this);

	for (size_t i = 0; i < row_; ++i)
	{
		for (size_t j = 0; j < col_; ++j)
		{
			result.getElement(i, j) = -this->getElement(i, j);
		}
	}

	return result;
}



const Matrix Matrix::operator*(const Matrix& other) const
{
	Matrix result(this->row_, other.col_);

	if (this->row_ != other.col_ || this->col_ != other.row_)
	{
		throw "Исключение: Размерности матриц не совпадают.";
	}

	for (size_t i = 0; i < row_; ++i)
	{
		for (size_t j_ = 0; j_ < other.col_; ++j_)
		{
			double value = 0;
			for (size_t j = 0; j < col_; ++j)
			{
				value += this->getElement(i, j) * other.getElement(j, j_);
			}
			result.getElement(i, j_) = value;
		}
	}

	return result;
}

const Matrix Matrix::operator*(const double arg) const noexcept
{
	Matrix result(*this);

	for (size_t i = 0; i < row_; ++i)
	{
		for (size_t j = 0; j < col_; ++j)
		{
			result.getElement(i, j) *= arg;
		}
	}

	return result;
}


void Matrix::operator*=(const double arg)
{
	*this = *this * arg;
}


Matrix Matrix::operator[](const size_t index) const
{
	return this->getCut(index, index + 1, 0, col_);
}

Matrix Matrix::operator()(const Matrix& other) const
{
	Matrix result(row_, 1);

	if (this->row_ != other.col_ || other.col_ != 1)
	{
		throw "Исключение: Размерности матриц не совпадают.";
	}

	for (size_t i = 0; i < row_; ++i)
	{
		for (size_t j = 0; j < col_; ++j)
		{
			result.getElement(i, j) *= other.getElement(j, 0);
		}
	}

	return result;
}

bool Matrix::operator==(const Matrix& other) const noexcept
{
	if (this->row_ != other.row_ || this->col_ != other.col_)
	{
		return false;
	}

	for (size_t i = 0; i < row_; ++i)
	{
		for (size_t j = 0; j < col_; ++j)
		{
			if (this->getElement(i, j) - other.getElement(i, j) > 0.001) return false;
		}
	}
	return true;
}


size_t Matrix::getQuantityRow(void) const noexcept
{
	return row_;
}

size_t Matrix::getQuantityCal(void) const noexcept
{
	return col_;
}

double& Matrix::getElement(size_t index_row, size_t index_col) const
{
	if (index_row >= row_ || index_col >= col_)
	{
		throw "Исключение: Обращение к не существующему элементу";
	}
	return ptr_[index_row * col_ + index_col];
}

Matrix Matrix::getCut(size_t begin_index_row, size_t end_index_row, size_t begin_index_col, size_t end_index_col) const
{
	Matrix cutting(end_index_row - begin_index_row, end_index_col - begin_index_col);

	for (size_t i = 0; i < cutting.getQuantityRow(); ++i)
	{
		for (size_t j = 0; j < cutting.getQuantityCal(); ++j)
		{
			cutting.getElement(i, j) = ptr_[(i + begin_index_row) * col_ + j + begin_index_col];
		}
	}

	return cutting;
}

double* Matrix::getData(void) const
{
	return ptr_;
}