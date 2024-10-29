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

Matrix Matrix::operator[](const size_t index) const
{
	return this->getCut(index, index + 1, 0, col_);
}

Matrix& Matrix::operator()(const std::function<double(double)>& func)
{
	for (size_t i = 0; i < row_; ++i)
	{
		for (size_t j = 0; j < col_; ++j)
		{
			getElement(i, j) = func(const_cast<const double&>(getElement(i, j)));
		}
	}
	return *this;
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

double& Matrix::getElement(int index_row, int index_col) const
{
	size_t relative_index_row = index_row < 0 ? row_ - 1 - index_row % row_ : index_row % row_;
	size_t relative_index_col = index_col < 0 ? col_ - 1 - index_col % col_ : index_col % col_;

	//if (index_row >= row_ || index_col >= col_)
	//{
	//	throw "Исключение: Обращение к не существующему элементу";
	//}

	return ptr_[relative_index_row * col_ + relative_index_col];
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