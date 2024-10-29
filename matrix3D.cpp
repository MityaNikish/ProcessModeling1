#include "matrix3D.h"
#include <limits>

namespace
{
	using enum Num;
}

Matrix3D::Matrix3D(const Vector3D& vec1, const Vector3D& vec2, const Vector3D& vec3)
{
	data_[0] = vec1;
	data_[1] = vec2;
	data_[2] = vec3;
}

Matrix3D Matrix3D::operator+(const Matrix3D& other) const
{
	return Matrix3D{ data_[0] + other.data_[0], data_[1] + other.data_[1], data_[2] + other.data_[2] };
}

Matrix3D Matrix3D::operator*(double value) const
{
	return Matrix3D{ data_[0] * value, data_[1] * value, data_[2] * value };
}

Matrix3D Matrix3D::operator/(double value) const
{
	return *this * (1 / value);
}

Vector3D Matrix3D::operator*(const Vector3D& vec) const
{
	return Vector3D{ (*this)[one] * vec , (*this)[two] * vec, (*this)[three] * vec };
}

Matrix3D Matrix3D::operator*(const Matrix3D& other) const
{
	Matrix3D result;

	Vector3D other_one = Vector3D{ other[one][one], other[two][one], other[three][one] };
	Vector3D other_two = Vector3D{ other[one][two], other[two][two], other[three][two] };
	Vector3D other_three = Vector3D{ other[one][three], other[two][three], other[three][three] };

	result[one][one] = (*this)[one] * other_one;
	result[one][two] = (*this)[one] * other_two;
	result[one][three] = (*this)[one] * other_three;

	result[two][one] = (*this)[two] * other_one;
	result[two][two] = (*this)[two] * other_two;
	result[two][three] = (*this)[two] * other_three;

	result[three][one] = (*this)[three] * other_one;
	result[three][two] = (*this)[three] * other_two;
	result[three][three] = (*this)[three] * other_three;

	return result;
}

Vector3D& Matrix3D::operator[](Num index)
{
	return data_[static_cast<size_t>(index)];
}

Vector3D Matrix3D::operator[](Num index) const
{
	return data_[static_cast<size_t>(index)];
}

double Matrix3D::getMaxElem()
{
	double max_element = DBL_MIN;
	for (size_t i = 0; i < 3; i++)
	{
		for (size_t j = 0; j < 3; j++)
		{
			double element = (*this)[static_cast<Num>(i)][static_cast<Num>(j)];
			if (element > max_element) max_element = element;
		}
	}
	return max_element;
}
