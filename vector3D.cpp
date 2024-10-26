#include "vector3D.h"

Vector3D::Vector3D()
{
	data_[0] = 0;
	data_[1] = 0;
	data_[2] = 0;
}

Vector3D::Vector3D(const double val1, const double val2, const double val3)
{
	data_[0] = val1;
	data_[1] = val2;
	data_[2] = val3;
}

Vector3D Vector3D::operator+(const Vector3D& other) const
{
	return Vector3D{ data_[0] + other.data_[0], data_[1] + other.data_[1], data_[2] + other.data_[2] };
}

Vector3D Vector3D::operator-(const Vector3D& other) const
{
	return *this + other * (-1.0);
}

Vector3D Vector3D::operator*(const double value) const
{
	return Vector3D{ data_[0] * value, data_[1] * value, data_[2] * value };
}

Vector3D Vector3D::operator/(const double value) const
{
	return *this * (1 / value);
}

double Vector3D::operator*(const Vector3D& other) const
{
	return data_[0] * other.data_[0] + data_[1] * other.data_[1] + data_[2] * other.data_[2];
}

double& Vector3D::operator[](const Num index)
{
	return data_[static_cast<size_t>(index)];
}

double Vector3D::operator[](const Num index) const
{
	return data_[static_cast<size_t>(index)];
}