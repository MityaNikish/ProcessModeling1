#pragma once
#include "vector3D.h"

class Matrix3D
{
	Vector3D data_[3];

public:
	Matrix3D() = default;
	Matrix3D(const Vector3D& vec1, const Vector3D& vec2, const Vector3D& vec3);

	Matrix3D operator+(const Matrix3D& other) const;
	Matrix3D operator*(double value) const;
	Matrix3D operator/(double value) const;
	Vector3D operator*(const Vector3D& vec) const;
	Matrix3D operator*(const Matrix3D& other) const;
	Vector3D& operator[](Num index);
	Vector3D operator[](Num index) const;
};