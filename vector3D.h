#pragma once

enum class Num
{
	one, two, three
};

class Vector3D
{
	double data_[3];

public:
	Vector3D();
	Vector3D(double val1, double val2, double val3);

	Vector3D operator+(const Vector3D& other) const;
	Vector3D operator-(const Vector3D& other) const;
	Vector3D operator*(double value) const;
	Vector3D operator/(double value) const;
	double operator*(const Vector3D& other) const;
	double& operator[](Num index);
	double operator[](Num index) const;
};