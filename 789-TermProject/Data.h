#pragma once

#include <cmath>
#include <iostream>
#include <assert.h>
#include <vector>
#include <iomanip>
using std::vector;

#include "Settings.h"

typedef unsigned int uint;
typedef long unsigned int luint;

class Mat3;
class Mat4;
class Vec3;
class Vec4;
class Quaternion;





typedef class Vec3 {
public:
	float x, y, z;
	Vec3();
	Vec3(const float n);
	Vec3(const float x1, const float y1, const float z1);

	Vec3 normalized();
	void normalize();

	float magnitude();

	float* toFloat3();


	// OPERATORS
	Vec3 operator/(const float& n) const;
	Vec3 operator/(const   int& n) const; // Scaler Div 
	Vec3 operator/(const  uint& n) const; // Scaler Div 

	friend Vec3 operator/(const float& n, const Vec3& vec);
	Vec3 operator^(const float& n) const;
	friend Vec3 operator^(const float& n, const Vec3& vec);
	Vec3 operator*(const float& n) const;
	friend Vec3 operator*(const float& n, const Vec3& vec);
	Vec3 operator+(const float& n) const;
	friend Vec3 operator+(const float& n, const Vec3& vec);
	Vec3 operator-(const float& n) const;
	friend Vec3 operator-(const float& n, const Vec3& vec);

	Vec3 operator*(const Vec3 rhs) const;
	Vec3 operator/(const Vec3 rhs) const;
	Vec3 operator-(const Vec3& rhs) const;
	Vec3 operator+(const Vec3& rhs) const;

	bool operator==(const Vec3& rhs) const;
	float& operator[](int idx);
	float operator[](int idx) const;
	friend std::ostream& operator<<(std::ostream& os, Vec3 rhs);
} vec3;

typedef class Vec4 {
public:
	float x, y, z, w;
	Vec4();
	Vec4(const float n);
	Vec4(Vec3 cpy, const float n);
	Vec4(const float n, Vec3 cpy);
	Vec4(const float x1, const float y1, const float z1, const float w1);

	Vec4 normalized();
	void normalize();

	float magnitude();
	float* toFloat4();

	// OPERATORS
	Vec4 operator/(const float& n) const;
	friend Vec4 operator/(const float& n, const Vec4& vec);
	Vec4 operator^(const float& n) const;
	friend Vec4 operator^(const float& n, const  Vec4& vec);
	Vec4 operator*(const float& n) const;
	friend Vec4 operator*(const float& n, const  Vec4& vec);
	Vec4 operator+(const float& n) const;
	friend Vec4 operator+(const float& n, const Vec4& vec);
	Vec4 operator-(const float& n) const;
	friend Vec4 operator-(const float& n, const Vec4& vec);

	Vec4 operator*(const Vec4& rhs) const;
	Vec4 operator/(const Vec4& rhs) const;
	Vec4 operator-(const Vec4& rhs) const;
	Vec4 operator+(const Vec4& rhs) const;

	float& operator[](int idx);
	float operator[](int idx) const;
	friend std::ostream& operator<<(std::ostream& os, Vec4 rhs);
} vec4;

typedef class Mat3 {
public:
	Vec3 x, y, z;

	Mat3();
	Mat3(const float n);
	Mat3(const float*& array);
	Mat3(float**& array);
	Mat3(vector<float>& array);
	Mat3(const Vec3& _x, const Vec3& _y, const Vec3& _z);


	friend std::ostream& operator<<(std::ostream& os, const Mat3& mat);
	Mat3 operator+(const float& n) const;
	friend Mat3 operator+(const float& n, const Mat3& mat);
	Mat3 operator-(const float& n) const;
	friend Mat3 operator-(const float& n, const Mat3& mat);
	Mat3 operator*(const float& n) const;
	friend Mat3 operator*(const float& n, const Mat3& mat);
	Mat3 operator/(const float& n) const;
	friend Mat3 operator/(const float& n, const Mat3& mat);


	Mat3 operator+(const Mat3& rhs) const;
	Mat3 operator-(const Mat3& rhs) const;
	Mat3 operator*(const Mat3& rhs) const;
	Vec3 operator*(const Vec3& vec) const;
	friend Vec3 operator*(const Vec3& vec, const Mat3& _mat);
	Mat3 operator/(const Mat3& rhs) const;
	Vec3& operator[](int idx);
	Vec3 operator[](int idx) const;

	float* toFloat3();
	float** toFloat4_2();


} mat3;

typedef class Mat4 {
public:
	Vec4 x, y, z, w;


	Mat4();
	Mat4(const float n);
	Mat4(float* array);
	Mat4(float** array);
	Mat4(vector<float>& array);
	Mat4(Mat3& x3, Vec3& col, Vec3& row, const float deep);
	Mat4(const Vec4& _x, const Vec4& _y, const Vec4& _z, const Vec4& _w);

	friend std::ostream& operator<<(std::ostream& os, const Mat4& mat);
	Mat4 operator+(const float& n) const;
	friend Mat4 operator+(const float& n, const Mat4& mat);
	Mat4 operator-(const float& n) const;
	friend Mat4 operator-(const float& n, const Mat4& mat);
	Mat4 operator*(const float& n) const;
	friend Mat4 operator*(const float& n, const Mat4& mat);
	Mat4 operator/(const float& n) const;
	friend Mat4 operator/(const float& n, const Mat4& mat);



	Mat4 operator+(const Mat4& rhs) const;
	Mat4 operator-(const Mat4& rhs) const;
	Mat4 operator*(const Mat4& rhs) const;
	Vec4 operator*(const Vec4& vec) const;
	friend Vec4 operator*(const Vec4& vec, const Mat4& _mat);
	Mat4 operator/(const Mat4& rhs) const;
	

	float* toFloat4();
	float** toFloat4_2();

	Vec4& operator[](int idx);
	Vec4 operator[](int idx) const;

	
} mat4;


class Data {
private://utility



public:
	Data() {}
   ~Data() {}

	static Mat3 outer(const Vec3& c, const Vec3& r);
	static Mat4 outer(const Vec4& c, const Vec4& r);

	static Vec3 cross(const Vec3& lhs, const Vec3& rhs);


	static float dot(const Vec3& lhs, const Vec3& rhs);
	static float dot(const Vec4& lhs, const Vec4& rhs);


	static Mat4 transpose(const Mat4& mat);
	static Mat3 transpose(const Mat3& mat);



	template <class T>
	static T adjoint(T mat, int size);

	template <class T>
	static float determinant(T mat, int n, int size);

	template <class T>
	static T getCofactor(T A, int p, int q, int n);

	
	static mat3 inverse( mat3& mat);
	static mat4 inverse( mat4& mat);

};
