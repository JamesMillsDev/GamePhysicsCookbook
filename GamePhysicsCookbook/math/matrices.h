#pragma once

#include "vectors.h"

typedef struct mat2
{
	union
	{
		struct
		{
			float m11, m12,
				  m21, m22;
		};

		float asArray[4];
	};

	inline float* operator[](int _index)
	{
		return &(asArray[_index * 2]);
	}

	inline mat2()
	{
		m11 = m22 = 1.0f;
		m12 = m21 = 0.0f;
	}

	inline mat2(float _f11, float _f12, float _f21, float _f22) 
	{
		m11 = _f11; m12 = _f12;
		m21 = _f21; m22 = _f22;
	}

} mat2;

typedef struct mat3 
{
	union 
	{
		struct 
		{
			float m11, m12, m13,
				  m21, m22, m23,
				  m31, m32, m33;
		};
		float asArray[9];
	};

	inline float* operator[](int _index) 
	{
		return &(asArray[_index * 3]);
	}

	inline mat3() 
	{
		m11 = m22 = m33 = 1.0f;
		m12 = m13 = m21 = 0.0f;
		m23 = m31 = m32 = 0.0f;
	}

	inline mat3(float _f11, float _f12, float _f13, float _f21, float _f22, float _f23, float _f31, float _f32, float _f33) 
	{
		m11 = _f11; m12 = _f12; m13 = _f13;
		m21 = _f21; m22 = _f22; m23 = _f23;
		m31 = _f31; m32 = _f32; m33 = _f33;
	}

} mat3;

typedef struct mat4 
{
	union 
	{
		struct 
		{
			float m11, m12, m13, m14,
				  m21, m22, m23, m24,
				  m31, m32, m33, m34,
				  m41, m42, m43, m44;
		};
		float asArray[16];
	};

	inline float* operator[](int _index) 
	{
		return &(asArray[_index * 4]);
	}

	inline mat4() 
	{
		m11 = m22 = m33 = m44 = 1.0f;
		m12 = m13 = m14 = m21 = 0.0f;
		m23 = m24 = m31 = m32 = 0.0f;
		m34 = m41 = m42 = m43 = 0.0f;
	}

	inline mat4(float _f11, float _f12, float _f13, float _f14, float _f21, float _f22, float _f23, float _f24, float _f31, float _f32, float _f33, float _f34, float _f41, float _f42, float _f43, float _f44) 
	{
		m11 = _f11; m12 = _f12; m13 = _f13; m14 = _f14;
		m21 = _f21; m22 = _f22; m23 = _f23; m24 = _f24;
		m31 = _f31; m32 = _f32; m33 = _f33; m34 = _f34;
		m41 = _f41; m42 = _f42; m43 = _f43; m44 = _f44;
	}

} mat4;

// ------------ Operators ------------ //

mat2 operator*(const mat2& _matrix, float _scalar);
mat3 operator*(const mat3& _matrix, float _scalar);
mat4 operator*(const mat4& _matrix, float _scalar);

mat2 operator*(const mat2& _matA, const mat2& _matB);
mat3 operator*(const mat3& _matA, const mat3& _matB);
mat4 operator*(const mat4& _matA, const mat4& _matB);

void Transpose(const float* _srcMat, float* _dstMat, int _srcRows, int _srcCols);
mat2 Transpose(const mat2& _matrix);
mat3 Transpose(const mat3& _matrix);
mat4 Transpose(const mat4& _matrix);

bool Multiply(float* _out, const float* _matA, int _aRows, int _aCols, const float* _matB, int _bRows, int _bCols);

float Determinant(const mat2& _matrix);
float Determinant(const mat3& _mat);
float Determinant(const mat4& _mat);

mat2 Cut(const mat3& _mat, int _row, int _col);
mat3 Cut(const mat4& mat, int row, int col);
mat2 Minor(const mat2& _mat);
mat3 Minor(const mat3& _mat);
mat4 Minor(const mat4& _mat);

void Cofactor(float* _out, const float* _minor, int _rows, int _cols);
mat3 Cofactor(const mat3& _mat);
mat2 Cofactor(const mat2& _mat);
mat4 Cofactor(const mat4& _mat);

mat2 Adjugate(const mat2& _mat);
mat3 Adjugate(const mat3& _mat);
mat4 Adjugate(const mat4& _mat);

mat2 Inverse(const mat2& _mat);
mat3 Inverse(const mat3& _mat);
mat4 Inverse(const mat4& _mat);

// ------------ Transformations ------------ //

mat4 Translation(float _x, float _y, float _z);
mat4 Translation(const vec3& _pos);
vec3 GetTranslation(const mat4& _mat);

mat4 Scale(float _x, float _y, float _z);
mat4 Scale(const vec3& _scale);
vec3 GetScale(const mat4& _mat);

mat4 Rotation(float _pitch, float _yaw, float _roll);
mat3 Rotation3x3(float _pitch, float _yaw, float _roll);

mat4 XRotation(float _angle);
mat3 XRotation3x3(float _angle);
mat4 YRotation(float _angle);
mat3 YRotation3x3(float _angle);
mat4 ZRotation(float _angle);
mat3 ZRotation3x3(float _angle);

mat4 AxisAngle(const vec3& _axis, float _angle);
mat3 AxisAngle3x3(const vec3& _axis, float _angle);

vec3 MultiplyPoint(const vec3& _vec, const mat4& _mat);
vec3 MultiplyVector(const vec3& _vec, const mat4& _mat);
vec3 MultiplyVector(const vec3& _vec, const mat3& _mat);