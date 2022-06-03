#include "matrices.h"

#include <cmath>
#include <cfloat>

#include "math.h"

mat2 operator*(const mat2& matrix, float scalar)
{
	mat2 result;
	for (int i = 0; i < 4; ++i)
	{
		result.asArray[i] = matrix.asArray[i] * scalar;
	}
	return result;
}

mat3 operator*(const mat3& _matrix, float _scalar) 
{
	mat3 result;
	for (int i = 0; i < 9; ++i)
	{
		result.asArray[i] = _matrix.asArray[i] * _scalar;
	}
	return result;
}

mat4 operator*(const mat4& _matrix, float _scalar)
{
	mat4 result;
	for (int i = 0; i < 16; ++i) 
	{
		result.asArray[i] = _matrix.asArray[i] * _scalar;
	}
	return result;
}

mat2 operator*(const mat2& _matA, const mat2& _matB) 
{
	mat2 res;
	Multiply(res.asArray, _matA.asArray, 2, 2, _matB.asArray, 2, 2);
	return res;
}

mat3 operator*(const mat3& _matA, const mat3& _matB) 
{
	mat3 res;
	Multiply(res.asArray, _matA.asArray, 3, 3, _matB.asArray, 3, 3);
	return res;
}

mat4 operator*(const mat4& _matA, const mat4& _matB) 
{
	mat4 res;
	Multiply(res.asArray, _matA.asArray, 4, 4, _matB.asArray, 4, 4);
	return res;
}

void Transpose(const float* _srcMat, float* _dstMat, int _srcRows, int _srcCols)
{
	for (int i = 0; i < _srcRows * _srcCols; i++)
	{
		int row = i / _srcRows;
		int col = i % _srcRows;
		_dstMat[i] = _srcMat[_srcCols * col + row];
	}
}

mat2 Transpose(const mat2& _matrix)
{
	mat2 result;
	Transpose(_matrix.asArray, result.asArray, 2, 2);
	return result;
}

mat3 Transpose(const mat3& _matrix)
{
	mat3 result;
	Transpose(_matrix.asArray, result.asArray, 3, 3);
	return result;
}

mat4 Transpose(const mat4& _matrix)
{
	mat4 result;
	Transpose(_matrix.asArray, result.asArray, 4, 4);
	return result;
}

bool Multiply(float* out, const float* matA, int aRows, int aCols, const float* matB, int bRows, int bCols) 
{
	if (aCols != bRows) 
	{
		return false;
	}

	for (int i = 0; i < aRows; ++i) 
	{
		for (int j = 0; j < bCols; ++j) 
		{
			out[bCols * i + j] = 0.0f;
			for (int k = 0; k < bRows; ++k) 
			{
				int a = aCols * i + k;
				int b = bCols * k + j;
				out[bCols * i + j] += matA[a] * matB[b];
			}
		}
	}
	return true;
}

float Determinant(const mat2& _matrix)
{
	return _matrix.m11 * _matrix.m22 -
		   _matrix.m12 * _matrix.m21;
}

float Determinant(const mat3& _mat) 
{
	float result = 0.0f;
	mat3 cofactor = Cofactor(_mat);
	for (int j = 0; j < 3; ++j) 
	{
		int index = 3 * 0 + j;
		result += _mat.asArray[index] * cofactor[0][j];
	}
	return result;
}

float Determinant(const mat4& _mat)
{
	float result = 0.0f;

	mat4 cofactor = Cofactor(_mat);
	for (int j = 0; j < 4; ++j) 
	{
		result += _mat.asArray[4 * 0 + j] * cofactor[0][j];
	}

	return result;
}

mat2 Cut(const mat3& _mat, int _row, int _col) 
{
	mat2 result;
	int index = 0;

	for (int i = 0; i < 3; ++i) 
	{
		for (int j = 0; j < 3; ++j) 
		{
			if (i == _row || j == _col)
				continue;

			int target = index++;
			int source = 3 * i + j;
			result.asArray[target] = _mat.asArray[source];
		}
	}

	return result;
}

mat3 Cut(const mat4& _mat, int _row, int _col)
{
	mat3 result;
	int index = 0;

	for (int i = 0; i < 4; ++i) 
	{
		for (int j = 0; j < 4; ++j) 
		{
			if (i == _row || j == _col)
				continue;

			int target = index++;
			int source = 4 * i + j;
			result.asArray[target] = _mat.asArray[source];
		}
	}

	return result;
}

mat2 Minor(const mat2& _mat) 
{
	return mat2(
		_mat.m22, _mat.m21,
		_mat.m12, _mat.m11
	);
}

mat3 Minor(const mat3& _mat) 
{
	mat3 result;

	for (int i = 0; i < 3; ++i) 
	{
		for (int j = 0; j < 3; ++j) 
		{
			result[i][j] = Determinant(Cut(_mat, i, j));
		}
	}

	return result;
}

mat4 Minor(const mat4& _mat)
{
	mat4 result;

	for (int i = 0; i < 4; ++i) 
	{
		for (int j = 0; j < 4; ++j) 
		{
			result[i][j] = Determinant(Cut(_mat, i, j));
		}
	}

	return result;
}

void Cofactor(float* _out, const float* _minor, int _rows, int _cols) 
{
	for (int i = 0; i < _rows; ++i) 
	{
		for (int j = 0; j < _cols; ++j) 
		{
			int t = _cols * j + i; // Target index
			int s = _cols * j + i; // Source index
			float sign = powf(-1.0f, i + j); // + or –
			_out[t] = _minor[s] * sign;
		}
	}
}

mat2 Cofactor(const mat2& _mat) 
{
	mat2 result;
	Cofactor(result.asArray, Minor(_mat).asArray, 2, 2);
	return result;
}

mat3 Cofactor(const mat3& _mat) 
{
	mat3 result;
	Cofactor(result.asArray, Minor(_mat).asArray, 3, 3);
	return result;
}

mat4 Cofactor(const mat4& _mat)
{
	mat4 result;
	Cofactor(result.asArray, Minor(_mat).asArray, 4, 4);
	return result;
}

mat2 Adjugate(const mat2& _mat)
{
	return Transpose(Cofactor(_mat));
}

mat3 Adjugate(const mat3& _mat)
{
	return Transpose(Cofactor(_mat));
}

mat4 Adjugate(const mat4& _mat)
{
	return Transpose(Cofactor(_mat));
}

mat2 Inverse(const mat2& _mat)
{
	float det = Determinant(_mat);
	if (CMP(det, 0.0f))
		return mat2();

	return Adjugate(_mat) * (1.0f / det);
}

mat3 Inverse(const mat3& _mat)
{
	float det = Determinant(_mat);
	if (CMP(det, 0.0f))
		return mat3();

	return Adjugate(_mat) * (1.0f / det);
}

mat4 Inverse(const mat4& _mat)
{
	float det = Determinant(_mat);
	if (CMP(det, 0.0f))
		return mat4();

	return Adjugate(_mat) * (1.0f / det);
}

mat4 Translation(float _x, float _y, float _z)
{
	return mat4(
		1.0f, 0.0f, 0.0f, 0.0f,
		0.0f, 1.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 1.0f, 0.0f,
		_x, _y, _z, 1.0f
	);
}

mat4 Translation(const vec3& _pos)
{
	return mat4(
		1.0f, 0.0f, 0.0f, 0.0f,
		0.0f, 1.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 1.0f, 0.0f,
		_pos.x, _pos.y, _pos.z, 1.0f
	);
}

vec3 GetTranslation(const mat4& _mat)
{
	return vec3(_mat.m41, _mat.m42, _mat.m43);
}

mat4 Scale(float _x, float _y, float _z)
{
	return mat4(
		_x, 0.0f, 0.0f, 0.0f,
		0.0f, _y, 0.0f, 0.0f,
		0.0f, 0.0f, _z, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	);
}

mat4 Scale(const vec3& _scale)
{
	return mat4(
		_scale.x, 0.0f, 0.0f, 0.0f,
		0.0f, _scale.y, 0.0f, 0.0f,
		0.0f, 0.0f, _scale.z, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	);
}

vec3 GetScale(const mat4& _mat)
{
	return vec3(_mat.m11, _mat.m22, _mat.m33);
}

mat4 Rotation(float _pitch, float _yaw, float _roll)
{
	return ZRotation(_roll) *
		XRotation(_pitch) *
		YRotation(_yaw);
}

mat3 Rotation3x3(float _pitch, float _yaw, float _roll)
{
	return ZRotation3x3(_roll) *
		XRotation3x3(_pitch) *
		YRotation3x3(_yaw);
}

mat4 XRotation(float _angle)
{
	_angle = DEG2RAD(_angle);
	return mat4(
		1.0f, 0.0f, 0.0f, 0.0f,
		0.0f, cosf(_angle), sinf(_angle), 0.0f,
		0.0f, -sinf(_angle), cos(_angle), 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	);
}

mat3 XRotation3x3(float _angle)
{
	_angle = DEG2RAD(_angle);
	return mat3(
		1.0f, 0.0f, 0.0f,
		0.0f, cosf(_angle), sinf(_angle),
		0.0f, -sinf(_angle), cos(_angle)
	);
}

mat4 YRotation(float _angle)
{
	_angle = DEG2RAD(_angle);
	return mat4(
		cosf(_angle), 0.0f, -sinf(_angle), 0.0f,
		0.0f, 1.0f, 0.0f, 0.0f,
		sinf(_angle), 0.0f, cosf(_angle), 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	);
}

mat3 YRotation3x3(float _angle)
{
	_angle = DEG2RAD(_angle);
	return mat3(
		cosf(_angle), 0.0f, -sinf(_angle),
		0.0f, 1.0f, 0.0f,
		sinf(_angle), 0.0f, cosf(_angle)
	);
}

mat4 ZRotation(float _angle)
{
	_angle = DEG2RAD(_angle);
	return mat4(
		cosf(_angle), sinf(_angle), 0.0f, 0.0f,
		-sinf(_angle), cosf(_angle), 0.0f, 0.0f,
		0.0f, 0.0f, 1.0f, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	);
}

mat3 ZRotation3x3(float _angle)
{
	_angle = DEG2RAD(_angle);
	return mat3(
		cosf(_angle), sinf(_angle), 0.0f,
		-sinf(_angle), cosf(_angle), 0.0f,
		0.0f, 0.0f, 1.0f
	);
}

mat4 AxisAngle(const vec3& _axis, float _angle)
{
	_angle = DEG2RAD(_angle);
	float c = cosf(_angle);
	float s = sinf(_angle);
	float t = 1.0f - cosf(_angle);

	float x = _axis.x;
	float y = _axis.y;
	float z = _axis.z;
	if (!CMP(MagnitudeSq(_axis), 1.0f)) 
	{
		float inv_len = 1.0f / Magnitude(_axis);
		x *= inv_len; // Normalize x
		y *= inv_len; // Normalize y
		z *= inv_len; // Normalize z
	} // x, y, and z are a normalized vector

	return mat4(
		t * (x * x) + c, t * x * y + s * z, t * x * z - s * y, 0.0f,
		t * x * y - s * z, t * (y * y) + c, t * y * z + s * x, 0.0f,
		t * x * z + s * y, t * y * z - s * x, t * (z * z) + c, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	);
}

mat3 AxisAngle3x3(const vec3& _axis, float _angle)
{
	_angle = DEG2RAD(_angle);
	float c = cosf(_angle);
	float s = sinf(_angle);
	float t = 1.0f - cosf(_angle);

	float x = _axis.x;
	float y = _axis.y;
	float z = _axis.z;
	if (!CMP(MagnitudeSq(_axis), 1.0f)) 
	{
		float inv_len = 1.0f / Magnitude(_axis);
		x *= inv_len;
		y *= inv_len;
		z *= inv_len;
	}

	return mat3(
		t * (x * x) + c, t * x * y + s * z, t * x * z - s * y,
		t * x * y - s * z, t * (y * y) + c, t * y * z + s * x,
		t * x * z + s * y, t * y * z - s * x, t * (z * z) + c
	);
}

vec3 MultiplyPoint(const vec3& _vec, const mat4& _mat)
{
	vec3 result;
	result.x = _vec.x * _mat.m11 + _vec.y * _mat.m21 +
			   _vec.z * _mat.m31 + 1.0f * _mat.m41;
	result.y = _vec.x * _mat.m12 + _vec.y * _mat.m22 +
			   _vec.z * _mat.m32 + 1.0f * _mat.m42;
	result.z = _vec.x * _mat.m13 + _vec.y * _mat.m23 +
			   _vec.z * _mat.m33 + 1.0f * _mat.m43;

	return result;
}

vec3 MultiplyVector(const vec3& _vec, const mat4& _mat)
{
	vec3 result;
	result.x = _vec.x * _mat.m11 + _vec.y * _mat.m21 +
			   _vec.z * _mat.m31 + 0.0f * _mat.m41;
	result.y = _vec.x * _mat.m12 + _vec.y * _mat.m22 +
			   _vec.z * _mat.m32 + 0.0f * _mat.m42;
	result.z = _vec.x * _mat.m13 + _vec.y * _mat.m23 +
			   _vec.z * _mat.m33 + 0.0f * _mat.m43;

	return result;
}

vec3 MultiplyVector(const vec3& _vec, const mat3& _mat)
{
	vec3 result;
	result.x = Dot(_vec, vec3(_mat.m11, _mat.m21, _mat.m31));
	result.y = Dot(_vec, vec3(_mat.m12, _mat.m22, _mat.m32));
	result.z = Dot(_vec, vec3(_mat.m13, _mat.m23, _mat.m33));

	return result;
}
