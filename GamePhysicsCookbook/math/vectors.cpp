#include "vectors.h"

#include <cmath>
#include <cfloat>

#include "math.h"

vec2 operator+(const vec2& _lhs, const vec2& _rhs)
{
    return { _lhs.x + _rhs.x, _lhs.y + _rhs.y };
}

vec3 operator+(const vec3& _lhs, const vec3& _rhs)
{
    return { _lhs.x + _rhs.x, _lhs.y + _rhs.y, _rhs.z + _rhs.z };
}

vec2 operator-(const vec2& _lhs, const vec2& _rhs)
{
    return { _lhs.x - _rhs.x, _lhs.y - _rhs.y };
}

vec3 operator-(const vec3& _lhs, const vec3& _rhs)
{
    return { _lhs.x - _rhs.x, _lhs.y - _rhs.y, _rhs.z - _rhs.z };
}

vec2 operator*(const vec2& _lhs, const vec2& _rhs)
{
    return { _lhs.x * _rhs.x, _lhs.y * _rhs.y };
}

vec3 operator*(const vec3& _lhs, const vec3& _rhs)
{
    return { _lhs.x * _rhs.x, _lhs.y * _rhs.y, _rhs.z * _rhs.z };
}

vec2 operator*(const vec2& _lhs, float _rhs)
{
    return { _lhs.x * _rhs, _lhs.y * _rhs };
}

vec3 operator*(const vec3& _lhs, float _rhs)
{
    return { _lhs.x * _rhs, _lhs.y * _rhs, _lhs.z * _rhs };
}

bool operator==(const vec2& _lhs, const vec2& _rhs)
{
    return CMP(_lhs.x, _rhs.x) && CMP(_lhs.y, _rhs.y);
}

bool operator==(const vec3& _lhs, const vec3& _rhs)
{
    return CMP(_lhs.x, _rhs.x) && CMP(_lhs.y, _rhs.y) && CMP(_lhs.z, _rhs.z);
}

bool operator!=(const vec2& _lhs, const vec2& _rhs)
{
    return !(_lhs == _rhs);
}

bool operator!=(const vec3& _lhs, const vec3& _rhs)
{
    return !(_lhs == _rhs);
}

float Dot(const vec2& _lhs, const vec2& _rhs) 
{
    return _lhs.x * _rhs.x + _lhs.y * _rhs.y;
}

float Dot(const vec3& _lhs, const vec3& _rhs) 
{
    return _lhs.x * _rhs.x + _lhs.y * _rhs.y + _lhs.z * _rhs.z;
}

float Magnitude(const vec2& _vector)
{
    return sqrtf(Dot(_vector, _vector));
}

float Magnitude(const vec3& _vector)
{
    return sqrtf(Dot(_vector, _vector));
}

float MagnitudeSq(const vec2& _vector)
{
    return Dot(_vector, _vector);
}

float MagnitudeSq(const vec3& _vector)
{
    return Dot(_vector, _vector);
}

float Distance(const vec2& _a, const vec2& _b)
{
    return Magnitude(_a - _b);
}

float Distance(const vec3& _a, const vec3& _b)
{
    return Magnitude(_a - _b);
}

void Normalize(vec2& _vector)
{
    _vector = _vector * (1.0f / Magnitude(_vector));
}

void Normalize(vec3& _vector)
{
    _vector = _vector * (1.0f / Magnitude(_vector));
}

vec2 Normalized(const vec2& _vector)
{
    return _vector * (1.0f / Magnitude(_vector));
}

vec3 Normalized(const vec3& _vector)
{
    return _vector * (1.0f / Magnitude(_vector));
}

vec3 Cross(const vec3& _lhs, const vec3& _rhs)
{
    vec3 result;
    result.x = _lhs.y * _rhs.z - _lhs.z * _rhs.y;
    result.y = _lhs.z * _rhs.x - _lhs.x * _rhs.z;
    result.z = _lhs.x * _rhs.y - _lhs.y * _rhs.x;
    return result;
}

float Angle(const vec2& _lhs, const vec2& _rhs)
{
    float m = sqrtf(MagnitudeSq(_lhs) * MagnitudeSq(_rhs));
    return acos(Dot(_lhs, _rhs) / m);
}

float Angle(const vec3& _lhs, const vec3& _rhs)
{
    float m = sqrtf(MagnitudeSq(_lhs) * MagnitudeSq(_rhs));
    return acos(Dot(_lhs, _rhs) / m);
}

vec2 Project(const vec2& _length, const vec2& _direction) 
{
    float dot = Dot(_length, _direction);
    float magSq = MagnitudeSq(_direction);
    return _direction * (dot / magSq);
}

vec3 Project(const vec3& _length, const vec3& _direction) 
{
    float dot = Dot(_length, _direction);
    float magSq = MagnitudeSq(_direction);
    return _direction * (dot / magSq);
}

vec2 Perpendicular(const vec2& _len, const vec2& _dir) 
{
    return _len - Project(_len, _dir);
}

vec3 Perpendicular(const vec3& _len, const vec3& _dir) 
{
    return _len - Project(_len, _dir);
}

vec2 Reflection(const vec2& _vec, const vec2& _normal)
{
    float d = Dot(_vec, _normal);
    return _vec - _normal * (d * 2.0f);
}

vec3 Reflection(const vec3& _vec, const vec3& _normal)
{
    float d = Dot(_vec, _normal);
    return _vec - _normal * (d * 2.0f);
}
