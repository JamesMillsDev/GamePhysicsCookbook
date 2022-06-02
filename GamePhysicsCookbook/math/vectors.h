#pragma once

#define ABSOLUTE(x, y) (fabsf(x-y) <= FLT_EPSILON)
#define RELATIVE(x, y) (fabsf(x-y) <= FLT_EPSILON * fmaxf(fabsf(x), fabsf(y)))

#define RAD2DEG(x) ((x) * 57.295754f)
#define DEG2RAD(x) ((x) * 0.0174533f)

typedef struct vec2
{
	union
	{
		struct
		{
			float x;
			float y;
		};

		float asArray[2];
	};

	float& operator[](int _index)
	{
		return asArray[_index];
	}

	vec2() : x(0.0f), y(0.0f) { }
	vec2(float _x, float _y) : x(_x), y(_y) { }

} vec2;

typedef struct vec3
{
	union
	{
		struct
		{
			float x;
			float y;
			float z;
		};

		float asArray[3];
	};

	float& operator[](int _index)
	{
		return asArray[_index];
	}

	vec3() : x(0.0f), y(0.0f), z(0.0f) { }
	vec3(float _x, float _y, float _z) : x(_x), y(_y), z(_z) { }

} vec3;

vec2 operator+(const vec2& _lhs, const vec2& _rhs);
vec3 operator+(const vec3& _lhs, const vec3& _rhs);
vec2 operator-(const vec2& _lhs, const vec2& _rhs);
vec3 operator-(const vec3& _lhs, const vec3& _rhs);
vec2 operator*(const vec2& _lhs, const vec2& _rhs);
vec3 operator*(const vec3& _lhs, const vec3& _rhs);
vec2 operator*(const vec2& _lhs, float _rhs);
vec3 operator*(const vec3& _lhs, float _rhs);
bool operator==(const vec2& _lhs, const vec2& _rhs);
bool operator==(const vec3& _lhs, const vec3& _rhs);
bool operator!=(const vec2& _lhs, const vec2& _rhs);
bool operator!=(const vec3& _lhs, const vec3& _rhs);

float Dot(const vec2& _lhs, const vec2& _rhs);
float Dot(const vec3& _lhs, const vec3& _rhs);

float Magnitude(const vec2& _vector);
float Magnitude(const vec3& _vector);

float MagnitudeSq(const vec2& _vector);
float MagnitudeSq(const vec3& _vector);

float Distance(const vec2& _a, const vec2& _b);
float Distance(const vec3& _a, const vec3& _b);

void Normalize(vec2& _vector);
void Normalize(vec3& _vector);

vec2 Normalized(const vec2& _vector);
vec3 Normalized(const vec3& _vector);

vec3 Cross(const vec3& _lhs, const vec3& _rhs);

float Angle(const vec2& _lhs, const vec2& _rhs);
float Angle(const vec3& _lhs, const vec3& _rhs);

vec2 Project(const vec2& _length, const vec2& _direction);
vec3 Project(const vec3& _length, const vec3& _direction);

vec2 Perpendicular(const vec2& _len, const vec2& _dir);
vec3 Perpendicular(const vec3& _len, const vec3& _dir);

vec2 Reflection(const vec2& _vec, const vec2& _normal);
vec3 Reflection(const vec3& _vec, const vec3& _normal);