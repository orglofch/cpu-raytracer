/*
* Copyright (c) 2015 Owen Glofcheski
*
* This software is provided 'as-is', without any express or implied
* warranty. In no event will the authors be held liable for any damages
* arising from the use of this software.
*
* Permission is granted to anyone to use this software for any purpose,
* including commercial applications, and to alter it and redistribute it
* freely, subject to the following restrictions:
*
*    1. The origin of this software must not be misrepresented; you must not
*    claim that you wrote the original software. If you use this software
*    in a product, an acknowledgment in the product documentation would be
*    appreciated but is not required.
*
*    2. Altered source versions must be plainly marked as such, and must not
*    be misrepresented as being the original software.
*
*    3. This notice may not be removed or altered from any source
*    distribution.
*/

#ifndef _COLOUR_HPP_
#define _COLOUR_HPP_

union Colour
{
	Colour() {
		r = g = b = 0;
	}
	Colour(float _r, float _g, float _b) {
		r = _r;
		g = _g;
		b = _b;
	}

	void operator +=(const Colour &other) {
		r += other.r;
		g += other.g;
		b += other.b;
	}

	void operator /=(const float s) {
		r /= s;
		g /= s;
		b /= s;
	}

	float &operator [] (int i) {
		return d[i];
	}

	float operator [] (int i) const {
		return d[i];
	}

	struct
	{
		float r, g, b;
	};
private:
	float d[3];
};

inline
Colour operator *(const Colour &a, const Colour &b) {
	return Colour(a.r * b.r, a.g * b.g, a.b * b.b);
}

inline
Colour operator *(const Colour &c, const float s) {
	return Colour(c.r * s, c.g * s, c.b * s);
}

inline
Colour operator *(const float s, const Colour &c) {
	return c * s;
}

inline
Colour operator +(const Colour &a, const Colour &b) {
	return Colour(a.r + b.r, a.g + b.g, a.b + b.b);
}

inline
Colour operator /(const Colour &c, const float s) {
	return Colour(c.r / s, c.g / s, c.b / s);
}

inline
Colour operator -(const Colour &a, const Colour &b) {
	return Colour(a.r - b.r, a.g - b.g, a.b - b.b);
}

#endif