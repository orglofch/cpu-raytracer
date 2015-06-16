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

	float &operator [] (int i) {
		return d[i];
	}

	float operator [] (int i) const {
		return d[i];
	}

	void operator +=(const Colour &other) {
		r += other.r;
		g += other.g;
		b += other.b;
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
Colour operator +(const Colour &a, const Colour &b) {
	return Colour(a.r + b.r, a.g + b.g, a.b + b.b); // TODO(orglofch): Possibly clamp
}

#endif