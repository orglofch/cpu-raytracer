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

#ifndef _UTIL_HPP_
#define _UTIL_HPP_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

#ifndef _DEBUG
#define _DEBUG 0
#endif

#define __FILENAME__ (strrchr(__FILE__, '\\') ? strrchr(__FILE__, '\\') + 1 : __FILE__)
#define LOG(fmt, ...) do { if (_DEBUG) { fprintf(stderr, "%s:%d: " fmt, __FILENAME__, __LINE__, __VA_ARGS__); } } while(0)

#define NOT_REACHED() assert(false)
#define NOT_NULL(data) assert(data)

inline
std::string ReadFile(const char *filename) {
	std::ifstream file(filename);
	if (!file.is_open()) {
		LOG("Unable to open %s\n", filename);
		return "";
	}

	std::stringstream file_data;
	file_data << file.rdbuf();
	file.close();

	return file_data.str();
}

inline
void clamp(float *val, float min, float max) {
	*val = std::min(max, std::max(min, *val));
}

inline
void clamp(int *val, int min, int max) {
	*val = std::min(max, std::max(min, *val));
}

inline
float Randf(float min, float max) {
	return (max - min) * rand() / RAND_MAX + min;
}

#endif