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

#include <GL/glew.h>
#include <GLFW\glfw3.h>

#include <boost\thread.hpp>
#include <cstdlib>
#include <ctime>
#include <limits>

#include "algebra.hpp"
#include "colour.hpp"
#include "image_util.hpp"
#include "gl_util.hpp"
#include "mesh_util.hpp"
#include "primitive_util.hpp"
#include "util.hpp"

using namespace std;

struct Light
{
	Point3 pos;
	Colour colour;
	Vector3 falloff;
};

struct SimulationState
{
	Point3 eye;
	Vector3 view;
	Vector3 up;

	float fov;

	Matrix4x4 pixel_to_world_matrix;

	Light light[2];
	Colour ambient;

	Sphere sphere[5];

	int max_aa_iterations;
	int max_reflection_iterations;

	int trace_thread_count;

	Texture image_texture;
};

void error_callback(int error, const char *description) {
	fputs(description, stderr);
}

void key_callback(GLFWwindow* window, int key, int scan_code, int action, int mods) {
	SimulationState *state = (SimulationState*)glfwGetWindowUserPointer(window);
	switch (key)
	{
		case GLFW_KEY_ESCAPE:
			glfwSetWindowShouldClose(window, GL_TRUE);
			break;
	}
}

void cursor_position_callback(GLFWwindow *window, double x_pos, double y_pos) {
	SimulationState *state = (SimulationState*)glfwGetWindowUserPointer(window);
}

void mouse_button_callback(GLFWwindow *window, int button, int action, int mods) {
	SimulationState *state = (SimulationState*)glfwGetWindowUserPointer(window);
}

void mouse_scroll_callback(GLFWwindow *window, double x_offset, double y_offset) {
	SimulationState *state = (SimulationState*)glfwGetWindowUserPointer(window);
}

bool FindIntersection(SimulationState &state, const Ray &ray, Intersection *intersection) {
	intersection->t = numeric_limits<double>::max();
	bool has_intersection = false;
	Intersection cur_intersection;
	for (int i = 0; i < 5; ++i) {
		if (SphereIntersect(state.sphere[i], ray, &cur_intersection) && cur_intersection.t < intersection->t) {
			*intersection = cur_intersection;
			has_intersection = true;
		}
	}
	return has_intersection;
}

void CreatePixelToWorldMatrix(const Point3 &eye,
							  const Vector3 &view,
							  const Vector3 &up,
							  const int fov,
							  const int width,
							  const int height,
							  Matrix4x4 *mat) {
	double d = view.length();
	double h = 2.0 * d * tan(toRad(fov) / 2.0); // TODO(orglofch): Test
	Matrix4x4 t1 = translation(Vector3(-width / 2.0, -height / 2.0, d));
	Matrix4x4 s2 = scaling(Vector3(-h / height, -h / height, 1.0));
	Matrix4x4 r3 = rotation(eye, view, up);
	Matrix4x4 t4 = translation(eye - Point3(0, 0, 0));
	*mat = t4 * r3 * s2 * t1;
}

Ray ReflectRay(const Ray &ray, const Point3 &pos, const Vector3 &normal) {
	Vector3 refl_dir = ray.dir - 2.0 * ray.dir.dot(normal) * normal;
	return Ray(pos, refl_dir);
}

Colour DirectLight(SimulationState &state, const Intersection &intersection, const Light &light) {
	Vector3 light_dir = light.pos - intersection.pos;
	double light_dist = light_dir.length();

	Ray light_ray(intersection.pos, light_dir);
	Intersection light_intersection;

	if (!FindIntersection(state, light_ray, &light_intersection) ||
		(light_intersection.pos - light.pos).length() >=
		(intersection.pos - light.pos).length()) {
		double attenuation = 1 / (light.falloff.x +
			light.falloff.y * light_dist +
			light.falloff.z * light_dist * light_dist);

		return light.colour * attenuation;
	} else {
		return Colour(0.0, 0.0, 0.0);
	}
}

bool IsBlack(const Colour &c) {
	return c.r == 0 && c.g == 0 && c.b == 0;
}

void CastRay(SimulationState &state, const Ray &ray, int iteration, Colour *colour);

void ColourForIntersection(SimulationState &state, const Intersection &intersection, const Ray &ray, int iteration, Colour *colour) {
	Colour col = state.ambient * intersection.material.diffuse;

	for (int i = 0; i < 2; ++i) {
		Vector3 light_dir = state.light[i].pos - intersection.pos;
		light_dir.normalize();

		Vector3 bisector = (state.eye - intersection.pos) + light_dir;
		bisector.normalize();

		double spec_angle = max(intersection.normal.dot(bisector), 0.0);
		double specular = pow(spec_angle, intersection.material.shininess);
		double lambertian = max(intersection.normal.dot(light_dir), 0.0);

		Colour light_colour = DirectLight(state, intersection, state.light[i]);

		if (lambertian > 0.0) {
			col += (lambertian * intersection.material.diffuse +
				specular * intersection.material.specular) * light_colour;
		} else {
			col += (lambertian * intersection.material.diffuse) * light_colour;
		}
	}
	if (!IsBlack(intersection.material.specular) && iteration <= state.max_reflection_iterations) {
		Ray refl_ray = ReflectRay(ray, intersection.pos, intersection.normal);
		Colour refl_colour;
		CastRay(state, refl_ray, iteration + 1, &refl_colour);
		col += intersection.material.specular * refl_colour;
	}
	*colour = col;
}

void CastRay(SimulationState &state, const Ray &ray, int iteration, Colour *colour) {
	Intersection intersection;
	if (FindIntersection(state, ray, &intersection)) {
		ColourForIntersection(state, intersection, ray, iteration, colour);
	} else {
		colour->r = 0.0f;
		colour->g = 0.0f;
		colour->b = 0.0f;
	}
}

Colour TraceArea(SimulationState &state, const Point3 &center, double size, int iteration) {
	Point3 pixel1 = center + Vector3(-size / 2, -size / 2, 0);
	Point3 pixel2 = center + Vector3(size / 2, -size / 2, 0);
	Point3 pixel3 = center + Vector3(-size / 2, size / 2, 0);
	Point3 pixel4 = center + Vector3(size / 2, size / 2, 0);
	Point3 pixel5 = center;

	Point3 p_world1 = state.pixel_to_world_matrix * pixel1;
	Point3 p_world2 = state.pixel_to_world_matrix * pixel2;
	Point3 p_world3 = state.pixel_to_world_matrix * pixel3;
	Point3 p_world4 = state.pixel_to_world_matrix * pixel4;
	Point3 p_world5 = state.pixel_to_world_matrix * pixel5;

	Vector3 dir1 = (p_world1 - state.eye);
	Vector3 dir2 = (p_world2 - state.eye);
	Vector3 dir3 = (p_world3 - state.eye);
	Vector3 dir4 = (p_world4 - state.eye);
	Vector3 dir5 = (p_world5 - state.eye);

	Ray ray1(state.eye, dir1);
	Ray ray2(state.eye, dir2);
	Ray ray3(state.eye, dir3);
	Ray ray4(state.eye, dir4);
	Ray ray5(state.eye, dir5);

	Colour colour1;
	CastRay(state, ray1, 1, &colour1);
	Colour colour2;
	CastRay(state, ray2, 1, &colour2);
	Colour colour3;
	CastRay(state, ray3, 1, &colour3);
	Colour colour4;
	CastRay(state, ray4, 1, &colour4);
	Colour colour5;
	CastRay(state, ray5, 1, &colour5);

	Colour diff_col1 = colour1 - colour5;
	Colour diff_col2 = colour2 - colour5;
	Colour diff_col3 = colour3 - colour5;
	Colour diff_col4 = colour4 - colour5;

	float diff1 = abs(diff_col1.r) + abs(diff_col1.g) + abs(diff_col1.b);
	float diff2 = abs(diff_col2.r) + abs(diff_col2.g) + abs(diff_col2.b);
	float diff3 = abs(diff_col3.r) + abs(diff_col3.g) + abs(diff_col3.b);
	float diff4 = abs(diff_col4.r) + abs(diff_col4.g) + abs(diff_col4.b);

	if (diff1 > 0.1 && iteration < state.max_aa_iterations) {
		colour1 = TraceArea(state, center + Vector3(-size / 4, -size / 4, 0), size / 2, iteration + 1);
	}
	if (diff2 > 0.1 && iteration < state.max_aa_iterations) {
		colour2 = TraceArea(state, center + Vector3(size / 4, -size / 4, 0), size / 2, iteration + 1);
	}
	if (diff3 > 0.1 && iteration < state.max_aa_iterations) {
		colour3 = TraceArea(state, center + Vector3(-size / 4, size / 4, 0), size / 2, iteration + 1);
	}
	if (diff4 > 0.1 && iteration < state.max_aa_iterations) {
		colour4 = TraceArea(state, center + Vector3(size / 4, size / 4, 0), size / 2, iteration + 1);
	}

	return (colour1 + colour2 + colour3 + colour4 + colour5) / 5;
}

void TracePixel(SimulationState &state, const Point3 &pixel, double size, Image *image) {
	int index = DataIndex(image->width, image->height, image->channels, pixel.x, pixel.y, 0);

	Colour colour = TraceArea(state, pixel, 0.5, 1);
	image->data[index + 0] = colour.r;
	image->data[index + 1] = colour.g;
	image->data[index + 2] = colour.b;
}

void TraceRect(SimulationState &state, const vector<Point3> &pixel_pool, Image *image) {
	for (const Point3 &p : pixel_pool) {
		TracePixel(state, p, 1, image); 
	}
}

void Trace(SimulationState &state, Image *image) {
	vector<vector<Point3> > thread_pixel_pools(state.trace_thread_count);
	for (unsigned int c = 0; c < image->width; ++c) {
		for (unsigned int r = 0; r < image->height; ++r) {
			int thread = rand() % state.trace_thread_count;
			// TODO(orglofch): Possibly remove 0.5 from pixel
			thread_pixel_pools[thread].push_back(Point3(c + 0.5, r + 0.5, 0.0));
		}
	}
	boost::thread_group trace_threads;
	for (int i = 0; i < state.trace_thread_count; ++i) {
		random_shuffle(thread_pixel_pools[i].begin(), thread_pixel_pools[i].end());
		boost::thread *worker_thread = new boost::thread(TraceRect, state, thread_pixel_pools[i], image);
		trace_threads.add_thread(worker_thread);
	}
	trace_threads.join_all();
}

void UpdateAndRender(SimulationState &state, Image &image, double dt) {
	glClear(GL_COLOR_BUFFER_BIT);

	glBindTexture(GL_TEXTURE_2D, state.image_texture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, image.width, image.height, 0, GL_RGB,
		GL_FLOAT, image.data);

	glBindTexture(GL_TEXTURE_2D, state.image_texture);
	glDrawRect(0, image.width, 0, image.height, -1);
}

int main(int argc, char **argv) {
	srand(time(NULL));

	string output_file = "output1.png";

	SimulationState state;
	state.eye = { 0.0, 0.0, 800.0 };
	state.view = { 0.0, 0.0, -1.0 };
	state.up = { 0.0, 1.0, 0.0 };
	state.fov = 50;
	state.ambient = { 0.3f, 0.3f, 0.3f };

	state.trace_thread_count = 8;

	state.max_aa_iterations = 1;
	state.max_reflection_iterations = 3;

	state.light[0].pos = { -100, 150, 400 };
	state.light[0].colour = { 0.9f, 0.9f, 0.9f };
	state.light[0].falloff = { 1, 0, 0 };
	state.light[1].pos = { 400, 100, 150 };
	state.light[1].colour = { 0.7f, 0.0f, 0.7f };
	state.light[1].falloff = { 1, 0, 0 };

	Material mat1;
	mat1.diffuse = { 0.7f, 1.0f, 0.7f };
	mat1.specular = { 0.5f, 0.7f, 0.5f };
	mat1.shininess = 25;

	Material mat2;
	mat2.diffuse = { 0.5f, 0.5f, 0.5f };
	mat2.specular = { 0.5f, 0.7f, 0.5f };
	mat2.shininess = 25;

	Material mat3;
	mat3.diffuse = { 1.0f, 0.6f, 0.1f };
	mat3.specular = { 0.5f, 0.7f, 0.5f };
	mat3.shininess = 25;

	state.sphere[0].pos = { 0, 0, -400 };
	state.sphere[0].radius = 100;
	state.sphere[0].material = mat1;
	state.sphere[1].pos = { 200, 50, -100 };
	state.sphere[1].radius = 150;
	state.sphere[1].material = mat1;
	state.sphere[2].pos = { 0, -1200, -500 };
	state.sphere[2].radius = 1000;
	state.sphere[2].material = mat2;
	state.sphere[3].pos = { -100, 25, -300 };
	state.sphere[3].radius = 50;
	state.sphere[3].material = mat3;
	state.sphere[4].pos = { 0, 100, -250 };
	state.sphere[4].radius = 25;
	state.sphere[4].material = mat1;

	Image image;
	image.width = 1024;
	image.height = 1024;
	image.channels = 3;

	// TODO(orglofch): Make implicit
	CreatePixelToWorldMatrix(state.eye, state.view, state.up, state.fov, 
		image.width, image.height, &state.pixel_to_world_matrix);

	GLFWwindow *window;
	if (!glfwInit()) {
		LOG("Failed to initialize glfw\n");
		exit(EXIT_FAILURE);
	}

	glfwSetErrorCallback(error_callback);

	window = glfwCreateWindow(image.width, image.height, "RayTracer", NULL, NULL);
	if (!window) {
		LOG("Failed to create window\n");
		glfwTerminate();
		exit(EXIT_FAILURE);
	}

	glfwSetWindowUserPointer(window, &state);

	glfwMakeContextCurrent(window);
	glfwSetKeyCallback(window, key_callback);
	glfwSetCursorPosCallback(window, cursor_position_callback);
	glfwSetMouseButtonCallback(window, mouse_button_callback);
	glfwSetScrollCallback(window, mouse_scroll_callback);
	glfwSwapInterval(1);

	if (glewInit() != GLEW_OK) {
		LOG("Failed to initialize glew\n");
		glfwTerminate();
		exit(EXIT_FAILURE);
	}

	glViewport(0, 0, image.width, image.height);
	glSetOrthographicProjection(0, image.width, 0, image.height, 0.1, 1000);

	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	glDisable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);
	glDisable(GL_BLEND);

	image.data = new float[image.width * image.height * 3]; // TODO(orglofch): Make implicit
	memset(image.data, 0, sizeof(float) * image.width * image.height * 3);
	glCreateTexture2D(&state.image_texture, image.width, image.height, 3, image.data);

	boost::thread worker_thread(Trace, state, &image);

	double time = glfwGetTime();
	while (!glfwWindowShouldClose(window)) {
		glfwPollEvents();

		UpdateAndRender(state, image, glfwGetTime() - time);

		glfwSwapBuffers(window);
	}

	worker_thread.join();

	WritePNG(image, output_file);

	glfwDestroyWindow(window);
	glfwTerminate();

	glDeleteTextures(1, &state.image_texture);

	exit(EXIT_SUCCESS);
}