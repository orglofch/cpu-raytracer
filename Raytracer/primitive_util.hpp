#ifndef _PRIMITIVE_UTIL_HPP_
#define _PRIMITIVE_UTIL_HPP_

#include <iostream>
#include <limits>
#include <vector>

#include "algebra.hpp"
#include "colour.hpp"

struct Material
{
	Material() : shininess(0) {}

	Colour diffuse;
	Colour specular;
	float shininess;
};

struct Intersection
{
	Point3 pos;
	Vector3 normal;
	double t;

	Material material;
};

struct Ray
{
	Ray(const Point3 &p, const Vector3 &d) : origin(p), dir(d) {}

	Point3 origin;
	Vector3 dir;
};

Point3 RayProjection(const Ray &ray, float t) {
	return ray.origin + t * (ray.dir);
}

struct Plane
{
	Point3 point;
	Vector3 normal;
};

inline
bool PlaneIntersect(const Plane &plane, const Ray &ray, Intersection *intersection) {
	float denom = plane.normal.dot(ray.dir);
	if (denom > EPSILON) {
		intersection->t = plane.normal.dot(plane.point - ray.origin) / denom;
		if (intersection->t > EPSILON) {
			return true;
		}
	}
	return false;
}

struct Primitive
{
	Material material;

	Matrix4x4 transform; // TODO(orglofch): Make implicitly chan / remove
	Matrix4x4 inv_transform;
};

struct Sphere : Primitive
{
	Sphere() : radius(0) {}
	Sphere(const Point3 &p, double rad) : pos(p), radius(rad) {}

	Point3 pos;
	double radius;
};

inline
bool SphereIntersect(const Sphere &sphere, const Ray &ray, Intersection *intersection) {
	double A = ray.dir.dot(ray.dir);
	Vector3 origin = ray.origin - sphere.pos;
	double B = 2 * (ray.dir.dot(origin));
	double C = origin.dot(origin) - sphere.radius * sphere.radius;

	float disc = B*B - 4 * A * C;
	if (disc < 0) {
		return false;
	}

	double roots[2];
	unsigned int num_roots = quadraticRoots(A, B, C, roots);

	if (num_roots == 1) {
		intersection->pos = RayProjection(ray, roots[0]);
		intersection->normal = intersection->pos - sphere.pos;
		intersection->normal.normalize();
		intersection->t = roots[0];
		intersection->material = sphere.material;
	} else if (num_roots == 2) {
		if (roots[0] > roots[1]) {
			double tmp = roots[0];
			roots[0] = roots[1];
			roots[1] = tmp;
		}
		if (roots[1] < EPSILON) {
			return false;
		}
		if (roots[0] < 0 && roots[1] > EPSILON) {
			intersection->pos = RayProjection(ray, roots[1]);
			intersection->normal = intersection->pos - sphere.pos;
			intersection->normal.normalize();
			intersection->t = roots[1];
			intersection->material = sphere.material;
			return true;
		} else if (roots[0] > EPSILON) {
			intersection->pos = RayProjection(ray, roots[0]);
			intersection->normal = intersection->pos - sphere.pos;
			intersection->normal.normalize();
			intersection->t = roots[0];
			intersection->material = sphere.material;
			return true;
		}
	}
	return false;
}

struct Polygon
{
	std::vector<Point3> vertices;
};

// TODO(orglofch): Possibly remove in favour of directly loading into polygons
/*inline
void PolygonsForMesh(const Mesh &mesh, std::vector<Polygon> *polygons) {
	for (const Face &face : mesh.faces) {
		Polygon polygon;

		for (int vertex : face) {
			polygon.vertices.push_back(mesh.transform * mesh.vertices[vertex]); // TODO(orglofch): Expensive, possibly remove
		}

		polygons->push_back(polygon);
	}
}*/

struct Triangle
{
	Point3 vertices[3];
};

inline
void TriangulatePolygons(const std::vector<Polygon> &polygons, std::vector<Triangle> *triangles) { // TODO(orglofch): Possibly optimize instead of connecting last vertex to all other vertices
	for (const Polygon &polygon : polygons) {
		int vertex_count = polygons.size();
		Point3 last_vertex = polygon.vertices[vertex_count - 1];
		for (int i = 0; i < vertex_count - 2; ++i) {
			Triangle triangle;

			triangle.vertices[0] = polygon.vertices[i];
			triangle.vertices[1] = polygon.vertices[i + 1];
			triangle.vertices[2] = last_vertex;

			triangles->push_back(triangle);
		}
	}
}

inline
bool TriangleIntersect(const Triangle &triangle, const Ray &ray, Intersection *intersection) {
	Point3 vertex_0 = triangle.vertices[0];

	Vector3 edge_1 = triangle.vertices[1] - vertex_0;
	Vector3 edge_2 = triangle.vertices[2] - vertex_0;

	Vector3 P = ray.dir.cross(edge_2);
	double det = edge_1.dot(P);
	if (det > -EPSILON && det < EPSILON) {
		return false;
	}
	double inv_det = 1.0f / det;

	Vector3 T = ray.origin - vertex_0;

	double u = T.dot(P) * inv_det;
	if (u < 0.0f || u > 1.0f) {
		return false;
	}

	Vector3 Q = T.cross(edge_1);

	double v = ray.dir.dot(Q) * inv_det;
	if (v < 0.0f || u + v > 1.0f) {
		return false;
	}

	double t = edge_2.dot(Q) * inv_det;
	if (t > EPSILON && t < intersection->t) {
		Point3 point = RayProjection(ray, t);
		Vector3 normal = edge_1.cross(edge_2);
		intersection->pos = point;
		intersection->normal = normal.dot(ray.dir) < 0 ? normal : -1 * normal;
		intersection->normal.normalize();
		intersection->t = t;
		// TODO(orglofch): intersection->material = mesh.material;
		return true;
	}
	return false;
}

#endif
