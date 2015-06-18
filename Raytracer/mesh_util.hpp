#ifndef _MESH_UTIL_HPP_
#define _MESH_UTIL_HPP_

#include "algebra.hpp"
#include "primitive_util.hpp"

typedef std::vector<int> Face;

struct Mesh : Primitive
{
	Mesh() {}
	Mesh(const std::vector<Point3> &v, const std::vector<Face> &f)
		: vertices(v), faces(f) {
	}

	std::vector<Point3> vertices;
	std::vector<Face> faces;
};

// TODO(orglofch): Possibly remove
inline
void TriangulateMesh(Mesh &mesh) { // TODO(orglofch): Possibly optimize instead of connecting last vertex to all other vertices
	std::vector<Face> new_faces;
	for (Face &face : mesh.faces) {
		int vertices = face.size();
		if (vertices > 3) {
			int last_vertex = face.at(vertices - 1);
			for (int i = 0; i < vertices - 2; ++i) {
				Face triangle;

				triangle.push_back(face[i]);
				triangle.push_back(face[i + 1]);
				triangle.push_back(last_vertex);

				new_faces.push_back(triangle);
			}
		}
	}
	mesh.faces.clear();
	mesh.faces.swap(new_faces);
}

// TODO(orglofch): Possibly remove
inline
bool MeshIntersect(const Mesh &mesh, const Ray &ray, Intersection *intersection) {
	intersection->t = std::numeric_limits<double>::max();
	bool has_intersection = false;
	for (const Face &tri : mesh.faces) {
		assert(tri.size() == 3);
		Point3 vertex_0 = mesh.vertices[tri[0]];

		Vector3 edge_1 = mesh.vertices[tri[1]] - vertex_0;
		Vector3 edge_2 = mesh.vertices[tri[2]] - vertex_0;

		Vector3 P = ray.dir.cross(edge_2);
		double det = edge_1.dot(P);
		if (det > -EPSILON && det < EPSILON) {
			continue;
		}
		double inv_det = 1.0f / det;

		Vector3 T = ray.origin - vertex_0;

		double u = T.dot(P) * inv_det;
		if (u < 0.0f || u > 1.0f) {
			continue;
		}

		Vector3 Q = T.cross(edge_1);

		double v = ray.dir.dot(Q) * inv_det;
		if (v < 0.0f || u + v > 1.0f) {
			continue;
		}

		double t = edge_2.dot(Q) * inv_det;
		if (t > 0 && t < intersection->t) {
			Point3 point = RayProjection(ray, t);
			Vector3 normal = edge_1.cross(edge_2);
			intersection->pos = point;
			intersection->normal = normal.dot(ray.dir) < 0 ? normal : -1 * normal;
			intersection->normal.normalize();
			intersection->t = t;
			intersection->material = mesh.material;
			has_intersection = true;
		}
	}
	return has_intersection;
}

inline
void LoadOBJ(const std::string &filename, Mesh *mesh) {
	std::vector<Point3> vertices;
	std::vector<Face> faces;

	std::ifstream ifs(filename);

	std::string lineHeader;
	while (ifs >> lineHeader) {
		if (lineHeader.compare("v") == 0) {
			Point3 vertex;
			ifs >> vertex[0] >> vertex[1] >> vertex[2];
			vertices.push_back(vertex);
		} else if (lineHeader.compare("vt") == 0) {
			// TODO(orglofch): Texture coordinate
		} else if (lineHeader.compare("vn") == 0) {
			// TODO(orglofch): Normal coordinate
		} else if (lineHeader.compare("f") == 0) {
			Face face;
			char line[256];
			ifs.getline(line, 256);
			std::istringstream iss(line);
			int vertex;
			while (iss) {
				iss >> vertex; // TODO(orglofch): Check if we have to read delimeter
				face.push_back(vertex - 1);
			}

			faces.push_back(face);
		}
	}
	mesh->vertices = vertices;
	mesh->faces = faces;
}

#endif