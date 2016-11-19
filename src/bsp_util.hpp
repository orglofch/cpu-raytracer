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

#ifndef _BSP_UTIL_HPP_
#define _BSP_UTIL_HPP_

#include <vector>

#include "primitive_util.hpp"

struct BSPTree
{
	Plane partition;

	std::vector<Polygon*> polygons;

	BSPTree *front, *back;
};

enum PartitionClassification
{
	PARTITION_CLASS_FRONT,
	PARTITION_CLASS_BACK,
	PARTITION_CLASS_DIVIDE,
	PARTITION_CLASS_INCIDENT,
};

inline
PartitionClassification ClassifyPolygon(const Polygon &polygon, const Plane &partition) {
	return PARTITION_CLASS_FRONT;
}

inline
void DividePolygon(const Polygon &polygon, Plane &plane, Polygon *front, Polygon *back) {
	// TODO(orglofch):
}

// TODO(orglofch): Possibly make this take a mesh and avoid the intermediate transfer to polygons
inline
void CreateBSP(const std::vector<Polygon> &polygons, BSPTree *tree) {
	const Polygon &p = polygons[0];
	tree->partition.point = p->vertices[0];
	tree->partition.normal = (p->vertices[1] - p->vertices[0]).cross(p->vertices[2] - p->vertices[0]);
	tree->polygons.push_back(p);
	tree->front = NULL;
	tree->back = NULL;

	std::vector<Polygon> front_list;
	std::vector<Polygon> back_list;

	for (const Polygon &polygon : polygons) {
		switch (ClassifyPolygon(polygon, tree->partition))
		{
			case PARTITION_CLASS_FRONT:
			{
				front_list.push_back(polygon);
				break;
			}
			case PARTITION_CLASS_BACK:
			{
				back_list.push_back(polygon);
				break;
			}
			case PARTITION_CLASS_DIVIDE:
			{
				Polygon front_piece, back_piece;
				DividePolygon(polygon, tree->partition, &front_piece, &back_piece);
				back_list.push_back(back_piece);
				front_list.push_back(front_piece);
				break;
			}
			case PARTITION_CLASS_INCIDENT:
			{
				tree->polygons.push_back(polygon);
				break;
			}
		}
	}
	
	if (!front_list.empty()) {
		tree->front = new BSPTree();
		CreateBSP(front_list, tree->front);
	}
	if (!back_list.empty()) {
		tree->back = new BSPTree();
		CreateBSP(back_list, tree->back);
	}
}

inline
void TriangulateBSPTree(BSPTree &tree) {
	// TODO(orglofch):
}

#endif