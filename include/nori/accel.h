/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include <nori/mesh.h>
#include <vector>

NORI_NAMESPACE_BEGIN

/**
 * \brief Octree node for acceleration data structure
 */
struct OctreeNode {
    /// Constructor for leaf node
    OctreeNode() : isLeaf(true) {
        for (int i = 0; i < 8; ++i)
            children[i] = nullptr;
    }
    
    /// Destructor
    ~OctreeNode() {
        if (!isLeaf) {
            for (int i = 0; i < 8; ++i)
                delete children[i];
        }
    }
    
    bool isLeaf;                      ///< True if this is a leaf node
    std::vector<uint32_t> triangles;  ///< Triangle indices (only for leaf nodes)
    OctreeNode* children[8];          ///< Child nodes (only for internal nodes)
};

/**
 * \brief Acceleration data structure for ray intersection queries
 *
 * Uses an octree to organize triangle data for efficient ray intersection.
 */
class Accel {
public:
    /// Destructor to clean up the octree
    ~Accel() { delete m_octree; }
    /**
     * \brief Register a triangle mesh for inclusion in the acceleration
     * data structure
     *
     * This function can only be used before \ref build() is called
     */
    void addMesh(Mesh *mesh);

    /// Build the acceleration data structure (octree)
    void build();

    /// Return an axis-aligned box that bounds the scene
    const BoundingBox3f &getBoundingBox() const { return m_bbox; }

    /**
     * \brief Intersect a ray against all triangles stored in the scene and
     * return detailed intersection information
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum extent
     *    information
     *
     * \param its
     *    A detailed intersection record, which will be filled by the
     *    intersection query
     *
     * \param shadowRay
     *    \c true if this is a shadow ray query, i.e. a query that only aims to
     *    find out whether the ray is blocked or not without returning detailed
     *    intersection information.
     *
     * \return \c true if an intersection was found
     */
    bool rayIntersect(const Ray3f &ray, Intersection &its, bool shadowRay) const;

private:
    Mesh         *m_mesh = nullptr; ///< Mesh (only a single one for now)
    BoundingBox3f m_bbox;           ///< Bounding box of the entire scene
    OctreeNode   *m_octree = nullptr; ///< Root of the octree
    
    /// Recursively build the octree
    OctreeNode* buildOctree(const BoundingBox3f& bbox, const std::vector<uint32_t>& triangles, int depth = 0);
    
    /// Recursively intersect ray with octree
    bool rayIntersectOctree(const OctreeNode* node, const BoundingBox3f& bbox, 
                           const Ray3f& ray, Intersection& its, bool shadowRay) const;
    
    /// Get the bounding box of a child octant
    BoundingBox3f getChildBBox(const BoundingBox3f& parentBBox, int childIndex) const;
};

NORI_NAMESPACE_END
