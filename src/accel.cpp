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

#include <nori/accel.h>
#include <Eigen/Geometry>
#include <algorithm>
#include <vector>

NORI_NAMESPACE_BEGIN

void Accel::addMesh(Mesh *mesh) {
    if (m_mesh)
        throw NoriException("Accel: only a single mesh is supported!");
    m_mesh = mesh;
    m_bbox = m_mesh->getBoundingBox();
}

void Accel::build() {
    if (!m_mesh) return;
    
    // Create list of all triangle indices
    std::vector<uint32_t> triangles;
    triangles.reserve(m_mesh->getTriangleCount());
    for (uint32_t i = 0; i < m_mesh->getTriangleCount(); ++i) {
        triangles.push_back(i);
    }
    
    // Build the octree
    m_octree = buildOctree(m_bbox, triangles);
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {
    if (!m_mesh || !m_octree) return false;
    
    Ray3f ray(ray_); /// Make a copy of the ray
    
    // Check if ray intersects the root bounding box
    if (!m_bbox.rayIntersect(ray)) {
        return false;
    }
    
    // Use octree for intersection
    return rayIntersectOctree(m_octree, m_bbox, ray, its, shadowRay);
}

OctreeNode* Accel::buildOctree(const BoundingBox3f& bbox, const std::vector<uint32_t>& triangles, int depth) {
    const int MAX_DEPTH = 10;  // Limit depth to avoid pathological cases
    const int MIN_TRIANGLES = 10;  // Create leaf if fewer than this many triangles
    
    // Create leaf node if stopping criteria are met
    if (triangles.empty()) {
        return nullptr;
    }
    
    if (triangles.size() <= MIN_TRIANGLES || depth >= MAX_DEPTH) {
        OctreeNode* leaf = new OctreeNode();
        leaf->isLeaf = true;
        leaf->triangles = triangles;
        return leaf;
    }
    
    // Create internal node
    OctreeNode* node = new OctreeNode();
    node->isLeaf = false;
    
    // Create 8 triangle lists for the 8 octants
    std::vector<uint32_t> childTriangles[8];
    
    // Classify triangles into octants
    for (uint32_t triIdx : triangles) {
        BoundingBox3f triBBox = m_mesh->getBoundingBox(triIdx);
        
        // Check which children this triangle overlaps
        for (int i = 0; i < 8; ++i) {
            BoundingBox3f childBBox = getChildBBox(bbox, i);
            if (triBBox.overlaps(childBBox)) {
                childTriangles[i].push_back(triIdx);
            }
        }
    }
    
    // Recursively build children
    for (int i = 0; i < 8; ++i) {
        if (!childTriangles[i].empty()) {
            BoundingBox3f childBBox = getChildBBox(bbox, i);
            node->children[i] = buildOctree(childBBox, childTriangles[i], depth + 1);
        }
    }
    
    return node;
}

BoundingBox3f Accel::getChildBBox(const BoundingBox3f& parentBBox, int childIndex) const {
    Point3f center = parentBBox.getCenter();
    Point3f min = parentBBox.min;
    Point3f max = parentBBox.max;
    
    // Child bounding box corners based on octant index
    Point3f childMin, childMax;
    
    childMin.x() = (childIndex & 1) ? center.x() : min.x();
    childMax.x() = (childIndex & 1) ? max.x() : center.x();
    
    childMin.y() = (childIndex & 2) ? center.y() : min.y();
    childMax.y() = (childIndex & 2) ? max.y() : center.y();
    
    childMin.z() = (childIndex & 4) ? center.z() : min.z();
    childMax.z() = (childIndex & 4) ? max.z() : center.z();
    
    return BoundingBox3f(childMin, childMax);
}

bool Accel::rayIntersectOctree(const OctreeNode* node, const BoundingBox3f& bbox, 
                              const Ray3f& ray, Intersection& its, bool shadowRay) const {
    if (!node) return false;
    
    // Check if ray intersects this bounding box
    if (!bbox.rayIntersect(ray)) {
        return false;
    }
    
    if (node->isLeaf) {
        // Leaf node: test all triangles
        bool foundIntersection = false;
        uint32_t closestTriangle = (uint32_t) -1;
        Ray3f rayLocal(ray);  // Local copy for updating maxt
        
        for (uint32_t triIdx : node->triangles) {
            float u, v, t;
            if (m_mesh->rayIntersect(triIdx, rayLocal, u, v, t)) {
                if (shadowRay) {
                    return true;  // Early exit for shadow rays
                }
                
                rayLocal.maxt = its.t = t;
                its.uv = Point2f(u, v);
                its.mesh = m_mesh;
                closestTriangle = triIdx;
                foundIntersection = true;
            }
        }
        
        if (foundIntersection) {
            // Fill in intersection details
            Vector3f bary;
            bary << 1-its.uv.sum(), its.uv;

            const MatrixXf &V  = m_mesh->getVertexPositions();
            const MatrixXf &N  = m_mesh->getVertexNormals();
            const MatrixXf &UV = m_mesh->getVertexTexCoords();
            const MatrixXu &F  = m_mesh->getIndices();

            uint32_t idx0 = F(0, closestTriangle), idx1 = F(1, closestTriangle), idx2 = F(2, closestTriangle);
            Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);

            its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;

            if (UV.size() > 0)
                its.uv = bary.x() * UV.col(idx0) + bary.y() * UV.col(idx1) + bary.z() * UV.col(idx2);

            its.geoFrame = Frame((p1-p0).cross(p2-p0).normalized());

            if (N.size() > 0) {
                its.shFrame = Frame(
                    (bary.x() * N.col(idx0) +
                     bary.y() * N.col(idx1) +
                     bary.z() * N.col(idx2)).normalized());
            } else {
                its.shFrame = its.geoFrame;
            }
        }
        
        return foundIntersection;
    } else {
        // Internal node: recursively test children in sorted order
        bool foundIntersection = false;
        Ray3f rayLocal(ray);
        
        // Structure to hold child information with distances
        struct ChildInfo {
            int index;
            float tNear;
            BoundingBox3f bbox;
            
            bool operator<(const ChildInfo& other) const {
                return tNear < other.tNear;
            }
        };
        
        std::vector<ChildInfo> childrenInfo;
        childrenInfo.reserve(8);
        
        // Collect all valid children with their intersection distances
        for (int i = 0; i < 8; ++i) {
            if (node->children[i]) {
                BoundingBox3f childBBox = getChildBBox(bbox, i);
                float tNear, tFar;
                
                if (childBBox.rayIntersect(ray, tNear, tFar)) {
                    ChildInfo info;
                    info.index = i;
                    info.tNear = tNear;
                    info.bbox = childBBox;
                    childrenInfo.push_back(info);
                }
            }
        }
        
        // Sort children by distance (closest first)
        std::sort(childrenInfo.begin(), childrenInfo.end());
        
        // Traverse children in sorted order with early termination
        for (const auto& childInfo : childrenInfo) {
            // Early termination: if we already found an intersection closer than
            // this child's bounding box, skip it
            if (foundIntersection && its.t < childInfo.tNear) {
                break;
            }
            
            Intersection childIts;
            if (rayIntersectOctree(node->children[childInfo.index], childInfo.bbox, rayLocal, childIts, shadowRay)) {
                if (shadowRay) {
                    return true;  // Early exit for shadow rays
                }
                
                if (!foundIntersection || childIts.t < its.t) {
                    its = childIts;
                    rayLocal.maxt = its.t;  // Update ray for closer intersections
                    foundIntersection = true;
                }
            }
        }
        
        return foundIntersection;
    }
}

NORI_NAMESPACE_END
