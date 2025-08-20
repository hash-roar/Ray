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
#include <vector>
#include <tbb/parallel_invoke.h>
#include <tbb/task_group.h>

NORI_NAMESPACE_BEGIN

bool Accel::s_useParallelConstruction = true;  // Default to parallel construction

void Accel::addMesh(Mesh *mesh) {
    m_meshes.push_back(mesh);
    
    // Update bounding box to include this mesh
    if (m_meshes.size() == 1) {
        m_bbox = mesh->getBoundingBox();
    } else {
        m_bbox.expandBy(mesh->getBoundingBox());
    }
}

void Accel::clear() {
    m_meshes.clear();
    delete m_octree;
    m_octree = nullptr;
    m_bbox = BoundingBox3f();
}

void Accel::build() {
    if (m_meshes.empty()) return;
    
    // Create list of all triangle references
    std::vector<TriangleRef> triangles;
    
    // Count total triangles to reserve space
    size_t totalTriangles = 0;
    for (size_t meshIdx = 0; meshIdx < m_meshes.size(); ++meshIdx) {
        totalTriangles += m_meshes[meshIdx]->getTriangleCount();
    }
    triangles.reserve(totalTriangles);
    
    // Add all triangles from all meshes
    for (size_t meshIdx = 0; meshIdx < m_meshes.size(); ++meshIdx) {
        Mesh* mesh = m_meshes[meshIdx];
        for (uint32_t triIdx = 0; triIdx < mesh->getTriangleCount(); ++triIdx) {
            triangles.emplace_back(static_cast<uint32_t>(meshIdx), triIdx);
        }
    }
    
    // Build the octree (choose parallel or serial based on flag)
    if (s_useParallelConstruction) {
        m_octree = buildOctreeParallel(m_bbox, triangles);
    } else {
        m_octree = buildOctree(m_bbox, triangles);
    }
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {
    if (m_meshes.empty() || !m_octree) return false;
    
    Ray3f ray(ray_); /// Make a copy of the ray
    
    // Check if ray intersects the root bounding box
    if (!m_bbox.rayIntersect(ray)) {
        return false;
    }
    
    // Use octree for intersection
    return rayIntersectOctree(m_octree, m_bbox, ray, its, shadowRay);
}

OctreeNode* Accel::buildOctree(const BoundingBox3f& bbox, const std::vector<TriangleRef>& triangles, int depth) {
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
    OctreeNode* node = new OctreeNode(bbox);
    node->isLeaf = false;
    
    // Create 8 triangle lists for the 8 octants
    std::vector<TriangleRef> childTriangles[8];
    
    // Classify triangles into octants using pre-computed bounding boxes
    for (const TriangleRef& triRef : triangles) {
        Mesh* mesh = m_meshes[triRef.meshIndex];
        BoundingBox3f triBBox = mesh->getBoundingBox(triRef.triangleIndex);
        
        // Check which children this triangle overlaps
        for (int i = 0; i < 8; ++i) {
            if (triBBox.overlaps(node->childBBoxes[i])) {
                childTriangles[i].push_back(triRef);
            }
        }
    }
    
    // Recursively build children
    for (int i = 0; i < 8; ++i) {
        if (!childTriangles[i].empty()) {
            node->children[i] = buildOctree(node->childBBoxes[i], childTriangles[i], depth + 1);
        }
    }
    
    return node;
}

OctreeNode* Accel::buildOctreeParallel(const BoundingBox3f& bbox, const std::vector<TriangleRef>& triangles, int depth) {
    const int MAX_DEPTH = 10;  // Limit depth to avoid pathological cases
    const int MIN_TRIANGLES = 10;  // Create leaf if fewer than this many triangles
    const int PARALLEL_THRESHOLD = 1000;  // Use parallel construction for subtrees with many triangles
    
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
    OctreeNode* node = new OctreeNode(bbox);
    node->isLeaf = false;
    
    // Create 8 triangle lists for the 8 octants
    std::vector<TriangleRef> childTriangles[8];
    
    // Classify triangles into octants using pre-computed bounding boxes
    for (const TriangleRef& triRef : triangles) {
        Mesh* mesh = m_meshes[triRef.meshIndex];
        BoundingBox3f triBBox = mesh->getBoundingBox(triRef.triangleIndex);
        
        // Check which children this triangle overlaps
        for (int i = 0; i < 8; ++i) {
            if (triBBox.overlaps(node->childBBoxes[i])) {
                childTriangles[i].push_back(triRef);
            }
        }
    }
    
    // Decide whether to build children in parallel or serial
    bool useParallel = (triangles.size() > PARALLEL_THRESHOLD && depth < 3);
    
    if (useParallel) {
        // Parallel construction using TBB task_group
        tbb::task_group tg;
        
        for (int i = 0; i < 8; ++i) {
            if (!childTriangles[i].empty()) {
                tg.run([=, &node]() {
                    node->children[i] = buildOctreeSerial(node->childBBoxes[i], childTriangles[i], depth + 1);
                });
            }
        }
        
        tg.wait();  // Wait for all child constructions to complete
    } else {
        // Serial construction for smaller subtrees
        for (int i = 0; i < 8; ++i) {
            if (!childTriangles[i].empty()) {
                node->children[i] = buildOctreeSerial(node->childBBoxes[i], childTriangles[i], depth + 1);
            }
        }
    }
    
    return node;
}

OctreeNode* Accel::buildOctreeSerial(const BoundingBox3f& bbox, const std::vector<TriangleRef>& triangles, int depth) {
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
    OctreeNode* node = new OctreeNode(bbox);
    node->isLeaf = false;
    
    // Create 8 triangle lists for the 8 octants
    std::vector<TriangleRef> childTriangles[8];
    
    // Classify triangles into octants using pre-computed bounding boxes
    for (const TriangleRef& triRef : triangles) {
        Mesh* mesh = m_meshes[triRef.meshIndex];
        BoundingBox3f triBBox = mesh->getBoundingBox(triRef.triangleIndex);
        
        // Check which children this triangle overlaps
        for (int i = 0; i < 8; ++i) {
            if (triBBox.overlaps(node->childBBoxes[i])) {
                childTriangles[i].push_back(triRef);
            }
        }
    }
    
    // Recursively build children (serial)
    for (int i = 0; i < 8; ++i) {
        if (!childTriangles[i].empty()) {
            node->children[i] = buildOctreeSerial(node->childBBoxes[i], childTriangles[i], depth + 1);
        }
    }
    
    return node;
}

bool Accel::rayIntersectOctree(const OctreeNode* node, const BoundingBox3f& bbox, 
                              const Ray3f& ray, Intersection& its, bool shadowRay) const {
    if (!node) return false;
    
    // Fast ray-box intersection test with early exit
    float tNear, tFar;
    if (!bbox.rayIntersect(ray, tNear, tFar)) {
        return false;
    }
    
    // Early termination: if we already have a closer intersection, skip this node
    if (tNear > ray.maxt) {
        return false;
    }
    
    if (node->isLeaf) {
        // Leaf node: test all triangles with optimized loop
        bool foundIntersection = false;
        uint32_t closestMeshIndex = (uint32_t) -1;
        uint32_t closestTriangle = (uint32_t) -1;
        float closestT = ray.maxt;
        
        const TriangleRef* triangles = node->triangles.data();
        const size_t numTriangles = node->triangles.size();
        
        for (size_t j = 0; j < numTriangles; ++j) {
            const TriangleRef& triRef = triangles[j];
            Mesh* mesh = m_meshes[triRef.meshIndex];
            float u, v, t;
            if (mesh->rayIntersect(triRef.triangleIndex, ray, u, v, t)) {
                if (shadowRay) {
                    return true;  // Early exit for shadow rays
                }
                
                if (t < closestT) {
                    closestT = t;
                    its.t = t;
                    its.uv = Point2f(u, v);
                    its.mesh = mesh;
                    closestMeshIndex = triRef.meshIndex;
                    closestTriangle = triRef.triangleIndex;
                    foundIntersection = true;
                }
            }
        }
        
        if (foundIntersection) {
            // Fill in intersection details
            Vector3f bary;
            bary << 1-its.uv.sum(), its.uv;

            Mesh* mesh = m_meshes[closestMeshIndex];
            const MatrixXf &V  = mesh->getVertexPositions();
            const MatrixXf &N  = mesh->getVertexNormals();
            const MatrixXf &UV = mesh->getVertexTexCoords();
            const MatrixXu &F  = mesh->getIndices();

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
        
        // Fixed-size arrays to avoid dynamic allocation
        struct ChildInfo {
            int index;
            float tNear;
        };
        
        ChildInfo childrenInfo[8];
        int validChildCount = 0;
        
        // Collect all valid children with their intersection distances
        for (int i = 0; i < 8; ++i) {
            if (node->children[i]) {
                float childTNear, childTFar;
                
                if (node->childBBoxes[i].rayIntersect(ray, childTNear, childTFar)) {
                    childrenInfo[validChildCount].index = i;
                    childrenInfo[validChildCount].tNear = childTNear;
                    validChildCount++;
                }
            }
        }
        
        // Simple insertion sort for small arrays (much faster than std::sort for <= 8 elements)
        for (int i = 1; i < validChildCount; ++i) {
            ChildInfo key = childrenInfo[i];
            int j = i - 1;
            while (j >= 0 && childrenInfo[j].tNear > key.tNear) {
                childrenInfo[j + 1] = childrenInfo[j];
                j--;
            }
            childrenInfo[j + 1] = key;
        }
        
        // Traverse children in sorted order with early termination
        for (int i = 0; i < validChildCount; ++i) {
            const ChildInfo& childInfo = childrenInfo[i];
            
            // Early termination: if we already found an intersection closer than
            // this child's bounding box, skip it
            if (foundIntersection && its.t < childInfo.tNear) {
                break;
            }
            
            Intersection childIts;
            if (rayIntersectOctree(node->children[childInfo.index], node->childBBoxes[childInfo.index], rayLocal, childIts, shadowRay)) {
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
