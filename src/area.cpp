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

#include <nori/emitter.h>
#include <nori/mesh.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Area light source emitter
 */
class AreaLight : public Emitter {
public:
    AreaLight(const PropertyList &props) {
        m_radiance = props.getColor("radiance", Color3f(1.0f));
    }

    /// Sample a point on the emitter surface
    Color3f sample(const Point3f &ref, const Point2f &sample, 
                   Point3f &p, Vector3f &n, float &pdf) const override {
        if (!m_mesh) {
            throw NoriException("AreaLight: No mesh associated with this emitter!");
        }
        
        // Sample a point on the mesh surface
        m_mesh->sampleSurface(sample, p, n, pdf);
        
        return m_radiance;
    }

    /// Evaluate the radiance emitted from point p in direction d
    Color3f eval(const Point3f &p, const Vector3f &d) const override {
        // Area lights emit uniformly in all directions above the surface
        // For now, we'll use a simple check - in practice, this should be more sophisticated
        // and should get the actual surface normal at point p
        return m_radiance;
    }

    /// Compute the probability density of sampling point p from reference point ref
    float pdf(const Point3f &ref, const Point3f &p) const override {
        if (!m_mesh) {
            return 0.0f;
        }
        
        // PDF is reciprocal of total surface area
        return 1.0f / m_mesh->getSurfaceArea();
    }

    void setParent(NoriObject *parent) override {
        if (parent->getClassType() == EMesh) {
            m_mesh = static_cast<const Mesh *>(parent);
        }
    }

    std::string toString() const override {
        return tfm::format(
            "AreaLight[\n"
            "  radiance = %s\n"
            "]",
            m_radiance.toString()
        );
    }

private:
    Color3f m_radiance;      ///< Radiance emitted by the area light
    const Mesh *m_mesh = nullptr;  ///< Associated mesh
};

NORI_REGISTER_CLASS(AreaLight, "area");

NORI_NAMESPACE_END
