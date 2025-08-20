#include "nori/common.h"
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

class SimpleIntegrator : public Integrator {
public:
    SimpleIntegrator(const PropertyList &props) {
        m_position = props.getPoint("position");
        m_energy = props.getColor("energy");
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        // Get intersection point
        Point3f x = its.p;
        
        // Vector from intersection point to light
        Vector3f lightDir = (m_position - x).normalized();
        float lightDistance = (m_position - x).norm();
        
        // Angle between light direction and surface normal
        float cosTheta = std::max(0.0f, its.shFrame.n.dot(lightDir));
        
        // Check visibility with shadow ray
        Ray3f shadowRay(x + its.shFrame.n * Epsilon, lightDir, Epsilon, lightDistance - Epsilon);
        bool visible = !scene->rayIntersect(shadowRay);
        
        if (!visible) {
            return Color3f(0.0f);
        }
        
        // Compute lighting: Φ / (4π * ||x-p||²) * max(0, cos θ) * V(x↔p)
        float falloff = 1.0f / (4.0f * M_PI * M_PI * lightDistance * lightDistance);
        Color3f result = m_energy * falloff * cosTheta;
        
        return result;
    }

    std::string toString() const {
        return tfm::format("SimpleIntegrator[\n"
                          "  position = %s,\n"
                          "  energy = %s\n"
                          "]", m_position.toString(), m_energy.toString());
    }

private:
    Point3f m_position;
    Color3f m_energy;
};

NORI_REGISTER_CLASS(SimpleIntegrator, "simple");
NORI_NAMESPACE_END
