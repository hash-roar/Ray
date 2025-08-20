#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

class AmbientOcclusionIntegrator : public Integrator {
public:
    AmbientOcclusionIntegrator(const PropertyList &props) {
        // No parameters needed for basic AO
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        // Get intersection point and create coordinate frame
        Point3f x = its.p;
        Frame frame(its.shFrame.n);
        
        // Sample a direction on the hemisphere with cosine weighting
        Point2f sample = sampler->next2D();
        Vector3f localDir = Warp::squareToCosineHemisphere(sample);
        
        // Transform from local coordinates to world coordinates
        Vector3f worldDir = frame.toWorld(localDir);
        
        // Create shadow ray in the sampled direction
        Ray3f shadowRay(x + its.shFrame.n * Epsilon, worldDir, Epsilon, std::numeric_limits<float>::infinity());
        
        // Check for occlusion
        bool occluded = scene->rayIntersect(shadowRay);
        
        // If not occluded, contribute to the lighting
        // The cosine term is already included in the cosine hemisphere sampling
        // so we just need to return the visibility
        if (!occluded) {
            return Color3f(1.0f);
        } else {
            return Color3f(0.0f);
        }
    }

    std::string toString() const {
        return "AmbientOcclusionIntegrator[]";
    }
};

NORI_REGISTER_CLASS(AmbientOcclusionIntegrator, "ao");
NORI_NAMESPACE_END
