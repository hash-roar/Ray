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

#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Direct illumination integrator using distribution ray tracing
 */
class WhittedIntegrator : public Integrator {
public:
    WhittedIntegrator(const PropertyList &props) {
        /* No parameters needed for this integrator */
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const override {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        Color3f result(0.0f);
        
        // Add direct emission if we hit an emitter
        if (its.mesh->isEmitter()) {
            Vector3f wo = -ray.d.normalized();
            result += its.mesh->getEmitter()->eval(its.p, wo);
        }

        // Sample direct illumination from all emitters
        const auto &emitters = scene->getEmitters();
        if (!emitters.empty()) {
            // Pick an emitter uniformly at random
            uint32_t emitterIndex = std::min((uint32_t)(sampler->next1D() * emitters.size()), 
                                           (uint32_t)emitters.size() - 1);
            const Emitter *emitter = emitters[emitterIndex];
            float emitterPdf = 1.0f / emitters.size();
            
            // Sample a point on the chosen emitter
            Point2f sample = sampler->next2D();
            Point3f lightPos;
            Vector3f lightNormal;
            float lightPdf;
            
            Color3f Le = emitter->sample(its.p, sample, lightPos, lightNormal, lightPdf);
            
            if (lightPdf > 0.0f) {
                // Compute direction from surface to light
                Vector3f lightDir = lightPos - its.p;
                float lightDistance = lightDir.norm();
                lightDir /= lightDistance; // normalize
                
                // Check visibility
                Ray3f shadowRay(its.p + its.shFrame.n * Epsilon, lightDir, 
                               Epsilon, lightDistance - Epsilon);
                bool visible = !scene->rayIntersect(shadowRay);
                
                if (visible) {
                    // Evaluate BSDF
                    BSDFQueryRecord bRec(its.toLocal(-ray.d), its.toLocal(lightDir), ESolidAngle);
                    Color3f fr = its.mesh->getBSDF()->eval(bRec);
                    
                    // Compute geometric term G(x ↔ y)
                    float cosTheta_x = std::max(0.0f, its.shFrame.n.dot(lightDir));
                    float cosTheta_y = std::max(0.0f, lightNormal.dot(-lightDir));
                    float G = visible ? (cosTheta_x * cosTheta_y) / (lightDistance * lightDistance) : 0.0f;
                    
                    // Add contribution: fr * G * Le / (lightPdf * emitterPdf)
                    result += fr * G * Le / (lightPdf * emitterPdf);
                }
            }
        }

        return result;
    }

    std::string toString() const override {
        return "WhittedIntegrator[]";
    }
};

NORI_REGISTER_CLASS(WhittedIntegrator, "whitted");

NORI_NAMESPACE_END
