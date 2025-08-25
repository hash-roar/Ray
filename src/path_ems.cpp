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
 * \brief Emitter sampling path tracer
 * 
 * This integrator implements a path tracer with explicit emitter sampling
 * to reduce variance compared to the naive material sampling approach.
 */
class PathEmsIntegrator : public Integrator {
public:
    PathEmsIntegrator(const PropertyList &props) {
        /* No parameters needed for this integrator */
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const override {
        Color3f radiance(0.0f);
        Color3f throughput(1.0f);
        Ray3f currentRay(ray);
        float eta = 1.0f;
        bool lastBounceSpecular = false;
        
        /* Iterative path tracing loop */
        for (int bounces = 0; ; ++bounces) {
            /* Find the surface that is visible in the requested direction */
            Intersection its;
            if (!scene->rayIntersect(currentRay, its)) {
                /* No intersection found - ray escaped to infinity */
                break;
            }

            /* Add emission if we hit an emitter */
            if (its.mesh->isEmitter()) {
                Vector3f wo = -currentRay.d.normalized();
                if (bounces == 0 || lastBounceSpecular) {
                    /* Direct hit from camera or after specular bounce - count emission */
                    radiance += throughput * its.mesh->getEmitter()->eval(its.p, wo);
                }
                /* For non-specular indirect hits, we skip adding emission to avoid double counting */
            }

            /* Get the BSDF at the intersection point */
            const BSDF *bsdf = its.mesh->getBSDF();

            /* Direct illumination via emitter sampling (only for diffuse materials) */
            const auto &emitters = scene->getEmitters();
            if (!emitters.empty() && bsdf->isDiffuse()) {
                /* Pick an emitter uniformly at random */
                uint32_t emitterIndex = std::min((uint32_t)(sampler->next1D() * emitters.size()), 
                                               (uint32_t)emitters.size() - 1);
                const Emitter *emitter = emitters[emitterIndex];
                float emitterPdf = 1.0f / emitters.size();
                
                /* Sample a point on the chosen emitter */
                Point2f sample = sampler->next2D();
                Point3f lightPos;
                Vector3f lightNormal;
                float lightPdf;
                
                Color3f Le = emitter->sample(its.p, sample, lightPos, lightNormal, lightPdf);
                
                if (lightPdf > 0.0f && !Le.isZero()) {
                    /* Compute direction from surface to light */
                    Vector3f lightDir = lightPos - its.p;
                    float lightDistance = lightDir.norm();
                    lightDir /= lightDistance; // normalize
                    
                    /* Check visibility */
                    Ray3f shadowRay(its.p + its.shFrame.n * Epsilon, lightDir, 
                                   Epsilon, lightDistance - Epsilon);
                    bool visible = !scene->rayIntersect(shadowRay);
                    
                    if (visible) {
                        /* Evaluate BSDF */
                        BSDFQueryRecord bRec(its.toLocal(-currentRay.d), its.toLocal(lightDir), ESolidAngle);
                        Color3f fr = bsdf->eval(bRec);
                        
                        if (!fr.isZero()) {
                            /* Compute geometric term */
                            float cosTheta_x = std::max(0.0f, its.shFrame.n.dot(lightDir));
                            float cosTheta_y = std::max(0.0f, lightNormal.dot(-lightDir));
                            float G = (cosTheta_x * cosTheta_y) / (lightDistance * lightDistance);
                            
                            /* Add direct illumination contribution */
                            radiance += throughput * fr * G * Le / (lightPdf * emitterPdf);
                        }
                    }
                }
            }

            /* Sample the BSDF to get a new direction for indirect illumination */
            BSDFQueryRecord bRec(its.toLocal(-currentRay.d));
            Point2f sample = sampler->next2D();
            Color3f f = bsdf->sample(bRec, sample);
            
            /* Check if the sampling was successful */
            if (f.isZero()) {
                break;
            }

            /* Update throughput */
            throughput *= f;
            
            /* Update eta (relative index of refraction) */
            eta *= bRec.eta;

            /* Track if this bounce was specular */
            lastBounceSpecular = !bsdf->isDiffuse();

            /* Russian Roulette termination (after at least 3 bounces) */
            if (bounces >= 3) {
                float maxComponent = throughput.maxCoeff();
                float continuationProb = std::min(maxComponent * eta * eta, 0.99f);
                
                if (sampler->next1D() > continuationProb) {
                    break;
                }
                
                /* Account for the continuation probability */
                throughput /= continuationProb;
            }

            /* Create the next ray */
            Vector3f nextDir = its.toWorld(bRec.wo);
            
            /* Choose appropriate surface offset */
            Vector3f offset = its.shFrame.n * Epsilon;
            if (nextDir.dot(its.shFrame.n) < 0) {
                /* Ray is going into the surface */
                offset = -offset;
            }
            /* Update current ray for next iteration */
            currentRay.o = its.p + offset;
            currentRay.d = nextDir;
            currentRay.mint = Epsilon;
            currentRay.maxt = INFINITY;
        }

        return radiance;
    }

    std::string toString() const override {
        return "PathEmsIntegrator[]";
    }
};

NORI_REGISTER_CLASS(PathEmsIntegrator, "path_ems");

NORI_NAMESPACE_END
