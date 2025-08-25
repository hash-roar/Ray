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
 * \brief Material sampling path tracer
 * 
 * This integrator implements a simple path tracer that relies on hitting
 * light sources by chance instead of explicit emitter sampling.
 */
class PathMatsIntegrator : public Integrator {
public:
    PathMatsIntegrator(const PropertyList &props) {
        /* No parameters needed for this integrator */
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const override {
        Color3f radiance(0.0f);
        Color3f throughput(1.0f);
        Ray3f currentRay(ray);
        float eta = 1.0f;
        
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
                radiance += throughput * its.mesh->getEmitter()->eval(its.p, wo);
            }

            /* Get the BSDF at the intersection point */
            const BSDF *bsdf = its.mesh->getBSDF();

            /* Sample the BSDF to get a new direction */
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
            currentRay = Ray3f(its.p + offset, nextDir);
        }

        return radiance;
    }

    std::string toString() const override {
        return "PathMatsIntegrator[]";
    }
};

NORI_REGISTER_CLASS(PathMatsIntegrator, "path_mats");

NORI_NAMESPACE_END
