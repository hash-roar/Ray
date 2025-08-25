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

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class Microfacet : public BSDF {
public:
    Microfacet(const PropertyList &propList) {
        /* RMS surface roughness */
        m_alpha = propList.getFloat("alpha", 0.1f);

        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);

        /* Albedo of the diffuse base material (a.k.a "kd") */
        m_kd = propList.getColor("kd", Color3f(0.5f));

        /* To ensure energy conservation, we must scale the 
           specular component by 1-kd. 

           While that is not a particularly realistic model of what 
           happens in reality, this will greatly simplify the 
           implementation. Please see the course staff if you're 
           interested in implementing a more realistic version 
           of this BRDF. */
        m_ks = 1 - m_kd.maxCoeff();
    }

    /// Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord &bRec) const {
        /* This is a BRDF, so no transmission is expected */
        if (Frame::cosTheta(bRec.wi) <= 0 || Frame::cosTheta(bRec.wo) <= 0)
            return Color3f(0.0f);

        /* Calculate the half-angle vector */
        Vector3f wh = (bRec.wi + bRec.wo).normalized();
        
        /* Compute the diffuse component */
        Color3f diffuse = m_kd / M_PI;
        
        /* Compute the specular component */
        float cosThetaI = Frame::cosTheta(bRec.wi);
        float cosThetaO = Frame::cosTheta(bRec.wo);
        float cosThetaH = Frame::cosTheta(wh);
        
        /* For microfacet BRDF, we need to be more lenient with grazing angles */
        /* Only reject if the configuration is clearly invalid */
        if (cosThetaH <= 0.0f) 
            return diffuse;
            
        /* Fresnel coefficient */
        float F = fresnel(wh.dot(bRec.wi), m_extIOR, m_intIOR);
        
        /* Beckmann microfacet distribution D(wh) */
        float tan2ThetaH = (1.0f - cosThetaH * cosThetaH) / (cosThetaH * cosThetaH);
        float D = std::exp(-tan2ThetaH / (m_alpha * m_alpha)) / 
                  (M_PI * m_alpha * m_alpha * std::pow(cosThetaH, 4.0f));

        /* Geometric term G using rational function approximation */
        auto G1 = [this](const Vector3f& w, const Vector3f& wh) -> float {
            float cosTheta = Frame::cosTheta(w);
            if (cosTheta <= 1e-6f) return 0.0f;
            
            float wh_dot_w = w.dot(wh);
            float chi_plus = (wh_dot_w > 0.0f) ? 1.0f : 0.0f;
            if (chi_plus == 0.0f) return 0.0f;
            
            float tanTheta = Frame::tanTheta(w);
            if (std::isinf(tanTheta) || tanTheta <= 0.0f) return 0.0f;
            
            float b = 1.0f / (m_alpha * tanTheta);
            if (b >= 1.6f) {
                return chi_plus;
            } else {
                return chi_plus * (3.535f * b + 2.181f * b * b) / 
                       (1.0f + 2.276f * b + 2.577f * b * b);
            }
        };
        
        float G = G1(bRec.wi, wh) * G1(bRec.wo, wh);
        
        /* Final specular term */
        Color3f specular = Color3f(m_ks * F * D * G / (4.0f * cosThetaI * cosThetaO));
        
        return diffuse + specular;
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
        /* This is a BRDF, so no transmission is expected */
        if (Frame::cosTheta(bRec.wi) <= 0 || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;

        /* Calculate the half-angle vector */
        Vector3f wh = (bRec.wi + bRec.wo).normalized();
        
        /* Probability density: weighted combination of diffuse and specular sampling */
        
        /* Diffuse component: cosine-weighted hemisphere sampling */
        float diffusePdf = Frame::cosTheta(bRec.wo) / M_PI;
        
        /* Specular component: importance sampling using Beckmann distribution */
        float beckmannPdf = Warp::squareToBeckmannPdf(wh, m_alpha);
        /* Convert from half-angle density to outgoing direction density using Jacobian */
        float jacobian = 1.0f / (4.0f * std::abs(wh.dot(bRec.wo)));
        float specularPdf = beckmannPdf * jacobian;
        
        /* Weighted combination based on ks */
        return (1.0f - m_ks) * diffusePdf + m_ks * specularPdf;
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const {
        if (Frame::cosTheta(bRec.wi) <= 0) 
            return Color3f(0.0f);

        bRec.measure = ESolidAngle;
        
        /* Make a copy of the sample for reuse */
        Point2f sample = _sample;
        
        /* Decide between diffuse and specular reflection based on ks */
        if (sample.x() < m_ks) {
            /* Sample specular component */
            /* Scale the sample to [0,1) for reuse */
            sample.x() /= m_ks;
            
            /* Sample the Beckmann distribution to get a microfacet normal */
            Vector3f wh = Warp::squareToBeckmann(sample, m_alpha);
            
            /* Reflect the incident direction about the microfacet normal */
            float wi_dot_wh = bRec.wi.dot(wh);
            if (wi_dot_wh <= 0.0f) {
                return Color3f(0.0f);
            }
            bRec.wo = 2.0f * wi_dot_wh * wh - bRec.wi;
            
            /* Check if the reflected direction is in the upper hemisphere */
            if (Frame::cosTheta(bRec.wo) <= 0.0f) {
                return Color3f(0.0f);
            }
        } else {
            /* Sample diffuse component */
            /* Scale the sample to [0,1) for reuse */
            sample.x() = (sample.x() - m_ks) / (1.0f - m_ks);
            
            /* Generate cosine-weighted direction on hemisphere */
            bRec.wo = Warp::squareToCosineHemisphere(sample);
        }
        
        /* Relative index of refraction: no change */
        bRec.eta = 1.0f;

        // Note: Once you have implemented the part that computes the scattered
        // direction, the last part of this function should simply return the
        // BRDF value divided by the solid angle density and multiplied by the
        // cosine factor from the reflection equation, i.e.
        return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
    }

    bool isDiffuse() const {
        /* While microfacet BRDFs are not perfectly diffuse, they can be
           handled by sampling techniques for diffuse/non-specular materials,
           hence we return true here */
        return false;
    }

    std::string toString() const {
        return tfm::format(
            "Microfacet[\n"
            "  alpha = %f,\n"
            "  intIOR = %f,\n"
            "  extIOR = %f,\n"
            "  kd = %s,\n"
            "  ks = %f\n"
            "]",
            m_alpha,
            m_intIOR,
            m_extIOR,
            m_kd.toString(),
            m_ks
        );
    }
private:
    float m_alpha;
    float m_intIOR, m_extIOR;
    float m_ks;
    Color3f m_kd;
};

NORI_REGISTER_CLASS(Microfacet, "microfacet");
NORI_NAMESPACE_END
