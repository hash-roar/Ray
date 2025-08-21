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

NORI_NAMESPACE_BEGIN

/// Ideal dielectric BSDF
class Dielectric : public BSDF {
public:
    Dielectric(const PropertyList &propList) {
        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);
    }

    Color3f eval(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return Color3f(0.0f);
    }

    float pdf(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return 0.0f;
    }

    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {
        float cosThetaI = Frame::cosTheta(bRec.wi);
        bool entering = cosThetaI > 0.0f;
        
        /* Determine the relative IOR and relevant indices */
        float etaI = entering ? m_extIOR : m_intIOR;
        float etaT = entering ? m_intIOR : m_extIOR;
        float eta = etaI / etaT;
        
        /* Make sure cosThetaI is positive for Fresnel computation */
        if (!entering) {
            cosThetaI = -cosThetaI;
        }
        
        /* Compute the Fresnel reflectance */
        float Fr = fresnel(cosThetaI, m_extIOR, m_intIOR);
        
        /* Compute sine squared of transmitted angle using Snell's law */
        float sin2ThetaT = eta * eta * (1.0f - cosThetaI * cosThetaI);
        
        /* Check for total internal reflection */
        bool totalInternalReflection = (sin2ThetaT >= 1.0f);
        
        if (sample.x() <= Fr || totalInternalReflection) {
            /* Sample reflection (either normal reflection or total internal reflection) */
            bRec.wo = Vector3f(
                -bRec.wi.x(),
                -bRec.wi.y(),
                bRec.wi.z()
            );
            bRec.measure = EDiscrete;
            bRec.eta = 1.0f;  /* No change in medium for reflection */
            
            return Color3f(1.0f);
        } else {
            /* Sample refraction */
            float cosThetaT = std::sqrt(1.0f - sin2ThetaT);
            
            /* Compute refracted direction */
            Vector3f wt = eta * (-bRec.wi) + (eta * cosThetaI - cosThetaT) * Vector3f(0, 0, 1);
            
            /* Handle the case when ray is coming from inside the material */
            if (!entering) {
                wt.z() = -wt.z();
            }
            
            bRec.wo = wt.normalized();
            bRec.measure = EDiscrete;
            bRec.eta = eta;
            
            /* The transmitted radiance is scaled by eta^2 due to solid angle compression */
            return Color3f(eta * eta);
        }
    }

    std::string toString() const {
        return tfm::format(
            "Dielectric[\n"
            "  intIOR = %f,\n"
            "  extIOR = %f\n"
            "]",
            m_intIOR, m_extIOR);
    }
private:
    float m_intIOR, m_extIOR;
};

NORI_REGISTER_CLASS(Dielectric, "dielectric");
NORI_NAMESPACE_END
