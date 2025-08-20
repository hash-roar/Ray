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

#pragma once

#include <nori/common.h>
#include <nori/sampler.h>
#include <nori/bitmap.h>

NORI_NAMESPACE_BEGIN

/// Mipmap pyramid for hierarchical sample warping
class HierarchicalSampler {
public:
    /// Constructor that builds the mipmap from an image file
    HierarchicalSampler(const std::string &filename);

    /// Sample a point according to the luminance distribution
    Point2f sample(const Point2f &sample) const;

    /// Evaluate the probability density at a given point
    float pdf(const Point2f &p) const;

private:
    std::vector<std::vector<float>> m_pyramid;  ///< Mipmap pyramid
    std::vector<Vector2i> m_sizes;              ///< Sizes at each level
    int m_levels;                               ///< Number of mipmap levels
    float m_normalization;                      ///< Normalization factor

    /// Build the mipmap pyramid from luminance values
    void buildPyramid(const std::vector<float> &luminance, int width, int height);

    /// Get luminance from a color value
    float getLuminance(const Color3f &color) const;
};

/// A collection of useful warping functions for importance sampling
class Warp {
public:
    /// Dummy warping function: takes uniformly distributed points in a square and just returns them
    static Point2f squareToUniformSquare(const Point2f &sample);

    /// Probability density of \ref squareToUniformSquare()
    static float squareToUniformSquarePdf(const Point2f &p);

    /// Sample a 2D tent distribution
    static Point2f squareToTent(const Point2f &sample);

    /// Probability density of \ref squareToTent()
    static float squareToTentPdf(const Point2f &p);

    /// Uniformly sample a vector on a 2D disk with radius 1, centered around the origin
    static Point2f squareToUniformDisk(const Point2f &sample);

    /// Probability density of \ref squareToUniformDisk()
    static float squareToUniformDiskPdf(const Point2f &p);

    /// Uniformly sample a vector on the unit sphere with respect to solid angles
    static Vector3f squareToUniformSphere(const Point2f &sample);

    /// Probability density of \ref squareToUniformSphere()
    static float squareToUniformSpherePdf(const Vector3f &v);

    /// Uniformly sample a vector on the unit hemisphere around the pole (0,0,1) with respect to solid angles
    static Vector3f squareToUniformHemisphere(const Point2f &sample);

    /// Probability density of \ref squareToUniformHemisphere()
    static float squareToUniformHemispherePdf(const Vector3f &v);

    /// Uniformly sample a vector on the unit hemisphere around the pole (0,0,1) with respect to projected solid angles
    static Vector3f squareToCosineHemisphere(const Point2f &sample);

    /// Probability density of \ref squareToCosineHemisphere()
    static float squareToCosineHemispherePdf(const Vector3f &v);

    /// Warp a uniformly distributed square sample to a Beckmann distribution * cosine for the given 'alpha' parameter
    static Vector3f squareToBeckmann(const Point2f &sample, float alpha);

    /// Probability density of \ref squareToBeckmann()
    static float squareToBeckmannPdf(const Vector3f &m, float alpha);

    /// Hierarchical sample warping using mipmapping for importance sampling
    static Point2f squareToHierarchical(const Point2f &sample, const std::string &filename);

    /// Probability density of \ref squareToHierarchical()
    static float squareToHierarchicalPdf(const Point2f &p, const std::string &filename);
};

NORI_NAMESPACE_END
