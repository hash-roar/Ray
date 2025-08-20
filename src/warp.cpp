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

#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>
#include <map>
#include <memory>

NORI_NAMESPACE_BEGIN

Point2f Warp::squareToUniformSquare(const Point2f &sample) {
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

Point2f Warp::squareToTent(const Point2f &sample) {
    // For tent distribution: p(t) = 1 - |t| for t in [-1, 1]
    // CDF: P(t) = (1 + t)^2 / 2 for t in [-1, 0]
    //      P(t) = 1 - (1 - t)^2 / 2 for t in [0, 1]
    // Inverse CDF for each dimension
    auto tentInverse = [](float xi) -> float {
        if (xi < 0.5f) {
            return std::sqrt(2.0f * xi) - 1.0f;
        } else {
            return 1.0f - std::sqrt(2.0f * (1.0f - xi));
        }
    };
    
    return Point2f(tentInverse(sample.x()), tentInverse(sample.y()));
}

float Warp::squareToTentPdf(const Point2f &p) {
    // Tent distribution: p(x,y) = p1(x) * p1(y)
    // where p1(t) = 1 - |t| for t in [-1, 1], 0 otherwise
    auto tentPdf1D = [](float t) -> float {
        if (t >= -1.0f && t <= 1.0f) {
            return 1.0f - std::abs(t);
        }
        return 0.0f;
    };
    
    return tentPdf1D(p.x()) * tentPdf1D(p.y());
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) {
    // Use polar coordinates: r = sqrt(xi1), theta = 2*pi*xi2
    float r = std::sqrt(sample.x());
    float theta = 2.0f * M_PI * sample.y();
    return Point2f(r * std::cos(theta), r * std::sin(theta));
}

float Warp::squareToUniformDiskPdf(const Point2f &p) {
    // Uniform distribution on unit disk has PDF = 1/pi for points inside
    if (p.x() * p.x() + p.y() * p.y() <= 1.0f) {
        return INV_PI;
    }
    return 0.0f;
}

Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
    // Uniform sampling on sphere using spherical coordinates
    // cos(theta) = 1 - 2*xi1 (uniform in cos(theta))
    // phi = 2*pi*xi2
    float cosTheta = 1.0f - 2.0f * sample.x();
    float sinTheta = std::sqrt(std::max(0.0f, 1.0f - cosTheta * cosTheta));
    float phi = 2.0f * M_PI * sample.y();
    
    return Vector3f(
        sinTheta * std::cos(phi),
        sinTheta * std::sin(phi),
        cosTheta
    );
}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {
    // Uniform distribution on unit sphere has PDF = 1/(4*pi)
    if (std::abs(v.norm() - 1.0f) < Epsilon) {
        return INV_FOURPI;
    }
    return 0.0f;
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
    // Uniform sampling on hemisphere: cos(theta) = xi1
    // phi = 2*pi*xi2
    float cosTheta = sample.x();
    float sinTheta = std::sqrt(std::max(0.0f, 1.0f - cosTheta * cosTheta));
    float phi = 2.0f * M_PI * sample.y();
    
    return Vector3f(
        sinTheta * std::cos(phi),
        sinTheta * std::sin(phi),
        cosTheta
    );
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
    // Uniform distribution on unit hemisphere has PDF = 1/(2*pi)
    // Check if vector is on unit hemisphere (z >= 0)
    if (v.z() >= 0.0f && std::abs(v.norm() - 1.0f) < Epsilon) {
        return INV_TWOPI;
    }
    return 0.0f;
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
    // Cosine-weighted hemisphere sampling
    // Use Malley's method: sample unit disk then project to hemisphere
    Point2f diskSample = squareToUniformDisk(sample);
    float r2 = diskSample.x() * diskSample.x() + diskSample.y() * diskSample.y();
    float z = std::sqrt(std::max(0.0f, 1.0f - r2));
    
    return Vector3f(diskSample.x(), diskSample.y(), z);
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
    // Cosine-weighted distribution: PDF = cos(theta)/pi
    // Check if vector is on unit hemisphere (z >= 0)
    if (v.z() >= 0.0f && std::abs(v.norm() - 1.0f) < Epsilon) {
        return Frame::cosTheta(v) * INV_PI;
    }
    return 0.0f;
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
    // Beckmann distribution sampling
    // phi = 2*pi*xi1 (uniform azimuth)
    float phi = 2.0f * M_PI * sample.x();
    
    // For theta, we need to invert the CDF of the longitudinal part
    // CDF: F(theta) = 1 - exp(-tan²(theta)/alpha²)
    // Inverse: theta = atan(sqrt(-alpha² * log(1 - xi2)))
    float tanTheta2 = -alpha * alpha * std::log(1.0f - sample.y());
    float tanTheta = std::sqrt(std::max(0.0f, tanTheta2));
    float cosTheta = 1.0f / std::sqrt(1.0f + tanTheta2);
    float sinTheta = tanTheta * cosTheta;
    
    return Vector3f(
        sinTheta * std::cos(phi),
        sinTheta * std::sin(phi),
        cosTheta
    );
}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
    // Beckmann distribution: D(theta,phi) = (1/2π) * (2*exp(-tan²θ/α²))/(α²*cos³θ)
    // Check if vector is on unit hemisphere (z >= 0)
    if (m.z() <= 0.0f || std::abs(m.norm() - 1.0f) > Epsilon) {
        return 0.0f;
    }
    
    float cosTheta = Frame::cosTheta(m);
    if (cosTheta <= 0.0f) return 0.0f;
    
    float tanTheta2 = Frame::tanTheta(m) * Frame::tanTheta(m);
    float cosTheta2 = cosTheta * cosTheta;
    float cosTheta3 = cosTheta2 * cosTheta;
    
    float numerator = 2.0f * std::exp(-tanTheta2 / (alpha * alpha));
    float denominator = alpha * alpha * cosTheta3;
    
    return INV_TWOPI * numerator / denominator;
}

// Global cache for hierarchical samplers to avoid reloading images
static std::map<std::string, std::shared_ptr<HierarchicalSampler>> g_hierarchicalSamplers;

Point2f Warp::squareToHierarchical(const Point2f &sample, const std::string &filename) {
    // Get or create the hierarchical sampler for this filename
    auto it = g_hierarchicalSamplers.find(filename);
    if (it == g_hierarchicalSamplers.end()) {
        auto sampler = std::make_shared<HierarchicalSampler>(filename);
        g_hierarchicalSamplers[filename] = sampler;
        return sampler->sample(sample);
    }
    return it->second->sample(sample);
}

float Warp::squareToHierarchicalPdf(const Point2f &p, const std::string &filename) {
    // Get or create the hierarchical sampler for this filename
    auto it = g_hierarchicalSamplers.find(filename);
    if (it == g_hierarchicalSamplers.end()) {
        auto sampler = std::make_shared<HierarchicalSampler>(filename);
        g_hierarchicalSamplers[filename] = sampler;
        return sampler->pdf(p);
    }
    return it->second->pdf(p);
}

// HierarchicalSampler implementation
HierarchicalSampler::HierarchicalSampler(const std::string &filename) {
    // Load the image
    Bitmap bitmap(filename);
    int width = bitmap.cols();
    int height = bitmap.rows();
    
    // Check that dimensions are powers of 2
    if ((width & (width - 1)) != 0 || (height & (height - 1)) != 0) {
        throw NoriException("Image dimensions must be powers of 2 for hierarchical sampling!");
    }
    
    // Convert to luminance
    std::vector<float> luminance(width * height);
    float totalLuminance = 0.0f;
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            float lum = getLuminance(bitmap(y, x));
            luminance[y * width + x] = lum;
            totalLuminance += lum;
        }
    }
    
    // Normalize
    m_normalization = totalLuminance;
    for (float &lum : luminance) {
        lum /= totalLuminance;
    }
    
    // Build the mipmap pyramid
    buildPyramid(luminance, width, height);
}

Point2f HierarchicalSampler::sample(const Point2f &sample) const {
    Point2f result = sample;
    int currentLevel = m_levels - 1; // Start from the top (smallest) level
    int x = 0, y = 0;
    
    // Traverse down the pyramid
    while (currentLevel >= 0) {
        int levelWidth = m_sizes[currentLevel].x();
        int levelHeight = m_sizes[currentLevel].y();
        
        if (levelWidth == 1 && levelHeight == 1) {
            // Base case: single pixel
            break;
        }
        
        // Current pixel coordinates at this level
        int currentX = std::min(x, levelWidth - 1);
        int currentY = std::min(y, levelHeight - 1);
        
        if (levelWidth > 1) {
            // Choose horizontal subdivision
            float leftSum = 0.0f, rightSum = 0.0f;
            
            // Calculate sum for left half
            for (int dy = 0; dy < std::min(2, levelHeight); ++dy) {
                int ty = std::min(currentY + dy, levelHeight - 1);
                leftSum += m_pyramid[currentLevel][ty * levelWidth + currentX];
                if (currentX + 1 < levelWidth) {
                    rightSum += m_pyramid[currentLevel][ty * levelWidth + (currentX + 1)];
                }
            }
            
            float totalSum = leftSum + rightSum;
            if (totalSum > 0.0f) {
                float leftProb = leftSum / totalSum;
                
                if (result.x() < leftProb) {
                    // Go left
                    result.x() = result.x() / leftProb;
                    x = currentX * 2;
                } else {
                    // Go right
                    result.x() = (result.x() - leftProb) / (1.0f - leftProb);
                    x = currentX * 2 + 1;
                }
            } else {
                x = currentX * 2;
            }
        }
        
        if (levelHeight > 1) {
            // Choose vertical subdivision
            float topSum = 0.0f, bottomSum = 0.0f;
            
            // Calculate sum for top half
            for (int dx = 0; dx < std::min(2, levelWidth); ++dx) {
                int tx = std::min(currentX + dx, levelWidth - 1);
                topSum += m_pyramid[currentLevel][currentY * levelWidth + tx];
                if (currentY + 1 < levelHeight) {
                    bottomSum += m_pyramid[currentLevel][(currentY + 1) * levelWidth + tx];
                }
            }
            
            float totalSum = topSum + bottomSum;
            if (totalSum > 0.0f) {
                float topProb = topSum / totalSum;
                
                if (result.y() < topProb) {
                    // Go top
                    result.y() = result.y() / topProb;
                    y = currentY * 2;
                } else {
                    // Go bottom
                    result.y() = (result.y() - topProb) / (1.0f - topProb);
                    y = currentY * 2 + 1;
                }
            } else {
                y = currentY * 2;
            }
        }
        
        currentLevel--;
    }
    
    // Convert discrete pixel coordinates to continuous [0,1) coordinates
    int baseWidth = m_sizes[0].x();
    int baseHeight = m_sizes[0].y();
    
    result.x() = (x + result.x()) / baseWidth;
    result.y() = (y + result.y()) / baseHeight;
    
    return result;
}

float HierarchicalSampler::pdf(const Point2f &p) const {
    if (p.x() < 0.0f || p.x() >= 1.0f || p.y() < 0.0f || p.y() >= 1.0f) {
        return 0.0f;
    }
    
    int baseWidth = m_sizes[0].x();
    int baseHeight = m_sizes[0].y();
    
    // Convert continuous coordinates to discrete pixel coordinates
    int x = std::min((int)(p.x() * baseWidth), baseWidth - 1);
    int y = std::min((int)(p.y() * baseHeight), baseHeight - 1);
    
    // Return normalized probability density
    float pixelProb = m_pyramid[0][y * baseWidth + x];
    return pixelProb * baseWidth * baseHeight; // Convert to density
}

void HierarchicalSampler::buildPyramid(const std::vector<float> &luminance, int width, int height) {
    m_levels = 0;
    int w = width, h = height;
    
    // Count levels
    while (w > 1 || h > 1) {
        m_levels++;
        w = std::max(1, w / 2);
        h = std::max(1, h / 2);
    }
    m_levels++; // Add the base level
    
    // Initialize pyramid
    m_pyramid.resize(m_levels);
    m_sizes.resize(m_levels);
    
    // Base level (level 0)
    m_pyramid[0] = luminance;
    m_sizes[0] = Vector2i(width, height);
    
    // Build higher levels
    for (int level = 1; level < m_levels; ++level) {
        int prevWidth = m_sizes[level - 1].x();
        int prevHeight = m_sizes[level - 1].y();
        int currWidth = std::max(1, prevWidth / 2);
        int currHeight = std::max(1, prevHeight / 2);
        
        m_sizes[level] = Vector2i(currWidth, currHeight);
        m_pyramid[level].resize(currWidth * currHeight);
        
        // Downsample from previous level
        for (int y = 0; y < currHeight; ++y) {
            for (int x = 0; x < currWidth; ++x) {
                float sum = 0.0f;
                
                // Sum 2x2 block from previous level
                for (int dy = 0; dy < 2; ++dy) {
                    for (int dx = 0; dx < 2; ++dx) {
                        int px = x * 2 + dx;
                        int py = y * 2 + dy;
                        
                        if (px < prevWidth && py < prevHeight) {
                            sum += m_pyramid[level - 1][py * prevWidth + px];
                        }
                    }
                }
                
                m_pyramid[level][y * currWidth + x] = sum;
            }
        }
    }
}

float HierarchicalSampler::getLuminance(const Color3f &color) const {
    // Standard luminance weights for RGB
    return 0.212671f * color.r() + 0.715160f * color.g() + 0.072169f * color.b();
}

NORI_NAMESPACE_END
