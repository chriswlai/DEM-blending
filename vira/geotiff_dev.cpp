#include "vira/vira.hpp"

#include <string>
#include <iostream>
#include <cmath>
#include <vector>

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/blocked_range2d.h"


using TFloat = double;
using TColor = vira::RGB;
using TMeshFloat = TFloat;

template <typename TMeshFloat>
const TMeshFloat INVALID = std::numeric_limits<TMeshFloat>::infinity();

template <typename T>
struct vec2 {
    T x, y;
};

float bilinearInterpolate(const vira::vec2<float>& point, vira::GeoTIFFImage& map) {
    int x0 = static_cast<int>(floor(point.x));
    int y0 = static_cast<int>(floor(point.y));
    int x1 = x0 + 1;
    int y1 = y0 + 1;

    float dx = point.x - x0;
    float dy = point.y - y0;

    float h00 = map(static_cast<uint32_t>(x0), static_cast<uint32_t>(y0));
    float h10 = map(static_cast<uint32_t>(x1), static_cast<uint32_t>(y0));
    float h01 = map(static_cast<uint32_t>(x0), static_cast<uint32_t>(y1));
    float h11 = map(static_cast<uint32_t>(x1), static_cast<uint32_t>(y1));

    float height = (1 - dx) * (1 - dy) * h00 +
                   dx * (1 - dy) * h10 +
                   (1 - dx) * dy * h01 +
                   dx * dy * h11;

    return height;
}

std::vector<std::vector<float>> createGaussianKernel(int size, float sigma) {
    std::vector<std::vector<float>> kernel(size, std::vector<float>(size));
    float sum = 0.0f;
    int halfSize = size / 2;
    float twoSigmaSquare = 2.0f * sigma * sigma;

    for (int x = -halfSize; x <= halfSize; ++x) {
        for (int y = -halfSize; y <= halfSize; ++y) {
            float r = std::sqrt(x*x + y*y);
            kernel[x + halfSize][y + halfSize] = std::exp(-(r*r) / twoSigmaSquare) / (M_PI * twoSigmaSquare);
            sum += kernel[x + halfSize][y + halfSize];
        }
    }

    // Normalize the kernel
    for (int x = 0; x < size; ++x) {
        for (int y = 0; y < size; ++y) {
            kernel[x][y] /= sum;
        }
    }

    return kernel;
}

void applyGaussianBlur(vira::Image<float>& smallHeightDeltas, int kernelSize, float sigma) {
    auto kernel = createGaussianKernel(kernelSize, sigma);
    vira::Image<float> blurredDeltas(smallHeightDeltas.resolution(), 0.0f);

    int halfSize = kernelSize / 2;
    for (uint32_t i = 0; i < smallHeightDeltas.resolution().x; ++i) {
        for (uint32_t j = 0; j < smallHeightDeltas.resolution().y; ++j) {
            float sum = 0.0f;
            for (int ki = -halfSize; ki <= halfSize; ++ki) {
                for (int kj = -halfSize; kj <= halfSize; ++kj) {
                    int ni = i + ki;
                    int nj = j + kj;
                    if (ni >= 0 && ni < smallHeightDeltas.resolution().x && nj >= 0 && nj < smallHeightDeltas.resolution().y) {
                        sum += smallHeightDeltas(ni, nj) * kernel[ki + halfSize][kj + halfSize];
                    }
                }
            }
            blurredDeltas(i, j) = sum;
        }
    }

    // Replace original deltas with blurred deltas
    smallHeightDeltas = blurredDeltas;
}

int main() {
    vira::GeoTIFFInterface geotiffInterface;

    // Paths to DEMs:
    std::string bigDEM_path = "REDACTED";
    std::string smallDEM_path = "REDACTED";
    bool appendRight = true;

    // std::string bigDEM_path = "REDACTED";
    // std::string smallDEM_path = "REDACTED";
    // bool appendRight = false;

    // Subsampling settings:
    size_t subsample_bigDEM = 1;
    size_t subsample_smallDEM = 25;

    // size_t subsample_bigDEM = 25;
    // size_t subsample_smallDEM = 200;

    // Load the bigDEM map:
    vira::GeoTIFFImage bigDEMHeightMap = geotiffInterface.load(bigDEM_path);
    if (appendRight) {
        bigDEMHeightMap.appendRight(bigDEMHeightMap);
    };
    bigDEMHeightMap.subsample(subsample_bigDEM);

    // Load smallDEM map:
    vira::GeoTIFFImage smallDEMHeightMap = geotiffInterface.load(smallDEM_path);
    smallDEMHeightMap.subsample(subsample_smallDEM);


    // ================================== //
    // === PROJECT bigDEM ONTO smallDEM === //
    // ================================== //

    // Project bigDEM vertices onto smallDEM:
    vira::Image<vira::vec2<float>> projectedMap = bigDEMHeightMap.projectOnto(smallDEMHeightMap);

    // Create index buffer:
    vira::IndexBuffer indexBuffer = vira::ImageUtils::imageToIndexBuffer(bigDEMHeightMap.getData());

    vira::Resolution resolution = smallDEMHeightMap.resolution();

    // Get the boundary pixels:
    std::vector<vira::vec2<uint32_t>> smallDEMperimeter = vira::ImageUtils::getValidPerimeter(smallDEMHeightMap.getData());

    const float N = 5.0;

    std::cout << "Removing overlap ... \n" << std::flush;
    std::vector<bool> triangleMarkedForDeletion(indexBuffer.size() / 3, false);

    // Looping over triangles (NOTICE: Parallelizing produces unknown behavior creating artifacts)
    for (size_t i = 0; i < indexBuffer.size() / 3; ++i) {
        std::array<uint32_t, 3> indices = { indexBuffer[3 * i], indexBuffer[3 * i + 1], indexBuffer[3 * i + 2] };
        bool allInside = true;
        for (uint32_t index : indices) {
            double x = static_cast<double>(projectedMap[index].x);
            double y = static_cast<double>(projectedMap[index].y);

            if ((smallDEMHeightMap.checkValidPixel(x, y)) && (!(x > N) || !(x < (resolution.x - 1 - N)) || !(y > N) || !(y < (resolution.y - 1 - N)))) {
                allInside = false;
                break;
            }

            if (!smallDEMHeightMap.checkValidPixel(x, y)) {
                allInside = false;
                break;
            }
        }
        if (allInside) {
            triangleMarkedForDeletion[i] = true;
        }
    }

    std::vector<int> vertexTriangleCount(projectedMap.size(), 0);
    std::vector<int> vertexTriangleDeletionCount(indexBuffer.size(), 0);

    for (size_t i = 0; i < indexBuffer.size(); i += 3) {
        std::array<uint32_t, 3> indices = { indexBuffer[i], indexBuffer[i + 1], indexBuffer[i + 2] };
        for (uint32_t index : indices) {
            vertexTriangleCount[index]++;
            if (triangleMarkedForDeletion[i / 3]) {
                vertexTriangleDeletionCount[index]++;
            }
        }
    }

    // vira::GeoTIFFConfig config2 = bigDEMHeightMap.getConfig();
    vira::Resolution resolution2 = bigDEMHeightMap.resolution();
    vira::Image<float> bigHeightDeltas(resolution2, 0.0f);

    float M = static_cast<float>(resolution2.x) / static_cast<float>(resolution.x);

    for (size_t i = 0; i < projectedMap.size(); ++i) {
        if (vertexTriangleCount[i] == vertexTriangleDeletionCount[i] && vertexTriangleCount[i] > 0) {
            bigDEMHeightMap[i] = std::numeric_limits<float>::infinity();
        }

        if (vertexTriangleCount[i] > 0 && vertexTriangleDeletionCount[i] > 1 && vertexTriangleCount[i] != vertexTriangleDeletionCount[i]) { 
            float innerHeightToChange = bigDEMHeightMap[i] + (bigHeightDeltas(projectedMap[i].x, projectedMap[i].y));
            float innerHeightToCompareTo = bilinearInterpolate(projectedMap[i], smallDEMHeightMap);

            if (i == 801016) {
                std::cout << "value\n";
            }
            float innerHeightToChange2 = bigDEMHeightMap[i+1];
            float innerHeightToCompareTo2 = bilinearInterpolate(projectedMap[i+1], smallDEMHeightMap);
        
            if (innerHeightToChange >= innerHeightToCompareTo) {
                float temp = bigHeightDeltas(projectedMap[i].x, projectedMap[i].y);
                bigHeightDeltas(projectedMap[i].x, projectedMap[i].y) = (innerHeightToCompareTo - innerHeightToChange) + temp;
                innerHeightToChange = innerHeightToCompareTo;
            }
                
            vira::vec2<float> direction = glm::normalize(projectedMap[i+1] - projectedMap[i]);
            float length = glm::length(projectedMap[i+1] - projectedMap[i]);
            float stepSize =  length / (M+1);
            float maximumDeviation = 0.f;
            float bigDEMheightDelta = (bigDEMHeightMap[i+1] - innerHeightToChange)/(M+1);

            for (size_t j = 1; j < M+1; ++j) {
                vira::vec2<float> samplePoint = projectedMap[i] + j*stepSize*direction;

                // Interpolate the height at the midpoint in smallDEM
                float samplePointHeightSmallDEM = bilinearInterpolate(samplePoint, smallDEMHeightMap);

                // Calculate sample point's height on bigDEM:
                float samplePointHeightBigDEM = innerHeightToChange + j*bigDEMheightDelta;    

                if (samplePointHeightBigDEM > samplePointHeightSmallDEM) {
                    float heightDifference = samplePointHeightBigDEM - samplePointHeightSmallDEM;
                    maximumDeviation = std::max(maximumDeviation, heightDifference);
                }
            }
            bigHeightDeltas(projectedMap[i].x, projectedMap[i].y) -= maximumDeviation;
            bigHeightDeltas(projectedMap[i+1].x, projectedMap[i+1].y) -= maximumDeviation;
        }
    }

    for (size_t i = 0; i < bigDEMHeightMap.size(); ++i) {
        if (std::isinf(projectedMap[i].x) || std::isinf(projectedMap[i].y) || std::isnan(projectedMap[i].x) || std::isnan(projectedMap[i].y)) {
            continue;
        }
        if (projectedMap[i].x < 0 || projectedMap[i].y < 0 || projectedMap[i].x >= resolution2.x || projectedMap[i].y >= resolution2.y) {
            continue;
        }
        bigDEMHeightMap[i] += bigHeightDeltas(projectedMap[i].x, projectedMap[i].y);
    }

    
    vira::Image<float> smallHeightDeltas(resolution, 0.0f);

    vira::Image<vira::vec2<float>> referencedMap = smallDEMHeightMap.projectOnto(bigDEMHeightMap);

    for (size_t i = 0; i < smallDEMperimeter.size()-1; ++i) {
        const vira::vec2<uint32_t>& perimeterPixel1 = smallDEMperimeter[i];
        const vira::vec2<uint32_t>& perimeterPixel2 = smallDEMperimeter[i+1];

        vira::vec2<float> projectedPixel1 = referencedMap(perimeterPixel1.x, perimeterPixel1.y);
        vira::vec2<float> projectedPixel2 = referencedMap(perimeterPixel2.x, perimeterPixel2.y);

        float outerHeightToChange1 = smallDEMHeightMap(perimeterPixel1.x, perimeterPixel1.y);
        float outerHeightToChange2 = smallDEMHeightMap(perimeterPixel2.x, perimeterPixel2.y);

        float outerHeightToCompareTo1 = bilinearInterpolate(projectedPixel1, bigDEMHeightMap);
        float outerHeightToCompareTo2 = bilinearInterpolate(projectedPixel2, bigDEMHeightMap);

        if (outerHeightToChange1 > outerHeightToCompareTo1) {
            float temp = smallHeightDeltas(perimeterPixel1.x, perimeterPixel1.y);
            smallHeightDeltas(perimeterPixel1.x, perimeterPixel1.y) = (outerHeightToCompareTo1 - outerHeightToChange1) + temp;
        }
        if (outerHeightToChange2 > outerHeightToCompareTo2) {
            float temp = smallHeightDeltas(perimeterPixel2.x, perimeterPixel2.y);
            smallHeightDeltas(perimeterPixel2.x, perimeterPixel2.y) = (outerHeightToCompareTo2 - outerHeightToChange2) + temp;
        }
    }

    vira::Image<float> fadeDeltas2 = smallHeightDeltas;
    vira::Image<uint8_t> counts2(fadeDeltas2.resolution(), 0);

    for (uint32_t i = 0; i < resolution.x; ++i) {
        for (uint32_t j = 0; j < resolution.y; ++j) {
            for (auto ii = std::max(0u, i-static_cast<uint32_t>(M)); ii < std::min(resolution.x, i+static_cast<uint32_t>(M)); ++ii) {
                for (auto jj = std::max(0u, j-static_cast<uint32_t>(M)); jj< std::min(resolution.y, j+static_cast<uint32_t>(M)); ++jj) {
                    float distance = std::sqrt(std::pow(static_cast<float>(i) - static_cast<float>(ii),2.f) + std::pow(static_cast<float>(j) - static_cast<float>(jj),2.f));
                    // fadeDeltas2(ii,jj) = smallHeightDeltas(i,j) * (M-distance);
                    fadeDeltas2(ii, jj) = std::min(0.0f, smallHeightDeltas(i,j) * (M - distance));
                    if (std::min(0.0f, smallHeightDeltas(i,j) * (M - distance)) > 0) {
                        std::cout << smallHeightDeltas(i,j) << " " << (M-distance) << " | ";
                    }
                    counts2(ii,jj)++;
                }
            }
        }
    }

    // Reset the original starting edges (the above process messes them up):
    for (vira::vec2<uint32_t>& edgeVertices : smallDEMperimeter) {
        fadeDeltas2(edgeVertices.x,edgeVertices.y) = smallHeightDeltas(edgeVertices.x,edgeVertices.y);
    }

    // Normalize based on the counts for each pixel:
    for (size_t i = 0; i < fadeDeltas2.size(); ++i) {
        fadeDeltas2[i] = fadeDeltas2[i]/static_cast<float>(counts2[i]);
    }

    // Apply a Gaussian blur to fadeDeltas2 before using it to adjust the smallDEMHeightMap
    int kernelSize = 5;
    float sigma = 1.0f;
    applyGaussianBlur(fadeDeltas2, kernelSize, sigma);

    for (size_t i = 0; i < smallDEMHeightMap.size(); ++i) {
        smallDEMHeightMap[i] += fadeDeltas2[i];
    }

    // Write to OBJ files:
    vira::OBJOptions objOptions;
    objOptions.writeUV = true;
    vira::OBJInterface<TMeshFloat> objInterface(objOptions);

    vira::GeoTIFF<TMeshFloat, TColor> bigDEM(bigDEMHeightMap);
    vira::Mesh<TMeshFloat, TColor> bigDEMMesh(bigDEM);
    bigDEMMesh.applyScale(.0001);
    objInterface.write("bigDEM.obj", bigDEMMesh);

    vira::GeoTIFF<TMeshFloat, TColor> smallDEM(smallDEMHeightMap);
    vira::Mesh<TMeshFloat, TColor> smallDEMMesh(smallDEM);
    smallDEMMesh.applyScale(.0001);
    objInterface.write("smallDEM.obj", smallDEMMesh);



    return 0;
};
