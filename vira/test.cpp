#include "vira/vira.hpp"

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/blocked_range2d.h"
#include <clipper2/clipper.h>
#include <nanoflann.hpp>
#include <triangle.h>

#include <string>

using TFloat = float;
using TColor = vira::RGB;
using TMeshFloat = float;

template <typename TMeshFloat>
const TMeshFloat INVALID = std::numeric_limits<TMeshFloat>::infinity();

namespace std {
    template <typename TMeshFloat>
    struct hash<vira::vec2<TMeshFloat>> {
        size_t operator()(const vira::vec2<TMeshFloat>& vec) const {
            size_t h1 = std::hash<TMeshFloat>()(vec.x);
            size_t h2 = std::hash<TMeshFloat>()(vec.y);
            return h1 ^ (h2 << 1);
        }
    };
}

struct PairHash {
    template <typename T1, typename T2>
    std::size_t operator()(const std::pair<T1, T2>& p) const {
        auto hash1 = std::hash<T1>{}(p.first);
        auto hash2 = std::hash<T2>{}(p.second);
        return hash1 ^ (hash2 << 1);
    }
};

// Write buffers directly to OBJ file:
void writeOBJ(const std::string& filepath, std::vector<vira::vec2<TMeshFloat>>& vb, vira::IndexBuffer& ib)
{
    // Be aware, this is not a feature complete OBJ writer!
    if (ib.size() % 3 != 0) {
        throw std::runtime_error("Invalid index buffer.  Must be divisible by 3!");
    }

    vira::makePath(filepath);
    std::ofstream fout;
    fout.open(filepath);
    fout << "OBJ file written by Vira\n";

    // Write the vertices:
    for (size_t i = 0; i < vb.size(); i++) {
        fout << "v " << vb[i].x << " " << vb[i].y << " " << 0 << "\n";
    };

    // Write the faces:
    for (size_t j = 0; j < ib.size() / 3; j++) {
        uint32_t f0 = ib[3 * j + 0];
        uint32_t f1 = ib[3 * j + 1];
        uint32_t f2 = ib[3 * j + 2];

        fout << "f " << f0 + 1 << " " << f1 + 1 << " " << f2 + 1 << "\n";
    };
};

void scaleVertices(std::vector<vira::vec2<TMeshFloat>>& vb, TMeshFloat scale)
{
    for (size_t i = 0; i < vb.size(); i++) {
        vb[i] = scale * vb[i];
    }
};

// TO BE DELETED --------------------------------------------------------------------------------
template <typename TMeshFloat>
void exportBoundingBoxToOBJ(const std::string& filename,
                            TMeshFloat smallMinX, TMeshFloat smallMinY,
                            TMeshFloat smallMaxX, TMeshFloat smallMaxY) {
    std::ofstream objFile(filename);
    if (!objFile.is_open()) {
        std::cerr << "Unable to open file: " << filename << std::endl;
        return;
    }

    // Define the four corner points of the bounding box
    std::vector<vira::vec2<TMeshFloat>> points = {
        {smallMinX, smallMinY},
        {smallMaxX, smallMinY},
        {smallMaxX, smallMaxY},
        {smallMinX, smallMaxY}
    };

    // Export points to OBJ format
    objFile << "# Bounding Box Points\n";
    for (const auto& point : points) {
        objFile << "v " << point.x << " " << point.y << " 0.0\n";
    }

    objFile << "# Bounding Box Edges\n";
    objFile << "l 1 2\n";
    objFile << "l 2 3\n";
    objFile << "l 3 4\n";
    objFile << "l 4 1\n";

    objFile.close();
    std::cout << "Bounding box exported to " << filename << std::endl;
}
// --------------------------------------------------------------------------------


template <typename TMeshFloat>
void removeOverlap(
    std::vector<vira::vec2<TMeshFloat>>& projectedVertices,
    vira::IndexBuffer& projectedIndices,
    std::unordered_set<vira::vec2<TMeshFloat>>& surroundingPoints,
    const Clipper2Lib::PathD& smallMeshEdges,
    Clipper2Lib::PathD& surroundingMesh)
{
    using namespace Clipper2Lib;
    surroundingPoints.clear();
    surroundingMesh.clear();

    double minX = std::numeric_limits<double>::max();
    double minY = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double maxY = std::numeric_limits<double>::lowest();
    for (const auto& point : smallMeshEdges) {
        if (point.x < minX) minX = point.x;
        if (point.y < minY) minY = point.y;
        if (point.x > maxX) maxX = point.x;
        if (point.y > maxY) maxY = point.y;
    }

    // // FOR LATER: Find the small mesh axis-aligned inner bounding box, automatically set all points inside to INVALID
    std::vector<bool> triangleMarkedForDeletion(projectedIndices.size() / 3, false);

    // Parallelizing produces unknown behavior creating artifacts
    for (size_t i = 0; i < projectedIndices.size() / 3; ++i) {
        std::array<uint32_t, 3> indices = {projectedIndices[3 * i], projectedIndices[3 * i + 1], projectedIndices[3 * i + 2]};
        bool allInside = true;
        
        for (uint32_t index : indices) {
            double x = static_cast<double>(projectedVertices[index].x);
            double y = static_cast<double>(projectedVertices[index].y);

            // FOR LATER: Preliminary check for the inner bounding box

            // Preliminary check for the outer bounding box
            if (x < minX || x > maxX || y < minY || y > maxY) {
                allInside = false;
                break;
            }

            // PointInPolygon check inside bounding box
            PointD point{x, y};
            if (PointInPolygon(point, smallMeshEdges) == PointInPolygonResult::IsOutside) {
                allInside = false;
                break;
            }
        }
        
        if (allInside) {
            triangleMarkedForDeletion[i] = true;
        }
    }

    std::vector<int> vertexTriangleCount(projectedVertices.size(), 0);
    std::vector<int> vertexTriangleDeletionCount(projectedVertices.size(), 0);

    for (size_t i = 0; i < projectedIndices.size(); i += 3) {
        std::array<uint32_t, 3> indices = {projectedIndices[i], projectedIndices[i + 1], projectedIndices[i + 2]};
        for (uint32_t index : indices) {
            vertexTriangleCount[index]++;
            if (triangleMarkedForDeletion[i / 3]) {
                vertexTriangleDeletionCount[index]++;
            }
        }
    }

    for (size_t i = 0; i < projectedVertices.size(); ++i) {
        if (vertexTriangleCount[i] == vertexTriangleDeletionCount[i] && vertexTriangleCount[i] > 0) {
            projectedVertices[i] = {::INVALID<TMeshFloat>, ::INVALID<TMeshFloat>};
        }
    }
}

// template <typename TMeshFloat>
// void removeOverlap(
//     std::vector<vira::vec2<TMeshFloat>>& projectedVertices,
//     vira::IndexBuffer& projectedIndices,
//     std::unordered_set<vira::vec2<TMeshFloat>>& surroundingPoints,
//     const Clipper2Lib::PathD& smallMeshEdges,
//     Clipper2Lib::PathD& surroundingMesh)
// {
//     using namespace Clipper2Lib;
//     surroundingPoints.clear();
//     surroundingMesh.clear();

//     double minX = std::numeric_limits<double>::max();
//     double minY = std::numeric_limits<double>::max();
//     double maxX = std::numeric_limits<double>::lowest();
//     double maxY = std::numeric_limits<double>::lowest();
//     for (const auto& point : smallMeshEdges) {
//         if (point.x < minX) minX = point.x;
//         if (point.y < minY) minY = point.y;
//         if (point.x > maxX) maxX = point.x;
//         if (point.y > maxY) maxY = point.y;
//     }

//     std::vector<bool> triangleMarkedForDeletion(projectedIndices.size() / 3, false);

//     for (size_t i = 0; i < projectedIndices.size() / 3; ++i) {
//         std::array<uint32_t, 3> indices = {projectedIndices[3 * i], projectedIndices[3 * i + 1], projectedIndices[3 * i + 2]};
//         bool allInside = true;

//         for (uint32_t index : indices) {
//             double x = static_cast<double>(projectedVertices[index].x);
//             double y = static_cast<double>(projectedVertices[index].y);

//             if (x < minX || x > maxX || y < minY || y > maxY) {
//                 allInside = false;
//                 break;
//             }

//             PointD point{x, y};
//             if (PointInPolygon(point, smallMeshEdges) == PointInPolygonResult::IsOutside) {
//                 allInside = false;
//                 break;
//             }
//         }

//         if (allInside) {
//             triangleMarkedForDeletion[i] = true;
//         }
//     }

//     std::vector<int> vertexTriangleCount(projectedVertices.size(), 0);
//     std::vector<int> vertexTriangleDeletionCount(projectedVertices.size(), 0);

//     for (size_t i = 0; i < projectedIndices.size(); i += 3) {
//         std::array<uint32_t, 3> indices = {projectedIndices[i], projectedIndices[i + 1], projectedIndices[i + 2]};
//         for (uint32_t index : indices) {
//             vertexTriangleCount[index]++;
//             if (triangleMarkedForDeletion[i / 3]) {
//                 vertexTriangleDeletionCount[index]++;
//             }
//         }
//     }

//     // Add points from triangles marked for deletion but not themselves deleted
//     for (size_t i = 0; i < projectedVertices.size(); ++i) {
//         if (vertexTriangleCount[i] > 0 && vertexTriangleDeletionCount[i] > 0 && vertexTriangleCount[i] != vertexTriangleDeletionCount[i]) {
//             surroundingPoints.insert(projectedVertices[i]);
//         }
//         if (vertexTriangleCount[i] == vertexTriangleDeletionCount[i] && vertexTriangleCount[i] > 0) {
//             projectedVertices[i] = {::INVALID<TMeshFloat>, ::INVALID<TMeshFloat>};
//         }
//     }

//     // Handle edge cases for corners and boundaries
//     for (const auto& point : surroundingPoints) {
//         // Additional checks for boundary conditions and precise corner handling can be done here
//     }
// }


// void adjustHeight()

// TO BE DELETED --------------------------------------------------------------------------------
template <typename TMeshFloat>
std::vector<vira::vec2<TMeshFloat>> generateBoundingPolygonPoints(const vira::vec2<int>& ldem60Resolution) {
    std::vector<vira::vec2<TMeshFloat>> boundingPolygon;

    // Bottom edge
    for (int x = 0; x < ldem60Resolution.x; ++x) {
        boundingPolygon.push_back({static_cast<TMeshFloat>(x), 0});
    }
    // Right edge
    for (int y = 1; y < ldem60Resolution.y; ++y) {  // Start from 1 to avoid repeating the corner
        boundingPolygon.push_back({static_cast<TMeshFloat>(ldem60Resolution.x - 1), static_cast<TMeshFloat>(y)});
    }
    // Top edge
    for (int x = ldem60Resolution.x - 2; x >= 0; --x) {  // Start from ldem60Resolution.x - 2 to avoid repeating the corner
        boundingPolygon.push_back({static_cast<TMeshFloat>(x), static_cast<TMeshFloat>(ldem60Resolution.y - 1)});
    }
    // Left edge
    for (int y = ldem60Resolution.y - 2; y > 0; --y) {  // Start from ldem60Resolution.y - 2 to avoid repeating the corner
        boundingPolygon.push_back({0, static_cast<TMeshFloat>(y)});
    }

    return boundingPolygon;
}

void printSmallMeshEdgeHeights(const std::vector<vira::vec2<TMeshFloat>>& smallMeshEdges, const vira::RasterMap& outerDEM) {
    for (const auto& point : smallMeshEdges) {
        // Assuming point.x and point.y are the coordinates and outerDEM stores the heights
        int x = static_cast<int>(point.x);
        int y = static_cast<int>(point.y);

        // Ensure the coordinates are within bounds of the DEM
        if (x >= 0 && x < outerDEM.getWidth() && y >= 0 && y < outerDEM.getHeight()) {
            TMeshFloat height = outerDEM.getValue(x, y);  // Use getValue method
            std::cout << "Height at (" << x << ", " << y << "): " << height << std::endl;
        } else {
            std::cout << "Coordinates (" << x << ", " << y << ") are out of bounds." << std::endl;
        }
    }
}

// --------------------------------------------------------------------------------

int main() {
    vira::GeoTIFFInterface geotiffInterface;

    // Load LDEM4:
    std::string ldem4path = "REDACTED";
    vira::RasterMap ldem4HeightMap = geotiffInterface.load(ldem4path);
    ldem4HeightMap.appendRight(ldem4HeightMap);
    vira::GeoTIFF<TMeshFloat, TColor> ldem4(ldem4HeightMap);

    // Load LDEM60S:
    std::string ldem60Spath = "REDACTED";
    vira::RasterMap ldem60HeightMap = geotiffInterface.load(ldem60Spath);
    vira::GeoTIFF<TMeshFloat, TColor> ldem60(ldem60HeightMap);
    
    // Sub-sample the DEMs for faster computations:
    ldem4.subsample(3);
    //ldem60.subsample(40);z
    ldem60.subsample(100);


    // ================================== //
    // === PROJECT LDEM4 ONTO LDEM60S === //
    // ================================== //
    // Obtain the resized maps:
    ldem4HeightMap = ldem4.getHeightMap();
    ldem60HeightMap = ldem60.getHeightMap();
    std::vector<vira::vec2<float>> edgeVertices = ldem60HeightMap.getPixelOutline();
    vira::IndexBuffer edgeIB;
    writeOBJ("perimeter.obj", edgeVertices, edgeIB);

    // Weird hack to only load the bottom half of the DEM (this skips the problematic final row of the image):
    vira::RasterConfig config = ldem4HeightMap.getConfig();
    config.yoff = 1;
    vira::Image<float> heights = ldem4HeightMap.getData();
    vira::Resolution newResolution = heights.getResolution();
    newResolution.y = newResolution.y - 1;
    vira::Image<float> newHeights(newResolution, 0.f);
    for (uint32_t i = 0; i < newResolution.x; i++) {
        for (uint32_t j = 0; j < newResolution.y; j++) {
            newHeights(i, j) = heights(i, j + 1);
        }
    }
    ldem4HeightMap = vira::RasterMap(newHeights, config);
    ldem4 = vira::GeoTIFF<TMeshFloat, TColor>(ldem4HeightMap);

    // Project LDEM4 vertices onto LDEM60S:
    std::vector<vira::vec2<float>> projectedVertices = ldem4HeightMap.projectOnto(ldem60HeightMap);
    vira::IndexBuffer projectedIndices = ldem4.makeIndexBuffer();


    // Generate Reference vertices:
    vira::Resolution ldem60Resolution = ldem60.getHeights().getResolution();
    std::vector<vira::vec2<float>> referenceVertices(ldem60Resolution.x * ldem60Resolution.y);
    size_t idx = 0;
    for (size_t i = 0; i < ldem60Resolution.x; i++) {
        for (size_t j = 0; j < ldem60Resolution.y; j++) {
            referenceVertices[idx] = vira::vec2<float>{i, j};
            idx++;
        }
    }
    vira::IndexBuffer referenceIndices = ldem60.makeIndexBuffer();

    // REMOVE LATER: Define the bounding polygon for LDEM60:
    std::vector<vira::vec2<TMeshFloat>> boundingPolygon = generateBoundingPolygonPoints<TMeshFloat>(ldem60Resolution);
    scaleVertices(boundingPolygon, .001);
    
    Clipper2Lib::PathD smallMeshEdges;
    for (const auto& edge : boundingPolygon) {
        smallMeshEdges.push_back({static_cast<double>(edge.x), static_cast<double>(edge.y)});
    }
    smallMeshEdges.push_back(smallMeshEdges.front());

    
    // =========================== //
    // === Your code goes here === //
    // =========================== //

    std::cout << "Removing overlap ... \n" << std::flush;
    // std::vector<vira::vec2<TMeshFloat>> surroundingPoints;
    std::unordered_set<vira::vec2<TMeshFloat>> surroundingPoints;
    Clipper2Lib::PathD surroundingMesh;
    scaleVertices(projectedVertices, .001);
    removeOverlap(projectedVertices, projectedIndices, surroundingPoints, smallMeshEdges, surroundingMesh);

    heights = ldem4HeightMap.getData();
    std::vector<float>& heightData = heights.getData();
    for (size_t i = 0; i < projectedVertices.size(); ++i){
        if (std::isinf(projectedVertices[i].x)) {
            heightData[i] = INVALID<float>;
        }
    }

    vira::RasterConfig config2 = ldem4HeightMap.getConfig();
    ldem4HeightMap = vira::RasterMap(heights, config2);
    ldem4 = vira::GeoTIFF<TMeshFloat, TColor>(ldem4HeightMap);

    printSmallMeshEdgeHeights(smallMeshEdges, ldem60HeightMap);

    // std::cout << "Calculating new infill triangles...\n" << std::flush;
    // uint32_t N = 10; // Placeholder (to be replaced)
    // std::vector<vira::vec2<TMeshFloat>> filledVertices(N);
    // vira::IndexBuffer filledIndexBuffer(3*N);

    // create_mesh(smallMeshEdges, surroundingMesh, filledVertices, filledIndexBuffer);

    // Write everything to to OBJ files:
    // scaleVertices(projectedVertices, .001);
    scaleVertices(referenceVertices, .001);
    writeOBJ("projected.obj", projectedVertices, projectedIndices);
    writeOBJ("reference.obj", referenceVertices, referenceIndices);
    // writeOBJ("filled.obj", filledVertices, filledIndexBuffer);

    vira::OBJOptions objOptions;
    objOptions.writeUV = true;
    vira::OBJInterface<TMeshFloat> objInterface(objOptions);

    vira::Mesh<TMeshFloat, TColor> ldem4Mesh(ldem4);
    vira::Mesh<TMeshFloat, TColor> ldem60Mesh(ldem60);


    vira::VertexBuffer<TMeshFloat, TColor> ldem4Vertices = ldem4Mesh.getVertexBuffer();
    vira::IndexBuffer ldem4Indices = ldem4Mesh.getIndexBuffer();
    // Apply removal to the actual mesh vertices here!
    ldem4Mesh = vira::Mesh<TMeshFloat, TColor>(ldem4Vertices, ldem4Indices);


    // Rescale so it loads in blender correctly:
    ldem4Mesh.applyScale(.0001);
    ldem60Mesh.applyScale(.0001);
    
    objInterface.write("ldem4.obj", ldem4Mesh);
    objInterface.write("ldem60S.obj", ldem60Mesh);

    return 0;
};