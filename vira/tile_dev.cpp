#include "vira/vira.hpp"

#include <string>
#include <iostream>

using TFloat = float;
using TColor = vira::RGB;
using TMeshFloat = TFloat;

int main() {
    vira::GeoTIFFInterface geotiffInterface;

    std::string path = "REDACTED";

    // Load the DEM:
    vira::GeoTIFFImage map = geotiffInterface.load(path);
    map.appendRight(map);

    vira::GeoTIFF<TMeshFloat, TColor> dem(map);

    // Create the tiles:
    std::vector<vira::GeoTIFF<TMeshFloat, TColor>> tiles = dem.tile(100000);

    // Load albedo map:
    vira::GeoTIFFImage ldamAlbedoMap = geotiffInterface.load("REDACTED");
    //ldamAlbedoMap = ldamAlbedoMap.subsample(10);

    // Convert albedo tif to PNG:
    vira::Image<float> albedo = ldamAlbedoMap.getData();
    vira::ImageInterface imageInterface;
    imageInterface.write("albedo.png", albedo);

    // Process each tile:
    for (size_t i = 0; i < tiles.size(); ++i) {
        tiles[i].UVMap(ldamAlbedoMap);
        tiles[i].subsample(3);
    }

    // Check recursive tiling:
    size_t index = 10;
    vira::GeoTIFF<TMeshFloat, TColor> tileToSplit = tiles[index];
    tiles.erase(tiles.begin() + index);
    std::vector<vira::GeoTIFF<TMeshFloat, TColor>> tiles2 = tileToSplit.tile(1000);
    for (size_t i = 0; i < tiles2.size(); ++i) {
        tiles2[i].subsample(2);
    }
    tiles.insert(tiles.end(), tiles2.begin(), tiles2.end());

    // Create the Meshes:
    vira::OBJOptions objOptions;
    objOptions.writeUV = true;
    vira::OBJInterface<TMeshFloat> objInterface(objOptions);
    for (size_t i = 0; i < tiles.size(); ++i) {
        vira::Mesh<TMeshFloat, TColor> mesh(tiles[i]);
        mesh.applyScale(.0001);
        objInterface.write("tile_" + std::to_string(i) + ".obj", mesh);
    }

    return 0;
};