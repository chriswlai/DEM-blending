# Working Mesh Merge without Grids

import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, Point
from shapely.ops import triangulate
import random
from scipy.spatial import Delaunay

def random_point_in_polygon(polygon):
    min_x, min_y, max_x, max_y = polygon.bounds
    while True:
        p = Point(random.uniform(min_x, max_x), random.uniform(min_y, max_y))
        if polygon.contains(p):
            return p
        
def create_dem(boundary_points, num_internal_points):
    boundary_polygon = Polygon(boundary_points)
    internal_points = [random_point_in_polygon(boundary_polygon) for _ in range(num_internal_points)]

    all_points = boundary_points[:-1] + [point.coords[0] for point in internal_points]
    return all_points, boundary_polygon

def scale_boundary_points(boundary_points, scale_factor):
    scaled_points = [(x * scale_factor, y * scale_factor) for x, y in boundary_points]
    return scaled_points

def rotate_boundary_points(boundary_points, angle, center):
    angle_rad = np.radians(angle)
    cos_angle = np.cos(angle_rad)
    sin_angle = np.sin(angle_rad)
    cx, cy = center
    rotated_points = []
    for x, y in boundary_points:
        x_new = cos_angle * (x - cx) - sin_angle * (y - cy) + cx
        y_new = sin_angle * (x - cx) + cos_angle * (y - cy) + cy
        rotated_points.append((x_new,y_new))
    return rotated_points

def translate_boundary_points(boundary_points, new_center):
    cx, cy = new_center
    ox, oy = Polygon(boundary_points).centroid.coords[0]
    translated_points = [(x + cx - ox, y + cy - oy) for x, y in boundary_points]
    return translated_points

def is_dem2_within_dem1(boundary_points_dem1, boundary_points_dem2):
    polygon_dem1 = Polygon(boundary_points_dem1)
    for point in boundary_points_dem2:
        if not polygon_dem1.contains(Point(point)):
            return False
    return True

initial_boundary_points_dem1 = [(0,0), (5,0), (5,5), (0,5), (0,0)]
initial_boundary_points_dem2 = [(1,1), (4,1), (4,4), (1,4), (1,1)]

scale_factor_dem1 = 2.0
scale_factor_dem2 = 0.5

rotation_angle_dem2 = random.uniform(0, 360)
random_position_dem2 = (random.uniform(0, 5), random.uniform(0, 5))

boundary_points_dem1 = scale_boundary_points(initial_boundary_points_dem1, scale_factor_dem1)
boundary_points_dem2 = scale_boundary_points(initial_boundary_points_dem2, scale_factor_dem2)

centroid_dem2 = Polygon(boundary_points_dem2).centroid.coords[0]
centroid_dem1 = Polygon(boundary_points_dem1).centroid.coords[0]

boundary_points_dem2 = rotate_boundary_points(boundary_points_dem2, rotation_angle_dem2, centroid_dem2)

# print("boundary", boundary_points_dem2)

translation_vector = (centroid_dem1[0] - centroid_dem2[0], centroid_dem1[1] - centroid_dem2[1])
boundary_points_dem2 = translate_boundary_points(boundary_points_dem2, translation_vector)

points_dem1, boundary_polygon_dem1 = create_dem(boundary_points_dem1, 20)
points_dem2, boundary_polygon_dem2 = create_dem(boundary_points_dem2, 5)

triangles_dem1 = triangulate(Polygon(points_dem1 + boundary_points_dem1[:-1]))

# filtered_triangles_dem1 = [tri for tri in triangles_dem1 if not boundary_polygon_dem2.intersects(tri)]
surrounding_points_dem1 = []
filtered_triangles_dem1 = []
for tri in triangles_dem1:
    if not boundary_polygon_dem2.intersects(tri):
        filtered_triangles_dem1.append(tri)
    else:
        for point in tri.exterior.coords:
            if not boundary_polygon_dem2.contains(Point(point)):
                surrounding_points_dem1.append(point)

# IMPROVE THIS PORTION TO FIND EXACTLY THE SURROUDING POINTS
# buffer_distance = 3
# surrounding_points_dem1 = [point for point in points_dem1 if boundary_polygon_dem2.buffer(buffer_distance).contains(Point(point))]

# print("boundary", boundary_points_dem2)

combined_points = boundary_points_dem2 + surrounding_points_dem1
print("combined", combined_points)

tri = Delaunay(np.array(combined_points))
triangles_between = [Polygon([combined_points[simplex[0]], combined_points[simplex[1]], combined_points[simplex[2]]]) for simplex in tri.simplices]

filtered_triangles_between = []
for tri_b in triangles_between:
    # if Polygon(surrounding_points_dem1).intersects(tri_b):
    #     triangles_between.remove(tri_b)
    # if not boundary_polygon_dem2.intersects(tri_b) and all(not tri_b.intersects(blue_triangle) for blue_triangle in triangles_dem1):
    #     filtered_triangles_between.append(tri_b)

    # First condition removes triangles formed outside DEM2 overlapping DEM1
    # Second condition removes triangles formed inside DEM2
    if Polygon(boundary_points_dem2).intersects(tri_b) and not tri_b.within(Polygon(boundary_points_dem2)):
        filtered_triangles_between.append(tri_b)
    
    # if not tri_b.within(Polygon(boundary_points_dem2)):
    #     if Polygon(surrounding_points_dem1).contains(tri_b):
    #         filtered_triangles_between.append(tri_b)




def plot_triangles(triangles, color, label):
    for tri in triangles:
        x, y = tri.exterior.xy
        plt.fill(x, y, alpha=0.5, color=color, label=label)

plt.figure()
plot_triangles(triangles_dem1, 'blue', 'DEM 1')
plot_triangles(triangulate(Polygon(points_dem2 + boundary_points_dem2)), 'green', 'DEM 2')
plt.title("Original DEMs")
plt.show()

plt.figure()
plot_triangles(filtered_triangles_dem1, 'blue', 'DEM 1')
x, y = zip(*boundary_points_dem1)
plt.plot(x, y, color='blue')
x, y = zip(*boundary_points_dem2)
plt.plot(x, y, color='green')
plt.title("Filtered DEM 1")
plt.show()

plt.figure()
plot_triangles(filtered_triangles_dem1, 'blue', 'DEM 1')
plot_triangles(triangulate(Polygon(points_dem2 + boundary_points_dem2)), 'green', 'DEM 2')
for triangle in filtered_triangles_between:
    x, y = triangle.exterior.xy
    plt.fill(x, y, alpha=0.5, color='red')
x, y = zip(*boundary_points_dem1)
plt.plot(x, y, color='blue')
x, y = zip(*boundary_points_dem2)
plt.plot(x, y, color='green')
plt.title("Generating filling triangles")
plt.show()