# Ladder Approach with Grids

import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, Point
from shapely.ops import triangulate
from scipy.spatial import Delaunay, ConvexHull

import random
import math

def create_grid(x_min, x_max, y_min, y_max, spacing):
    x_coords = np.arange(x_min, x_max + spacing, spacing)
    y_coords = np.arange(y_min, y_max + spacing, spacing)
    grid_points = [(x, y) for x in x_coords for y in y_coords]

    edge_points = set()
    for x in x_coords:
        edge_points.add((x, y_min))
        edge_points.add((x, y_max))
    for y in y_coords:
        edge_points.add((x_min, y))
        edge_points.add((x_max, y))

    return grid_points, list(edge_points)

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

    # Has issue when calculating for border points
    # ox, oy = Polygon(boundary_points).centroid.coords[0]
    ox, oy = 5, 5
    translated_points = [(x + cx - ox, y + cy - oy) for x, y in boundary_points]
    return translated_points

def is_dem2_within_dem1(boundary_points_dem1, boundary_points_dem2):
    polygon_dem1 = Polygon(boundary_points_dem1)
    for point in boundary_points_dem2:
        if not polygon_dem1.contains(Point(point)):
            return False
    return True

x_min, x_max, y_min, y_max = 0, 10, 0, 10
corners_dem2 = [(0,0), (10,0), (10,10), (0,10), (0,0)]

spacing_dem1 = 0.3
spacing_dem2 = 0.25

scale_factor_dem1 = 3.0
scale_factor_dem2 = 1.0

grid_points_dem1, edge_dem1 = create_grid(x_min, x_max, y_min, y_max, spacing_dem1)
grid_points_dem2, edge_dem2 = create_grid(x_min, x_max, y_min, y_max, spacing_dem2)

# print(edge_dem2)

boundary_points_dem1 = scale_boundary_points(grid_points_dem1, scale_factor_dem1)
boundary_points_dem2 = scale_boundary_points(grid_points_dem2, scale_factor_dem2)
corners_dem2 = scale_boundary_points(corners_dem2, scale_factor_dem2)
edge_dem2 = scale_boundary_points(edge_dem2, scale_factor_dem2)

# rotation_angle_dem2 = 40
# new_center_dem2 = (15,15)
rotation_angle_dem2 = random.uniform(0, 360)
new_center_dem2 = (random.uniform(10, 20), random.uniform(10, 20))

boundary_points_dem2 = translate_boundary_points(boundary_points_dem2, new_center_dem2)
boundary_points_dem2 = rotate_boundary_points(boundary_points_dem2, rotation_angle_dem2, new_center_dem2)

corners_dem2 = translate_boundary_points(corners_dem2, new_center_dem2)
corners_dem2 = rotate_boundary_points(corners_dem2, rotation_angle_dem2, new_center_dem2)

edge_dem2 = translate_boundary_points(edge_dem2, new_center_dem2)
edge_dem2 = rotate_boundary_points(edge_dem2, rotation_angle_dem2, new_center_dem2)

dem1 = Delaunay(boundary_points_dem1)
dem2 = Delaunay(boundary_points_dem2)

tri_dem1 = [Polygon([boundary_points_dem1[simplex[0]], boundary_points_dem1[simplex[1]], boundary_points_dem1[simplex[2]]]) for simplex in dem1.simplices]
tri_dem2 = [Polygon([boundary_points_dem2[simplex[0]], boundary_points_dem2[simplex[1]], boundary_points_dem2[simplex[2]]]) for simplex in dem2.simplices]

# print(corners_dem2)

surrounding_points_dem1 = set()
filtered_triangles_dem1 = []
for tri in tri_dem1:
    if not Polygon(corners_dem2).intersects(tri):
        filtered_triangles_dem1.append(tri)
    else:
        for point in tri.exterior.coords:
            if not Polygon(corners_dem2).contains(Point(point)):
                surrounding_points_dem1.add(point)
surrounding_points_dem1 = list(surrounding_points_dem1)

def distance(point1, point2):
    return math.sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2)

def order_points(points):
    if not points:
        return []
    
    unvisited = points[:]
    ordered_points = []

    current_point = random.choice(unvisited)
    ordered_points.append(current_point)
    unvisited.remove(current_point)

    while unvisited:
        nearest_point = min(unvisited, key=lambda p:distance(current_point, p))
        ordered_points.append(nearest_point)
        unvisited.remove(nearest_point)
        current_point = nearest_point

    return ordered_points


def create_mesh(edge_dem2, surrounding_points_dem1):

    print("edge_dem2 len", len(edge_dem2))
    print("surrounding_points_dem1 len", len(surrounding_points_dem1))

    edge_dem2 = order_points(edge_dem2)
    surrounding_points_dem1 = order_points(surrounding_points_dem1)

    print("edge_dem2 len", len(edge_dem2))
    print("surrounding_points_dem1 len", len(surrounding_points_dem1))
    print(edge_dem2)

    if len(edge_dem2) >= len(surrounding_points_dem1):
        main_edge = edge_dem2
        sec_edge = surrounding_points_dem1
    else:
        main_edge = surrounding_points_dem1
        sec_edge = edge_dem2

    mesh = []

    main_len = len(main_edge)
    sec_len = len(sec_edge)

    print(main_len, sec_len)

    used_vertices = set()
    sec_edge_copy = sec_edge[:]

    for i in range(main_len):
        current_vertex = main_edge[i]
        prev_vertex = main_edge[i-1]

        closest_vertex = min(sec_edge, key=lambda v: np.linalg.norm(np.array(v) - np.array(current_vertex)))
        # closest_vertex = min((v for v in sec_edge_copy if v not in used_vertices), key=lambda v: np.linalg.norm(np.array(v) - np.array(current_vertex)), default=None)

        if i > 0:
            used_vertices.add(closest_vertex)
            mesh.append((prev_vertex, current_vertex, closest_vertex))

        # if closest_vertex is not None:
        #     used_vertices.add(closest_vertex)
        #     mesh.append((prev_vertex, current_vertex, closest_vertex))

    for vertex in sec_edge_copy:
        if vertex not in used_vertices:
            print("UNUSED: ", vertex)
            closest_main_vertex = min(main_edge, key=lambda v: np.linalg.norm(np.array(v) - np.array(vertex)))
            print(closest_main_vertex)
            # mesh.append((prev_vertex, current_vertex, closest_vertex))

    return mesh

mesh = create_mesh(edge_dem2, surrounding_points_dem1)

def plot_mesh(mesh, color, label):
    for tri in mesh:
        x = [point[0] for point in tri] + [tri[0][0]]
        y = [point[1] for point in tri] + [tri[0][1]]
        plt.fill(x, y, alpha=0.5, color=color, label=label)

def plot_triangles(triangles, color, label):
    for tri in triangles:
        x, y = tri.exterior.xy
        plt.fill(x, y, alpha=0.5, color=color, label=label)

plt.figure()
plt.scatter(*zip(*boundary_points_dem1), color='blue', label='DEM 1', s=10)
plt.scatter(*zip(*boundary_points_dem2), color='green', label='DEM 2', s=10)
plot_triangles(tri_dem1, 'blue', 'DEM1')
plot_triangles(tri_dem2, 'green', 'DEM2')
plt.title("Unfiltered DEMs")
plt.gca().set_aspect('equal', adjustable='box')
# plt.show()

plt.figure()
plt.scatter(*zip(*surrounding_points_dem1), color='blue', label='DEM 1', s=10)
plt.scatter(*zip(*boundary_points_dem2), color='green', label='DEM 2', s=10)
plot_triangles(filtered_triangles_dem1, 'blue', 'DEM 1')
plot_triangles(tri_dem2, 'green', 'DEM2')
plt.gca().set_aspect('equal', adjustable='box')
plt.title("Filtered DEM 1")
# plt.show()

plt.figure()
plt.scatter(*zip(*surrounding_points_dem1), color='blue', label='DEM 1', s=10)
plt.scatter(*zip(*boundary_points_dem2), color='green', label='DEM 2', s=10)
plot_triangles(filtered_triangles_dem1, 'blue', 'DEM 1')
plot_triangles(tri_dem2, 'green', 'DEM2')
plot_mesh(mesh, 'red', 'fill')
plt.gca().set_aspect('equal', adjustable='box')
plt.title("Generating filling triangles")
plt.show()