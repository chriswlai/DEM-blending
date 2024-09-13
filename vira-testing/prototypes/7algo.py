# Polar transformation ladder without Delaunay

import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, Point, LineString
from scipy.spatial import Delaunay
from collections import defaultdict
from itertools import combinations
import scipy

import math
import random

def create_circle(radius, num_points, center):
    theta = np.linspace(0, 2 * np.pi, num_points)
    cx, cy = center
    circle = [(cx + radius * np.cos(t), cy + radius * np.sin(t)) for t in theta]
    return circle

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
        rotated_points.append((x_new, y_new))
    return rotated_points

def translate_boundary_points(boundary_points, new_center):
    cx, cy = new_center
    ox, oy = 5, 5
    translated_points = [(x + cx - ox, y + cy - oy) for x, y in boundary_points]
    return translated_points

def is_valid_triangle(triangle, existing_triangles):
    pass

def cartesian_to_polar(points, center):
    cx, cy = center
    polar = []
    for x, y in points:
        dx = x - cx
        dy = y - cy
        r = np.sqrt(dx**2 + dy**2)
        theta = np.arctan2(dy,dx)
        polar.append((r, theta))
    return polar

def polar_in_cartesian(points, center):
    # cx, cy = center
    cartesian = []
    for r, theta in points:
        # x = r * np.cos(theta) + cx
        # y = r * np.sin(theta) + cy
        x = theta / math.pi * 180
        y = r
        # print((x, y))
        cartesian.append((x, y))
    return cartesian

def create_mesh(edge_dem2, surrounding_points_dem1):
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
    
    i, j = 0, 0
    while i < main_len - 1 and j < sec_len:
        current_vertex = main_edge[i]
        next_vertex = main_edge[i+1]
        sec_vertex = sec_edge[j]

        valid_points = []
        for k in range(i + 1, main_len):
            point = main_edge[k]
            dist_to_sec = np.linalg.norm(np.array(point) - np.array(sec_vertex))
            dist_to_curr = np.linalg.norm(np.array(point) - np.array(current_vertex))
            dist_to_next = np.linalg.norm(np.array(point) - np.array(next_vertex))
            if dist_to_sec < dist_to_curr and dist_to_sec < dist_to_next:
                valid_points.append(point)

        for point in valid_points:
            mid_point = ((sec_vertex[0] + point[0]) / 2, (sec_vertex[1] + point[1]) / 2)
            triangle = [current_vertex, next_vertex, mid_point]

            if not any(Polygon(triangle).intersects(Polygon(tri)) for tri in mesh):
                mesh.append(triangle)

        j += 1

        if j == sec_len and i < main_len - 2:
            mid_point = ((main_edge[i][0] + main_edge[i+1][0])/2, (main_edge[i][1]+main_edge[i+1][1])/2)
            triangle = [current_vertex, next_vertex, mid_point]

            if not any(Polygon(triangle).intersects(Polygon(tri)) for tri in mesh):
                mesh.append(triangle)

        i += 1

    return mesh

    

x_min, x_max, y_min, y_max = 0, 10, 0, 10
corners_dem2 = [(0, 0), (10, 0), (10, 10), (0, 10), (0, 0)]

spacing_dem1 = 0.3
spacing_dem2 = 0.25

scale_factor_dem1 = 3.0
scale_factor_dem2 = 1.0

grid_points_dem1, edge_dem1 = create_grid(x_min, x_max, y_min, y_max, spacing_dem1)
grid_points_dem2, edge_dem2 = create_grid(x_min, x_max, y_min, y_max, spacing_dem2)

boundary_points_dem1 = scale_boundary_points(grid_points_dem1, scale_factor_dem1)
boundary_points_dem2 = scale_boundary_points(grid_points_dem2, scale_factor_dem2)
corners_dem2 = scale_boundary_points(corners_dem2, scale_factor_dem2)
edge_dem2 = scale_boundary_points(edge_dem2, scale_factor_dem2)

rotation_angle_dem2 = random.uniform(0, 360)
new_center_dem2 = (random.uniform(10, 20), random.uniform(10, 20))

print("rotation: ", rotation_angle_dem2)
print("new center: ", new_center_dem2)

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

surrounding_points_dem1 = []
filtered_triangles_dem1 = []
for tri in tri_dem1:
    if not Polygon(corners_dem2).intersects(tri):
        filtered_triangles_dem1.append(tri)
    else:
        for point in tri.exterior.coords:
            if not Polygon(corners_dem2).contains(Point(point)):
                surrounding_points_dem1.append(point)

corner_tri = Delaunay(corners_dem2)
corner_tri = [Polygon([corners_dem2[simplex[0]], corners_dem2[simplex[1]], corners_dem2[simplex[2]]]) for simplex in corner_tri.simplices]

circle = create_circle(5, 100, new_center_dem2)

edge_dem2_polar = cartesian_to_polar(edge_dem2, new_center_dem2)
surrounding_points_dem1_polar = cartesian_to_polar(surrounding_points_dem1, new_center_dem2)
circle_polar = cartesian_to_polar(circle, new_center_dem2)

# print(circle_polar)

edge_dem2_cartesian = polar_in_cartesian(edge_dem2_polar, new_center_dem2)
surrounding_points_dem1_cartesian = polar_in_cartesian(surrounding_points_dem1_polar, new_center_dem2)
circle_cartesian = polar_in_cartesian(circle_polar, new_center_dem2)

edge_dem2_cartesian = list(set(edge_dem2_cartesian))
surrounding_points_dem1_cartesian = list(set(surrounding_points_dem1_cartesian))

edge_dem2_cartesian.sort()
surrounding_points_dem1_cartesian.sort()

mesh = create_mesh(edge_dem2_cartesian, surrounding_points_dem1_cartesian)

def plot_mesh(mesh, color, label):
    for tri in mesh:
        x = [point[0] for point in tri] + [tri[0][0]]
        y = [point[1] for point in tri] + [tri[0][1]]
        plt.fill(x, y, alpha=0.5, color=color, label=label)

def plot_triangles(triangles, color, label):
    for tri in triangles:
        x, y = tri.exterior.xy
        plt.fill(x, y, alpha=0.5, color=color, label=label)

def plot_polar(points, color, label):
    r = [point[0] for point in points]
    theta = [point[1] for point in points]
    plt.polar(theta, r, 'o', color=color, label=label)

def plot_cartesian(points, color, label):
    # x = [point[0] for point in points]
    # y = [point[1] for point in points]
    x, y = zip(*points)
    plt.scatter(x, y, color=color, label=label, s=10)

plt.figure()
plt.scatter(*zip(*boundary_points_dem1), color='blue', label='DEM 1', s=10)
plt.scatter(*zip(*boundary_points_dem2), color='green', label='DEM 2', s=10)
plt.scatter(*zip(*circle), color="red", label='circle', s=10)
plot_triangles(tri_dem1, 'blue', 'DEM1')
plot_triangles(tri_dem2, 'green', 'DEM2')
plt.title("Unfiltered DEMs")
plt.gca().set_aspect('equal', adjustable='box')

plt.figure()
plot_polar(edge_dem2_polar, 'blue', 'DEM 1')
plot_polar(surrounding_points_dem1_polar, 'green', 'DEM 2')
plot_polar(circle_polar, 'red', 'circle')

plt.figure()
plot_cartesian(edge_dem2_cartesian, 'blue', 'DEM 2 edge')
plot_cartesian(surrounding_points_dem1_cartesian, 'green', 'DEM 1 border')
plot_cartesian(circle_cartesian, 'red', 'circle')
plot_mesh(mesh, 'red', 'fill')
plt.show()


# plt.figure()
# plt.scatter(*zip(*surrounding_points_dem1), color='blue', label='DEM 1', s=10)
# plt.scatter(*zip(*boundary_points_dem2), color='green', label='DEM 2', s=10)
# plot_triangles(filtered_triangles_dem1, 'blue', 'DEM 1')
# plot_triangles(tri_dem2, 'green', 'DEM2')
# plt.gca().set_aspect('equal', adjustable='box')
# plt.title("Filtered DEM 1")

# plt.figure()
# plt.scatter(*zip(*surrounding_points_dem1), color='blue', label='DEM 1', s=10)
# plt.scatter(*zip(*boundary_points_dem2), color='green', label='DEM 2', s=10)
# plot_triangles(filtered_triangles_dem1, 'blue', 'DEM 1')
# plot_triangles(tri_dem2, 'green', 'DEM2')
# plot_mesh(mesh, 'red', 'fill')
# plt.gca().set_aspect('equal', adjustable='box')
# plt.title("Generating filling triangles")
# plt.show()