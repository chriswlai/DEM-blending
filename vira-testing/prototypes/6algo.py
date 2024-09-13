# Polar transformation ladder with Delaunay

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
    cartesian = []
    for r, theta in points:
        x = theta # / math.pi * 180
        y = r
        cartesian.append((x, y))
    return cartesian

def polar_to_cartesian(points, center):
    cx, cy = center
    cartesian = []
    for theta, r in points:
        theta = theta # / 180 * math.pi
        x = r * np.cos(theta) + cx
        y = r * np.sin(theta) + cy
        cartesian.append((x, y))
    return cartesian

def distance(point1, point2):
    return math.sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2)

def order_points(points):
    if not points:
        return []

    points.sort()

    for i in range(len(points)-2):
        if distance(points[i], points[i+2]) < distance(points[i], points[i+1]):
            temp = points[i+1]
            points[i+1] = points[i+2]
            points[i+2] = temp
            i += 2
        
    return points

    # unvisited = points[:]
    # ordered_points = []

    # current_point = min(unvisited)
    # print(current_point)
    # ordered_points.append(current_point)
    # unvisited.remove(current_point)

    # while unvisited:
    #     nearest_point = min(unvisited, key=lambda p:distance(current_point, p))
    #     ordered_points.append(nearest_point)
    #     unvisited.remove(nearest_point)
    #     current_point = nearest_point

    # return ordered_points

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

    def mesh_add(group, mesh):
        group = np.array(group)
        tri = Delaunay(group)
        for simplex in tri.simplices:
            vertices = group[simplex]
            in_edge = [tuple(v) in main_edge for v in vertices]
            in_surrounding = [tuple(v) in sec_edge for v in vertices]

            if any(in_edge) and any(in_surrounding) and not any(Polygon(vertices).intersection(Polygon(tri)).area > 0 for tri in mesh):
                mesh.append(vertices)

    group = []
    i, j = 0, 0
    group.append(sec_edge[j])
    while i < main_len - 1 and j < sec_len - 1:    
        group.append(sec_edge[j+1])
        while main_edge[i][0] < sec_edge[j+1][0]:
            # print(i, j+1)
            group.append(main_edge[i])
            i += 1
        
        n_group = []
        n_group.append(main_edge[i-1])
        n_group.append(sec_edge[j+1])
        j += 1

        # print(i, j, group)
        print(sec_edge[j], polar_to_cartesian([sec_edge[j]], new_center_dem2))
        if len(group) >= 3:
            mesh_add(group, mesh)
            group = n_group
        else:
            group.append(n_group)

    # Add ending
    group = []
    if i < main_len - 1:
        print("1")
        group.append(sec_edge[j])
        while i < main_len + 1:
            group.append(main_edge[i-1])
            i += 1
        mesh_add(group, mesh)
    elif j < sec_len - 1:
        print("2")
        group.append(main_edge[i])
        while j < sec_len + 1:
            group.append(sec_edge[j-1])
            j += 1
        while i < main_len + 1:
            group.append(main_edge[i-1])
            i += 1
        mesh_add(group, mesh)
    else:
        print("3")
        group.append(sec_edge[j])
        while i < main_len + 1:
            group.append(main_edge[i-1])
            i += 1
        mesh_add(group, mesh)
    
    print(group)
    print(polar_to_cartesian(group, new_center_dem2))

    # print(main_edge[0])
    # print(main_edge[-1])

    # print(polar_to_cartesian([main_edge[0]], new_center_dem2))
    # print(polar_to_cartesian([main_edge[-1]], new_center_dem2))

    # print(sec_edge[0])   
    # print(sec_edge[-1])

    # print(polar_to_cartesian([sec_edge[0]], new_center_dem2))
    # print(polar_to_cartesian([sec_edge[-1]], new_center_dem2))

    # Add beginning + end connection
    group = []
    group.append(main_edge[0])
    group.append(sec_edge[0])
    group.append(main_edge[-1])
    group.append(sec_edge[-1])

    print(polar_to_cartesian(group, new_center_dem2))

    group = np.array(group)
    tri = Delaunay(group)
    for simplex in tri.simplices:
        vertices = group[simplex]
        in_edge = [tuple(v) in main_edge for v in vertices]
        in_surrounding = [tuple(v) in sec_edge for v in vertices]

        if any(in_edge) and any(in_surrounding):
            mesh.append(vertices)

    mesh_add(group, mesh)

    # ----------------------------------------------------------------
            
    # i, j = 0, 0
    # while i < main_len - 1 and j < sec_len:
    #     group = [main_edge[i], main_edge[i+1], sec_edge[j]]
    #     mesh_add(group)

    #     if j + 1 < sec_len and main_edge[i+1][0] > sec_edge[j+1][0]:
    #         j += 1
    #     else:
    #         i += 1

    # i, j = 0, 0
    # curr_group = []
    # while i < main_len and j < sec_len:
    #     if main_edge[i][0] < sec_edge[j][0]:
    #         curr_group.append(main_edge[i])
    #         i += 1
    #     else:
    #         curr_group.append(sec_edge[j])
    #         j += 1

    #     if any(p in curr_group for p in main_edge) and any(p in curr_group for p in sec_edge):
    #         if len(curr_group) >= 3:
    #             mesh_add(curr_group)
    #         curr_group = []

    #     if curr_group:
    #         if len(curr_group) >= 3:
    #             mesh_add(curr_group)
        
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

# rotation_angle_dem2 = random.uniform(0, 360)
# new_center_dem2 = (random.uniform(10, 20), random.uniform(10, 20))

rotation_angle_dem2 =  285.2477427260016
new_center_dem2 =  (19.978115414473777, 11.174187721358601)

# rotation_angle_dem2 =  173.82664253393395
# new_center_dem2 =  (12.014897836877648, 11.910980116356068)

# rotation_angle_dem2 = 38.65513902838576
# new_center_dem2 = (14.11462082602359, 18.769192623796997)

# rotation_angle_dem2 = 35.822667233819544
# new_center_dem2 = (12.431269454420525, 11.601026779884446)

# rotation:  253.30829483134886
# new center:  (19.315218281703135, 17.88599911328245)

# rotation:  163.83589011761663
# new center:  (13.9839067338282, 17.08852541999605)


print("rotation_angle_dem2 = ", rotation_angle_dem2)
print("new_center_dem2 = ", new_center_dem2)

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

# mesh = create_mesh(order_points(edge_dem2_cartesian), order_points(surrounding_points_dem1_cartesian))

# print(mesh)

cartesian_mesh = []
for triangle in mesh:
    cartesian_triangle = polar_to_cartesian(triangle, new_center_dem2)
    cartesian_mesh.append(cartesian_triangle)

# print(cartesian_mesh)

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

# plt.figure()
# plt.scatter(*zip(*boundary_points_dem1), color='blue', label='DEM 1', s=10)
# plt.scatter(*zip(*boundary_points_dem2), color='green', label='DEM 2', s=10)
# # plt.scatter(*zip(*circle), color="red", label='circle', s=10)
# plot_triangles(tri_dem1, 'blue', 'DEM1')
# plot_triangles(tri_dem2, 'green', 'DEM2')
# plt.title("Unfiltered DEMs")
# plt.gca().set_aspect('equal', adjustable='box')

# plt.figure()
# plot_polar(edge_dem2_polar, 'blue', 'DEM 1')
# plot_polar(surrounding_points_dem1_polar, 'green', 'DEM 2')
# # plot_polar(circle_polar, 'red', 'circle')

plt.figure()
plot_cartesian(edge_dem2_cartesian, 'green', 'DEM 2 edge')
plot_cartesian(surrounding_points_dem1_cartesian, 'blue', 'DEM 1 border')
# plot_cartesian(circle_cartesian, 'red', 'circle')
plot_mesh(mesh, 'red', 'fill')

# plt.figure()
# plt.scatter(*zip(*surrounding_points_dem1), color='blue', label='DEM 1', s=10)
# plt.scatter(*zip(*boundary_points_dem2), color='green', label='DEM 2', s=10)
# plot_triangles(filtered_triangles_dem1, 'blue', 'DEM 1')
# plot_triangles(tri_dem2, 'green', 'DEM2')
# plt.gca().set_aspect('equal', adjustable='box')
# plt.title("Filtered DEM 1")

plt.figure()
plt.scatter(*zip(*surrounding_points_dem1), color='blue', label='DEM 1', s=10)
plt.scatter(*zip(*boundary_points_dem2), color='green', label='DEM 2', s=10)
plot_triangles(filtered_triangles_dem1, 'blue', 'DEM 1')
plot_triangles(tri_dem2, 'green', 'DEM2')
plot_mesh(cartesian_mesh, 'red', 'fill')
plt.gca().set_aspect('equal', adjustable='box')
plt.title("Generating filling triangles")
plt.show()