import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, Point
from shapely.ops import triangulate
from scipy.spatial import Delaunay

import random

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

# print(corners_dem2)

surrounding_points_dem1 = []
filtered_triangles_dem1 = []
for tri in tri_dem1:
    if not Polygon(corners_dem2).intersects(tri):
        filtered_triangles_dem1.append(tri)
    else:
        for point in tri.exterior.coords:
            if not Polygon(corners_dem2).contains(Point(point)):
                surrounding_points_dem1.append(point)

def point_in_tri_exclusive(pt, v1, v2, v3):
    def sign(p1, p2, p3):
        return (p1[0]-p3[0]) * (p2[1]-p3[1]) - (p2[0]-p3[0]) * (p1[1]-p3[1])
    
    d1 = sign(pt, v1, v2)
    d2 = sign(pt, v2, v3)
    d3 = sign(pt, v3, v1)

    return (d1 > 0 and d2 > 0 and d3 > 0) or (d1 < 0 and d2 < 0 and d3 < 0)

def points_in_tri(tri, combined_points):
    v1, v2, v3 = tri
    for point in combined_points:
        if point_in_tri_exclusive(point, v1, v2, v3):
            return True
    return False

def do_lines_intersect(p1, q1, p2, q2):
    def orientation(p, q, r):
        val = (q[1] - p[1]) * (r[0] - q[0]) - (q[0] - p[0]) * (r[1] - q[1])

        if val == 0:
            return 0 # collinear
        elif val > 0:
            return 1 # clockwise
        else:
            return 2 # counterclockwise
        
    # def on_segment(p, q, r):
    #     if (q[0] <= max(p[0], r[0]) and q[0] >= min(p[0], r[0]) and q[1] <= max(p[1], r[1]) and q[1] >= min(p[1], r[1])):
    #         return True
    #     return False
    
    o1 = orientation(p1, q1, q2)
    o2 = orientation(p1, q1, q2)
    o3 = orientation(p2, q2, q1)
    o4 = orientation(p2, q2, q1)

    if o1 == 0 and o2 == 0 and o3 == 0 and o4 == 0:
        if p1 == p2 or p1 == q2 or q1 == p2 or q1 == q2:
            return False
        else:
            return True
        
    return True

    # if o1 != o2 and o3 != o4:
    #     return True
    
    # if o1 == 0 and on_segment(p1, p2, q1):
    #     return True
    
    # if o2 == 0 and on_segment(p1, q2, q1):
    #     return True
    
    # if o3 == 0 and on_segment(p2, p1, q2):
    #     return True
    
    # if o4 == 0 and on_segment(p2, q1, q2):
    #     return True
    
    # return False

def does_triangle_overlap(exist_tri, new_tri):
    # exist_tri = exist_tri.exterior.coords

    # exist_tri = [(exist_tri[i], exist_tri[i+1], exist_tri[0]) for i in range(1, len(exist_tri-1))]

    new_edges = [
        (new_tri[0], new_tri[1]),
        (new_tri[1], new_tri[2]),
        (new_tri[2], new_tri[0])
    ]

    for tri in exist_tri: 
        existing_edges = [
            (tri.exterior.coords[0], tri.exterior.coords[1]),
            (tri.exterior.coords[1], tri.exterior.coords[2]),
            (tri.exterior.coords[2], tri.exterior.coords[0])
        ]
            # existing_edges = [
            #     (tri[0], tri[1]),
            #     (tri[1], tri[2]),
            #     (tri[2], tri[0])
            # ]
        for new_edge in new_edges:
            for existing_edge in existing_edges:
                if do_lines_intersect(new_edge[0], new_edge[1], existing_edge[0], existing_edge[1]):
                    return True
    return False

def create_mesh(edge_dem2, surrounding_points_dem1):
    combined_points = np.array(edge_dem2 + surrounding_points_dem1)
    # print("combined", combined_points)

    delaunay = Delaunay(combined_points)

    valid_triangles = []
    for simplex in delaunay.simplices:
        vertices = combined_points[simplex]
        in_edge = [tuple(v) in edge_dem2 for v in vertices]
        in_surrounding = [tuple(v) in surrounding_points_dem1 for v in vertices]

        if any(in_edge) and any(in_surrounding):
            print("checkpoint 1")
            if not points_in_tri(vertices, combined_points):
            # if not does_triangle_overlap(filtered_triangles_dem1, vertices):
                print("checkpoint 2")
                valid_triangles.append(vertices)

    return valid_triangles

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