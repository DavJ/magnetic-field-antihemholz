from stl import mesh
import math
import numpy as np
from field_calculator import AntiHemholtz

generator_radius = 1
offset_radius = 0.05
generator_offset = 0.5
offset_offset = 0.025
steps_vertical = 100
steps_horizontal = 720
number_of_vertices = steps_vertical * steps_horizontal


generator = AntiHemholtz(1, generator_radius, generator_offset)

theta0 = math.atan((generator_radius + offset_radius)/(generator_offset - offset_offset))
r = generator_radius + offset_radius
z = generator_offset
dtheta = (math.pi/2 - theta0) / steps_vertical
dfi = 2 * math.pi / steps_horizontal

vertices = []
faces = []
i, j = 0, 0
for fi in np.linspace(0, 2*math.pi, steps_horizontal):
    for theta in np.linspace(theta0, math.pi/2, steps_vertical):
        r1, z1 = r, z
        nv_r, nv_z = generator.normal_vector(r, z)
        r, z = r + nv_r * 2. * math.pi * dtheta, z + nv_z * 2. * math.pi * dtheta

        vert1 = [r1 * math.cos(fi), r1 * math.sin(fi), z1]
        vert2 = [r * math.cos(fi), r * math.sin(fi), z]
        vertices.append(vert1)
        vert_index = len(vert1) - 1

        face1 = [vert_index, vert_index + 1, (vert_index + steps_vertical) % number_of_vertices]
        face2 = [vert_index + 1, (vert_index + 1 + steps_vertical) % number_of_vertices, (vert_index + steps_vertical) % number_of_vertices]
        faces.append(face1)
        faces.append(face2)
    vertices.append(vert2)

np_faces = np.array(faces)
np_vertices = np.array(vertices)

# Create the mesh
ufo = mesh.Mesh(np.zeros(np_faces.shape[0], dtype=mesh.Mesh.dtype))
for i, f in enumerate(faces):
    for j in range(3):
        ufo.vectors[i][j] = np_vertices[f[j],:]

# Write the mesh to file "cube.stl"
ufo.save('ufo.stl')