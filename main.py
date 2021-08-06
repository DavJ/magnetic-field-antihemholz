from stl import mesh
import math
import numpy as np
from field_calculator import AntiHemholtz

generator_radius = 1
generator_offset = 0.5
steps_vertical = 100
steps_horizontal = 720


generator = AntiHemholtz(1, generator_radius, generator_offset)

theta0 = math.asin(generator_radius/generator_offset)
r = generator_radius
z = generator_offset
dtheta = (math.pi/2 - theta0) / steps_vertical
dfi = 2 * math.pi / steps_horizontal

vertices = []
faces = []

for fi in np.linspace(0, 2*math.pi, steps_horizontal):
    for theta in np.linspace(theta0, math.pi/2, steps_vertical):
        r1, z1 = r, z
        r, z = r, z + generator.normal_vector(r, z) * 2 * dtheta

        vert1 = [r1 * math.cos(fi), r1 * math.sin(fi), z1]
        vert2 = [r * math.cos(fi), r * math.sin(fi), z]
        vert3 = [r1 * math.cos(fi + dfi), r1 * math.sin(fi + dfi), z1]
        vert4 = [r * math.cos(fi + dfi), r * math.sin(fi + dfi), z]
        vertices.append([vert1, vert2, vert3, vert4])

        face1 = [vert1, vert2, vert3]
        face2 = [vert2, vert4, vert3]
        faces.append([face1, face2])

# Create the mesh
ufo = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
for i, f in enumerate(faces):
    for j in range(3):
        ufo.vectors[i][j] = vertices[f[j],:]

# Write the mesh to file "cube.stl"
ufo.save('ufo.stl')