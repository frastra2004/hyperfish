from gen_ellipse import ellipse
import meshio
import matplotlib.pyplot as plt

import argparse

def parse_args():
    p = argparse.ArgumentParser(description="Run fenics job from GUI")
    p.add_argument("dome_b", type=float, help="Dome b (m)")
    p.add_argument("dome_a_b", type=float, help="Dome a/b")
    p.add_argument("iris_b", type=float, help="Iris b (m)")
    p.add_argument("iris_a_b", type=float, help="Iris a/b")
    p.add_argument("bore_rad", type=float, help="Bore radius (m)")
    p.add_argument("eq_radius", type=float, help="Equator radius (m)")
    p.add_argument("wall_angle", type=float, help="Wall angle (deg)")
    p.add_argument("length", type=float, help="Length (m)")
    return p.parse_args()

args = parse_args()

ellll= ellipse(L=args.length,R_eq=args.eq_radius,R_bore=args.bore_rad,Wall_angle=args.wall_angle,Dome_b=args.dome_b, Dome_a_over_b=args.dome_a_b, Iris_b=args.iris_b, Iris_a_over_b=args.iris_a_b )
ellll.generate()


# read with meshio
msh_filename2= '/Users/straniero/Documents/Dphil/hyperfish/ell_shape_new.msh'
mesh = meshio.read(msh_filename2)

# Extract points and triangular cells
points = mesh.points[:, :2]   # (x,y) but represents (z,r)
cells = None
for cell_block in mesh.cells:
    if cell_block.type == "triangle":
        cells = cell_block.data
        break


if cells is None:
    raise RuntimeError("No triangular cells found in mesh.")

print("Nodes:", points.shape[0], "Triangles:", cells.shape[0])

plt.triplot(points[:,0], points[:,1], cells, lw=0.5)
plt.gca().set_aspect('equal')
plt.xlabel('z (m)')
plt.ylabel('r (m)')
plt.title('2D axisymmetric mesh (z,r)')
plt.show()
