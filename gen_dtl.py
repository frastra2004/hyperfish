import gmsh
import os
import numpy as np


import pandas as pd

class DTL():
    def __init__(self, L, R_eq, R_bore, angle, g, d):
        self.L = L
        self.R_eq = R_eq
        self.R_bore = R_bore
        self.angle = (angle*np.pi)/180
        self.g = g
        self.d = d

    def generate(self):
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 1)
        gmsh.model.add("dtl_cav")
        p1 = gmsh.model.geo.addPoint(0, 0, 0)
        p2 = gmsh.model.geo.addPoint(self.L, 0, 0)
        p3 = gmsh.model.geo.addPoint(self.L, self.R_bore, 0)
        p4 = gmsh.model.geo.addPoint(self.g, self.R_bore, 0)
        p5 = gmsh.model.geo.addPoint(self.g + np.tan(self.angle)*(self.d-self.R_bore), self.d, 0)
        p6 = gmsh.model.geo.addPoint(self.L, self.d, 0)
        p7 = gmsh.model.geo.addPoint(self.L, self.R_eq, 0)
        p8 = gmsh.model.geo.addPoint(0, self.R_eq, 0)

        l1 = gmsh.model.geo.addLine(p1,p2)
        l2 = gmsh.model.geo.addLine(p2,p3)
        l3 = gmsh.model.geo.addLine(p3,p4)
        l4 = gmsh.model.geo.addLine(p4,p5)
        l5 = gmsh.model.geo.addLine(p5,p6)
        l6 = gmsh.model.geo.addLine(p6,p7)
        l7 = gmsh.model.geo.addLine(p7,p8)
        l8 = gmsh.model.geo.addLine(p8,p1)

        cl1 = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4, l5, l6, l7, l8])
        s1 = gmsh.model.geo.addPlaneSurface([cl1])

        gmsh.model.geo.synchronize()
        surf_tag = gmsh.model.addPhysicalGroup(2, [s1])
        gmsh.model.setPhysicalName(2, surf_tag, "domain")
        line_tags = [l2, l3, l4, l5, l6, l7, l8]
        bnd_tag = gmsh.model.addPhysicalGroup(1, line_tags)
        gmsh.model.setPhysicalName(1, bnd_tag, "boundary")

        mesh_size = 0.01*(self.L+self.R_eq) # decrease for finer mesh
        gmsh.model.mesh.setSize(gmsh.model.getEntities(0), mesh_size)
        #print('hi')
        gmsh.model.mesh.generate(2)
        #print('hey')

        #msh_filename = "ell_shape_new.msh"
        out = os.path.expanduser("~/Documents/DPhil/hyperfish/dtl_shape_new.msh")
        try:
            gmsh.write(out)
            print(f"\nSaved mesh to {out}")
        except Exception as e:
            print("\nError writing mesh:", e)
            # print gmsh last error if available
            try:
                print("Gmsh last error:", gmsh.logger.getLastError())
            except Exception:
                pass

        print(f"Saved mesh to {out}")

        # finalize gmsh
        gmsh.finalize()


        