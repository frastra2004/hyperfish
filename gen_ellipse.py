import gmsh
import os
import numpy as np
from mpi4py import MPI

import pandas as pd

class ellipse():
    def __init__(self, L, R_eq, R_bore, Wall_angle, Dome_b, Dome_a_over_b, Iris_b, Iris_a_over_b):
        self.L = L
        self.R_eq = R_eq
        self.R_bore = R_bore
        self.Wall_angle = (Wall_angle*np.pi)/180
        self.Dome_b = Dome_b
        self.Dome_a_over_b = Dome_a_over_b
        self.Iris_b = Iris_b
        self.Iris_a_over_b = Iris_a_over_b

        #iris ellipse parameters
        self.iris_center_y =self.R_bore + self.Iris_b
        self.iris_center_x = 0

        #dome ellipse parameters
        self.dome_center_y = self.L - self.Dome_b
        self.dome_center_x = self.L/2



    def generate(self):
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 1)
        gmsh.model.add("ell_cav")
        p1 = gmsh.model.geo.addPoint(0, 0, 0)
        p2 = gmsh.model.geo.addPoint(0, self.R_bore, 0)
        p3 = gmsh.model.geo.addPoint(self.Iris_a_over_b * self.Iris_b*np.cos(np.arctan(-np.tan(self.Wall_angle) / self.Iris_a_over_b)) + self.iris_center_x, self.iris_center_y+ self.Iris_b*np.sin(np.arctan(-np.tan(self.Wall_angle) / self.Iris_a_over_b)), 0)
        p4 = gmsh.model.geo.addPoint(self.dome_center_x - self.Dome_a_over_b * self.Dome_b*np.cos(np.arctan(-np.tan(self.Wall_angle) / self.Dome_a_over_b)), self.dome_center_y + self.Dome_b*np.sin(np.arctan(np.tan(self.Wall_angle) / self.Dome_a_over_b)), 0)
        p5 = gmsh.model.geo.addPoint(self.dome_center_x, self.dome_center_y+self.Dome_b, 0)
        p6 = gmsh.model.geo.addPoint(self.dome_center_x + self.Dome_a_over_b * self.Dome_b*np.cos(np.arctan(-np.tan(self.Wall_angle) / self.Dome_a_over_b)), self.dome_center_y + self.Dome_b*np.sin(np.arctan(np.tan(self.Wall_angle) / self.Dome_a_over_b)), 0)
        p7 = gmsh.model.geo.addPoint(-self.Iris_a_over_b * self.Iris_b*np.cos(np.arctan(-np.tan(self.Wall_angle) / self.Iris_a_over_b)) + self.iris_center_x + self.L, self.iris_center_y+ self.Iris_b*np.sin(np.arctan(-np.tan(self.Wall_angle) / self.Iris_a_over_b)), 0)
        p8 = gmsh.model.geo.addPoint(self.L, self.R_bore, 0)
        p9 = gmsh.model.geo.addPoint(self.L, 0, 0)
        
        

        if self.Iris_a_over_b<1.0:  # major axis is vertical
            p10 = gmsh.model.geo.addPoint(self.iris_center_x, self.iris_center_y, 0) # Left iris
            p11 = gmsh.model.geo.addPoint(self.iris_center_x, self.iris_center_y+0.1*self.Iris_b, 0) # defines point on maj axis

            p14 = gmsh.model.geo.addPoint(self.iris_center_x+self.L, self.iris_center_y, 0) # right iris
            p15= gmsh.model.geo.addPoint(self.iris_center_x+self.L, self.iris_center_y+0.1*self.Iris_b, 0) # defines point on min axis

        else:                                                                                       # major axis is horizontal
            p10 = gmsh.model.geo.addPoint(self.iris_center_x, self.iris_center_y, 0)
            p11 = gmsh.model.geo.addPoint(self.iris_center_x+0.1*self.Iris_b*self.Iris_a_over_b, self.iris_center_y, 0) # defines point on maj axis
            p14 = gmsh.model.geo.addPoint(self.iris_center_x+self.L, self.iris_center_y, 0) # right iris
            p15= gmsh.model.geo.addPoint(self.iris_center_x+self.L+0.1*self.Iris_b*self.Iris_a_over_b, self.iris_center_y, 0) # defines point on min axis

        
        if self.Dome_a_over_b<1.0:
            p12 = gmsh.model.geo.addPoint(self.dome_center_x, self.dome_center_y, 0)
            p13 = gmsh.model.geo.addPoint(self.dome_center_x, self.dome_center_y + 0.1*self.Dome_b, 0)
        else:
            p12 = gmsh.model.geo.addPoint(self.dome_center_x, self.dome_center_y, 0)
            p13 = gmsh.model.geo.addPoint(self.dome_center_x+ 0.1*self.Dome_b*self.Dome_a_over_b, self.dome_center_y, 0)

        l1 = gmsh.model.geo.addLine(p1,p2)
        l2 = gmsh.model.geo.addEllipseArc(p2, p10, p11, p3)
        l3 = gmsh.model.geo.addLine(p3, p4)
        l4 = gmsh.model.geo.addEllipseArc(p4, p12, p13, p5)
        l5 = gmsh.model.geo.addEllipseArc(p5, p12, p13, p6)
        l6 = gmsh.model.geo.addLine(p6, p7)
        l7 = gmsh.model.geo.addEllipseArc(p7, p14, p15, p8)
        l8 = gmsh.model.geo.addLine(p8,p9)
        l9 = gmsh.model.geo.addLine(p9,p1)

        cl1 = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4, l5, l6, l7, l8,l9])
        s1 = gmsh.model.geo.addPlaneSurface([cl1])

        gmsh.model.geo.synchronize()
        surf_tag = gmsh.model.addPhysicalGroup(2, [s1])
        gmsh.model.setPhysicalName(2, surf_tag, "domain")
        line_tags = [l1, l2, l3, l4, l5, l6, l7, l8]
        bnd_tag = gmsh.model.addPhysicalGroup(1, line_tags)
        gmsh.model.setPhysicalName(1, bnd_tag, "boundary")

        mesh_size = 0.02*self.L # decrease for finer mesh
        gmsh.model.mesh.setSize(gmsh.model.getEntities(0), mesh_size)
        #print('hi')
        gmsh.model.mesh.generate(2)
        #print('hey')

        #msh_filename = "ell_shape_new.msh"
        out = os.path.expanduser("~/Documents/DPhil/hyperfish/ell_shape_new.msh")
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
                



        
#ellll= ellipse(L=0.1153035,R_eq=0.1092505,R_bore=0.035,Wall_angle=14,Dome_b=0.016, Dome_a_over_b=1.2, Iris_b=0.07602, Iris_a_over_b=0.5 )
#ellll.generate()



