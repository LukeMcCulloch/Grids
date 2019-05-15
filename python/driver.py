# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 23:31:29 2015

@author: lukemcculloch

This Python version was built to debug a C language version 
of the same Delaunay mesher using Lawson's algorithm.

 /* Alma's description (also good wih my code):
 
     ELEMENT NODE INDEXING CONVENTIONS
     ---------------------------------
     triangle index notation:
                       2---1        tri[t][0]   t is the number of the triangle
                        \ /         tri[t][1]   t ranges from 0 to nt-1
                         0          tri[t][2]

    Access node numbers as follows:
    n0 = tri[t][0];
    n1 = tri[t][1];
    n2 = tri[t][2];
    Retrieve coordinates as follows:
    x0 = x[n0];
    y0 = y[n0];
    x1 = x[n1];
    y1 = y[n1];
    x2 = x[n2];
    y2 = y[n2];

     ELEMENT NEIGHBOR INDEXING CONVENTIONS
     ---------------------------------
                         1
                       2---1        nbrs[t][0]  is neighbor on side with nodes 0 & 1
                      2 \ / 0       nbrs[t][1]  is neighbor on side with nodes 1 & 2
                         0          nbrs[t][2]  is neighbor on side with nodes 2 & 0

  */
"""
import numpy as np
from trimesh import trimesh

from random import randint, seed

class mesh_driver(object):
    """Pick a file to test with!
    
    Warning, anything but btest.pts is hideously slow 
    in this version!
    """
    def __init__(self,
                 #filename = r'test_files/30P30N.pts'):
                 #filename = r'test_files/naca0012.pts'):
                 #filename = r'test_files/small.pts'):
                 filename = r'test_files/btest.pts'):
                 #filename = r'test_files/btest2.pts'):
        self.read_input(filename)
        self.get_pts()
        self.get_boundary()
        return
        
    def read_input(self, filename = r'btest.pts'):
        f = open(filename)
        self.lines = f.readlines()
        f.close()
        return 
        
    def get_pts(self):
        self.nn = int(self.lines[1].split()[0])
        self.ipts = 2
        self.epts = 2+self.nn
        self.pts = np.zeros((self.nn,2),float)
        for i, el in enumerate(range(self.ipts,self.epts)):
            self.pts[i,:] = [
                            float(el) for el in (
                                    self.lines[el].split()[0:2])
                            ]
        # or could match random input
        #        seed(4)
        #        xs = [randint(1, 98) for x in range(self.epts)]
        #        ys = [randint(1, 98) for x in range(self.epts)]
        #        for i, el in enumerate(
        #                range(self.ipts,self.epts)):
        #            self.pts[i,:] = xs[i],ys[i]
        return
        
    def get_boundary(self):
        self.nb = int(self.lines[self.epts+1].split()[0])
        self.bs = []
        p = self.epts+2
        self.nbs = []
        for i in range(self.nb):
            nseg = int(self.lines[p+1].split()[0])
            self.nbs.append(nseg)
            for j in range(nseg):
                segment = [int(el) for el in self.lines[p+2+j].split()[0:2]]
                segment.insert(0,i)
                self.bs.append(segment)
            p = p+ 3 + nseg -1
        self.bs = np.asarray(self.bs) 
        return

        
    
        

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    print "\n===================================================="
    print "\n     Driver for 2d unstructured meshing module      "
    print "\n===================================================="
    
    mesh_driver1 = mesh_driver()
    #mesh_driver1.read_input()
    #mesh_driver1.get_pts()
    mesh_driver1.get_boundary()
    mesh = trimesh(mesh_driver1)
    mesh.verbose = False
    mesh.veryverbose = False
    mesh.insert_points()
    self = mesh
    #mesh.plot_mesh()
    
    print "\n===================================================="
    print "\n     Starting Boundary Reconstruction   "
    print "\n===================================================="
    #"""
    mesh.dopics = False
    mesh.verbose = False
    mesh.veryverbose = False
    mesh.boundary_reconstruction(
            savesequence=False)
    mesh.plot_mesh()
    #"""