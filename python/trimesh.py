# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 23:31:29 2015

@author: lukemcculloch

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
import matplotlib.pyplot as plt
from make_neighbors import make_nbrs
from search import search, tlm_search
from circle_test import circle_test
from circle_test import My_circle_test
from circle_test import inCircleRobust

import time

eps = 1.e-12


def cross2d(u2d, v2d):
    u = np.asarray([u2d[0],u2d[1],0.])
    v = np.asarray([v2d[0],v2d[1],0.])
    return np.array((
        u[1]*v[2] - u[2]*v[1],
        u[2]*v[0] - u[0]*v[2],
        u[0]*v[1] - u[1]*v[0]
        ))
    
    
def cross(u, v):
    return np.array((
        u[1]*v[2] - u[2]*v[1],
        u[2]*v[0] - u[0]*v[2],
        u[0]*v[1] - u[1]*v[0]
        ))
    
    
class triangle(object):
    def __init__(self, p1,p2,p3):
        self.make_nodes(p1,p2,p3)
        self.make_edges(p1,p2,p3)
        return
    def make_nodes(self, p1,p2,p3):
        self.nodes = {}
        self.nodes[0] = p1
        self.nodes[1] = p2
        self.nodes[2] = p3
        return
    def make_edges(self, p1,p2,p3):
        self.edge = {}
        self.edge[1] = np.asarray([p2,p1])
        self.edge[2] = np.asarray([p3,p2])
        self.edge[3] = np.asarray([p1,p3])
        return
        


       
class trimesh(object):
    
    def __init__(self, mesh):
        self.dopics         = False
        self.verbose        = False
        self.veryverbose    = False
        self.printtime      = False
        self.figcount       = 0
        self.piccount       = 0
        self.seed           = 0
        self.pts            = mesh.pts
        self.nn             = mesh.nn
        self.bs             = mesh.bs
        self.nb             = mesh.nb
        self.nbs            = mesh.nbs
        self.tdim           = self.nn*3
        self.tri = np.asarray([[0,0,0] for el in range(mesh.nn*3)])
        self.circles = {}
        self.nt             = 0
        self.nbr            = []
        self.set_limits()
        self.create_superpoints()
        self.nbr = [[0,0,0] for el in range(self.tdim)]
        make_nbrs(self.nn, self.nt, self.tri, self.nbr)
        return
        
    def set_limits(self):
        self.xmin=10e20
        self.ymin=10e20
        self.xmax=10e-20
        self.ymax=10e-20
#        for n in range(self.mesh.nn):
#            self.xmin = min(self.xmin, self.mesh.pts[0,n])
#            self.xmax = min(self.xmax, self.mesh.pts[0,n])
#            self.ymin = min(self.ymin, self.mesh.pts[1,n])
#            self.ymax = min(self.ymax, self.mesh.pts[1,n])
        self.xmin = 1.2*min(self.xmin, np.min(self.pts[:,0]) )
        self.xmax = 1.2*max(self.xmax, np.max(self.pts[:,0]) )
        self.ymin = 1.2*min(self.ymin, np.min(self.pts[:,1]) )
        self.ymax = 1.2*max(self.ymax, np.max(self.pts[:,1]) )
        return
    
    def create_superpoints(self):
        self.supernodes = [0,1,2,3]
        self.tempx = np.zeros((self.nn+4),float)
        self.tempy = np.zeros((self.nn+4),float)
        for n in range(self.nn):
            self.tempx[n+4] = self.pts[n,0]
            self.tempy[n+4] = self.pts[n,1]
        self.tempx[0]=self.xmin 
        self.tempx[1]=self.xmax
        self.tempx[2]=self.xmax 
        self.tempx[3]=self.xmin
        self.tempy[0]=self.ymin 
        self.tempy[1]=self.ymin 
        self.tempy[2]=self.ymax 
        self.tempy[3]=self.ymax
        
        self.tri[0] = [0,1,2]
        self.tri[1] = [0,2,3]
        
        self.nt = 2
        return
        
        
    def insert_points(self):
        """
        n=4
        
        nn=self.nn
        ntri = self.nt
        tri = self.tri
        nbr = self.nbr
        pts = self.pts
        xt = self.tempx[n]
        yt = self.tempy[n]
        t=0
        i=0
        #"""
        self.set_plot_pts()
        #for n in range(4,4):
        for n in range(4, self.nn + 4): #self.nn): # 
            self.insert_point(n)
        return    
    
    def delete_node(self,index):
        return
    def delete_tri(self,index):
        return
    
    def insert_point(self, n):
        print 'n = {}'.format(n)
        if self.verbose:
            print ''
            print '----------------------'
            print 'TOP OF THE INSERT LOOP'
            print 'n = {}'.format(n)
        if self.nt>self.nn*3:
            self.nt=-1
            #print'break loop'
            return
        #else:
        self.seed = 0
        
        ##
        ## Optimize list
        ##
        self.Optimize = []

        self.xt = self.tempx[n]
        self.yt = self.tempy[n]
        
        #        if self.verbose:
        #            print 'self.seed = {}'.format(self.seed)
        t0 = time.time()
        self.seed = search(self.seed, self.xt, self.yt, 
                           self.tempx,self.tempy, 
                           self.tri,self.nbr)
        tSolve = time.time() - t0
        if self.printtime: print("tlm_search took {:.5f} seconds.".format(tSolve))
        
        if self.verbose:
            print 'the                    '
            print '    search             '
            print '           has         '
            print '               returned'
        if self.veryverbose:
            print 'returned seed = {}'.format(self.seed)
        
        #/ save nodes and neighbors of old triangle
        p0=self.tri[self.seed][0]
        p1=self.tri[self.seed][1]
        p2=self.tri[self.seed][2]  
        
        if self.veryverbose:
            self.plot_tris([self.tri[self.seed]], 
                           [self.xt,self.yt])
            self.plot_mesh()
        
        #neighbor triangles
        #        nbrtemp0=self.nbr[self.seed][0]
        #        nbrtemp1=self.nbr[self.seed][1]
        #        nbrtemp2=self.nbr[self.seed][2]
        
        #new point
        p3 = n
        
        ##// Subdivide Triangle
        ##// Add 4 points for 3 triangles to the optimize list

        # print "\nReassign 1st Triangle"
        self.tri[self.seed][0]=p0
        self.tri[self.seed][1]=p1
        self.tri[self.seed][2]=p3
        self.Optimize.append(self.seed)
        
        # print "\nReassign 2nd Triangle"
        self.tri[self.nt][0]=p1
        self.tri[self.nt][1]=p2
        self.tri[self.nt][2]=p3
        self.Optimize.append(self.nt)
        #we have added a triangle
        self.nt += 1
        
        # print "\nReassign 3rd Triangle"
        self.tri[self.nt][0]=p2
        self.tri[self.nt][1]=p0
        self.tri[self.nt][2]=p3
        self.Optimize.append(self.nt)
        #we have added a triangle
        self.nt += 1
        
        # run make neighbors - this will be sllllooooowww
        t0 = time.time()
        make_nbrs(self.nn+4, self.tdim, self.tri, self.nbr)
        tSolve = time.time() - t0
        if self.printtime: print("make_nbrs took {:.5f} seconds.".format(tSolve))
        
        self.nopt = 3 #len(self.nopt) # num of tris to optimize

        ## optimize:
        
        
        t0 = time.time()
        self.compute_optimized_nodes()
        tSolve = time.time() - t0
        if self.printtime: print("Optimize took {:.5f} seconds.".format(tSolve))
        
        if self.dopics:
            self.plot_mesh()
        return
    
    def br(self):
        flipcount = 0
        #"""node to tri hash table"""
        nhash = []
        for n in range(self.nn+4):
            nhash.append([])
            
        for t in range(self.nt):
            for i in range(3):
                n = self.tri[t][i]
                nhash[n].append(t)
        return
        
    
    
    
    
    
    def boundary_reconstruction(self):
        #flipcount = 0
        """
        nhash ==> node : tri   hashtable
            
            {node:tri}
            
            #len(nhash) == len(self.pts)+4
            
            nhash[0] = [0, 1, 5]
            node 0 is used in tris 0,1,5
            
            
            self.bs[0,0,1]
            self.bs[0,1,2]
            self.bs[0,2,3]
            
            bs[0] is the 1st segment of side 0 & goes from node 0+4 to node 1+4
            bs[1] is the 2nd segment of side 0 & goes from node 1+4 to node 2+4
            
            
        
        """
        nhash = [] 
        for n in range(self.nn+4):
            nhash.append([])
            
        for t in range(self.nt):
            for i in range(3):
                n = self.tri[t][i]
                nhash[n].append(t)
                
        #        for i in range(self.nb): #loop boundaries
        #            for j in range(self.nbs[i]): #loop segments
        #                b1=self.bs[i][j][0]+4
        #                b2=self.bs[i][j][1]+4
        ## oops, messed up bs indexes.  
        ## lets use a bit more python
        ## instead of fixing the indexes:
        #chgs = True
        #while chgs:
        #
        # bs=self.bs[0]
        #
        """
        z=4
        currenttriangle=15
        
        m=1
        n=4
        
        """
        success = False
        for seg,bs in enumerate(self.bs):
            print '----------------------------'
            print 'start of the boundaries loop'
            
            
            #print '\n\n ck New Boundary segment = ',[el+4 for el in bs]
            #nhash ==> node : tri   hashtable
            b1 = bs[1]+4
            b2 = bs[2]+4
            print "b1 = {}, b2 = {}".format(b1,b2)
            #numtriangles = len(nhash[b1]) # number of tri's using node b1
            # for triangle using node b1,
            #currenttriangle <= node b1 : tri num nhash[b1][z]
            
            """  ## debug:
            import numpy as np
            import matplotlib.pyplot as plt
            from make_neighbors import make_nbrs
            from search import search, tlm_search
            from circle_test import circle_test
            from circle_test import My_circle_test
            from circle_test import inCircleRobust
            
            nhash = [] 
            for n in range(self.nn+4):
                nhash.append([])
            
            
            for t in range(self.nt):
                for i in range(3):
                    n = self.tri[t][i]
                    nhash[n].append(t)
            
            
            
            bs=self.bs[1]
            
            b1 = bs[1]+4
            b2 = bs[2]+4
            
            
            n= 5
            m= 1
            z=3
            currenttriangle =  18
            b1= 5
            b2= 6
            self.currenttriangle = currenttriangle
            """
            
#            nhash = [] 
#            for n in range(self.nn+4):
#                nhash.append([])
#            
#            
#            for t in range(self.nt):
#                for i in range(3):
#                    n = self.tri[t][i]
#                    nhash[n].append(t)
            
            #for z, currenttriangle in enumerate(nhash[b1]):
            while len(nhash[b1])>0:
                currenttriangle =  nhash[b1].pop()
                #print 'currenttriangle = ',currenttriangle
                #z=currenttriangle
                # self.currenttriangle = currenttriangle
                # m = 0
                self.currenttriangle = currenttriangle
                #print ' currenttriangle = ',currenttriangle
                for m in range(3): #loop nodes of current triangle
                    #print m #index of current node
                    n = self.tri[currenttriangle][m]  #print ' current node = ',n
                    print '-----------------------------------'
                    print ' The current test triangle ={} , tri node m = {}, n = {} '.format(
                            currenttriangle,m,n)
                    print 'trying tri ',n
                    if b1==n: #if ray origin
                        print 'found a ray origin!',n
                        #for nbr_this in self.nbr[currenttriangle]:
                        #searchlist = [currenttriangle]+self.nbr[currenttriangle]
                        #for nbr_this in self.nbr[currenttriangle]:
                        #                        if currenttriangle == 16 and b1 == 5:
                        #                            print 'n=',n
                        #                            print 'm=',m
                        #                            print 'currenttriangle = ',currenttriangle
                        #                            print 'b1=',b1
                        #                            print 'b2=',b2
                        
                        ## get the test point - the unshared m  x,y coords                       
                        #self.xt=self.tempx[self.tri[m][q]]
                        #self.yt=self.tempy[self.tri[m][q]]
                        #xb4    =self.tempx[self.tri[nbr_this][q]]
                        #yb4    =self.tempy[self.tri[nbr_this][q]]
                        
                        bnode0 = m
                        bnode1=(m+1)%3
                        bnode2=(m+2)%3
                        
                        
                        
                        
                        # self.tri[self.currenttriangle] <=> nodes of curr triangle
                        #x,y of points on tri
                        xb1=self.tempx[self.tri[self.currenttriangle][bnode0] ]
                        xb2=self.tempx[self.tri[self.currenttriangle][bnode1]]
                        xb3=self.tempx[self.tri[self.currenttriangle][bnode2]]
                        #
                        xb4 = self.tempx[b2]
                        
                        
                        yb1=self.tempy[self.tri[self.currenttriangle][bnode0] ]
                        yb2=self.tempy[self.tri[self.currenttriangle][bnode1]]
                        yb3=self.tempy[self.tri[self.currenttriangle][bnode2]]
                        #
                        yb4 = self.tempy[b2]
                        
                        v1x=xb2-xb1
                        v1y=yb2-yb1
                        
                        v2x=xb3-xb1
                        v2y=yb3-yb1
                        
                        bxray=xb4-xb1
                        byray=yb4-yb1
                        
                        #get normals facing in
                        """
                        when both normals are +
                        you have trapped the ray
                        """
                        vec1 = np.asarray([xb1,yb1])
                        vec2 = np.asarray([xb2,yb2])
                        vec3 = np.asarray([xb3,yb3])
                        
                        #https://stackoverflow.com/questions/1243614/
                        #               how-do-i-calculate-the-normal-vector-of-a-line-segment
                        
                        vmag1 = np.linalg.norm(vec2-vec1)
                        xnormal1 = v1x/vmag1
                        ynormal1 = v1y/vmag1
                        #
                        dp1 = np.dot(np.asarray([ynormal1 , -xnormal1]),
                                     np.asarray([bxray,byray]))
                        
                        
                        vmag2 = np.linalg.norm(vec3-vec1)
                        xnormal2 = v2x/vmag2
                        ynormal2 = v2y/vmag2
                        #
                        dp2 = np.dot(np.asarray([ynormal2 , -xnormal2]),
                                     np.asarray([bxray,byray]))
                        #print 'dp1 = ',dp1
                        #print 'dp2 = ',dp2
                        
                        vt1 = vec2-vec1
                        vt2 = vec3-vec1
                        
                        ct1 = -cross2d(np.asarray([bxray,byray]),vt1)[2]
                        ct2 = cross2d(np.asarray([bxray,byray]),vt2)[2]
                        
#                        if currenttriangle==15:
#                            print("\nCheck area cur tri=15")
#                            print 'v1y = ',v1y
#                            print("bnode0 = {}".format(bnode0))
#                            print("bnode1 = {}".format(bnode1))
#                            print("bnode2 = {}".format(bnode2))
#                            print "vmag1 = {}, vmag2 = {}".format(vmag1, vmag2)
#                            print "xnormal1 = {}, ynormal1 = {}".format(xnormal1, ynormal1)
#                            print "byray = {}, bxray = {}".format(byray, bxray)
#                            print("dp1 = P{}, dp2 = {}".format(dp1,dp2))
                        
                        if currenttriangle==18 or \
                            currenttriangle==16 or \
                            currenttriangle==15 or \
                            currenttriangle==10 :
                            print("\nCheck area cur tri={}".format(currenttriangle))
                            print 'v1y = ',v1y
                            print("bnode0 = {}".format(bnode0))
                            print("bnode1 = {}".format(bnode1))
                            print("bnode2 = {}".format(bnode2))
                            print "vmag1 = {}, vmag2 = {}".format(vmag1, vmag2)
                            print "xnormal1 = {}, ynormal1 = {}".format(xnormal1, ynormal1)
                            print "byray = {}, bxray = {}".format(byray, bxray)
                            print("dp1 = P{}, dp2 = {}".format(dp1,dp2))
                        
                        #column 29
                        #if dp1>0 and dp2>0.:
                        #if dp1<0 and dp2<0.:
                        #if dp1>0 and dp2<0.:
                        if dp1<0. and dp2>0.:
                        #if ct1>0. and ct2>0.:
                            print 'trapped the ray!'
                            print '   tri = ',currenttriangle
                            raynbr = self.nbr[currenttriangle][bnode1]
                            if raynbr == -1:
                                pass #print 'DONOT flip :: -1 boundary edge ', 
                            if raynbr != -1: #-1 denotes a real edge in nbrs
                                #boundary flip needs to flip! 
                                #print 'flip tri boundary edge ',(currenttriangle,'--',raynbr)
                                #print 'curr tri = ', currenttriangle
                                #print 'neighbor tri = ', raynbr
                                #print ''
                                # (that is what makes it different)
                                testflip = np.zeros(6,int)
                                success = self.boundary_flip(
                                                        currenttriangle, 
                                                        raynbr, 
                                                        testflip)
                                #self.edge_flip(currenttriangle, 
                                #               raynbr, 
                                #               testflip)
                                #flipcount +=1
                                if success:
                                    print 'edge flipped!, curtri was ',currenttriangle
                                    print '\n'
                                    make_nbrs(self.nn+4,
                                              self.tdim, 
                                              self.tri, 
                                              self.nbr)
                                    
                                    
                                    nhash = []
                                    for n in range(self.nn+4):
                                        nhash.append([])
                                        
                                    for t in range(self.nt):
                                        for i in range(3):
                                            n = self.tri[t][i]
                                            nhash[n].append(t)
                                
                                    if self.dopics and success: 
                                        #print 'successful flip'
                                        self.plot_mesh(plotboundary=True,
                                                                    show=True)
                                    break
                                #elif success:
                                #    #print 'successful flip'
                                #else:
                                #    print 'failed flip'
                                #break
                            else:
                                #print 'did nothing to tri currenttriangle = ',currenttriangle
                                success = False
                        
                    elif b2==n: # found the node on the boundary.
                        #print 'boundary is on line'
                        #print 'curr tri = ', currenttriangle,'neighbor tri = ', raynbr
                        break
                    else:
                        pass
                
                #if success:
                #    break
            #if success:
            #    break
                

                            
        return
        
    
    
    def compute_optimized_nodes(self):
        ##"""Grab the ith triangle and put it on j"""
        #for i in range(self.nopt):
            #j = self.Optimize[i]
            #            if self.verbose:
            #                print 'i = ',i
            #                print 'Optimize list = ',self.Optimize
            #                print '  Optimization of triangle {}'.format(j)
        
        while self.Optimize:
            j = self.Optimize.pop()
            if self.verbose:
                print 'Optimize list = ',self.Optimize
                print '  Optimization of triangle {}'.format(j)
                
            #""" check seed, run circletest on seed triangle 
            # against seed triangles 3 neighbors"""
            #for m in self.nbr[j]:
            for k in range(3):
                #print 'new k, k = ',k
                m = self.nbr[j][k]
                if self.verbose:
                    print 'Checking triangle {}'.format(m)
                    
                    
                if m==-1:
                    if self.veryverbose:
                        print 'empty neighbor {}'.format(m)
                        print 'testing continue, k = ',k
                    continue 
                elif m!=-1:
                    #"""
                    #    If those neighbor's outnode points 
                    #        pass the circle test 
                    #    then leave them alone
                    #"""
                    if self.veryverbose:
                        print 'm = ',m 
                        
                    #'find this triangles
                    # node in triangle m that is '
                    # not in triangle j
                    # triangle j is the seed triangle
                    q=self.outnode(j,m)
                    if self.veryverbose:
                        print 'neighbor = ',m
                        print 'tri = ',j 
                        print 'testing outnode tri[{}][{}] = {}'.format(m,q,self.tri[m][q])
                    
                    # get the circle test points
                    x1=self.tempx[self.tri[j][0]]
                    x2=self.tempx[self.tri[j][1]]
                    x3=self.tempx[self.tri[j][2]]  
                    y1=self.tempy[self.tri[j][0]]
                    y2=self.tempy[self.tri[j][1]]
                    y3=self.tempy[self.tri[j][2]]
                    
                    ## get the test point - the unshared m  x,y coords                       
                    self.xt=self.tempx[self.tri[m][q]]
                    self.yt=self.tempy[self.tri[m][q]]
 
                    if self.veryverbose:
                        print 'calling the circle test'
                    flag=circle_test([x1,y1],[x2,y2],[x3,y3],
                                     [self.xt,self.yt])
                    #flag = My_circle_test([x1,y1],[x2,y2],[x3,y3],
                    #                      [self.xt,self.yt])
                    #flag = inCircleRobust(
                    #        [x1,y1],[x2,y2],[x3,y3],
                    #        p=np.asarray([self.xt,self.yt]))
                    
                    areaflag = False
                    if  flag:
                        if self.veryverbose:
                            print 'failed circle test'
                            print 'do edge flip'
                        testflip = np.zeros((6),int)
                        areaflag = self.edge_flip(j,m,testflip)
                    else:
                        if self.veryverbose:
                            print 'triangle {} passed the circle test'.format(j)
                            print 'testing continue, k = ',k
                        #continue
                    if areaflag:
                        if self.verbose:
                            print 'Flip Connectivity Accepted'
                            print 'Area Edge Flip Passed'.format(j)
                        make_nbrs(self.nn+4, self.tdim, self.tri, self.nbr)
                        self.Optimize.append(j)
                        self.Optimize.append(m)
                        self.nopt = len(self.Optimize)
                        continue
                    #elif areaflag is False:
                    else:
                        if self.veryverbose:
                            print 'Area flip failed'
        return
        

    
    def outnode(self, m, j):
        """
            find the node in triangle j 
            not shared with triangle m
        """
        #print 'outnode(m = {}, j={})'.format(m,j)
        #print ' currenttriangle = ',m
        #print 'nbr tri = ',j
        #ans = -1
        test1 = self.tri[j][0]
        test2 = self.tri[j][1]
        test3 = self.tri[j][2]
        
        if (test1 not in self.tri[m]):
            #ans=test1
            ans = 0
            #print 'outnode is correct in test1'
        elif (test2 not in self.tri[m]):
            #ans=test2
            ans = 1
            #print 'outnode is correct test2'
        else:
            #ans=test3
            ans = 2
            assert(test3 not in self.tri[m]),'bad answer outnode at outnode 2'
#            if (test3 not in self.tri[m]):
#                #print 'outnode is correct test3'
#                pass
#            else:
#                #print 'bad answer outnode is 2'
        if self.veryverbose:
            print 'find the node in tri[{}] = {}'.format(j,self.tri[j])
            print 'that is not in tri[{}] = {}'.format(m,self.tri[m])
            print 'found node = {} in tri[{}]'.format(ans, self.tri[j][ans])
        
        
        #print 'outnode of nbr = ',self.tri[j][ans]
        return ans
        
        
    def boundary_flip(self, j, m, testflip):
        """
        j = current triangle
        m = triangle to flip with
        """
        aa=self.outnode(j,m)  # index of node in tri[m] not in tri[j]
        bb=self.outnode(m,j)  # index of node in tri[j] not in tri[m]
        for i in range(3):
            testflip[i]     = self.tri[j][i]
            testflip[i+3]   = self.tri[m][i]
            
            ##testflip is not used!!
            
        #        p0=b
        #        p1=(b+1)%4
        #        p2=(b+2)%4
        #        p3=(b+3)%4
        
        switch = {0:[0,1,2],
                  1:[1,2,0],
                  2:[2,0,1]}
        case = switch[bb]
        a=case[0]
        b=case[1]
        d=case[2]
        
        #switch = {0:0,1:1,2:2}
        c = aa#switch[aa]
        #c = switch[aa]
        #c = bb
        
        a=self.tri[j][a]
        b=self.tri[j][b]
        d=self.tri[j][d]        
        c=self.tri[m][c]
        
        area1 = ((self.tempx[c]-self.tempx[b])*
                 (self.tempy[a]-self.tempy[b])-
                 (self.tempx[a]-self.tempx[b])*
                 (self.tempy[c]-self.tempy[b]))
                 
        area2 = ((self.tempx[a]-self.tempx[d])*
                 (self.tempy[c]-self.tempy[d])-
                 (self.tempx[c]-self.tempx[d])*
                 (self.tempy[a]-self.tempy[d]))
        
        #print 'area1 = ', area1
        #print 'area2 = ', area2
        if ( (area1>0.0) and (area2>0.0)):
            #print 'flipping boundary!, tri = {}, raynbr = {}'.format(j,m)
            self.tri[m][0]=c
            self.tri[m][1]=d
            self.tri[m][2]=a
        
            self.tri[j][0]=a
            self.tri[j][1]=b
            self.tri[j][2]=c
            return True
        
        else:
            #print 'cannot flip boundary'
            #print 'area1 = ',area1
            #print 'area2 = ',area2
            return False

    def edge_flip(self, j, m, testflip):
        aa=self.outnode(j,m)  # index of node in tri[m] not in tri[j]
        if self.verbose:
            print 'The node in m that is not in j'
            print ' is aa = {}, tri[{}][{}] = {}'.format(aa, m, aa,self.tri[m][aa])
        bb=self.outnode(m,j) # index of node in tri[j] not in tri[m]
        if self.verbose:
            print 'The node in j that is not in m'
            print ' is bb = {}, tri[{}][{}] = {}'.format(bb, j, bb,self.tri[j][bb])
        
        
        if self.tri[m][aa] == self.tri[j][bb]:
            print 'Total'
            print '    Logic'
            print '        Failure'
            print 'OutNodes Found Equal to Each Other'
            assert(False),'OutNodes Found Equal to Each Other'
        else:
            if self.veryverbose:
                print 'nodes ok' 
        #populate the testflip with tri and m
        for i in range(3):
            testflip[i]     = self.tri[j][i]
            testflip[i+3]   = self.tri[m][i]
            
            ##testflip is not used!!
        
        #index the seed nodes
#        p0 = b
#        p1 = (b+1)%4
#        p2 = (b+2)%4
#        p3 = (b+3)%4
        
        switch = {0:[0,1,2],
                  1:[1,2,0],
                  2:[2,0,1]}
        case = switch[bb]
        a=case[0]
        b=case[1]
        d=case[2]
        
        switch = {0:0,1:1,2:2}
        c = aa#switch[aa]
        
        a=self.tri[j][a]
        b=self.tri[j][b]
        d=self.tri[j][d]        
        c=self.tri[m][c]
        
        area1 = ((self.tempx[c]-self.tempx[b])*
                 (self.tempy[a]-self.tempy[b])-
                 (self.tempx[a]-self.tempx[b])*
                 (self.tempy[c]-self.tempy[b]))
                 
        area2 = ((self.tempx[a]-self.tempx[d])*
                 (self.tempy[c]-self.tempy[d])-
                 (self.tempx[c]-self.tempx[d])*
                 (self.tempy[a]-self.tempy[d]))
        
        #// Now do the Circle Test
        #// on the new J with outlier m
        x1=self.tempx[a]
        x2=self.tempx[b]
        x3=self.tempx[c]
    	  
        y1=self.tempy[a]
        y2=self.tempy[b]
        y3=self.tempy[c]
        #//then get the test point - unshared m
        xt=self.tempx[d]
        yt=self.tempy[d]
        
        cflagj = circle_test([x1,y1],[x2,y2],[x3,y3],
                             [xt,yt])
        #cflagj = inCircleRobust(
        #                    [x1,y1],[x2,y2],[x3,y3],
        #                    p=np.asarray([self.xt,self.yt]))
        #cflagj = My_circle_test([x1,y1],[x2,y2],[x3,y3],
        #                        [self.xt,self.yt])
        
        #// Now do the Circle Test
        #// on the new m with outlier j
        x1=self.tempx[c]
        x2=self.tempx[d]
        x3=self.tempx[a]
        
        y1=self.tempy[c]
        y2=self.tempy[d]
        y3=self.tempy[a]
        #//then get the test point - unshared j
        xt=self.tempx[b]
        yt=self.tempy[b]
        
        cflagm = circle_test([x1,y1],[x2,y2],[x3,y3],
                             [xt,yt])
        #cflagm = inCircleRobust(
        #                    [x1,y1],[x2,y2],[x3,y3],
        #                    p=np.asarray([self.xt,self.yt]))
        #cflagm = My_circle_test([x1,y1],[x2,y2],[x3,y3],
        #                        [self.xt,self.yt])
        
        areaflag = False
        if self.veryverbose:
            print 'area1 = ',area1
            print 'area2 = ',area2
            print 'ctestj = ',cflagj
            print 'ctestm = ',cflagm
        if (area1>0. and area2>0. and (not cflagj) and (not cflagm) ):
            areaflag=True
            if self.verbose:
                print 'Area Flip Passed'
            
            #self.plot_tris([self.tri[m],self.tri[j]])
            self.tri[m][0]=c
            self.tri[m][1]=d
            self.tri[m][2]=a
        
            self.tri[j][0]=a
            self.tri[j][1]=b
            self.tri[j][2]=c
        return areaflag
        

        
        
    def ccw(self, tri):
        """Tests whether the triangle is ccw"""
        A = tri[0]
        B = tri[1]
        C = tri[2]
        return (B[0] - A[0]) * (C[1] - A[1]) > (B[1] - A[1]) * (C[0] - A[0])
        
          
          
          

    def plot_tris(self, tris, 
                  pt=None, which = None,
                  plotboundary=False,
                  show=False):
        pts = self.plot_pts
        for i,el in enumerate(tris):
            #print el
            xpts = np.asarray([ pts[el[0],0],pts[el[1],0],pts[el[2],0], pts[el[0],0] ])
            ypts = np.asarray([ pts[el[0],1],pts[el[1],1],pts[el[2],1], pts[el[0],1] ])
            plt.plot(xpts,ypts )
            avgx =  xpts[:-1].sum()/len(xpts[:-1])
            avgy = ypts[:-1].sum()/len(ypts[:-1])
            #if self.verbose:
            plt.plot(pts[el[0],0],
                     pts[el[0],1],marker = 'o', color = 'black')
            plt.annotate(str(el[0]), 
                         xy=(pts[el[0],0],pts[el[0],1]))
            plt.plot(pts[el[1],0],
                     pts[el[1],1],marker = 'o', color = 'black' )
            plt.annotate(str(el[1]), 
                         xy=(pts[el[1],0],pts[el[1],1]))
            plt.plot(pts[el[2],0],
                     pts[el[2],1],marker = 'o', color = 'black' )
            plt.annotate(str(el[2]), 
                         xy=(pts[el[2],0],pts[el[2],1]))
            plt.plot(avgx,
                     avgy)
            plt.annotate('{}'.format(i), xy=(avgx,avgy))
        if which is not None:
            for i,el in enumerate(tris):
                #print el
                xpts = np.asarray([ pts[el[0],0],pts[el[1],0],pts[el[2],0], pts[el[0],0] ])
                ypts = np.asarray([ pts[el[0],1],pts[el[1],1],pts[el[2],1], pts[el[0],1] ])
                plt.plot(xpts,ypts )
                avgx =  xpts.sum()/len(xpts)
                avgy = ypts.sum()/len(ypts)
                #if self.verbose:
                plt.plot(pts[el[0],0],pts[el[0],1],marker = 'o', color = 'black'  )
                plt.annotate(str(el[0]), xy=(pts[el[0],0],pts[el[0],1]))#, xytext=(-1,-1))
                plt.plot(pts[el[1],0],pts[el[1],1],marker = 'o', color = 'black'  )
                plt.annotate(str(el[1]), xy=(pts[el[1],0],pts[el[1],1]))
                plt.plot(pts[el[2],0],pts[el[2],1],marker = 'o', color = 'black'  )
                plt.annotate(str(el[2]), xy=(pts[el[2],0],pts[el[2],1]))
                plt.plot(avgx,avgy)
                plt.annotate('T{}'.format(which), xy=(avgx,avgy))
                which +=1
        if pt is not None:
            plt.plot(pt[0],pt[1], marker = 'o', color = 'red' )
        
        if plotboundary:
            #todo: self.plot_boundary()
            for el in self.bs:
                p1 = self.pts[el[1]]
                p2 = self.pts[el[2]]
                plt.plot([p1[0],p2[0]],
                         [p1[1],p2[1]])
        if show: plt.show()
        #plt.savefig('tris/'+str(self.figcount)+'.png', format='png')
        self.figcount += 1
        return 
        
    def plot_trial(self):
        for el in self.tri:
            x1 = self.tempx[el[0]]
            y1 = self.tempy[el[0]]
            x2 = self.tempx[el[1]]
            y2 = self.tempy[el[1]]
            x3 = self.tempx[el[2]]
            y3 = self.tempy[el[2]]
            plt.plot([x1,x2,x3,x1],[y1,y2,y3,y1])
            plt.show()
        return 
        
    def plot_boundary(self):
        for el in self.bs:
            p1 = self.pts[el[1]]
            p2 = self.pts[el[2]]
            #print p1,p2
            plt.plot([p1[0],p2[0]],
                     [p1[1],p2[1]])
        #plt.show()
        return 

    def plot_mesh(self,
                  plotboundary=False,
                  show=False):
        self.plot_tris(self.tri[0:self.nt],
                       plotboundary=plotboundary,
                       show=show)
        plt.axis('equal')
        #plt.savefig('mesh/'+str(self.piccount)+'.png', format='png')
        self.piccount += 1
        return
    
    def set_plot_pts(self):
        xpts = self.tempx
        ypts = self.tempy
        pts = []
        for x,y in zip(xpts,ypts):
            pts.append([x,y])
        self.plot_pts = np.asarray(pts)
        return
        
    def print_neighbors(self):
        for t in range(self.nt):
            print 'element {} has neighbors {},{},{}'.format(t, 
                                                            self.nbr[t][0],
                                                            self.nbr[t][1],
                                                            self.nbr[t][2])
    def print_tri_nodes(self):
        """not really for use - since it shows duplicates"""
        for tri in self.tri:
            for node in tri:
                print self.tempx[node],self.tempy[node]
        return
    
    def test_search(self,pt):
        """
        seed = search(self.seed, pt[0], pt[1], 
                           self.tempx,self.tempy, 
                           self.tri,self.nbr)
        #"""
        seed = tlm_search(self.seed, pt[0], pt[1], 
                           self.tempx,self.tempy, 
                           self.tri,self.nbr)
        
        save = self.verbose
        self.verbose = True
        self.plot_tris(self.tri[seed:seed+1], pt, seed)
        #plt.plot(pt[0],pt[1],'ro' )
        self.verbose = save
        return seed
        
                                                            
if __name__ == "__main__":
    print "\n===================================================="
    print "\n     Driver for 2d unstructured meshing module      "
    print "\n===================================================="
    
    #fi = open(trimesh.dat)
    #lines = fi.readlines()
    from driver import mesh_driver
    #mesh_driver1 = mesh_driver(filename = r'test_files/test.pts')
    mesh_driver1 = mesh_driver(filename = r'test_files/btest.pts')
    #mesh_driver1 = mesh_driver(filename = r'test_files/small.pts')
    #mesh_driver1.read_input()
    #mesh_driver1.get_pts()
    #mesh_driver1.get_boundary()
    mesh = trimesh(mesh_driver1)
    mesh.insert_points()
    self = mesh
    self.plot_mesh()
    self.test_search([-.0,.3])
    self.test_search([-.2,.3])
    self.test_search([-.4,.3])