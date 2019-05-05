# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 23:22:04 2015

@author: lukemcculloch
"""
import numpy as np

def ccw(A, B, C):
    """Tests whether the turn 
        formed by A, B, and C is ccw"""
    return (B[0] - A[0]) * (C[1] - A[1]) > (B[1] - A[1]) * (C[0] - A[0])
    
    
def My_circle_test(p1,p2,p3,p4):
    """
        is test point 4
        within the circle defined 
        by three points on the circumference?
        
        Return True is p4 is in the circle
        
        Return Flase if p4 is not in the circle
    """
    #print 'p1 = ',p1
    #print 'p2 = ',p2
    #print 'p3 = ',p3
    #print 'p4 = ',p4
    """
    p1 =  [-0.59999999999999998, -0.59999999999999998]
    p2 =  [0.59999999999999998, -0.59999999999999998]
    p3 =  [-0.5, -0.5]
    p4 =  [-0.59999999999999998, -0.59999999999999998]
    
    p1 = [x1,y1]
    p2 = [x2,y2]
    p3 = [x3,y3]
    p4 = [self.xt, self.yt]
    """
    
    # check that p1,p2,p3 are ccw
    ccw_check = ccw(p1,p2,p3)
    
    if ccw_check == True:
        #print 'points '
        #print '      are' 
        #print '          ccw'
        ssq1 = np.dot(p1,p1)
        ssq2 = np.dot(p2,p2)
        ssq3 = np.dot(p3,p3)
        
        a = np.linalg.det(np.asarray([[p1[0],p1[1],1.],
                                      [p2[0],p2[1],1.],
                                      [p3[0],p3[1],1.]]))
        
        bx = np.linalg.det(np.asarray([[ssq1,p1[1],1.],
                                       [ssq2,p2[1],1.],
                                       [ssq3,p3[1],1.]]))
        
        by = np.linalg.det(np.asarray([[ssq1,p1[0],p1[1]],
                                       [ssq2,p2[0],p2[1]],
                                       [ssq3,p3[0],p3[1]]]))
                           
        c = np.linalg.det(np.asarray([[ssq1,p1[0],p1[1]],
                                      [ssq2,p2[0],p2[1]],
                                      [ssq3,p3[0],p3[1]]]))
        
        xo =  bx/(2.*a)
        yo = -by/(2.*a)
        
        radius = np.sqrt((bx*bx)+(by*by)+(4.0*a*c))/(2.0*abs(a))
        
        if abs(a)<10.e-13:# or radius != radius:
            #print 'colinear points'
            xo=(p1[0]+p2[0]+p3[0])/3.0
            yo=(p1[1]+p2[1]+p3[1])/3.0
            radius=10.e20
        
        #print 'a={}, p4 = {}'.format([xo,yo], p4)
        dist = np.linalg.norm(np.asarray([xo,yo])-np.asarray(p4))
        #print 'radius = {}, dist = {}'.format(radius, dist)
        
        if dist >= radius:
            flag=False
        elif dist<radius:
            flag=True
        else:
            flag = None
        return flag
    else:
        #        print 'Points '
        #        print '       Are'
        #        print '           NOT'
        #        print '              CCW!!'
        return True
        
        
        



def circle_test_nah( p1,p2,p3,p4):#x1,  y1,  x2,  y2,  x3,  y3,  x4,  y4):
    """ use Cramer's rule to solve the system of eqns
         to figure out vars, 
         check out the Wikipedia entry on Cramer's rule
    """
    x1 = p1[0]
    y1 = p1[1]
    x2 = p2[0]
    y2 = p2[1]
    x3 = p3[0]
    y3 = p3[1]
    x4 = p3[0]
    y4 = p3[1]
    flag = 0
    meep = False
    
    d1 = x1*x1 + y1*y1
    d2 = x2*x2 + y2*y2
    d3 = x3*x3 + y3*y3
    
    a = x2-x1; b = y2-y1; c = x3-x1; d = y3-y1;
    e = 0.5*(d2-d1); f = 0.5*(d3-d1); 
    
    den = (a*d - b*c)
    if(abs(den) < 1e-12):
        den = 1e-20
        meep = True #//colinear
    px = (e*d - b*f)/den
    py = (a*f - e*c)/den
    
    r = (x1-px)*(x1-px) + (y1-py)*(y1-py) #i will test for distance, so no need to sqrt it
	
    dist = (px-x4)*(px-x4) + (py-y4)*(py-y4) #the distance between 4th point and the center

    diff = abs(dist-r)
    
    if(dist<r):
        flag = True #1
    if(diff < 1e-12):
        flag = False #0 #co-circular points; dont flip, you will have to flip again
    if(meep):
        flag = True #1
    
    return flag






        
def circle_test(a,b,c,d):
    ccw_check = ccw(a,b,c)
    #if ccw_check == True:
    #assert(ccw_check),'ERROR, coliniear nodes'
    #print 'ccw ok'
    dpx2 = d[0]**2
    dpy2 = d[1]**2
    det = np.linalg.det(np.asarray([[a[0]-d[0], a[1]-d[1], (a[0]**2-dpx2)+(a[1]**2-dpy2)],
                                    [b[0]-d[0], b[1]-d[1], (b[0]**2-dpx2)+(b[1]**2-dpy2)],
                                    [c[0]-d[0], c[1]-d[1], (c[0]**2-dpx2)+(c[1]**2-dpy2)]]))
    if (det>0.):
        return True
    elif det<0.:
        return False
    else:
        print 'colinear'
        return True
    


#def dot(a, b):
#    return a.x*b.x + a.y*b.y + a.z*b.z
#
#def cross(a, b):
#    return Point( a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x )

def circumcenter(self, tri):
    """Compute circumcenter and circumradius of a triangle in 2D.
    Uses an extension of the method described here:
    http://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html
    """
    pts = np.asarray([self.coords[v] for v in tri])
    pts2 = np.dot(pts, pts.T)
    A = np.bmat([[2 * pts2, [[1],
                             [1],
                             [1]]],
                  [[[1, 1, 1, 0]]]])

    b = np.hstack((np.sum(pts * pts, axis=1), [1]))
    x = np.linalg.solve(A, b)
    bary_coords = x[:-1]
    center = np.dot(bary_coords, pts)

    # radius = np.linalg.norm(pts[0] - center) # euclidean distance
    radius = np.sum(np.square(pts[0] - center))  # squared distance
    return (center, radius)

def inCircleFast(self, tri, p):
    """Check if point p is inside of precomputed circumcircle of tri.
    """
    center, radius = self.circles[tri]
    return np.sum(np.square(center - p)) <= radius


def inCircleRobust(#tri, 
                   a,b,c,
                   p):
    """Check if point p is inside of circumcircle around the triangle tri.
    This is a robust predicate, slower than compare distance to centers
    ref: http://www.cs.cmu.edu/~quake/robust.html
    """
    tri = [a,b,c]
    #m1 = np.asarray([self.coords[v] - p for v in tri])
    m1 = np.asarray([np.asarray(v) - p for v in tri])
    m2 = np.sum(np.square(m1), axis=1).reshape((3, 1))
    m = np.hstack((m1, m2))    # The 3x3 matrix to check
    return np.linalg.det(m) <= 0


#def IsInCircumcircleOf(self, T):
#
#    a = T.v[0] - T.v[2]
#    b = T.v[1] - T.v[2]
#
#    # Ref: https://en.wikipedia.org/wiki/Circumscribed_circle#Circumcircle_equations
#    z = cross(a,b)
#    p0 = cross(dot(a,a)*b-dot(b,b)*a, z)*(0.5/dot(z,z)) + T.v[2]
#
#    r2 = 0.25*dot(a, a)*dot(b,b)*dot(a-b, a-b)/dot(z, z)
#
#    #print "IsInC"
#    #print self, p0
#    #print sqrt(r2), "\n"
#    #print 
#
#    return dot(self-p0, self-p0) <= r2