# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 23:31:29 2015

@author: lukemcculloch
"""
import numpy as np

def search(seed, xt, yt, x, y, tri, nbr ):
    """
    seed    = self.seed
    xt,yt   = the point we want to find
    x,y     = the poins that make up triangle element nodes
    
    seed    = self.seed
    pts     = [.2,.3]
    xt      = pt[0]
    yt      = pt[1]
    x       = self.tempx
    y       = self.tempy    
    pts     = self.pts
    tri     = self.tri
    nbr     = self.nbr
    """

    count = 0
    test = 1.
    while test>0:
        tr = seed
        #nodes
        a = tri[seed][0]
        b = tri[seed][1]
        c = tri[seed][2]
        #edge vectors
        ax = x[b]-x[a]
        bx = x[c]-x[b]
        cx = x[a]-x[c]
        ay = y[b]-y[a]
        by = y[c]-y[b]
        cy = y[a]-y[c]
        #normalized edge vectors
        nrml = np.sqrt(ax*ax + ay*ay)#np.linalg.norm([ax,ay])
        ax = ax/nrml
        ay = ay/nrml   
        nrml = np.sqrt(bx*bx + by*by)#np.linalg.norm([bx,by])
        bx = bx/nrml
        by = by/nrml
        nrml = np.sqrt(cx*cx + cy*cy)#np.linalg.norm([cx,cy])
        cx = cx/nrml
        cy = cy/nrml
        #vector from edge midpts to new point
        tvax = xt - (x[b]+x[a])/2.
        tvay = yt - (y[b]+y[a])/2.
        tvbx = xt - (x[c]+x[b])/2.
        tvby = yt - (y[c]+y[b])/2.
        tvcx = xt - (x[a]+x[c])/2.
        tvcy = yt - (y[a]+y[c])/2.
        #normalized vector from edge mid points to new point
        nrml = np.sqrt(tvax*tvax+tvay*tvay)
        tvax = tvax/nrml
        tvay = tvay/nrml
        nrml = np.sqrt(tvbx*tvbx+tvby*tvby)
        tvbx = tvbx/nrml
        tvby = tvby/nrml
        nrml = np.sqrt(tvcx*tvcx+tvcy*tvcy)
        tvcx = tvcx/nrml
        tvcy = tvcy/nrml
        # dot normal vector from edge midpt to new point 
        #    with 
        # normalized "flipped"  edge vector, ie the normal vector to the edge
        dota = ay*tvax-ax*tvay
        dotb = by*tvbx-bx*tvby
        dotc = cy*tvcx-cx*tvcy
        #        dota = ax*tvax-ay*tvay
        #        dotb = bx*tvbx-by*tvby
        #        dotc = cx*tvcx-cy*tvcy
        if(dota>dotb):
            if(dota>dotc):
                test = dota;
                i = 0;
                seed = nbr[seed][0];
                #print 'seed = ',seed
            else:
                test = dotc;
                i=2;
                seed = nbr[seed][2];
                #print 'seed = ',seed
                
        else:
            if (dotb>dotc):
                test = dotb;
                i=1;
                seed = nbr[seed][1];
                #print 'seed = ',seed
            else:
                test = dotc;
                i=2;
                seed = nbr[seed][2];
                #print 'seed = ',seed
    #print 'tri = ',tr
    return tr
    


def sumsquares(x,y):
    return  x*x + y*y

def distances(x1,x2):
    return x2-x1    
    
def magnitudes(x0,y0,x1,y1):
    return np.sqrt(   (x0-x1)**2 + (y0-y1)**2   )
    
def tlm_search(seed, xt, yt, x, y, tri, nbr ):
    """/*find the triangle that contains the pt given*/
            
    seed    = self.seed
    pts     = [.2,.3]
    xt      = pt[0]
    yt      = pt[1]
    x       = self.tempx
    y       = self.tempy    
    pts     = self.pts
    tri     = self.tri
    nbr     = self.nbr
    """
#    maxold=0.0
#    maxnew=0.0
#    flag=0.0
    loop=1.
    
    while loop>0.:
        t=seed
        #        inflag=0
        #        outflag=0
        #count=0
        results = []
        for i in range(3):
            #count +=1
            n0 = tri[t][i]
            n1 = tri[t][(i+1)%3]

            vec = np.asarray([x[n1], y[n1]]) - np.asarray([x[n0], y[n0]])
            mag = np.linalg.norm(vec)
            
            xnormal = vec[0]/mag#( x[n1] - x[n0] )/mag
            ynormal = vec[1]/mag#( y[n1] - y[n0] )/mag
            
            midx = 0.5*( x[n1] + x[n0] )
            midy = 0.5*( y[n1] + y[n0] )
            
            xtv = xt - midx
            ytv = yt - midy
            
            tv_mag = np.linalg.norm([xtv,ytv])
            xtv_norm = xtv/tv_mag
            ytv_norm = ytv/tv_mag 
            
            results.append(
                np.dot( np.asarray([ynormal,-xnormal]) , 
                                  np.asarray([xtv_norm,ytv_norm]) )
                                  )

        if(results[0]> results[1]):
            if(results[0]>results[2]):
                loop = results[0]
                i = 0
                seed = nbr[seed][0]
            else:
                loop = results[2]
                i=2
                seed = nbr[seed][2]
                
        else:
            if (results[1]>results[2]):
                loop = results[1]
                i=1
                seed = nbr[seed][1]
            else:
                loop = results[2]
                i=2
                seed = nbr[seed][2] 
                
#            flag = np.dot(np.asarray([ynormal,xnormal]) , np.asarray([xtv_norm,ytv_norm]) )
            
#            if flag <= 0.0:
#                inflag +=1
#                if inflag > 2:
#                    seed = t
#                    loop = 0
#            elif flag > 0.0:
#                maxold = flag
#                maxnew = max(maxnew, flag)
#                outflag += 1
#                if maxnew >= maxold:
#                    t=nbr[t][i]
#                if count > 2:
#                    continue
    return seed

if __name__ == '__main__':
    print 'test to be implemented'
    bdim = 80
    
