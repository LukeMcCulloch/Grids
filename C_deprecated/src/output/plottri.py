import numpy as np
import matplotlib.pyplot as plt

from random import randint, seed

class MeshViewer(object):
    """Read the tri.mesh
    file output by the trimesher 
    and plot the triangles
    """
    def __init__(self, filename = r'tri.mesh',
                 verbose=False):
        self.cur_ln         = 0
        self.figcount       = 0
        self.piccount       = 0
        self.verbose = verbose
        print 'reading input'
        self.read_input(filename)
        print 'get pts'
        self.get_pts()
        print 'get blocks'
        self.get_blocks()
        print 'get tris'
        self.get_tris()
        #self.get_boundary()
        
        
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
            self.cur_ln = el
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
    
    def get_blocks(self):
        self.cur_ln +=1 
        cur_line = self.lines[self.cur_ln]
        lsplit = cur_line.split()
        if lsplit[0] == '#Number' and \
            lsplit[2] == 'blocks':
            self.cur_ln +=1
        self.nblocks = int( 
                    self.lines[self.cur_ln].split()[0] )
        return
    
    def get_tris(self):
        self.cur_ln +=1
        cur_line = self.lines[self.cur_ln]
        lsplit = cur_line.split()
        if lsplit[0] == '#Number' and \
            lsplit[2] == 'triangular':
            self.cur_ln +=1
        self.ntris = int( self.lines[self.cur_ln].split()[0] )
        self.nt = self.ntris
        
        self.tri = np.zeros((self.ntris+4,3),int)
        
        self.set_limits()
        self.create_superpoints()
        
        for i, el in enumerate(range(0,self.ntris)):
            self.cur_ln += 1
            self.tri[i,:] = [
                            int(el) for el in (
                                    self.lines[self.cur_ln].split()[0:3])
                            ]
        return
        
        
    def get_boundary(self):
        """TODO: fix this for tri.mesh output
        """
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
    
    
    
    def set_limits(self):
        """
        to match plotting for the moment
        fixme: remove supernodes from trimesh results
        and thus from the need to handle in plotting here.
        """
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
        """
        to match plotting for the moment
        fixme: remove supernodes from trimesh results
        and thus from the need to handle in plotting here.
        """
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
        
        self.nt += 2
        return
    
          
          

    
          

    def plot_tris(self, tris, 
                  pt=None, which = None,
                  plotboundary=False,
                  show=False):
        """
        tris = self.tri 
        pt=None
        which = None
        plotboundary=False
        show=False
        """
        print 'plot tris'
        self.set_plot_pts()
        pts = self.plot_pts
        mmark = " "
        if self.verbose:
            mmark = 'o'
        for i,el in enumerate(tris):
            #print el
            elx = el[0]-1
            ely = el[1]-1
            elz = el[2]-1
            xpts = np.asarray([ pts[elx,0],pts[ely,0],pts[elz,0], pts[elx,0] ])
            ypts = np.asarray([ pts[elx,1],pts[ely,1],pts[elz,1], pts[elx,1] ])
            plt.plot(xpts,ypts )
            avgx =  xpts[:-1].sum()/len(xpts[:-1])
            avgy = ypts[:-1].sum()/len(ypts[:-1])
            plt.plot(pts[elx,0],
                     pts[elx,1],marker = mmark, color = 'black')
            plt.plot(pts[ely,0],
                     pts[ely,1],marker = mmark, color = 'black' )
            plt.plot(pts[elz,0],
                     pts[elz,1],marker = mmark, color = 'black' )
            plt.plot(avgx,
                     avgy, color = 'black' )
            if self.verbose:
                plt.annotate(str(elx), 
                         xy=(pts[elx,0],pts[elx,1]))
                plt.annotate(str(ely), 
                         xy=(pts[ely,0],pts[ely,1]))
                plt.annotate(str(ely), 
                         xy=(pts[ely,0],pts[ely,1]))
                plt.annotate(str(elz), 
                         xy=(pts[elz,0],pts[elz,1]))
                plt.annotate('{}'.format(i), xy=(avgx,avgy))
        if which is not None:
            for i,el in enumerate(tris):
                #print el
                elx = el[0]-1
                ely = el[1]-1
                elz = el[2]-1
                xpts = np.asarray([pts[elx,0],
                                   pts[ely,0],
                                   pts[elz,0], 
                                   pts[elx,0] ])
                ypts = np.asarray([pts[elx,1],
                                   pts[ely,1],
                                   pts[elz,1], 
                                   pts[elx,1] ])
                plt.plot(xpts,ypts )
                avgx =  xpts.sum()/len(xpts)
                avgy = ypts.sum()/len(ypts)
                #if self.verbose:
                plt.plot(pts[elx,0],
                         pts[elx,1],marker = 'o', color = 'black'  )
                plt.annotate(str(elx), 
                             xy=(pts[elx,0],
                                 pts[elx,1]))#, xytext=(-1,-1))
                plt.plot(pts[ely,0],
                         pts[ely,1],marker = 'o', color = 'black'  )
                plt.annotate(str(ely), 
                             xy=(pts[ely,0],
                                 pts[ely,1]))
                plt.plot(pts[elz,0],
                         pts[elz,1],marker = 'o', color = 'black'  )
                plt.annotate(str(elz),
                             xy=(pts[elz,0],
                                 pts[elz,1]))
                plt.plot(avgx,avgy)
                if self.verbose:
                    plt.annotate('T{}'.format(which), xy=(avgx,avgy))
                which +=1
        if pt is not None:
            plt.plot(pt[0],pt[1], marker = 'o', color = 'black' )
        
        if plotboundary:
            print 'plot boundary'
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
        """
        plotboundary=False
        show=Flase
        """
        print 'plot mesh'
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
        
    
if __name__ == """__main__""":
    self = MeshViewer()
    self.plot_mesh()