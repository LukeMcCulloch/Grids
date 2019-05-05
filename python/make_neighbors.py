# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 23:22:04 2015

@author: lukemcculloch
"""

"""
Notes on the conversion from C:

The break statement, like in C, 
breaks out of the smallest enclosing 
for or while loop.

The continue statement, also borrowed from C, 
continues with the next iteration of the loop:
"""

def make_nbrs(nn, ntri, tri, nbr):
    #node to tri hash table
    """
        e.g.
            node 0 : is in triangles 0,1
            node 1 : is in triangles 0
            node 2 : is in triangles 0,1
            node 2 : is in triangles 1
    """
    nhash = {}
    for n in range(nn):
        nhash[n] = []#list()
        
    for t in range(ntri):
        for i in range(3):
            n = tri[t][i]
            nhash[n].append(t)
    
    #make neighbors
    for t in range(ntri):
        for i in range(3):
            nbr[t][i] = -1
            n0 = tri[t][i]
            n1 = tri[t][(i+1)%3]
            #for j in range(len(nhash[n0])):
            for el in nhash[n0]:
                #a = nhash[n0][j]
                a = el
                if a==t:
                    continue 
                if a in nhash[n1]:
                    nbr[t][i] = a
                    break
                
    return