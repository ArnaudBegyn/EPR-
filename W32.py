#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 14:44:15 2023

@author: arnaud
"""



def sigma(x,y):
    """produit symplectique sur F_2^N"""
    N = len(x) // 2
    result = 0
    for k in range(N):
        result ^= x[k]*y[N+k] ^ x[N+k]*y[k]
        #  * prioritaire sur ^ (=XOR)
    return result

def somme(x,y):
    """somme de deux vecteurs de F_2^N"""
    return [ x[k]^y[k] for k in range(len(x)) ]
    
def droites():
    """donne la liste de toutes les droites de F_2^4
    sous la forme u , v, u+v donnes en liste
    dans l'ordre lexicographique"""
    liste = []
    F = [ [i, j, k, l] for i in range(2)\
         for j in range(2) for k in range(2)\
             for l in range(2)\
             if (i, j, k, l) != (0, 0, 0, 0)]
    for u in F :
        for v in F :
            if v != u :
                w = somme(u, v)
                d = sorted([u, v , w]) 
                # ordre lexicographique
                if d not in liste :
                    liste.append(d)
    return liste
                  
def EstIsotrope(V):
    """teste si V est totalement isotrope"""
    for x in V:
        for y in V:
            if sigma(x, y) == 1 :
                return False
    return True
    
def droites_isotropes():
    """donne toutes les droites totalement isotropes
    de F_2^N"""
    liste = []
    D = droites()
    for d in D :
        if EstIsotrope(d) :
            liste.append(d)
    return liste

def eq(t, u) :
    a, b, c, d = t
    x, y, z, t = u
    return a*x^b*y^c*z^d*t

def plans():
    """donne la liste de touts les plans de F_2^4"""
    F = [ [i, j, k, l] for i in range(2)\
         for j in range(2) for k in range(2) for l in range(2)\
             if (i, j, k, l) != (0, 0, 0, 0)]
    liste = []
    for a in range(2) :
        for b in range(2) :
            for c in range(2) :
                for d in range(2) :
                    t = (a, b, c, d)
                    if t != (0, 0, 0, 0) :
                        P = []
                        for u in F :    
                            if eq(t, u) == 0 :
                                P.append(u)
                        liste.append(sorted(P))
    return liste

def plans_isotropes():
    """donne tous les plans totalement isotropes de F_2^4"""
    liste = []
    P = plans()
    for p in P :
        if EstIsotrope(p) :
            liste.append(p)
    return liste

                    
                        
    
    
    