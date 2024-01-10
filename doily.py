#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 14:44:15 2023

@author: arnaud
"""

I = [[1,0],[0,1]]
X = [[0,1],[1,0]]
Y = [[0,-1j],[1j,0]]
Z = [[1,0],[0,-1]]

def ProdTens(A,B):
    """produit tensoriel" des matrices A et B"""
    nA = len(A)
    pA = len(A[0])
    nB = len(B)
    pB = len(B[0])
    p = pA*pB
    n = nA*nB
    result = [ [ 0 for j in range(p)] for i in range(n)]
    for i in range(n) :
        for j in range(p) :
            result[i][j] = A[i//nA][j//pA] * B[i%nB][j%pB]
    return result

# phase = -1  1  -i  i

def sigma(x,y):
    """produit symplectique d'observables de N-qubits en binaire"""
    N = len(x) // 2
    result = 0
    for k in range(N):
        result ^= x[k]*y[N+k] ^ x[N+k]*y[k]
        #  * prioritaire sur ^
    return result
    
def commute(x,y):
    """ renvoie True si x et y commutent et False si x et y anticommutent """
    return sigma(x,y) == 0

def bin1(S):
    """code les matrices de Pauli en binaire"""
    match S:
        case 'I' :
            return [0,0]
        case 'X' :
            return [0,1]
        case 'Z' :
            return [1,0]
        case 'Y' :
            return [1,1]

def binN(S):
    """code un observable de N-qubit en binaire"""
    abscisses = []
    ordonnees = []
    for lettre in S:
        code = bin1(lettre)
        abscisses.append(code[0]) 
        ordonnees.append(code[1])
    return abscisses+ordonnees
        
# print(binN('XYZ')==[0,1,1,1,1,0])
# print(binN('ZIX')==[1,0,0,0,0,1])
# print(binN('YYY')==[1,1,1,1,1,1])

def bin2Pauli(L):
    """decode un observable de N-qubit code en binaire"""
    N = len(L)//2
    result = ''
    for k in range(N) :
        code = [L[k],L[k+N]]
        match code:
            case [0,0] :
                result += 'I' 
            case [0,1] :
                result += 'X' 
            case [1,0] :
                result += 'Z' 
            case [1,1] :
                result += 'Y' 
    return result

# print(commute(binN('XYZ'),binN('ZIX')))

def prodMatBin(O1,O2):
    """ produit matriciel de deux observables de Pauli sous forme binaire:
    on perd la phase"""
    return [ O1[k]^O2[k] for k in range(len(O1))]

# print(bin2Pauli(prodMatBin([0,1,1,1], [1,1,1,0]))) # XY . YZ

def isAline(O1,O2,O3): 
    """ teste si les 3 observables N-qubit codes en binaire 
    sont sur la meme ligne"""
    N = len(O1) // 2
    product = prodMatBin(prodMatBin(O1,O2),O3)
    return commute(O1,O2) and commute(O1,O3) and commute(O2,O3) and product == [0]*2*N

#print(isAline(binN('XYZ'),binN('ZIX'),binN('YYY')))

def isAOvoid(O1,O2,O3,O4,O5): 
    """ teste si les 5 observables N-qubit codes en binaire 
    forment un ovoide"""
    N = len(O1) // 2
    product = prodMatBin(prodMatBin(prodMatBin(prodMatBin(O1,O2),O3),O4),O5)
    return not commute(O1,O2) and not commute(O1,O3) and not commute(O2,O3) \
        and not commute(O1,O4) and not commute(O1,O5)\
            and not commute(O2,O4)  and not commute(O2,O5)\
                and not commute(O3,O4) and not commute(O3,O5)\
                    and not commute(O4,O5) and product == [0]*2*N

#print(isAOvoid(binN('IX'),binN('IZ'),binN('XY'),binN('ZY'),binN('YY')))
#print(isAOvoid(binN('IX'),binN('ZZ'),binN('IY'),binN('YZ'),binN('XZ')))

def isARoot(O1,O2,O3,O4,O5,c): 
    """ teste si les 5 observables N-qubit codes en binaire 
    forment une racine avec c"""
    s = sigma(O1,c)+sigma(O2,c)+sigma(O3,c)+sigma(O4,c)+sigma(O5,c)    
    return isAOvoid(O1, O2, O3, O4, O5) and s == 2

#print(isARoot(binN('IX'),binN('IZ'),binN('XY'),binN('ZY'),binN('YY'),binN('XI')))

#print(sigma(binN('IX'),binN('XI')), sigma(binN('IZ'),binN('XI')), \
#sigma(binN('XY'),binN('XI')), sigma(binN('ZY'),binN('XI')), \
#sigma(binN('YY'),binN('XI')))

#print(sigma(binN('IX'),binN('XX')), sigma(binN('IZ'),binN('XX')), \
#sigma(binN('XY'),binN('XX')), sigma(binN('ZY'),binN('XX')), \
#sigma(binN('YY'),binN('XX')))
 
#print(isARoot(binN('IX'),binN('ZZ'),binN('IY'),binN('YZ'),binN('XZ'),binN('IX')))

def OvoidSorted(L):
    """trie les 5 ovoides donnes dans l'ordre croissant"""
    sortie = []
    for k in range(6) :
        O = L[k] # k-ieme ovoid
        O = [ binN(q) for q in O ]
        O = sorted(O)
        #code = []
        #for q in O:
        #    code = code + q
        #sortie.append(code)
        sortie.append(O)
    #return sortie 
    sortie=sorted(sortie)
    sortie = [ [ bin2Pauli(q) for q in O] for O in sortie ]
    return sortie

O1 = ['ZZZ','YYY','YYX','XZZ','YIZ'] 
O2 = ['ZXY','ZZI','ZXX','YII','XZZ'] 
O3 = ['ZXY','IYY','XXY','IIZ','YYX'] 
O4 = ['IYX','ZZI','IYY','XZI','YIZ'] 
O5 = ['IYX','ZXX','XXX','IIZ','YYY'] 
O6 = ['XZI','XXY','XXX','YII','ZZZ'] 
L = OvoidSorted([O1,O2,O3,O4,O5,O6])

"""
[['XXX', 'IIZ', 'IYX', 'ZXX', 'YYY'],   c = XZZ
 ['XXX', 'XXY', 'XZI', 'YII', 'ZZZ'],   c = XZZ
 ['IIZ', 'XXY', 'IYY', 'ZXY', 'YYX'],   c = XZZ
 ['IYX', 'XZI', 'IYY', 'YIZ', 'ZZI'],   c = XZZ
 ['XZZ', 'ZXX', 'YII', 'ZXY', 'ZZI'],   c = IYY
 ['XZZ', 'YIZ', 'YYX', 'ZZZ', 'YYY']]   c = XXX
"""


def GenDoily(O,c):
    sortie = { 'o1' : O[0], 'o2' : O[1], 'o3' : O[2], 'o4' : O[3], 'o5' : O[4], 'c' : c }
    sortie['c.o1']= bin2Pauli(prodMatBin(binN(c),binN(sortie['o1'])))
    sortie['c.o2']= bin2Pauli(prodMatBin(binN(c),binN(sortie['o2'])))
    sortie['c.o3']= bin2Pauli(prodMatBin(binN(c),binN(sortie['o3'])))
    sortie['(c.o1).o4']= bin2Pauli(prodMatBin(binN(sortie['c.o1']),binN(sortie['o4'])))
    sortie['(c.o2).o4']= bin2Pauli(prodMatBin(binN(sortie['c.o2']),binN(sortie['o4'])))
    sortie['(c.o3).o4']= bin2Pauli(prodMatBin(binN(sortie['c.o3']),binN(sortie['o4'])))
    sortie['(c.o1).o5']= bin2Pauli(prodMatBin(binN(sortie['c.o1']),binN(sortie['o5'])))
    sortie['(c.o2).o5']= bin2Pauli(prodMatBin(binN(sortie['c.o2']),binN(sortie['o5'])))
    sortie['(c.o3).o5']= bin2Pauli(prodMatBin(binN(sortie['c.o3']),binN(sortie['o5'])))
    return sortie
    return set(sortie.values())

#print(sorted(GenDoily(L[0],'XZZ').values()) ==  sorted(GenDoily(L[1],'XZZ').values()))
#print(sorted(GenDoily(L[0],'XZZ').values()) ==  sorted(GenDoily(L[2],'XZZ').values()))
#print(sorted(GenDoily(L[0],'XZZ').values()) ==  sorted(GenDoily(L[3],'XZZ').values()))
#print(sorted(GenDoily(L[0],'XZZ').values()) ==  sorted(GenDoily(L[4],'IYY').values()))
#print(sorted(GenDoily(L[0],'XZZ').values()) ==  sorted(GenDoily(L[4],'XXX').values()))

print(GenDoily(L[0],'XZZ') ==  GenDoily(L[1],'XZZ'))
print(GenDoily(L[0],'XZZ') ==  GenDoily(L[2],'XZZ'))
print(GenDoily(L[0],'XZZ') ==  GenDoily(L[3],'XZZ'))
print(GenDoily(L[0],'XZZ') ==  GenDoily(L[4],'IYY'))
print(GenDoily(L[0],'XZZ') ==  GenDoily(L[5],'XXX'))