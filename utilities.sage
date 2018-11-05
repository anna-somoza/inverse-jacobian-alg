#*****************************************************************************
# Copyright (C) 2018 Anna Somoza <anna.somoza.henares@gmail.com>
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#*****************************************************************************

#### SAGE UTILITIES

#TODO:
#   Add documentation

##Functions to transform an anti-hermitian matrix:

def transformation(E,X,i,j,k,lst=None):
    E = E.with_added_multiple_of_row(i,j,k);
    E = E.with_added_multiple_of_column(i,j,k.conjugate());
    X = X.with_added_multiple_of_row(i,j,k);
    if lst != None:
        lst.append(['T',i,j,k])
        return E, X, lst
    return E, X
def swap(E,X,i,j,lst=None):
    E.swap_rows(i,j);
    E.swap_columns(i,j);
    X.swap_rows(i,j);
    if lst != None:
        lst.append(['S',i,j])
        return E, X, lst
    return E, X
def rescale(E,X,i,k,lst=None):
    E = E.with_rescaled_row(i,k)
    E = E.with_rescaled_col(i,k.conjugate())
    X = X.with_rescaled_row(i,k)
    if lst != None:
        lst.append(['R',i,k])
        return E, X, lst
    return E,X
def chain_of_trans(E,X,lst):
    for L in lst:
        if L[0] == 'T':
            E,X = transformation(E,X,L[1],L[2],L[3])
        elif L[0] == 'S':
            E,X = swap(E,X,L[1],L[2])
        elif L[0] == 'R':
            E,X = rescale(E,X,L[1],L[2])
        else:
            raise TypeError('Char should be T, S or R')
    return E, X


##Functions to embed
def embed(x, z, k=None):
    l = fnc(x)
    xC = sum([l[i]*z^i for i in range(4)])
    if k!=None:
        assert abs(xC.n(100) - K.embeddings(ComplexField(100))[k](x)) < 1e-3
    return xC.expand()
    

##Functions in relation to binary quintics
#Reference:
#Théorie des formes binaires (by Francesco Faà di Bruno, 1825-1888)
def canonical_form(P):
    assert P.is_homogeneous(), 'It is not a binary quintic'
    (x, y) = P.parent().gens()
    CC = P.base_ring()
    R.<X> = PolynomialRing(CC)
    
    C, M = (P.coefficients(), P.monomials())
    V = [0 for i in range(6)]
    for i in range(6):
        try:
            k = M.index(y^i*x^(5-i))
            V[i] = C[k]/binomial(5,i)
        except ValueError:
            V[i] = 0
    N = matrix(R,4,[[1, X, X^2, X^3], V[:4], V[1:5], V[2:6]])
    try:
    	K.<g> = N.determinant().splitting_field()
    	(alpha, beta, gamma) = [a[0] for a in N.determinant().base_extend(K).roots()]
    	L = matrix(K, 3, [[z^i for z in [alpha, beta, gamma]] for i in range(3)])
    	l = vector(K,3,V[:3])
    except NotImplementedError:
    	(alpha, beta, gamma) = N.determinant().complex_roots()
    	L = matrix(CC, 3, [[z^i for z in [alpha, beta, gamma]] for i in range(3)])
    	l = vector(CC,3,V[:3])
    (p,q,r) = L\l
    (l,m,n) = (p/(gamma-beta)^5, q/(alpha-gamma)^5, r/(alpha-beta)^5)
    return l,m,n
    
def invariants((l,m,n)):
    I4 = (l*m + l*n + m*n)^2 - 4*l*m*n^2
    I8 = l^2*m^2*n^2*(l*m-l*n-m*n)
    I12 = l^4*m^4*n^4
    I18 = 4*l^5*m^5*n^5*(l-m)*(l+n)*(m+n)
    return I4, I8, I12, I18
