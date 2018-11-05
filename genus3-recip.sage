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

#TODO:
#   Add documentation

import warnings
warnings.simplefilter("ignore", UserWarning)
load('utilities.sage')
load('period_matrices.sage')
load('thetafunctions.sage')

from recip import *
from recip.period_matrices import _symplectic_basis, _big_period_matrix
Parallelism().set() #Uses as many parallel threads as possible

def mult_k(k, g2, g3, g4):
    assert k.is_square(), 'k must be a square'
    return (g2*k, g3*sqrt(k^3), g4*k^2)
        

def find_gs(j1, j2):
    g2 = j1
    g3 = j1^2
    g4 = j1^2*j2
    for g in [g2,g3,g4]:
        (g2, g3, g4) = mult_k(g.denom()*g.denom().squarefree_part(),g2,g3,g4)
        if all([g2 in ZZ, g3 in ZZ, g4 in ZZ]):
            break
    for fac in g2.factor():
        if fac[1] > 1:
            while all([a in ZZ for a in mult_k(fac[0]^(-2),g2,g3,g4)]):
                (g2, g3, g4) = mult_k(fac[0]^(-2),g2,g3,g4)
    return g2, g3, g4

def reduce_periodmatrix(Z): 
    bas = list(Z.basis())
    xi = Z.xi()
    CMtype = Z.CM_type()
    K = Z.CM_field()
    symplectic_basis = _symplectic_basis(bas, xi, conjugate, double_check=True)
    red_basis = reduce_Siegel_epsilon(vector(symplectic_basis),
        [phi.post_compose(K.embedding()) for phi in CMtype])[0].list()
    return PeriodMatrix_CM(CM_type = CMtype, basis = red_basis, xi = xi)
    
def rational_representation(Z,el):
    CMtype = Z._CM_type
    basis = Z._basis
    g = len(basis)/2
    K = CMtype.codomain()
    PI = _big_period_matrix(CMtype,basis)
    BigPi  = block_matrix(2,1,[PI, PI.conjugate()], subdivide=False)
    BigRho = block_diagonal_matrix(diagonal_matrix([phi(el) for phi in CMtype]),
        diagonal_matrix([phi(el) for phi in CMtype]).conjugate(), subdivide=False)
    return (BigPi.inverse()*BigRho*BigPi)
    
def compatible_CMtype(K,condition,z5,pw):
    lst = K.CM_types(K)
    for idx in range(len(lst)):
        Phi = lst[idx]
        Phiz5 = [pw(phi(z5)) for phi in Phi]
        Phiz5.sort()
        if Phiz5 == [pw(z5^cd) for cd in condition]:
            return Phi
    


#INPUT:
#   Totally real cubic field
poly = x^3 - x^2 -2*x + 1
#   Precision
prec = 600

Js = [[],[]]
g = 3

print "Precision = ",prec, '\n'

R.<x> = PolynomialRing(QQ)
K0.<t0> = NumberField(R(poly))
k.<z> = CyclotomicField(3)
[KK] = K0.composite_fields(k, 't')
K = CM_Field(KK.optimized_representation()[0])
roots = K.roots_of_unity()
zeta3 = roots[[el.complex_embedding(100).n() for el in  K.roots_of_unity()].index(exp(2*pi.n(200)*I/3).n())]
pw = zeta3.coordinates_in_terms_of_powers()

CC = ComplexField(prec)
K.embed(CC)
RCC.<X,Z> = PolynomialRing(CC,2)
Psi = K.complex_embeddings(200)[0]

CMtype = compatible_CMtype(K, [2,1,1], zeta3, pw)
print 'CM type on zeta_3'
for el in CMtype:
    print pw(el(zeta3))
print '\n'

print 'Field: '+str(K.polynomial())+'\n'+'Class number of K: '+str(K.class_number())+'\n'+'Clas number of K_+: '+str(K.real_field().class_number())+'\n'

for Omega_nonred in K.period_matrices_iter(CMtype, reduced=False):
    print Omega_nonred._ideal, Omega_nonred._xi
    Om = Omega_nonred.complex_matrix();
    print 'Eigenvalues of Im(Omega): '+str([vap.n() for vap in Om.apply_map(lambda i: i.imag_part()).eigenvalues()])+'\n'
    Omega = reduce_periodmatrix(reduce_periodmatrix(Omega_nonred))
    N = rational_representation(Omega, zeta3).transpose()
    Om = Omega.complex_matrix();
    print 'Eigenvalues of Im(Omega) after reducing: '+str([vap.n() for vap in Om.apply_map(imag_part).eigenvalues()])+'\n'

    ###Computing characteristics###
    A = N[:g,:g]
    B = N[:g,g:]
    C = N[g:,:g]
    D = N[g:,g:]

    ###Riemann theta constant###
    N2 = matrix(FiniteField(2), identity_matrix(2*g) - N.transpose().inverse())
    v = vector(FiniteField(2),(C*(D.transpose())).diagonal()+(A*(B.transpose())).diagonal())
    Delta = 1/2*vector(QQ,N2.solve_right(-v).list());
    print 'Riemann theta constant:\n'+str(Delta)+'\n'

    ###(1-rho)-torsion points###
    Ng = matrix(FiniteField(3), identity_matrix(2*g) - N.transpose().inverse())
    W = Ng.right_kernel().list()
    TorsionPoints = [1/3*vector(QQ, t.list()) for t in W]
    assert len(TorsionPoints) == 27, 'Too many or too few taus'

    ###Torsion points in the zero locus of theta###
    Zeros = [];
    for val in sorted(list(Theta_char_Pari([(Delta + TP, Om, 50) for TP in TorsionPoints]))):
        if abs(val[1]) < 1e-2:
            i = TorsionPoints.index(val[0][0][0] - Delta)
            if i != 0:
                Zeros.append(i)
    assert len(Zeros) == 14, 'Too many or too few zeros'
    print 'Zeros done'

    ###Choice of D_i###
    def are_li(x1,y1,z1):
        return rank(matrix(GF(3),[3*TorsionPoints[x1], 3*TorsionPoints[y1], 3*TorsionPoints[z1]])) == 3

    def find_ds():
        for a, b, c in Combinations(Zeros,3):
            forth = ((-TorsionPoints[a]-TorsionPoints[b]-TorsionPoints[c]).apply_map(lambda i : i-numerator(i)//denominator(i)))
            LI = are_li(a,b,c);
            if forth in TorsionZeros and LI:
                d = Zeros[TorsionZeros.index(forth)]
                assert are_li(d,b,c), 'Not li family'
                assert are_li(a,d,c), 'Not li family'
                assert are_li(a,b,d), 'Not li family'
                assert all([ el in ZZ for el in (TorsionPoints[a]+TorsionPoints[b]+TorsionPoints[c]+forth)]), 'sum not zero'
                for d0, d1, dlam, dmu in Permutations([a,b,c,d]):
                    D0   = TorsionPoints[d0]
                    D1   = TorsionPoints[d1]
                    Dlam = TorsionPoints[dlam]
                    Dmu  = TorsionPoints[dmu]
                    Div1 = 2*Dlam + D1
                    Div2 = 2*D1 + Dlam
                    Div3 = 2*Dmu + D1
                    Div4 = 2*D1 + Dmu
                    L = [Div1 - D0, Div2 - D0, Div3 - D0, Div4 - D0]
                    if all([l.apply_map(lambda i : i-numerator(i)//denominator(i)) not in TorsionZeros for l in L]):
                        return d0,d1,dlam,dmu
        raise TypeError('unsuccessful')

    TorsionZeros = [TorsionPoints[aux] for aux in Zeros]

    (D0,D1,Dlam, Dmu) = (TorsionPoints[d] for d in find_ds())

    Div1 = 2*Dlam + D1
    Div2 = 2*D1 + Dlam
    Div3 = 2*Dmu + D1
    Div4 = 2*D1 + Dmu
    
    DS = [Div1 - D0 - Delta, Div2 - D0 - Delta, Div3 - D0 - Delta, Div4 - D0 - Delta]
    TC = [-1 for i in range(8)]

    for val in sorted(list(Theta_char_Pari([(ch, Om, prec, True) for ch in DS]))):
            i = DS.index(val[0][0][0])
            if TC[i] == -1:
                TC[i] = val[1]


    LAM = (TC[0]/TC[1])^3*exp(6*pi*I*((DS[1]-DS[0])[:3]*D0[3:]+(-Delta + 2*D1 + Dlam)[:3]*(2*Delta - 3*D1 - 3*Dlam)[3:])).n(prec)
    MU = (TC[2]/TC[3])^3*exp(6*pi*I*((DS[2]-DS[0])[:3]*D0[3:]+(-Delta + 2*D1 + Dmu)[:3]*(2*Delta - 3*D1 - 3*Dmu)[3:])).n(prec)

    C.<X> = PolynomialRing(ComplexField(prec))
    P = X*(X-1)*(X-LAM)*(X-MU);
    P2=P.subs(X = X-P[3]/4);
    j1 = P2[1]^2/P2[2]^3;
    j2 = P2[0]/P2[2]^2;

    Js[0].append(j1)
    Js[1].append(j2)
    
if K.is_isomorphic(CyclotomicField(9)):
    print '(g2, g3, g4) = (0, -1, 0)'
else:     
    n = len(Js[0])
    CC = ComplexField(prec)
    Coeffs_Hj = [[CC((-1)^k*sum([prod(y) for y in Combinations(Js[i],k)])) for k in range(n+1)] for i in range(2)]
    HJs = map(lambda i : flatten(i.factor())[0:2*len(i.factor()):2], [sum([Coeffs_Hj[i][j].algdep(1).change_ring(QQ).roots()[0][0]*x^(n-j) for j in range(n+1)]) for i in range(2)])
    for case in range(len(HJs[0])):
        if HJs[0][case].degree() == 1:
            (j1, j2) = (HJs[0][case].roots()[0][0], HJs[1][case].roots()[0][0])
            (g2,g3,g4) = find_gs(j1,j2)
            print 'Case #'+str(case)+':\n (g2, g3, g4) = ('+', '.join([str(g2.factor()),str(g3.factor()),str(g4.factor())])+')'
        else:
            KJs = NumberField(HJs[0][case],'a')
            assert KJs.is_isomorphic(NumberField(HJs[1][case],'b'))
            if HJs[0][case].change_ring(K).is_irreducible():
                Fs = K.extension(HJs[0][case],'b').absolute_field('b').optimized_representation()[0].subfields(3)
                F.<g> = Fs[[fld[0].is_isomorphic(KJs) for fld in Fs].index(True)][0].change_names()
            else:
                F.<g> = K.subfields(3)[0][0].change_names()
            Js_alg = [HJs[i][case].change_ring(F).roots() for i in range(2)]
            (j1, j2) = (Js_alg[0][0][0], Js_alg[1][0][0])
            print 'Case #'+str(case)+':\n (j1, j2) = ('+str(j1)+', '+str(j2)+') where '+str(g.minpoly().subs(x = var('g')))
