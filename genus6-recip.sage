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
#   List of cubic fields
data = [x^3 - 3*x - 1,
        x^3 - x^2 -2*x + 1,
        x^3 - x^2 - 4*x - 1,
        x^3 - 12*x - 14]
#   Precision
prec = 100


results = []
for poly in data: 
    g = 6

    print "Precision = ",prec, '\n'

    R.<x> = PolynomialRing(QQ)
    K0.<t0> = NumberField(R(poly))
    k.<z> = CyclotomicField(5)
    [KK] = K0.composite_fields(k, 't')
    K = CM_Field(KK.optimized_representation()[0])
    roots = K.roots_of_unity()
    zeta5 = roots[[el.complex_embedding(100).n() for el in  K.roots_of_unity()].index(exp(2*pi.n(200)*I/5).n())]
    pw = zeta5.coordinates_in_terms_of_powers()
    
    CC = ComplexField(prec)
    K.embed(CC)
    RCC.<X,Z> = PolynomialRing(CC,2)
    Psi = K.complex_embeddings(200)[0]

    CMtype = compatible_CMtype(K, [3,2,2,1,1,1], zeta5, pw)
    print 'CM type on zeta_5'
    for el in CMtype:
        print pw(el(zeta5))
    print '\n'

    print 'Field: '+str(K.polynomial())+'\n'+'Class number of K: '+str(K.class_number())+'\n'+'Clas number of K_+: '+str(K.real_field().class_number())+'\n'

    for Omega_nonred in K.period_matrices_iter(CMtype, reduced=False):
        print Omega_nonred._ideal, Omega_nonred._xi
        Om = Omega_nonred.complex_matrix();
        print 'Eigenvalues of Im(Omega): '+str([vap.n() for vap in Om.apply_map(lambda i: i.imag_part()).eigenvalues()])+'\n'
        Omega = reduce_periodmatrix(reduce_periodmatrix(Omega_nonred))
        N = rational_representation(Omega, zeta5).transpose()
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
        Ng = matrix(FiniteField(5), identity_matrix(2*g) - N.transpose().inverse())
        W = Ng.right_kernel().list()
        TorsionPoints = [1/5*vector(QQ, t.list()) for t in W]
        assert len(TorsionPoints) == 125, 'Too many or too few taus'

        ###Torsion points in the zero locus of theta###
        Zeros = [];
        for val in sorted(list(Theta_char_Pari([(Delta + TP, Om, 25) for TP in TorsionPoints]))):
            if abs(val[1]) < 1e-2:
                i = TorsionPoints.index(val[0][0][0] - Delta)
                if i != 0:
                    Zeros.append(i)
        assert len(Zeros) == 100, 'Too many or too few zeros'
        print 'Zeros done'

        ###Choice of D_i###
        def are_li(x1,y1,z1):
            return rank(matrix(GF(5),[5*TorsionPoints[x1], 5*TorsionPoints[y1], 5*TorsionPoints[z1]])) == 3

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
                        Div1 = 3*Dlam + 2*D1 + Dmu
                        Div2 = 3*D1 + 2*Dlam + Dmu
                        Div3 = 3*Dmu + 2*D1 + Dlam
                        Div4 = 3*D1 + 2*Dmu + Dlam
                        L = [Div1, Div1 - D0, Div2, Div2 - D0, Div3, Div3 - D0, Div4, Div4 - D0]
                        if all([l.apply_map(lambda i : i-numerator(i)//denominator(i)) not in TorsionZeros for l in L]):
                            return d0,d1,dlam,dmu
            raise TypeError('unsuccessful')

        TorsionZeros = [TorsionPoints[aux] for aux in Zeros]

        (D0,D1,Dlam, Dmu) = (TorsionPoints[d] for d in find_ds())

        Div1 = 3*Dlam + 2*D1 + Dmu
        Div2 = 3*D1 + 2*Dlam + Dmu
        Div3 = 3*Dmu + 2*D1 + Dlam
        Div4 = 3*D1 + 2*Dmu + Dlam

        DS = [Div1 - D0 - Delta, Div1 - Delta, Div2 - D0 - Delta, Div2 - Delta, Div3 - D0 - Delta, Div3 - Delta, Div4 - D0 - Delta, Div4 - Delta]
        TC = [-1 for i in range(8)]

        for val in sorted(list(Theta_char_Pari([(ch, Om, prec, True) for ch in DS]))):
                i = DS.index(val[0][0][0])
                if TC[i] == -1:
                    TC[i] = val[1]

        lambda_apr = CC(exp(10*pi.n(prec)*I*(DS[0] - DS[2])[:g]*D0[g:]))*(TC[0]/TC[1]*TC[3]/TC[2])^5
        mu_apr = CC(exp(10*pi.n(prec)*I*(DS[4] - DS[6])[:g]*D0[g:]))*(TC[4]/TC[5]*TC[7]/TC[6])^5

        f = Z*X*(X-Z)*(X-CC(lambda_apr)*Z)*(X-CC(mu_apr)*Z)
        print poly, '\n', lambda_apr, mu_apr, '\n', invariants(canonical_form(RCC(f)))
        results.append([lambda_apr, mu_apr])
        print '--------------\n'
