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
#   Generalize to any genus g

@parallel
def Theta_char(cc,Omega,prec):
    Pi = pi.n(prec)
    eig = min((Omega.apply_map(lambda aux: aux.imag_part())).eigenvalues())
    dig = floor(RR(prec * log(2,10)))
    bound = ceil(RR(1./2 + sqrt(1./4 + dig*ln(10)/(pi.n(prec)*eig))))
    
    delta = vector(cc[:g]) 
    eps   = vector(cc[g:])
    Z = eps+Omega*delta
    
    s = 0
    for ml in cartesian_product_iterator([range(-bound,bound+1) for i in range(g)]):
        m  = vector(ml)
        s += exp(Pi*I*(2*m*Z + m*Omega*m)).n(prec)
    return exp(Pi*I*(2*delta*eps+ delta*Omega*delta)).n(prec)*s
       	
@parallel
def Theta_char_Pari(char, Z, prec, partial=False):
    if partial:
        TC = list(Partial_Theta_char_Pari([(char, Z, prec, comb) for comb in cartesian_product([["+","-"]]*3)]))
        return ComplexField(prec)(sum([val[1] for val in TC]))
    else:
        eig = min((Z.apply_map(lambda aux: aux.imag_part())).eigenvalues())
        
        dig = floor(RR(prec * log(2,10)))
        bound = ceil(RR(1./2 + sqrt(1./4 + dig*ln(10)/(pi.n(prec)*eig)))) 
        
        idx = 0
        fnc = ''
        exe = ''
        var_list = ''
        val_list = ''
        for i in range(g):
            for j in range(i,g):
                fnc += 'Z%s = %s;'%(idx,Z[i][j])
                var_list += 'z%s,'%idx
                val_list += 'Z%s,'%idx
                idx += 1
        for i in range(g):
            exe += 'D%s = %s;'%(i,char[i])
            exe += 'E%s = %s;'%(i,char[i+g])
            var_list += 'd%s,'%i
            val_list += 'D%s,'%i
            var_list += 'e%s,'%i
            val_list += 'E%s,'%i
        fnc += 'bound = %s;'%bound
        fnc += 'dig = %s;'%dig
        fnc += 'default(realprecision,dig); \
        thetaval('+var_list+'bd)= \
        {s=0; t=0; \
        cutoff = -dig*log(10.0)-2*log(bd);'
        for idx in range(g):
            fnc += 'for(i%s = -bd, bd, '%idx
        if g == 3:
            fnc += 't = Pi*I*(z0*i0^2 + 2*z1*i0*i1 + z3*i1^2 + 2*z2*i0*i2 + 2*z4*i1*i2 + z5*i2^2 + 2*z0*d0*i0 + 2*z1*d1*i0 + 2*z2*d2*i0 + 2*z1*d0*i1 + 2*z3*d1*i1 + 2*d2*z4*i1 + 2*z2*d0*i2 + 2*d1*z4*i2 + 2*d2*z5*i2 + z0*d0^2 + 2*z1*d0*d1 + z3*d1^2 + 2*z2*d0*d2 + 2*d1*d2*z4 + d2^2*z5 + 2*e0*i0 + 2*e1*i1 + 2*e2*i2 + 2*d0*e0 + 2*d1*e1 + 2*d2*e2);\
        if(real(t)>cutoff,s = s + exp(t))'+')'*g+';  return(s);}'
        if g == 6:
            fnc += 't = Pi*I*(i0^2*z0 + 2*i0*d0*z0 + d0^2*z0 + 2*i0*i1*z1 + 2*i1*d0*z1 + 2*i0*d1*z1 + 2*d0*d1*z1 + 2*i0*i2*z2 + 2*i2*d0*z2 + 2*i0*d2*z2 + 2*d0*d2*z2 + 2*i0*i3*z3 + 2*i3*d0*z3 + 2*i0*d3*z3 + 2*d0*d3*z3 + 2*i0*i4*z4 + 2*i4*d0*z4 + 2*i0*d4*z4 + 2*d0*d4*z4 + 2*i0*i5*z5 + 2*i5*d0*z5 + 2*i0*d5*z5 + 2*d0*d5*z5 + i1^2*z6 + 2*i1*d1*z6 + d1^2*z6 + 2*i1*i2*z7 + 2*i2*d1*z7 + 2*i1*d2*z7 + 2*d1*d2*z7 + 2*i1*i3*z8 + 2*i3*d1*z8 + 2*i1*d3*z8 + 2*d1*d3*z8 + 2*i1*i4*z9 + 2*i4*d1*z9 + 2*i1*d4*z9 + 2*d1*d4*z9 + 2*i1*i5*z10 + 2*i5*d1*z10 + 2*i1*d5*z10 + 2*d1*d5*z10 + i2^2*z11 + 2*i2*d2*z11 + d2^2*z11 + 2*i2*i3*z12 + 2*i3*d2*z12 + 2*i2*d3*z12 + 2*d2*d3*z12 + 2*i2*i4*z13 + 2*i4*d2*z13 + 2*i2*d4*z13 + 2*d2*d4*z13 + 2*i2*i5*z14 + 2*i5*d2*z14 + 2*i2*d5*z14 + 2*d2*d5*z14 + i3^2*z15 + 2*i3*d3*z15 + d3^2*z15 + 2*i3*i4*z16 + 2*i4*d3*z16 + 2*i3*d4*z16 + 2*d3*d4*z16 + 2*i3*i5*z17 + 2*i5*d3*z17 + 2*i3*d5*z17 + 2*d3*d5*z17 + i4^2*z18 + 2*i4*d4*z18 + d4^2*z18 + 2*i4*i5*z19 + 2*i5*d4*z19 + 2*i4*d5*z19 + 2*d4*d5*z19 + i5^2*z20 + 2*i5*d5*z20 + d5^2*z20 + 2*i0*e0 + 2*d0*e0 + 2*i1*e1 + 2*d1*e1 + 2*i2*e2 + 2*d2*e2 + 2*i3*e3 + 2*d3*e3 + 2*i4*e4 + 2*d4*e4 + 2*i5*e5 + 2*d5*e5);\
        if(real(t)>cutoff,s = s + exp(t))'+')'*g+';  return(s);}'
        gp(fnc)
        exe += 'TC = thetaval('+val_list+'bound);'
        gp(exe)
        theta = gp.eval('TC')
        return ComplexField(prec)(theta)
    
@parallel
def Partial_Theta_char_Pari(char,Z, prec, domain):
    eig = min((Z.apply_map(lambda aux: aux.imag_part())).eigenvalues())

    dig = floor(RR(prec * log(2,10)))
    bound = ceil(RR(1./2 + sqrt(1./4 + dig*ln(10)/(pi.n(prec)*eig)))) 

    idx = 0
    fnc = ''
    exe = ''
    var_list = ''
    val_list = ''
    for i in range(g):
        for j in range(i,g):
            fnc += 'Z%s = %s;'%(idx,Z[i][j])
            var_list += 'z%s,'%idx
            val_list += 'Z%s,'%idx
            idx += 1
    for i in range(g):
        exe += 'D%s = %s;'%(i,char[i])
        exe += 'E%s = %s;'%(i,char[i+g])
        var_list += 'd%s,'%i
        val_list += 'D%s,'%i
        var_list += 'e%s,'%i
        val_list += 'E%s,'%i
    fnc += 'bound = %s;'%bound
    for i in range(g):
        if i < 3:
            if domain[i] == '+':
                exe += 'B%s0 = %s; B%s1 = %s;'%(i,0,i,bound)
            else:
                exe += 'B%s0 = %s; B%s1 = %s;'%(i,-bound,i,-1)
        else:
            exe += 'B%s0 = %s; B%s1 = %s;'%(i,-bound,i,bound)
        var_list += 'b%s0, b%s1,'%(i,i)
        val_list += 'B%s0, B%s1,'%(i,i)
    fnc += 'dig = %s;'%dig
    fnc += 'default(realprecision,dig); \
    thetaval('+var_list+'bd)= \
    {s=0; t=0; \
    cutoff = -dig*log(10.0)-2*log(bd);'
    for idx in range(g):
        fnc += 'for(i%s = b%s0, b%s1, '%(idx, idx, idx)
    if g == 3:
        fnc += 't = Pi*I*(z0*i0^2 + 2*z1*i0*i1 + z3*i1^2 + 2*z2*i0*i2 + 2*z4*i1*i2 + z5*i2^2 + 2*z0*d0*i0 + 2*z1*d1*i0 + 2*z2*d2*i0 + 2*z1*d0*i1 + 2*z3*d1*i1 + 2*d2*z4*i1 + 2*z2*d0*i2 + 2*d1*z4*i2 + 2*d2*z5*i2 + z0*d0^2 + 2*z1*d0*d1 + z3*d1^2 + 2*z2*d0*d2 + 2*d1*d2*z4 + d2^2*z5 + 2*e0*i0 + 2*e1*i1 + 2*e2*i2 + 2*d0*e0 + 2*d1*e1 + 2*d2*e2);\
    if(real(t)>cutoff,s = s + exp(t))'+')'*g+';  return(s);}'
    if g == 6:
        fnc += 't = Pi*I*(i0^2*z0 + 2*i0*d0*z0 + d0^2*z0 + 2*i0*i1*z1 + 2*i1*d0*z1 + 2*i0*d1*z1 + 2*d0*d1*z1 + 2*i0*i2*z2 + 2*i2*d0*z2 + 2*i0*d2*z2 + 2*d0*d2*z2 + 2*i0*i3*z3 + 2*i3*d0*z3 + 2*i0*d3*z3 + 2*d0*d3*z3 + 2*i0*i4*z4 + 2*i4*d0*z4 + 2*i0*d4*z4 + 2*d0*d4*z4 + 2*i0*i5*z5 + 2*i5*d0*z5 + 2*i0*d5*z5 + 2*d0*d5*z5 + i1^2*z6 + 2*i1*d1*z6 + d1^2*z6 + 2*i1*i2*z7 + 2*i2*d1*z7 + 2*i1*d2*z7 + 2*d1*d2*z7 + 2*i1*i3*z8 + 2*i3*d1*z8 + 2*i1*d3*z8 + 2*d1*d3*z8 + 2*i1*i4*z9 + 2*i4*d1*z9 + 2*i1*d4*z9 + 2*d1*d4*z9 + 2*i1*i5*z10 + 2*i5*d1*z10 + 2*i1*d5*z10 + 2*d1*d5*z10 + i2^2*z11 + 2*i2*d2*z11 + d2^2*z11 + 2*i2*i3*z12 + 2*i3*d2*z12 + 2*i2*d3*z12 + 2*d2*d3*z12 + 2*i2*i4*z13 + 2*i4*d2*z13 + 2*i2*d4*z13 + 2*d2*d4*z13 + 2*i2*i5*z14 + 2*i5*d2*z14 + 2*i2*d5*z14 + 2*d2*d5*z14 + i3^2*z15 + 2*i3*d3*z15 + d3^2*z15 + 2*i3*i4*z16 + 2*i4*d3*z16 + 2*i3*d4*z16 + 2*d3*d4*z16 + 2*i3*i5*z17 + 2*i5*d3*z17 + 2*i3*d5*z17 + 2*d3*d5*z17 + i4^2*z18 + 2*i4*d4*z18 + d4^2*z18 + 2*i4*i5*z19 + 2*i5*d4*z19 + 2*i4*d5*z19 + 2*d4*d5*z19 + i5^2*z20 + 2*i5*d5*z20 + d5^2*z20 + 2*i0*e0 + 2*d0*e0 + 2*i1*e1 + 2*d1*e1 + 2*i2*e2 + 2*d2*e2 + 2*i3*e3 + 2*d3*e3 + 2*i4*e4 + 2*d4*e4 + 2*i5*e5 + 2*d5*e5);\
    if(real(t)>cutoff,s = s + exp(t))'+')'*g+';  return(s);}'
    gp(fnc)
    exe += 'TC = thetaval('+val_list+'bound);'
    gp(exe)
    theta = gp.eval('TC')
    return ComplexField(prec)(theta)
