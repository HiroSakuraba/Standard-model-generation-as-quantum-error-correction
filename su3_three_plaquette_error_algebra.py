"""Exact 27D SU(3) plaquette error algebra via a finite-field rank witness.

All matrices lie in Q(sqrt(3)); entries have denominators powers of 3.
Specializing sqrt(3)->10 in F_97 is valid because 10^2=3 mod 97.
Rank over F_97 LOWER-BOUNDS rank over Q(sqrt(3)).
Matching explicit commutant lower bounds proves the exact dimensions.
This is NOT an approximate floating-point commutant calculation.
"""
import sympy as sp
from su3_three_plaquette_chain import operators, charge_conjugation

P=97
ROOT=10
assert ROOT*ROOT%P==3


def finite_field_entry(value):
    value=sp.expand(value)
    b=value.coeff(sp.sqrt(3))
    a=sp.expand(value-b*sp.sqrt(3))
    assert a.is_Rational and b.is_Rational,(value,a,b)
    def rational(x):
        return (int(x.p)*pow(int(x.q),-1,P))%P
    return (rational(a)+ROOT*rational(b))%P


def rank_of_commutator_equations(generators):
    """Sparse exact Gaussian elimination on [X,M]=0 over F_97.

    Variables are X_(a,b). A rank k minor surviving modulo P is
    nonzero in characteristic zero. All arithmetic here is exact.
    """
    dim=generators[0].rows
    basis={}
    for generator in generators:
        assert generator.shape==(dim,dim)
        columns=[[(k,finite_field_entry(generator[k,j])) for k in range(dim)
                  if generator[k,j]!=0] for j in range(dim)]
        rows=[[(k,finite_field_entry(generator[i,k])) for k in range(dim)
               if generator[i,k]!=0] for i in range(dim)]
        for i in range(dim):
            for j in range(dim):
                eq={}
                for k,v in columns[j]:
                    pos=i*dim+k
                    eq[pos]=(eq.get(pos,0)+v)%P
                for k,v in rows[i]:
                    pos=k*dim+j
                    eq[pos]=(eq.get(pos,0)-v)%P
                eq={pos:coeff for pos,coeff in eq.items() if coeff}
                while eq:
                    leading=min(eq)
                    prior=basis.get(leading)
                    if prior is None:
                        inv=pow(eq[leading],-1,P)
                        basis[leading]={k:(v*inv)%P for k,v in eq.items()}
                        break
                    factor=eq[leading]
                    for k,v in prior.items():
                        residue=(eq.get(k,0)-factor*v)%P
                        if residue:eq[k]=residue
                        else:eq.pop(k,None)
    return len(basis)


def test():
    states,D,W,B=operators(3)
    C=charge_conjugation(states)
    I=sp.eye(27)
    assert C!=I and C*C==I
    assert all(C*g==g*C for g in [D]+B)
    real_rank=rank_of_commutator_equations([D]+B)
    assert real_rank==27**2-2,(real_rank,'expected modular upper bound commutant 2')
    # I and C are independent exact rational commutant members; therefore
    # characteristic-zero nullity >=2. Modular rank gives nullity <=2.
    # Hence exact characteristic-zero commutant dimension =2.
    even=(I+C)/2;odd=(I-C)/2
    assert (even.rank(),odd.rank())==(14,13)
    full_rank=rank_of_commutator_equations([D]+W+[w.T for w in W])
    assert full_rank==27**2-1
    # I gives exact nullity >=1, modular rank gives <=1.
    print('FINITELY VERIFIED exact commutant real D+B0+B1+B2 = 2 (mod 97 rank 727; explicit I,C)')
    print('Hermitian error algebra = M_14(C) + M_13(C), dimension=365, no noiseless qubit')
    print('FINITELY VERIFIED exact commutant complex D+W_i+W_i^dagger = 1 (mod 97 rank 728; explicit I)')
    print('Full directed-loop error algebra = M_27(C), dimension=729')
    print('ALL EXACT MODULAR CERTIFICATE TESTS PASS')


if __name__=='__main__':
    test()
