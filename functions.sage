# discQ(D) : integer -> boolean
# return true if D is a discriminant of an order, false otherwise
def discQ(D):
    return ((D%4==0 or D%4==1) and not(is_square(D)))

# fund_discQ(D) : integer -> boolean
# return true if D is a discriminant of a maximal quadratic order, false otherwise
def fund_discQ(D):
    return ((D%4==0 and is_squarefree(D/4)) or (D%4==1 and is_squarefree(D))) and not(is_square(D))

# genus(D) : discriminant -> integer
# return the order of the genus group,
# or equivalently, the order of the 2-torsion subgroup of the narrow class group
def genus(D):
    r=len(prime_factors(odd_part(D)))
    if D%4==1 or D%16==4:
        return 2^(r-1)
    if D%16==8 or D%16==12 or D%32==16:
        return 2^r
    if D%32==0:
        return 2^(r+1)

# genus_slow(D) : discriminant -> integer
# return the order of the genus group,
# or equivalently, the order of the 2-torsion subgroup of the narrow class group
def genus_slow(D):
    pos_divs=divisors(D)
    neg_divs=[-d for d in pos_divs]
    splits={1}
    for d in (pos_divs+neg_divs):
        if fund_discQ(d) and discQ(D/d):
            splits.add(min(squarefree_part(d),squarefree_part(D/d)))
    return len(splits)

# wide_genus_slow(D) : discriminant -> integer
# returns the order of the 2-torsion subgroup of the wide class group
def wide_genus_slow(D):
    pos_divs=divisors(D)
    splits={1}
    for d in pos_divs:
        if fund_discQ(d) and discQ(D/d):
            splits.add(min(squarefree_part(d),squarefree_part(D/d)))
    return len(splits)

# fd(D) : discriminant -> list
# returns [f,D0] where D0 is a fundamental discriminant and D=f^2*D0
def fd(D):
    s=squarefree_part(D)
    if (s%4==0 or s%4==1):
        return [(D/s)^(1/2),s]
    else:
        return [(D/(4*s))^(1/2),4*s]

# integerQ(n) : real -> boolean
# return true n is an integer, false otherwise
def integerQ(n):
    return n==floor(n)

# negpell(D) : discriminant -> boolean
# return true if the order of discriminant D contains a unit of negative norm, false otherwise
def negpell(D):
    K=QQ[(D+D^(1/2))/2]
    G=K.unit_group()
    u=G.gens_values()[1]
    uconj=u.trace()-u
    n=1
    while true:
        if integerQ((u^n-uconj^n)/D^(1/2)):
            return norm(u^n)==-1
        n+=1

# wide_genus(D) : discriminant -> integer
# returns the order of the 2-torsion subgroup of the wide class group
def wide_genus(D):
    # To do: Write a faster implementation using the factorization of D rather than iterating through the divisiors
    return wide_genus_slow(D)

# split_ntwQ(D) : discriminant -> boolean
# returns True if the exact sequence 1 -> ker -> Cl^+(O_D) -> Cl(O_D) -> 1 is split, and False otherwise
def split_ntwQ(D):
    ng=genus(D)
    wg=wide_genus(D)
    if negpell(D):
        return True
    elif ng==wg:
        return False
    else:
        return True