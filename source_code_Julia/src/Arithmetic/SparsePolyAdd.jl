"""
SparsePolyAdd: Sparse polynomial addition

- Julia version:
- Author: Mars Semenova
- Date: 2024-06-12

# Examples

```jldoctest
julia>
```

NEXT STEPS (TODO)
- MutableArithmetics pkg (gives interface for modifying objs in place)
 = should be making SMP (a mathematical obj) compatible w its interface
     - also MultivariatePolynomials interface
 - pass Base.+ or equivalent instead of passing an anon func
"""

import .Utils._truncate!

"""
Helper function to allow addition functionality to be extended for subtraction
by setting the operation between coefficients via an arg instead of negating the second
polynomial for better performance.

@param a - First polynomial operand. Must have the same coefficient type and
	number of variables as b.
@param b - Second polynomial operand. Must have the same coefficient type and
	number of variables as a.
@param op - Operation to apply to the coefficients of the passed polynomials
	when their terms' exps are the same.
@return Resulting SparsePolynomial{C,N}.
"""
function _mergeWithCoefOp(a::SparsePolynomial{C,N}, b::SparsePolynomial{C,N}, op::Function) where {C,N}
	_assertCanonical(a)
	_assertCanonical(b)

	if iszero(b)
        return SparsePolynomial(a)
	elseif iszero(a)
		return op(1, 1) == 0 ? negate(b) : SparsePolynomial(b)
    end

	asize = size(a)
    bsize = size(b)

    # init a new poly
	c = allocPoly(C, a.vars, asize + bsize)

	# add terms
    i = 1
    j = 1
    k = 1
    while (i <= asize && j <= bsize)
        if a.terms[i].exp < b.terms[j].exp # a < b
            c.terms[k] = Term{C,N}(b.terms[j])
			if op(1,1) == 0
				c.terms[k].coef = negateCoef(c.terms[k].coef)
			end
            k += 1
            j += 1
        elseif a.terms[i].exp == b.terms[j].exp # a = b
            newCoef = op(a.terms[i].coef, b.terms[j].coef)
            if newCoef != 0
    			cexp = a.terms[i].exp
                c.terms[k] = Term{C,N}(newCoef, a.terms[i].exp)
                k += 1
			end
            i += 1
            j += 1
        else # a > b
            c.terms[k] = Term{C,N}(a.terms[i])
            k += 1
            i += 1
        end
    end

	while i <= asize
		c.terms[k] = Term{C,N}(a.terms[i])
		k += 1
		i += 1
	end

	while j <= bsize
		c.terms[k] = Term{C,N}(b.terms[j])
		if op(1,1) == 0
			c.terms[k].coef = negateCoef(c.terms[k].coef)
		end
		k += 1
		j += 1
	end

	if k == 1
		return zero(SparsePolynomial{C,N})
	end

	_truncate!(c.terms, k-1)

	return c
end

"""
Function which performs the addition of 2 polynomials.

@param a - First polynomial operand. Must have the same coefficient type and
	number of variables as b.
@param b - Second polynomial operand. Must have the same coefficient type and
	number of variables as a.
@return Resulting SparsePolynomial{C,N}.
"""
function addPolynomials!(a::SparsePolynomial{C,N}, b::SparsePolynomial{C,N}) where {C,N}
	# TODO
	# TODO: replace implementation of addPolynomials w copy of inplace
end

function addPolynomials(a::SparsePolynomial{C,N}, b::SparsePolynomial{C,N}) where {C,N}
    return _mergeWithCoefOp(a, b, ((acoef,bcoef) -> acoef + bcoef))
end

"""
Function which performs the subtraction of 2 polynomials.

@param a - First polynomial operand. Must have the same coefficient type and
	number of variables as b.
@param b - Second polynomial operand. Must have the same coefficient type and
	number of variables as a.
@return Resulting SparsePolynomial{C,N}.
"""
function subPolynomials!(a::SparsePolynomial{C,N}, b::SparsePolynomial{C,N}) where {C,N}
	# TODO
	# TODO: replace implementation of subPolynomials w copy of inplace
end

function subPolynomials(a::SparsePolynomial{C,N}, b::SparsePolynomial{C,N}) where {C,N}
	return _mergeWithCoefOp(a, b, ((acoef,bcoef) -> acoef - bcoef))
end

# Base defs
import Base.+, Base.-
+(a::SparsePolynomial{C,N}, b::SparsePolynomial{C,N}) where {C,N} = addPolynomials(a,b)
-(a::SparsePolynomial{C,N}, b::SparsePolynomial{C,N}) where {C,N} = subPolynomials(a,b)