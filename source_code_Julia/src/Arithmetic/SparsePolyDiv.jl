#=
SparsePolyDiv: Sparse polynomial division

- Julia version:
- Author: Mars Semenova
- Date: 2024-06-28

# Examples

```jldoctest
julia>
```

NEXT STEPS (TODO)
-
=#

import .Utils._truncate!

"""
Function which divides a polynomial by a scalar.

@param p - Polynomial.
@param s - Scalar. Must be the same type as p's coefficients.
@return Resulting SparsePolynomial{C,N}.
"""
function divPolyByScalar!(p::SparsePolynomial{C,N}, s::C) where {C,N}
	# TODO
    # TODO: replace implementation of divPolyByScalar w copy of inplace
end

function divPolyByScalar(p::SparsePolynomial{C,N}, s::C) where {C,N}
	return _iterWithScalarOp(p, s, ((c,x) -> c/x))
end

"""
Function which tests whether a coefficient is 0. Uses Julia's
built-in eps() function as the margin of error for floats.

@param c - Coefficient of type C.
@return Whether c is 0 or is approximately 0.
"""
function isCoefZero(c::C)::Bool where {C}
	return c == C(0)
end

function isCoefZero(c::C)::Bool where {C<:AbstractFloat}
	return isapprox(c, C(0))
end

"""
Function which divides a polynomial by a single term.

@param p - SparsePolynomial to divide. Must have the same coefficient type and
	number of variables as t.
@param t - Term to divide p by. Must have the same coefficient type and
	number of variables as p.
@return The resulting quotient and remainder, in that order.
"""
function divideBySingleTerm!(p::SparsePolynomial{C,N}, t::Term{C,N}) where {C,N}
	# TODO
	# TODO: replace implementation of divideBySingleTerm w copy of inplace
end

function divideBySingleTerm(p::SparsePolynomial{C,N}, t::Term{C,N}) where {C,N}
	_assertCanonical(p)

	# TODO: add iszero(t) check when it's implemented

	if iszero(p)
		return zero(p), zero(p)
    end

	if iszero(t.exp) # div by scalar
		return divPolyByScalar(p, t.coef), zero(p)
	end

	if length(p) == 1 && p.terms[1] == t
		return SparsePolynomial([Term{C,N}(C(1))], p.vars), zero(p)
	end

	maxSize = length(p) + 1
	q = allocPoly(C, p.vars, maxSize)
	r = allocPoly(C, p.vars, maxSize)

	curP = 1
	curQ = 1
	curR = 1
	while curP <= length(p)
		if isDivisible(p.terms[curP].exp, t.exp)
			coef = p.terms[curP].coef / t.coef
			exp = p.terms[curP].exp - t.exp
			if !isCoefZero(coef)
				q.terms[curQ] = Term{C,N}(coef, exp)
				curQ += 1
			end
		else
			r.terms[curR] = Term{C,N}(p.terms[curP])
			curR += 1
		end
		curP += 1
	end

	if curQ == 1
		q = zero(p)
	else
		_truncate!(q.terms, curQ-1)
	end

	if curR == 1
		r = zero(p)
	else
		_truncate!(r.terms, curR-1)
	end

	return q, r
end

# TODO: make inplace
function _reallocTerms(terms::Vector{Term{C,N}}, newSize::Int) where {C,N}
	newTerms = Vector{Term{C,N}}(undef, newSize)

	for x = 1:length(terms)
		newTerms[x] = Term{C,N}(terms[x]) # TODO: need to cpy?
	end

	return newTerms
end

"""
Helper method which computes the coefficient of the next term of the quotient or remainder.

@param h - ProductHeap used to computer quotient*b.
@param q - Quotient polynomial.
@param b - Divisor polynomial.
@return Coefficient of type C of the next term of the quotient or remainder,
"""
function _divisionGetNextTerm(h::ProductHeap, q::SparsePolynomial{C,N}, b::SparsePolynomial{C,N})::C where {C,N}
	if h.heapSize == 0
		throw(ArgumentError("Heap is empty (unexpected program state)."))
	end

	coef = C(0)
	lastB = h.lastB
	insertChain = nothing
	maxTerm = nothing
	nextMaxTerm = nothing
	nextDegs = _prodHeapPeek(h)
	maxDegs = ExpVec(nextDegs)

	while nextDegs == maxDegs
		insertChain = maxTerm
		maxTerm = _prodHeapRemove!(h)

		while !isnothing(maxTerm)
			coef += q.terms[maxTerm.aidx].coef * b.terms[maxTerm.bidx].coef

			nextMaxTerm = maxTerm.next
			if maxTerm.bidx != lastB
				maxTerm.bidx += 1
				maxTerm.next = insertChain
				insertChain = maxTerm
			end

			maxTerm = nextMaxTerm
		end

		maxTerm = insertChain
		nextDegs = _prodHeapPeek(h)
	end

	# insert all successors back to heap
	while !isnothing(maxTerm)
		nextMaxTerm = maxTerm.next
		maxTerm.next = nothing
		_prodHeapInsert!(h, maxTerm, q.terms[maxTerm.aidx].exp + b.terms[maxTerm.bidx].exp)
		maxTerm = nextMaxTerm
	end

	return coef
end

"""
Function which performs the division of 2 polynomials.

@param a - First polynomial operand. Must have the same coefficient type and
	number of variables as b.
@param b - Second polynomial operand. Must have the same coefficient type and
	number of variables as a.
@return Resulting SparsePolynomial{C,N}.
"""
function dividePolynomials!(a::SparsePolynomial{C,N}, b::SparsePolynomial{C,N}) where {C,N}
	# TODO
	# TODO: replace implementation of dividePolynomials w copy of inplace
end

"""
Convert a number to 0 if it is very close.
If abs(c) < tolerance, then return 0, otherwise return c.

@param c - the number to evaluate
@param tolerance - the maximum value which should be considered as non-zero.
@return 0 if abs(c) < tolerance, otherwise c.
"""
function _fuzzyZero(c::C, tolerance=zero(C)::C) where {C<:Number}
	if abs(c) < tolerance
		return zero(C)
	else
		return c
	end
end

"""
Remove coefficients from a polynomial if they are very close.
For each coefficient c, if abs(c) < tolerance, then remove that term.

@param a - the polynomial to evaluate
@param tolerance - the maximum value which should be considered as non-zero.
@return a polynomial equal to a but with coefficients less than tolerance removed.
"""
function _fuzzyZero(a::SparsePolynomial{C,N}, tolerance=zero(C)::C) where {C<:Number,N}
	if (iszero(tolerance))
		return a
	end

	fuzzyTerms = Vector{Term{C,N}}(undef, length(a))
	idx = 1
	for term in a.terms
		c = _fuzzyZero(term.coef)
		if !iszero(c)
			fuzzyTerms[idx] = Term(c, term.exp)
			idx += 1
		end
	end

	if (idx == 1)
		return zero(a)
	else
		resize!(fuzzyTerms, idx)
		return SparsePolynomial(fuzzyTerms, a.vars)
	end
end


function dividePolynomials(a::SparsePolynomial{C,N}, b::SparsePolynomial{C,N}) where {C,N}
	_assertCanonical(a)
	_assertCanonical(b)

	if iszero(b) # assert not div by 0 or ret 0, 0 if diving 0 poly
        throw(DivideError())
	elseif iszero(a)
		return zero(a), zero(a)
    end



	if a == b
		return SparsePolynomial([Term{C,N}(C(1))], a.vars), zero(a)
	end

	if length(b) == 1
		return divideBySingleTerm(a, b.terms[1])
	end

	maxQSize = length(a) < 5 ? 5 : length(a)
	maxRSize = maxQSize
	leadB = b.terms[1]
	q = allocPoly(C, a.vars, maxQSize)
	r = allocPoly(C, a.vars, maxRSize)
	leadAInd = 1
	curR = 1

	# init q w lt(a)/lt(b)
	while leadAInd <= length(a) && !isDivisible(a.terms[leadAInd].exp, leadB.exp)
		r.terms[curR] = Term{C,N}(a.terms[leadAInd]) # a.terms >= lt(b)
		curR += 1

		if curR > maxRSize
			maxRSize <<= 1
			r.terms = _reallocTerms(r.terms, maxRSize)
		end

		leadAInd += 1
	end

	if leadAInd > length(a) # no div needed
		return zero(a), SparsePolynomial(a)
	end

	curQ = 1
	coef = a.terms[leadAInd].coef / leadB.coef
	exp = a.terms[leadAInd].exp - leadB.exp
	q.terms[curQ] = Term{C,N}(coef, exp) # lt(a)/lt(b)
	curA = leadAInd + 1

	# init prodHeap
	h = _prodHeapInit(q, b)
	_prodHeapGrow!(h, maxQSize)
	h.lastB = length(b)
	_prodHeapInsert!(h, q, curQ, b, 2)
	curQB = _prodHeapPeek(h)
	curQ += 1

	# mult between q (quotient) + b (divisor)
	comp = 0
	coef = 0
	exp = 0
	while curQB != nothing || curA <= length(a)
		if curA > length(a)
			if curQB == nothing
				break
			end
			comp = 1
		elseif curQB == nothing
			comp = -1
		else
			comp = cmp(curQB, a.terms[curA].exp)
		end

		# compute coef + exp of a - qb
		if comp > 0
			exp = curQB
			coef = _divisionGetNextTerm(h, q, b)

			if isCoefZero(coef)
				curQB = _prodHeapPeek(h)
				continue
			else
				coef = negateCoef(coef)
			end
		elseif comp == 0
			exp = curQB
			coef = _divisionGetNextTerm(h, q, b)

			if isCoefZero(coef) # coefs cancelled out
				curQB = _prodHeapPeek(h)
				continue
			else
				coef = a.terms[curA].coef - coef
				curA += 1
				if isCoefZero(coef)
					curQB = _prodHeapPeek(h)
					continue
				end
			end
		else
			exp = a.terms[curA].exp
			coef = a.terms[curA].coef
			curA += 1
		end

		if isDivisible(exp, leadB.exp) # new q term
			exp -= leadB.exp
			coef /= leadB.coef
			q.terms[curQ] = Term{C,N}(coef, exp)

			if (curQ + 1 > maxQSize)
				maxQSize <<= 1
				q.terms = _reallocTerms(q.terms, maxQSize)
				_prodHeapGrow!(h, maxQSize)
			end

			_prodHeapInsert!(h, q, curQ, b, 2)
			curQ += 1
		else # new r term
			r.terms[curR] = Term{C,N}(coef, exp)
			curR += 1

			if curR > maxRSize
				maxRSize <<= 1
				r.terms = _reallocTerms(r.terms, maxRSize)
			end
		end

		curQB = _prodHeapPeek(h)
	end

	if curQ - 1 < maxQSize
		_truncate!(q.terms, curQ - 1)
	end

	if curR == 1
		r = zero(a)
	else
		if curR - 1 < maxRSize
			_truncate!(r.terms, curR - 1)
		end
	end

	return q, r
end

# Base defs
import Base./
/(a::SparsePolynomial{C,N}, b::SparsePolynomial{C,N}) where {C,N} = dividePolynomials(a, b)
/(p::SparsePolynomial{C,N}, t::Term{C,N}) where {C,N} = divideBySingleTerm(p, t)
/(p::SparsePolynomial{C,N}, s::C) where {C,N} = divPolyByScalar(p, s)
