#######################
# Definitions for a Term of a polynomial: a coefficient and monomial pair.
#
# Author: Alex Brandt
# Date: 2024-06-17
#
#######################

"""
A single term within a sparse polynomial containing a coefficient monomial pair.

`C` is the coefficient type which must be a `Number`.
`N` is the number of variables in the monomial.

See also [`ExpVec`](@ref).

"""
mutable struct TermT{C<:Number, N, T <: Unsigned}
	coef::C
	exp::ExpVec{T,N}
end

function TermT(c::C, e::ExpVec{T,N}) where {C,N,T}
	return TermT(c, e)
end

const Term{C,N} = TermT{C,N,DEFAULT_EXP_T} where {C,N}
Term(c::C, e::ExpVec{DEFAULT_EXP_T,N}) where {C,N} = TermT(c,e)


"""
    Term{C,N}(t::Term{C,N})

Construct a term with the same parameters as t.

Used for ease of programming as a copy mechanism for terms.
"""
Term(t::Term{C,N}) where {C<:Number,N} = Term{C,N}(t.coef,t.exp)
Term{C,N}(t::Term{C,N}) where {C<:Number,N} = Term{C,N}(t.coef, t.exp)


"""
    Term{C,N}()

Construct the constant term with 0 as coefficient.

Used for ease of programming and initialization.
Terms with 0 as coefficient should not be stored in a SparsePolynomial.
"""
Term{C,N}() where {C<:Number,N} = Term{C,N}(zero(C),zero(ExpVec{DEFAULT_EXP_T,N}));

"""
    Term{C,N}(c::C)

Construct the constant term with c as coefficient.

Used for ease of programming and initialization.
Terms with 0 as coefficient should not be stored in a SparsePolynomial.
"""
Term{C,N}(c::C) where {C<:Number,N} = Term{C,N}(c,zero(ExpVec{DEFAULT_EXP_T,N}));


"""
	isless(term1, term2)

Returns true if term1 is less than term2 using a lexicographical term ordering.

Terms with equal monomials are compared based on coefficients.

Equivalent to `term1 < term2`.
"""
function Base.isless(t1::Term{C,N}, t2::Term{C,N}) where {C<:Number, N}
	if t1.exp == t2.exp
		return t1.coef < t2.coef
	end

	return Base.isless(t1.exp, t2.exp);
end

"""
	isequal(term1, term2)

Returns true if term1 exactly equals term2.

In particular, terms which are mathematically equal, but have a different
variable ordering are *not* considered equal.

Equivalent to `term1 == term2`.
"""
function Base.isequal(t1::Term{C,N}, t2::Term{C,N}) where {C<:Number, N}
	return t1.coef == t2.coef && t1.exp == t2.exp;
end

function Base.isapprox(t1::Term{C,N}, t2::Term{C,N}; kwds...) where {C<:Number, N}
	return t1.exp == t2.exp && isapprox(t1.coef, t2.coef; kwds...)
end

import Base.==, Base.<;
==(t1::Term{C,N}, t2::Term{C,N}) where {C<:Number, N} = Base.isequal(t1,t2);
<(t1::Term{C,N}, t2::Term{C,N}) where {C<:Number, N} = Base.isless(t1,t2);

import Base.iszero, Base.isone, Base.zero, Base.one

function Base.zero(t::Term{C,N}) where {C,N}
	return Term{C,N}()
end

function Base.zero(tt::Type{Term{C,N}}) where {C,N}
	return Term{C,N}()
end

function Base.one(t::Term{C,N}) where {C,N}
	return Term{C,N}(one(C))
end

function Base.one(tt::Type{Term{C,N}}) where {C,N}
	return Term{C,N}(one(C))
end

function Base.iszero(t::Term{C,N}) where {C,N}
	return iszero(t.coef)
end

function Base.isone(t::Term{C,N}) where {C,N}
	return iszero(t.exp) && isone(t.coef)
end

function isconstant(t::Term{C,N}) where {C,N}
	return iszero(t.exp)
end


"""
Implementation of what should be displayed for a Term object.
"""
function Base.show(io::IO, t::Term{C,N}) where {C,N}
	compact = get(io, :compact, true)
	if compact
		_printTerm(io, t)
	else
		_printTerm(io, t) # TODO: diff format
	end
end

function _printTerm(io::IO, term::Term{C,N}, endStr::String="") where {C<:Number, N}
	vars = Vector{Symbol}(undef,N)
	for i in 1:N
		vars[i] = Symbol("_Z" * string(i))
	end
	_printTerm(io, term, vars, endStr)
end

function _printTerm(io::IO, term::Term{C,N}, vars::Vector{Symbol}, endStr="") where {C<:Number, N}
	if isone(abs(term.coef)) && term.coef != abs(term.coef)
		print(io, "-")
	end
	if iszero(term.exp) || !isone(abs(term.coef))
		print(io, term.coef)
	end
	printMonomial(io, term.exp, vars)
	print(io, endStr)
end

function _printTermAbs(io::IO, term::Term{C,N}, vars::Vector) where {C<:Number, N}
	if iszero(term.exp) || !isone(abs(term.coef))
		print(io, abs(term.coef))
	end
	printMonomial(io, term.exp, vars);
end

