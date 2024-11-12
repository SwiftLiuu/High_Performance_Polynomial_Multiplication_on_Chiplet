#######################
# Definitions for the SparsePolynomial data structure.
# Operations on this data structure are defined elsewhere,
# see src/SparsePolynomials.jl.
#
# Author: Alex Brandt
# Date: 2024-06-17
#
#######################

import .SparsePolynomials.Utils._expandRaggedArray!, .SparsePolynomials.Utils._isValidParam

"""
A sparse multivariate polynomial.

`C` is the coefficient type which must be a `Number`.
`N` is the number of variables in the monomial.

See also [`Term`](@ref).
"""
mutable struct SparsePolynomial{C<:Number,N}
	nvar::Int #TODO not really needed... use length(vars), update constructors
	vars::Vector{Symbol}
	terms::Vector{Term{C,N}} #There are significant performance improvements using Term{C,N} rather than Term{C} where {N}!!
end

"""
An alias to SparsePolynomial.

`SMP` is an alias to `SparsePolynomial`.
"""
SMP{C,N} = SparsePolynomial{C,N} where {C<:Number,N}



##############
# CONSTRUCTORS
##############

"""
Construct a `SparsePolynomial` from a `Vector` of `Term`s.
"""
function SparsePolynomial(terms::Vector{Term{C,N}}) where {C,N}
	nVars = N;
	vars = Vector{Symbol}(undef,nVars);
	for i in 1:nVars
		vars[i] = Symbol("_Z" * string(i))
	end
	return SparsePolynomial{C,N}(N, vars, terms)
end

"""
Construct a `SparsePolynomial` from a `Vector` of `Term`s with given variables.
"""
function SparsePolynomial(terms::Vector{Term{C,N}}, vars::Vector{Symbol}) where {C,N}
	return SparsePolynomial{C,N}(N, vars, terms)
end

"""
Construct a `SparsePolynomial` from two parallel arrays.

`coefs` is a list of coefficients
`exps` is a list of lists of integers where each inner list represents the exponents of a multivariate monomial.
"""
function SparsePolynomial(coefs::Vector{C}, exps::Matrix{T}) where {C,T<:Integer}

	nExps, nVars = size(exps);
	vars = Vector{Symbol}(undef,nVars);
	for i in 1:nVars
		vars[i] = Symbol("_Z" * string(i))
	end
	return SparsePolynomial(coefs, exps, vars)
end

#TODO integrate this into above doc.
function SparsePolynomial(coefs::Vector{C}, exps::Matrix{T}, vars::Vector{String}) where {C<:Number, T<:Integer}

	varSyms = map(Symbol, vars)
	return SparsePolynomial(coefs, exps, varSyms);
end

#TODO integrate this into above doc.
function SparsePolynomial(coefs::Vector{C}, exps::Matrix{T}, vars::Vector{Symbol}) where {C<:Number, T<:Integer}

	nExps, nVars = size(exps)
	nTerms = min(length(coefs), nExps)
	terms = Vector{Term{C,nVars}}(undef, nTerms)
	for i in 1:nTerms
		terms[i] = Term(coefs[i],makeExpVec(exps[i,:]))
	end

	_canonicalizePoly!(terms)
	return SparsePolynomial(nVars, vars, terms);
end

"""
    SparsePolynomial(p::SparsePolynomial{C,N})

Construct a polynomial with the same parameters as p.

Used for ease of programming as a copy mechanism for SparsePolynomials.
"""
function SparsePolynomial(p::SparsePolynomial{C,N}) where {C<:Number,N}
	terms = Vector{Term{C,N}}(undef, length(p))
	for x in 1:length(p)
		terms[x] = Term(p.terms[x])
	end
	return SparsePolynomial(p.nvar, p.vars, terms)
end

"""
    SparsePolynomial(::Type{C1}, p::SparsePolynomial{C2,N})

Convert a polynomial with coefficients of type C2 to a polynomial coefficients of type C1.

Returns the newly constructed polynomial.
"""
function SparsePolynomial(::Type{C1}, p::SparsePolynomial{C2,N}) where {C1<:Number,C2<:Number,N}
	nt = length(p)
	terms = Vector{Term{C1,N}}(undef, nt)
	for i in 1:nt
		terms[i] = Term(C1(p.terms[i].coef), p.terms[i].exp)
	end
	return SparsePolynomial(terms, p.vars)
end

#############
# ZERO / ONE
#############

function Base.iszero(x::SparsePolynomial{C,N}) where {C,N}
	return size(x) == 0 || iszero(x.terms[1])
end

function Base.isone(x::SparsePolynomial{C,N}) where {C,N}
	return size(x) == 1 && isone(x.terms[1])
end

function isconstant(x::SparsePolynomial{C,N}) where {C,N}
	return size(x) == 1 && iszero(x.terms[1].exp)
end

function Base.zero(x::SparsePolynomial{C,N}) where {C,N}
	if size(x) > 0
		t = Term(zero(C), zero(x.terms[1].exp))
	else
		t = Term(zero(C), zero(ExpVec{DEFAULT_EXP_T, x.nvar}))
	end
	return SMP(x.nvar, x.vars, [t])
end

function Base.zero(x::Type{SparsePolynomial{C,N}}) where {C,N}
	c = zero(C);
	e = zero(ExpVec{DEFAULT_EXP_T, N})
    t = Term(c, e);
	return SMP([t])
end

function Base.one(x::SparsePolynomial{C,N}) where {C,N}
	if size(x) > 0
		t = Term(one(C), zero(x.terms[1].exp))
	else
		t = Term(one(C), zero(ExpVec{DEFAULT_EXP_T,N}))
	end
	return SMP(x.nvar, x.vars, [t])
end

function Base.one(x::Type{SparsePolynomial{C,N}}) where {C,N}
	c = one(C)
	e = zero(ExpVec{DEFAULT_EXP_T, N})
	return SMP([Term(c,e)]);
end

function Base.isequal(x::SparsePolynomial{C,N}, y::SparsePolynomial{C,N}) where {C,N}
	if x.vars == y.vars
		return x.terms == y.terms
	end

	return false;
end
==(x::SparsePolynomial{C,N}, y::SparsePolynomial{C,N}) where {C,N} = Base.isequal(x,y);

function Base.isapprox(x::SparsePolynomial{C,N}, y::SparsePolynomial{C,N}; kwds...) where {C,N}
	if x.vars == y.vars && length(x.terms) == length(y.terms)
		ret = true
		for i in 1:length(x.terms)
			ret = ret && isapprox(x.terms[i], y.terms[i]; kwds...)
		end
		return ret
	end
	return false
end

Base.length(x::SparsePolynomial{C,N}) where {C,N} = length(x.terms)

Base.size(x::SparsePolynomial{C,N}) where {C,N} = length(x.terms)

"""
A low-level method to append a new term to the end of a polynomial's list of terms.

This method does not sort terms and does not condense like-terms.
The term to be appended must be less than every term currently existing in `p`.
If it is not, `p` is left in an undefined state.
"""
function Base.push!(p::SparsePolynomial{C,N}, x::Term{C,N}...) where {C,N}
	if (iszero(p))
		#needed so that we don't push new non-zero terms on top of a zero-term.
		empty!(p.terms)
	end

	push!(p.terms, x...)
end

function Base.resize!(p::SparsePolynomial{C,N}) where {C,N}
	#TODO!!
end

function _condenseSortedTerms!(terms::Vector{Term{C,N}}) where {C<:Number,N}
	insertIdx = 1;
	for comapreIdx in 2:length(terms)
		if terms[insertIdx].exp == terms[comapreIdx].exp
			terms[insertIdx].coef += terms[comapreIdx].coef;
		elseif (comapreIdx - insertIdx > 1)
			if terms[insertIdx].coef != C(0)
				insertIdx += 1;
			end
			terms[insertIdx].exp = terms[comapreIdx].exp
			terms[insertIdx].coef = terms[comapreIdx].coef #TODO, a faster move/swap? depending on C?
		elseif terms[insertIdx].coef != C(0)
			insertIdx += 1;
		end
	end

	if terms[insertIdx].coef == C(0)
		insertIdx -= 1;
	end

	resize!(terms, insertIdx)

end

function _canonicalizePoly!(terms::Vector{Term{C,N}}) where {C<:Number, N}
	sort!(terms, rev=true)
	_condenseSortedTerms!(terms)
end

"""
Get a list of the partial degrees for each variable in a `SparsePolynomial`.
"""
function partialDegrees(p::SparsePolynomial{C,N}) where {C,N}

	ret = zeros(Int, p.nvar)
	for term in p.terms
		for i in 1:p.nvar
			pdeg = getExp(term.exp, i)
			if pdeg > ret[i]
				ret[i] = pdeg
			end
		end
	end

	return ret;
end

function printPoly(io::IO, poly::SMP{C,N}, endStr="", debug::Bool=false) where {C<:Number,N}
	if length(poly) == 0
		print(io, "0")
	else
		if isassigned(poly.terms, 1)
			_printTerm(io, poly.terms[1], poly.vars)
		else
			if debug
				print(io, "undef ")
			end
		end
		for i in 2:length(poly)
			if isassigned(poly.terms, i)
				if poly.terms[i].coef > 0
					print(io, " + ");
				else
					print(io, " - ");
				end
				_printTermAbs(io, poly.terms[i], poly.vars);
			else
				if debug
					print(io, "undef ")
				end
			end
		end
	end
	print(io, endStr)
end

function _showPoly(io::IO, poly::SMP{C,N}) where {C,N}
	print(io, typeof(poly));
	print(io, ":\n  ");
	printPoly(io, poly, "");
end


function Base.show(io::IO, poly::SMP{C,N}) where {C,N}
	compact = get(io, :compact, true)
	if compact
		printPoly(io, poly)
	else
		_showPoly(io, poly);
	end
end

function Base.show(io::IO, mime::MIME"text/plain", poly::SMP{C,N}) where {C,N}
	compact = get(io, :compact, false)
	if compact
		printPoly(io, poly)
	else
		_showPoly(io, poly);
	end
end

function Base.parse(::Type{SMP{C,N}}, polyStr::AbstractString) where {C,N}
	return parsePoly(C, polyStr)
end

function _parseMonomial(monomial::AbstractString)
	if length(monomial) == 0
		return []
	elseif !occursin("^", monomial) && !occursin("*", monomial)
		return [Pair(monomial, 1)];
	end

	outList = []
	regExp1 = Regex("^([^\\*]+?)\\^([0-9]+)\\*?");
	regExp2 = Regex("^([^\\*]+?)\\*?");
	while true
		m = match(regExp1, monomial)
		if !isnothing(m)
            var,exp = m.captures[1:2]
		else
			m = match(regExp2, monomial)
			if isnothing(m) break; end
			var = m.captures[1]
			exp = "1"
		end

		push!(outList, Pair(var, parse(Int, exp)))

		monomial = monomial[length(m.match)+1:end]
	end

	return outList

end

function parsePoly(coefType::Type{C}, polyStr::AbstractString) where {C}
	termsAndArith = split(polyStr, " ");
	isNeg = false;
	varDict = Dict();
	varList = String[];

	coefList = C[]
	expListOfLists = Vector{Int64}[];
	for token in termsAndArith
		token = strip(token)
		if occursin("+", token)
			isNeg = false;
		elseif occursin("-", token)
			isNeg = true;
		end

		m = match(Regex("([0-9]*)\\*?([^0-9\\-+\\*][^\\-+]*)"), token);
		if !isnothing(m)
            num, mon = m.captures[1:2]
		else
			m = match(Regex("([0-9]+)"), token)
			if isnothing(m) continue; end
			num = m.captures[1]
			mon = ""
		end

		if length(num) > 0
			number = parse(coefType, num);
		else
			number = one(C)
		end
		number = isNeg ? C(-1) * number : number;
		isNeg = false;
		push!(coefList, number);

		varExpList = _parseMonomial(mon);
		for pair in varExpList
			var, exp = pair.first, pair.second;
			if !(var in keys(varDict))
				varDict[var] = length(keys(varDict)) + 1
				push!(varList, var);
			end
		end

		expList = [0 for _ in 1:length(keys(varDict))];
		for pair in varExpList
			var, exp = pair.first, pair.second;
			expList[varDict[var]] = exp;
		end
		push!(expListOfLists, expList)
	end

	_expandRaggedArray!(expListOfLists)
	M = stack(expListOfLists, dims=1);

	return SparsePolynomial(coefList, M, varList);
end


function _concatenatePolys(polys::Vector{SMP{C,N}}) where {C,N}
	finalSize = 0
	for p in polys
		finalSize += length(p)
	end

	terms = Vector{Term{C,N}}(undef, finalSize);
	termIdx = 1;
	for p in polys
		for i in eachindex(p.terms)
			terms[termIdx] = p.terms[i]
			termIdx += 1
		end
	end

	return SMP{C,N}(N, polys[1].vars, terms);
end

"""
Function which allocates a polynomial of a specified type, number of variables, and variable symbols
with a specified number of undefined terms.

@param coefType - Coefficients' type.
@param vars - Vector{Symbol} of variables. Must be of a valid
	length (1-10, inclusive).
@param alloc - Number of terms to allocate.
@return Allocated SparsePolynomial{C,N}.
"""
function allocPoly(coefType::Type{C}, vars::Vector{Symbol}, alloc::Int) where {C}
	nvars = length(vars)
	if _isValidParam(nvars, R_NVARS) && alloc > 0
		alloc = alloc > 6000000 ? 6000000 : alloc
		return SparsePolynomial(nvars, vars, Vector{Term{C,nvars}}(undef, alloc))
	end
end


function _getExpVecT(poly::SMP{C,N}) where {C,N}
	if (length(poly.terms) == 0)
		return DEFAULT_EXP_T
	else
		return typeof(poly.terms[1]).parameters[3]
	end
end

"""
Helper function to assert that a polynomial is canonical.

@param p - SparsePolynomial to check.
@throw AssertionError if p has undefined terms.
"""
function _assertCanonical(p::SparsePolynomial{C,N}) where {C,N}
	for x = 1:length(p)
		@assert isassigned(p.terms, x) "terms[$x] undefined!"
	end
	# TODO: add more checks?
end
