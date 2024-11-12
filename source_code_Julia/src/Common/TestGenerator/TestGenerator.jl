"""
# module TestGenerator

- Julia version:
- Author: Mars Semenova
- Date: 2024-05-27

# Examples

```jldoctest
julia>
```

NEXT STEPS (TODO)
- implement custom coefficient generation
	- user passes in a function which generates coefficients + there
	are basic ones if no preference
- for the struct GeneratedPolyInfo, make coefBound of type C and modify
all dependent code
- verify _getNextDegrees() works as intended + dense & univariate behaviour is preserved
- if multiple return types possible for local variables declare the expected one
if possible
- implement a variation of buildRandomPoly which works like ...FromMax methods in BPAS
- clear info before adding new poly info bc since it's passed by ref it may have old info
+ will break generateTests() as it indexes into info
- deal w case where invalid coefType passed
- refactor info::GeneratedPolyInfo > polyInfo
- instead of returning a pair return a tuple w simpler syntax (i.e. just return x, y)
- see if it's poss to refactor =nothing params into keyword params so that there is not need to
put a lot of nothings if including only later optional params
"""

module TestGenerator

using ..SparsePolynomials, ..Utils

import ..Utils._isValidParam

import Random

export buildRandomPoly, generateTests, giveMeAPolynomial, giveMeAUniPolynomial, giveMeADensePolynomial, giveMeAMultiPolynomial, GeneratedPolyInfo

"""
Struct used in providing info about parameters used in generating a polynomial.
"""
mutable struct GeneratedPolyInfo
	nvars::Int
	nterms::Int
	sparsity::Int
	includeNeg::Bool
	coefBound # TODO: make ::C
	seed
end

"""
Helper function used to print GeneratedPolyInfo objects.
"""
function _printGeneratedPolyInfo(io::IO, polyInfo::GeneratedPolyInfo, endStr="\n")
	print(io, typeof(polyInfo))
	print(io, ":\n  ");
	print(io, "nvars: ", string(polyInfo.nvars), endStr)
	print(io, "  nterms: ", string(polyInfo.nterms), endStr)
	print(io, "  sparsity: ", string(polyInfo.sparsity), endStr)
	print(io, "  includeNeg: ", string(polyInfo.includeNeg), endStr)
	print(io, "  coefBound: ", string(polyInfo.coefBound), endStr)
	print(io, "  seed: ", string(polyInfo.seed))
end

"""
Implementation of what should be displayed for a GeneratedPolyInfo object.
"""
function Base.show(io::IO, polyInfo::GeneratedPolyInfo)
	compact = get(io, :compact, true)
	if compact
		_printGeneratedPolyInfo(io, polyInfo)
	else
		_printGeneratedPolyInfo(io, polyInfo, ", ")
	end
end

"""
Get the next degree_t for a polynomial given the previous term's degree_t
and a "step" to take. In the univariate case this is prev+step.
In the multivariate case we consider an integer of radix maxUniDeg with
coefficients described our degrees_t. We step such that the returned value
is equal to prev + step in this radix maxUniDeg representation.
e.g: prev = [1,2,7], step = 5, maxUniDeg = 10. Then next is [1,3,2];

@param prev - Array of previous degrees.
@param step - Step to take.
@param maxUniDegree - Max degree which a variable can be.
@return Next degrees.
"""
function _getNextDegrees(prev::Array{Int, 1}, step::Int, maxUniDeg::Int)::Array{Int, 1}
	weightedDegree = 0
	weight = 0
	nvars = length(prev)

	for x = 1:nvars
		weight = maxUniDeg ^ (nvars - x)
		weightedDegree += weight * prev[x]
	end
	weightedDegree += step

	nextDegs = Array{Int}(undef, nvars)
	for x = 1:nvars
		weight = maxUniDeg ^ (nvars - x)
		nextDeg = weightedDegree
		if weight != 0
			nextDeg = trunc(Int, (weightedDegree / weight))
		end
		if nextDeg < 0 # TODO: y happening here but not in c :'), rem after ensured no neg
			nextDeg = 0
		end
		weightedDegree -= weight * nextDeg
		nextDegs[x] = nextDeg
	end
	return nextDegs
end

"""
Helper function used to generate or validate a seed.

@param setSeed - Either a Boolean set to T indicating that a seed needs
	to be generated or a seed which is to be used in random generation.
	If passed seed is invalid setSeed is set to false and no seeding is
	done.
@return A pair with a Boolean indicating whether seeding is to be done
	and a seed.
"""
function _getSeed(setSeed)::Tuple{Bool, Int}
	# verify param
	if (setSeed isa Bool && !setSeed) || (!(setSeed isa Int) && !(setSeed isa Bool)) || (setSeed isa Int && setSeed <= 0)
		return false, -1
	end

	# generate or set seed
	seed = (setSeed isa Bool) ? abs(rand(Int)) : setSeed
	Random.seed!(seed)

	return true, seed
end

"""
Helper function which returns a validated coefBound.

@param coefBound - Bound to validate.
@return Valid coefBound.
"""
function _getCoefBound(coefBound::C)::C where{C}
	if coefBound <= C(0)
		coefBound = typemax(C)
		if isinf(coefBound)
			coefBound = prevfloat(coefBound)
		end
	end

	return coefBound
end

function _getCoefBound(coefBound::C)::C where{C<:BigInt}
	if coefBound <= C(0)
		coefBound = BIG_INT_MAX
	end

	return coefBound
end

"""
Helper function to validate a range passed to buildRandomPoly().

@param range - Range to validate.
@param defRange - Associated default range for the parameter.
@return A valid Pair representing the range to be used with rand().
"""
function _validRange(range, defRange::Pair)::Pair
	min = defRange.first
	max = defRange.second
    return (range isa UnitRange{<:Integer} && range.start >= min && range.start <= range.stop) ? Pair(range.start, range.stop) : defRange # TODO: set upper limit?
end

"""
Helper function to verify includeNeg and negate coefficients.
"""
function _doIncludeNeg(c::C, includeNeg=nothing) where {C}
	if !(includeNeg isa Bool)
		includeNeg = rand(Bool)
	end

	if includeNeg # 50/50 chance of negative if includeNeg = T
		isNeg = rand(Bool)
		c = isNeg ? negateCoef(c) : c
	end

	return c
end

"""
Helper function to gen a random coef based on type and coefBound.

@param coefType - Coefficient's type.
@param includeNeg (optional) - Boolean indicating whether to allow negative numbers
	in the polynomials. If not provided value will be generated randomly.
@param coefBound (optional) - Value with the same data type as the polynomial's coefficients
	which sets the upper bound for randomly generating coefficients. It is not inclusive of the
	value. If the value is not provided or is invalid the data type's default maximum value will
	be used, or a preset maximum for BigInt.
@return Randomly generated coef.
"""
function _genRandCoef(coefType::Type{C}, includeNeg=nothing, coefBound::C=C(0))::C where{C}
	coefBound = _getCoefBound(coefBound)
	coef = C(0)

	while coef == C(0)
		coef = C(rand(1:(coefBound-1)))
		coef = _doIncludeNeg(coef, includeNeg)

		if coef >= coefBound # to account for potential Unsigned overflow
			coef %= coefBound
		end
	end

	return coef
end

function _genRandCoef(coefType::Type{C}, includeNeg=nothing, coefBound::C=C(0))::C where{C<:AbstractFloat}
	coefBound = _getCoefBound(coefBound)

	coef = C(0)
	while coef == C(0)
		coef = rand(C)
		coef *= C(10) ^ rand(0:(length(string(coef)) >= 2 ? length(string(coef))-2 : 0))
		coef %= coefBound
	end

	coef = _doIncludeNeg(coef, includeNeg)

	return coef
end

"""
Function which generates a random polynomial based on passed parameters.

@param coefType - Coefficients' type.
@param nvars (optional) - Number of variables in the polynomial. Must either
	be a valid Int (1-10, inclusive) or a valid Int range between which random
	values will be generated. If not provided, value will be generated randomly
	using a default range.
@param nterms (optional) - Number of terms in the polynomial. Must either
	be a valid Int (0-TODO, inclusive) or a valid Int range between
	which random values will be generated. If not provided or invalid,
	value will be generated randomly using a default range.
@param sparsity (optional) - Sparsity of the polynomial. Must either
	be a valid Int (2-100, inclusive) or a valid Int range between which
	random values will be generated. If not provided or invalid, value
	will be generated randomly using a default range.
@param includeNeg (optional) - Boolean indicating whether to allow negative numbers
	in the polynomials. If not provided value will be generated randomly.
@param coefBound (optional) - Value with the same data type as the polynomial's coefficients
	which sets the upper bound for randomly generating coefficients. It is not inclusive of the
	value. If the value is not provided or is invalid the data type's default maximum value will
	be used, or a preset maximum for BigInt.
@param info (optional) - A GeneratedPolyInfo[] into which a polynomial's parameters
	(nvars, nterms, sparsity, includeNeg, coefBound, seed) for generating it is stored.
@param setSeed (optional) - Either a boolean indicating whether to generate a seed for random
	generation or a positive Int seed. No seeding by default or if seed is invalid. If this
	parameter is true or a seed is passed a Pair with the polynomial and the seed will be returned.
@return A SparsePolynomial object or a Pair with a SparsePolynomial object and the seed used to
	generate it if seed != F or is valid.
"""
function buildRandomPoly(coefType::Type{C}, nvars=nothing, nterms=nothing, sparsity=nothing, includeNeg=nothing, coefBound::C=C(0), info=nothing, setSeed=false) where {C}
	# validate params
	nvarsCustom = _isValidParam(nvars, R_NVARS)
	nvarsRange = _validRange(nvars, R_NVARS)

	ntermsCustom = _isValidParam(nterms, R_NTERMS)
	ntermsRange = _validRange(nterms, R_NTERMS)

	sparsityCustom = _isValidParam(sparsity, R_SPARSITY)
	sparsityRange = _validRange(sparsity, R_SPARSITY)

	coefBound = _getCoefBound(coefBound)
	setInfo = (info isa Array{GeneratedPolyInfo, 1})

	# seeding
	setSeed, seed = _getSeed(setSeed)

	# generate params if needed
	if !nvarsCustom
		nvars = rand(nvarsRange.first:nvarsRange.second)
	end
	if !ntermsCustom
		nterms = rand(ntermsRange.first:ntermsRange.second)
	end
	if !sparsityCustom
		sparsity = rand(sparsityRange.first:sparsityRange.second)
	end
	if !(includeNeg isa Bool)
		includeNeg = rand(Bool)
	end

	# update info
	if setInfo
		polyInfo = GeneratedPolyInfo(nvars, nterms, sparsity, includeNeg, coefBound, (setSeed ? seed : false))
		push!(info, polyInfo)
	end

	# if nterms = 0 return zero poly
	if nterms == 0
		if setSeed
    		return Pair(zero(SparsePolynomial{C,nvars}), seed)
		end
		return zero(SparsePolynomial{C,nvars})
	end

	# gen poly
	coefList = Array{C}(undef, nterms)
	expListOfLists = Array{Array{Int}}(undef, nterms)
	maxTotalDeg = sparsity * nterms
	maxUniDeg = trunc(Int, ceil(maxTotalDeg ^ (1.0 / nvars)))
	degs = fill(0, nvars)
	degs[nvars] = (sparsity == 2 || rand(Bool)) ? 0 : 1 # if dense always include the const + otherwise 50/50

	# generate each term
	for x = nterms:-1:1
		# gen coefficient
		coef = _genRandCoef(coefType, includeNeg, coefBound)

		# append term
		expListOfLists[x] = degs
		coefList[x] = coef

		# gen next monomial
		step = 0
		while step == 0 # calculate step based on sparsity
			step = abs(rand(Int)) % sparsity
			step = step < sparsity / 2 ? step + trunc(Int, (sparsity / 2)) : step
		end

		degs = _getNextDegrees(degs, step, maxUniDeg)
	end

	M = stack(expListOfLists, dims=1)

	if setSeed
		return Pair(SparsePolynomial(coefList, M), seed) # TODO: add N
	end
    return SparsePolynomial(coefList, M) # TODO: add N
end

"""
Function which generates a test set of polynomials based on passed parameters.

@param ntests - Number of tests. Must be an Int.
@param coefType - Coefficients' type.
@param nvars (optional) - Number of variables in the polynomial. Must either
	be a valid Int (1-10, inclusive), which will remain constant for all the tests, or
	a valid Int range between which random values will be generated. If not
	provided, value will be generated randomly using a default range.
@param nterms (optional) - Number of terms in the polynomial. Must either
	be a valid Int (0-TODO, inclusive), which will remain constant for all the tests, or
	a valid Int range between which random values will be generated. If not
	provided or invalid, value will be generated randomly using a default range.
@param sparsity (optional) - Sparsity of the polynomial. Must either
	be a valid Int (2-100, inclusive), which will remain constant for all the tests, or
	a valid Int range between which random values will be generated. If not
	provided or invalid, value will be generated randomly using a default range.
@param includeNeg (optional) - Boolean indicating whether to allow negative numbers
	in the polynomials. If not provided value will be generated randomly.
@param coefBound (optional) - Value with the same data type as the polynomial's coefficients
	which sets the upper bound for randomly generating coefficients. It is not inclusive of the
	value. If the value is not provided or is invalid the data type's default maximum value will
	be used, or a preset maximum for BigInt.
@param info (optional) - A GeneratedPolyInfo[] into which each test's parameters
	(nvars, nterms, sparsity, includeNeg, coefBound, seed) for generating the polynomial is stored.
@param setSeed (optional) - Either a boolean indicating whether to generate a seed for random
	generation or a positive Int seed. No seeding by default or if seed is invalid. If this
	parameter is true or a seed is passed the seed will be stored in the generated polynomials'
	GeneratedPolyInfo[] objects if polyInfo has been passed and returned in a Pair
	with the array of polynomials.
@return An array of generated polynomials or a Pair with the same array and the seed
	used to generate the set of polynomials if setSeed = T or a valid seed was provided.
"""
function generateTests(ntests::Int, coefType::Type{C}, nvars=nothing, nterms=nothing, sparsity=nothing, includeNeg=nothing, coefBound::C=0, info=nothing, setSeed=false) where {C}
	if ntests <= 0 # assert ntests > 0
    	return SparsePolynomial[]
	end

	# seeding
	seed = -1
	if !(setSeed isa Bool && !setSeed)
		setSeed, seed = _getSeed(setSeed)
	end

	# generate tests
	polyList = Array{SparsePolynomial}(undef, ntests)
	for x = 1:ntests
		poly = buildRandomPoly(coefType, nvars, nterms, sparsity, includeNeg, coefBound, info)
		polyList[x] = (poly isa Pair) ? poly.first : poly
		if setSeed
			info[x].seed = seed
		end
	end

	if setSeed
		return Pair(polyList, seed)
	end
	return polyList
end

"""
Function which generates a random test set of polynomials.

@param ntests - Number of tests. Must be an int.
@param coefType - Coefficients' type.
@param info (optional) - A GeneratedPolyInfo[] into which each test's parameters
	(nvars, nterms, sparsity, includeNeg, coefBound, seed) for generating the polynomial is stored.
@param setSeed (optional) - Either a boolean indicating whether to generate a seed for random
	generation or a positive Int seed. No seeding by default or if seed is invalid. If this
	parameter is true or a seed is passed the seed will be stored in the generated polynomials'
	GeneratedPolyInfo[] objects if polyInfo has been passed and returned in a Pair
	with the array of polynomials.
@return An array of generated polynomials or a Pair with the same array and the seed
	used to generate the set of polynomials if setSeed = T or a valid seed was provided.
"""
function generateTests(ntests::Int, coefType::Type{C}, info=nothing, setSeed=false) where {C}
	return generateTests(ntests, coefType, nothing, nothing, nothing, nothing, C(0), info, setSeed)
end

"""
Function which returns a random polynomial.

@param coefType - Coefficients' type.
@param info (optional) - A GeneratedPolyInfo[] into which each test's parameters
	(nvars, nterms, sparsity, includeNeg, coefBound, seed) for generating the polynomial is stored.
@param setSeed (optional) - Either a boolean indicating whether to generate a seed for random
	generation or a positive Int seed. No seeding by default or if seed is invalid. If this
	parameter is true or a seed is passed a Pair with the polynomial and the seed will be returned.
@return A SparsePolynomial object or a Pair with a SparsePolynomial object and the seed used to
	generate it if seed != F or is valid.
"""
function giveMeAPolynomial(coefType::Type{C}, info=nothing, setSeed=false) where {C}
	return buildRandomPoly(coefType, nothing, nothing, nothing, nothing, C(0), info, setSeed)
end

"""
Function which returns a random univariate polynomial.

@param coefType - Coefficients' type.
@param info (optional) - A GeneratedPolyInfo[] into which each test's parameters
	(nvars, nterms, sparsity, includeNeg, coefBound, seed) for generating the polynomial is stored.
@param setSeed (optional) - Either a boolean indicating whether to generate a seed for random
	generation or a positive Int seed. No seeding by default or if seed is invalid. If this
	parameter is true or a seed is passed a Pair with the polynomial and the seed will be returned.
@return A SparsePolynomial object or a Pair with a SparsePolynomial object and the seed used to
	generate it if seed != F or is valid.
"""
function giveMeAUniPolynomial(coefType::Type{C}, info=nothing, setSeed=false) where {C}
	return buildRandomPoly(coefType, 1, nothing, nothing, nothing, C(0), info, setSeed)
end

"""
Function which returns a random dense polynomial.

@param coefType - Coefficients' type.
@param info (optional) - A GeneratedPolyInfo[] into which each test's parameters
	(nvars, nterms, sparsity, includeNeg, coefBound, seed) for generating the polynomial is stored.
@param setSeed (optional) - Either a boolean indicating whether to generate a seed for random
	generation or a positive Int seed. No seeding by default or if seed is invalid. If this
	parameter is true or a seed is passed a Pair with the polynomial and the seed will be returned.
@return A SparsePolynomial object or a Pair with a SparsePolynomial object and the seed used to
	generate it if seed != F or is valid.
"""
function giveMeADensePolynomial(coefType::Type{C}, info=nothing, setSeed=false) where {C}
	return buildRandomPoly(coefType, nothing, nothing, 2, nothing, C(0), info, setSeed)
end

"""
Function which returns a random multivariate polynomial.

@param coefType - Coefficients' type.
@param info (optional) - A GeneratedPolyInfo[] into which each test's parameters
	(nvars, nterms, sparsity, includeNeg, coefBound, seed) for generating the polynomial is stored.
@param setSeed (optional) - Either a boolean indicating whether to generate a seed for random
	generation or a positive Int seed. No seeding by default or if seed is invalid. If this
	parameter is true or a seed is passed a Pair with the polynomial and the seed will be returned.
@return A SparsePolynomial object or a Pair with a SparsePolynomial object and the seed used to
	generate it if seed != F or is valid.
"""
function giveMeAMultiPolynomial(coefType::Type{C}, info=nothing, setSeed=false) where {C}
	return buildRandomPoly(coefType, rand(2:R_NVARS.second), nothing, nothing, nothing, C(0), info, setSeed)
end
end #end module