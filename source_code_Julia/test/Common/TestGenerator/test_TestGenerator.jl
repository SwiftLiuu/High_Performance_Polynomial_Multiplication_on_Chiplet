#=
test_TestGenerator: Unit tests for TestGenerator.jl

- Julia version: 
- Author: Mars Semenova
- Date: 2024-05-30
=#

#include("../../../src/SparsePolynomials.jl")
#include("../../TestUtils/TestUtils.jl")

#using .SparsePolynomials

import .SparsePolynomials.Utils._isValidParam, .SparsePolynomials.TestGenerator._getCoefBound

"""
GENERAL TESTING NOTES
- ensure that the changes in SparsePolynomials.jl printing methods aren't buggy
    - ensure that if the leading term's coef = -1 -[vars] is printed instead
    of -1[vars]
    - ensure printing 1 term works bc refactored printPoly
- double check that there are no duplicate includes bc there will be issues

FUTURE WORK (TODO)
- write dedicated tests for GeneratedPolyInfo
- write dedicated tests for _getNextDegrees() (e.g. ensure univariate +
dense behaviour preserved) + add checks to existing tests that ensure
that getNextDegrees() works properly
- write dedicated tests for giveMeAPolynomial()
- write dedicated tests for _getCoefBound() to ensure max included
- make onFail a macro to make unit testing cleaner
"""

"""
Helper function to verify sparsity is correct.
"""
function _verifySparsity(p::SparsePolynomial{C,N}, sparsity::Int)::Bool where {C,N}
    # TODO
    return true
end

"""
Helper function to verify that a poly was constructed correctly.
"""
function _verifyPoly(p::SparsePolynomial{C,N}, info::GeneratedPolyInfo) where {C,N}
    # correct types
    _onFail(@test p isa SparsePolynomial{C,info.nvars}) do
       _logInfo(p, info)
    end
    _onFail(@test p.terms isa Vector{Term{C,info.nvars}}) do
       _logInfo(p, info)
    end

    # correct values gened
    _onFail(@test (_isValidParam(length(p.vars), R_NVARS) && _isValidParam(length(p.terms), R_NTERMS) && _isValidParam(info.sparsity, R_SPARSITY))) do
       _logInfo(p, info)
    end
    validCoefs = true
    for x = 1:length(p)
        if p.terms[x].coef >= info.coefBound
            validCoefs = false
            continue
        end
    end
    _onFail(@test validCoefs) do
       _logInfo(p, info)
    end
    if !info.includeNeg
        nonNegCoefs = true
        for x = 1:length(p)
            if p.terms[x].coef < 0
                nonNegCoefs = false
                continue
            end
        end
        _onFail(@test nonNegCoefs) do
           _logInfo(p, info)
        end
    end
    _onFail(@test _verifySparsity(p, info.sparsity)) do
       _logInfo(p, info)
    end

    # info set properly
    _onFail(@test info.nvars == length(p.vars) && (info.nterms == length(p.terms) || iszero(p) && info.nterms == 0)) do
       _logInfo(p, info, length(p.vars), length(p.terms))
    end
end

@testset verbose = true "Unit tests for Common/TestGenerator/TestGenerator.jl" begin
    funcs = Dict("buildRandomPoly()" => buildRandomPoly, "generateTests()" => generateTests) # TODO implement iter over funcs

    for f in keys(funcs)
        # unit tests for buildRandomPoly() + generateTests()
        @testset "$f tests" begin
            # test buildRandomPoly() + generateTests() functionality
            @testset "Test $f functionality" begin
                for C in coefTypes
                    # only coefType provided - should gen a rand poly
                    info = GeneratedPolyInfo[]
                    ntests = abs(rand(1:10))
                    res = funcs[f] == buildRandomPoly ? giveMeAPolynomial(C, info, true) : funcs[f](ntests, C, info, true)
                    p = res.first

                    if funcs[f] == buildRandomPoly
                        _verifyPoly(p, info[1])
                    else
                        for x = 1:ntests
                            _verifyPoly(p[x], info[x])
                        end
                    end

                    # setSeed = T > setSeed = seed (Int) - should generate a seed + ret it + be able to replicate the gened poly
                    seed = res.second
                    info2 = GeneratedPolyInfo[]
                    res2 = funcs[f] == buildRandomPoly ? giveMeAPolynomial(C, info2, seed) : funcs[f](ntests, C, info2, seed)
                    p2 = res2.first

                    if funcs[f] == buildRandomPoly
                        _onFail(@test info[1].seed == seed) do
                           _logInfo(p, info[1])
                        end
                        _onFail(@test res2 == res) do
                           _logInfo(p, info[1], p2, info2[1])
                        end
                    else
                        for x = 1:ntests
                            _onFail(@test info[x].seed == seed) do
                               _logInfo(p[x], info[x])
                            end
                        end
                        _onFail(@test res2 == res) do
                           _logInfo(p, info, p2, info2)
                        end
                    end

                    # set coefBound - should gen polys (rand or custom) w coefs < coefBound
                    info2 = GeneratedPolyInfo[]
                    res2 = funcs[f] == buildRandomPoly ? funcs[f](C, nothing, nothing, nothing, nothing, C(100), info2, true) : funcs[f](ntests, C, nothing, nothing, nothing, nothing, C(100), info2, true)
                    p2 = res2.first

                    if funcs[f] == buildRandomPoly
                        _onFail(@test info[1].coefBound == _getCoefBound(C(0))) do
                           _logInfo(p, info[1])
                        end
                        _verifyPoly(p2, info2[1])
                        _onFail(@test info2[1].coefBound == C(100)) do
                           _logInfo(p2, info2[1])
                        end
                    else
                        for x = 1:ntests
                            _onFail(@test info[x].coefBound == _getCoefBound(C(0))) do
                               _logInfo(p[x], info[x])
                            end
                            _verifyPoly(p2[x], info2[x])
                            _onFail(@test info2[x].coefBound == C(100)) do
                               _logInfo(p2[x], info2[x])
                            end
                        end
                    end

                    # all const params set - should gen a poly according to params
                    info = GeneratedPolyInfo[]
                    res = funcs[f] == buildRandomPoly ? funcs[f](C, 2, 3, 4, false, C(100), info, true) : funcs[f](ntests, C, 2, 3, 4, false, C(100), info, true)
                    p = res.first

                    if funcs[f] == buildRandomPoly
                        _verifyPoly(p, info[1])
                    else
                        for x = 1:ntests
                            _verifyPoly(p[x], info[x])
                        end
                    end

                    # some const params set - should gen a set w said params const + others rand
                    info = GeneratedPolyInfo[]
                    res = funcs[f] == buildRandomPoly ? funcs[f](C, 2, nothing, 4, nothing, C(100), info, true) : funcs[f](ntests, C, 2, nothing, 4, nothing, C(100), info, true)
                    p = res.first

                    if funcs[f] == buildRandomPoly
                        _verifyPoly(p, info[1])
                    else
                        for x = 1:ntests
                            _verifyPoly(p[x], info[x])
                        end
                    end

                    # params are ranges - should gen polys within provided ranges
                    info = GeneratedPolyInfo[]
                    res = funcs[f] == buildRandomPoly ? funcs[f](C, 2:6, 3:10, 2:7, true, C(100), info, true) : funcs[f](ntests, C, 2:6, 3:10, 2:7, true, C(100), info, true)
                    p = res.first

                    if funcs[f] == buildRandomPoly
                        _verifyPoly(p, info[1])
                        _onFail(@test length(p.vars) >= 2 && length(p.vars) <= 6) do
                           _logInfo(p, info[1])
                        end
                        _onFail(@test length(p.terms) >= 3 && length(p.terms) <= 10) do
                           _logInfo(p, info[1])
                        end
                        _onFail(@test info[1].sparsity >= 2 && info[1].sparsity <= 7) do
                           _logInfo(p, info[1])
                        end
                    else
                        for x = 1:ntests
                            _verifyPoly(p[x], info[x])
                            _onFail(@test length(p[x].vars) >= 2 && length(p[x].vars) <= 6) do
                               _logInfo(p[x], info[x])
                            end
                            _onFail(@test length(p[x].terms) >= 3 && length(p[x].terms) <= 10) do
                               _logInfo(p[x], info[x])
                            end
                            _onFail(@test info[x].sparsity >= 2 && info[x].sparsity <= 7) do
                               _logInfo(p[x], info[x])
                            end
                        end
                    end

                    # info arr GeneratedPolyInfo[] passed - passed ref should contain correct info abt gened poly
                    info = GeneratedPolyInfo[]
                    res = funcs[f] == buildRandomPoly ? funcs[f](C, 2, 3, 4, true, C(100), info, true) : funcs[f](ntests, C, 2, 3, 4, true, C(100), info, true)
                    p = res.first
                    seed = res.second

                    if funcs[f] == buildRandomPoly
                        _onFail(@test info[1].nvars == 2 && info[1].nterms == 3 && info[1].sparsity == 4 && (info[1].includeNeg || (C<:Unsigned && !info[1].includeNeg)) && info[1].coefBound == C(100) && info[1].seed == seed) do
                           _logInfo(p, info[1])
                        end
                        _verifyPoly(p, info[1])
                    else
                        for x = 1:ntests
                            _onFail(@test info[x].nvars == 2 && info[x].nterms == 3 && info[x].sparsity == 4 && (info[x].includeNeg || (C<:Unsigned && !info[x].includeNeg)) && info[x].coefBound == C(100) && info[x].seed == seed) do
                               _logInfo(p[x], info[x])
                            end
                            _verifyPoly(p[x], info[x])
                        end
                    end

					# nterms == 0 - should ret 0 poly + passed ref should contain correct info abt gened poly
					info = GeneratedPolyInfo[]
                    res = funcs[f] == buildRandomPoly ? funcs[f](C, 2, 0, 4, true, C(100), info, true) : funcs[f](ntests, C, 2, 0, 4, true, C(100), info, true)
                    p = res.first

					if funcs[f] == buildRandomPoly
                        _onFail(@test iszero(p)) do
                           _logInfo(p, info[1])
                        end
                        _verifyPoly(p, info[1])
                    else
                        for x = 1:ntests
                            _onFail(@test iszero(p[x])) do
                               _logInfo(p[x], info[x])
                            end
                            _verifyPoly(p[x], info[x])
                        end
                    end

                    # ntests = 0 > return empty SparsePolynomial[]
                    if funcs[f] == generateTests
                        ps = generateTests(0, C)
                        _onFail(@test ps == SparsePolynomial[]) do
                           _logInfo(ps)
                        end
                    end
                end
            end

            # test buildRandomPoly() + generateTests() error trapping
            @testset "Test $f error trapping" begin
                # invalid coefType - TODO
                # TODO

                for C in coefTypes
                    # invalid nvars or nvars range - should def to rand nvars + def nvars range


                    # invalid nterms or nterms range - should def to rand nterms + def nterms range


                    # invalid sparsity or sparsity range - should def to rand sparsity + def sparsity range


                    # invalid includeNeg - should def to rand includeNeg


                    # includeNeg = true for Unsigned types = override to false


                    # invalid coefBound - should def to def range for C


                    # invalid info - should behave as if nothing was passed


                    # invalid seed - should behave as if nothing was passed


                    # invalid types passed to untyped params

                end
            end
        end
    end

    # unit tests for _getCoefBound()
    @testset "_getCoefBound() tests" begin
        # test _getCoefBound() functionality
        @testset "Test _getCoefBound() functionality" begin
            for C in coefTypes
                coefBound = _getCoefBound(C(10))

                _onFail(@test coefBound isa C && coefBound == C(10)) do
                   _logInfo(C, coefBound)
                end
            end
        end

         @testset "Test _getCoefBound() error trapping" begin
            for C in coefTypes
                # bound = 0 - ret max of type
                coefBound = _getCoefBound(C(0))
                maxBound = 0
                if C == BigInt
                    maxBound = BIG_INT_MAX
                else
			        maxBound = typemax(C)
                    if isinf(maxBound)
                        maxBound = prevfloat(maxBound)
                    end
                end

                _onFail(@test coefBound isa C && coefBound == maxBound) do
                   _logInfo(C, coefBound, maxBound)
                end

                # bound neg - ret max of type
                if !(C <: Unsigned)
                    coefBound = _getCoefBound(C(-1))

                    _onFail(@test coefBound isa C && coefBound == maxBound) do
                       _logInfo(C, coefBound, maxBound)
                    end
                end
            end
         end
    end

    # unit tests for giveMeAUniPolynomial()
    @testset "giveMeAUniPolynomial() tests" begin
        # test giveMeAUniPolynomial() functionality
        for C in coefTypes
            info = GeneratedPolyInfo[]
            res = giveMeAUniPolynomial(C, info, true)
            p = res.first

            _verifyPoly(p, info[1])
            _onFail(@test info[1].nvars == 1 && length(p.vars) == 1) do
               _logInfo(p, info[1])
            end
        end
    end

    # unit tests for giveMeADensePolynomial()
    @testset "giveMeADensePolynomial() tests" begin
        # test giveMeADensePolynomial() functionality
        for C in coefTypes
            info = GeneratedPolyInfo[]
            res = giveMeADensePolynomial(C, info, true)
            p = res.first

            _verifyPoly(p, info[1])
            _onFail(@test _verifySparsity(p, 2)) do
               _logInfo(p, info[1])
            end
        end
    end

    # unit tests for giveMeAMultiPolynomial()
    @testset "giveMeAMultiPolynomial() tests" begin
        # test giveMeAMultiPolynomial() functionality
        for C in coefTypes
            info = GeneratedPolyInfo[]
            res = giveMeAMultiPolynomial(C, info, true)
            p = res.first

            _verifyPoly(p, info[1])
            _onFail(@test info[1].nvars > 1 && length(p.vars) > 1) do
               _logInfo(p, info[1])
            end
        end
    end
end