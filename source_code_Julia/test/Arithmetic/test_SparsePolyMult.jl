#include("../../src/SparsePolynomials.jl")
#include("../TestUtils/TestUtils.jl")

#using .SparsePolynomials

import .SparsePolynomials.TestGenerator._genRandCoef

"""
GENERAL TESTING NOTES
-

FUTURE WORK (TODO)
-
"""

function _verifyScalarMult(p::SparsePolynomial{C,N}, s::C, np::SparsePolynomial{C,N})::Bool where {C,N}
    for x = 1:length(np)
        if np.terms[x].coef != p.terms[x].coef * s
            return false
        end
    end
    return true
end

function _verifyUnsignedNegate(p::SparsePolynomial{C,N}, np::SparsePolynomial{C,N})::Bool where {C<:Unsigned,N}
    for x = 1:length(np)
        if np.terms[x].coef != typemax(C) - p.terms[x].coef + 1
            return false
        end
    end
    return true
end

@testset verbose = true "Unit tests for SparsePolyMult.jl" begin

    # TODO

    # unit tests for multPolyByScalar()
    @testset "multPolyByScalar() tests" begin
        # test multPolyByScalar() functionality
        @testset "Test multPolyByScalar() functionality" begin
            for C in coefTypes
                info = GeneratedPolyInfo[]
                res = giveMeAPolynomial(C, info, true)
                p = res.first
                s = _genRandCoef(C)
                np = s * p

                _onFail(@test _verifyScalarMult(p, s, np)) do
                   _logInfo(info[1], s, p, np)
                end

                # s = 0 - ret zero poly
                _onFail(@test iszero(p * C(0))) do
                    println("zero s passed to multPolyByScalar() doesn't work as intended")
                end

                # p = zero poly - ret zero poly
                _onFail(@test iszero(zero(p) * _genRandCoef(C))) do
                    println("zero poly passed to multPolyByScalar() doesn't work as intended")
                end
            end
        end

        # test multPolyByScalar() error trapping
        @testset "Test multPolyByScalar() error trapping" begin
            # p has undef terms - AssertionError thrown
            p = allocPoly(Int, [Symbol("x"),Symbol("y"),Symbol("z")], 10)
            _onFail(@test_throws AssertionError p * _genRandCoef(Int)) do
                println("poly w undef terms passed to multPolyByScalar() doesn't throw AssertionError")
            end
        end
    end

    # unit tests for negate()
    @testset "negate() tests" begin
        # test negate() functionality
        @testset "Test negate() functionality" begin
            for C in coefTypes
                info = GeneratedPolyInfo[]
                res = giveMeAPolynomial(C, info, true)
                p = res.first

                if !(C<:Unsigned)
                    np = negate(p)
                    _onFail(@test _verifyScalarMult(p, C(-1), np)) do
                       _logInfo(info[1], p, np)
                    end
                end

                # p = zero poly - ret zero poly
                p = zero(p)
                _onFail(@test iszero(negate(p))) do
                    println("zero poly passed to negate() doesn't work as intended")
                end
            end
        end

        # test negate() error trapping
        @testset "Test negate() error trapping" begin
            # negate Unsigned types - ret overflowed coefs
            for C in coefTypes
                if C <: Unsigned
                    info = GeneratedPolyInfo[]
                    res = giveMeAPolynomial(C, info, true)
                    p = res.first

                    np = negate(p)
                    _onFail(@test _verifyUnsignedNegate(p, np)) do
                       _logInfo(info[1], typemax(C), p, np)
                    end
                end
            end

            # p has undef terms - AssertionError thrown
            p = allocPoly(Int, [Symbol("x"),Symbol("y"),Symbol("z")], 10)
            _onFail(@test_throws AssertionError negate(p)) do
                println("poly w undef terms passed to negate() doesn't throw AssertionError")
            end
        end
    end
end
