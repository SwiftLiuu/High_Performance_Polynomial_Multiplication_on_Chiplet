#=
test_SparsePolyAdd: Unit tests for SparsePolyAdd.jl
- Julia version:
- Author: Mars Semenova
- Date: 2024-06-12
=#
#include("../../src/SparsePolynomials.jl")
#include("../TestUtils/TestUtils.jl")

#using .SparsePolynomials

import .SparsePolynomials._canonicalizePoly!

"""
GENERAL TESTING NOTES
-

FUTURE WORK (TODO)
-
"""

"""
Helper method which tests the implementation of addition with a basic
addition alg.
"""
function _verifyAdd(a::SparsePolynomial{C,N}, b::SparsePolynomial{C,N}, res::SparsePolynomial{C,N})::Bool where {C,N}
    terms = vcat(a.terms, b.terms)
    SparsePolynomials._canonicalizePoly!(terms)
    poly = SparsePolynomial(terms, a.vars)
    return poly == res
end

@testset verbose = true "Unit tests for SparsePolyAdd.jl" begin
    funcs = Dict("addPolynomials()" => addPolynomials, "subPolynomials()" => subPolynomials) # TODO implement iter over funcs

    for f in keys(funcs)
        # unit tests for addPolynomials() + subPolynomials()
        @testset "$f tests" begin
            # test addPolynomials() + subPolynomials() functionality
            @testset "Test $f functionality" begin
                for C in coefTypes
                    info1 = GeneratedPolyInfo[]
                    res1 = buildRandomPoly(C, 1:5, 1:50, 2:5, true, C(100), info1, true)
                    p1 = res1.first
                    info2 = GeneratedPolyInfo[]
                    res2 = buildRandomPoly(C, info1[1].nvars, 1:50, 2:5, true, C(100), info2, true)
                    p2 = res2.first
                    np = funcs[f](p1, p2)

                    # basic test
                    valid = funcs[f] == addPolynomials ? _verifyAdd(p1, p2, np) : _verifyAdd(p1, negate(p2), np)
                    _onFail(@test valid) do
                        println(f)
                       _logInfo(info1[1], p1, info2[1], p2, np)
                    end

                    # a = zero poly - ret b (add) or ret -b (sub)
                    np = funcs[f](zero(p1), p2)
                    valid = funcs[f] == addPolynomials ? np == p2 : np == negate(p2)
                    _onFail(@test valid) do
                        println(f)
                       _logInfo(info2[1], p2, np)
                    end

                    # b = zero poly - ret a
                    np = funcs[f](p1, zero(p2))
                    _onFail(@test np == p1) do
                        println(f)
                       _logInfo(info1[1], p1, np)
                    end

                    # a + -a or a - a - ret 0 poly
                    np = funcs[f] == addPolynomials ? p1 + negate(p1) : p1 - p1
                    _onFail(@test np == zero(p1)) do
                        println(f)
                       _logInfo(info1[1], p1, np)
                    end

                    # a + b = zero polys = ret 0 poly
                    np = funcs[f](zero(p1), zero(p2))
                    _onFail(@test np == zero(p1)) do
                       println("$f called w zero polys doesn't work as intended")
                    end
                end
            end

            # test addPolynomials() + subPolynomials() error trapping
            @testset "Test $f error trapping" begin
                # p1 +/or p2 has undef terms - AssertionError thrown
                p = buildRandomPoly(Int, 3, 10)
                p1 = allocPoly(Int, [Symbol("x"),Symbol("y"),Symbol("z")], 10)
                p2 = allocPoly(Int, [Symbol("x"),Symbol("y"),Symbol("z")], 10)

                _onFail(@test_throws AssertionError funcs[f](p1, p2)) do
                    println("polys w undef terms passed to $f doesn't throw AssertionError")
                end

                _onFail(@test_throws AssertionError funcs[f](p1, p)) do
                    println("polys w undef terms passed to $f doesn't throw AssertionError")
                end

                _onFail(@test_throws AssertionError funcs[f](p, p2)) do
                    println("polys w undef terms passed to $f doesn't throw AssertionError")
                end
            end
       end
   end
end
