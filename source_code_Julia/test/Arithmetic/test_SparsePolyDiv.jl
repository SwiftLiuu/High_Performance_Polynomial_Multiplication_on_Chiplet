#=
test_SparsePolyDiv: Unit tests for SparsePolyDiv.jl

- Julia version:
- Author: Mars Semenova, Alex Brandt
- Date: 2024-06-28
=#
include("../../src/SparsePolynomials.jl")
include("../TestUtils/TestUtils.jl")

using .SparsePolynomials

import .SparsePolynomials._multiplySMP, .SparsePolynomials.TestGenerator._genRandCoef

"""
GENERAL TESTING NOTES
-

FUTURE WORK (TODO)
- ensure that all cases covered since random gen not feasible without extra work (e.g. importing
    div method from C)
"""


@testset verbose = true "Sanity tests for SparsePolyDiv.jl" begin
    @testset "divPolyByScalar() tests" begin
        p = parsePoly(Float32, "223x^2*y + 12x*y + y^3 + 12y")
        s = 124.0f0
        m = s*p
        q = divPolyByScalar(m, s)
        @test isapprox(p,q)

        for C in floatCoefTypes
            info = GeneratedPolyInfo[]
            p = buildRandomPoly(C, 1:5, 1:30, 2:5, true, C(100), info, true).first
            s = _genRandCoef(C, true, C(100))
            mult = s * p
            q = divPolyByScalar(mult, s)

            #_onFail(@test true) do
            _onFail(@test isapprox(q, p)) do
                println(mult)
                println(s)
                println(p)
                println(q)
                # _logInfo(mult, s, q, p, info[1], C, eps(C))
            end

            # info1 = GeneratedPolyInfo[]
            # info2 = GeneratedPolyInfo[]
            # p1 = buildRandomPoly(C, 1:5, 1:10, 2:5, true, C(100), info1, true).first
            # p2 = buildRandomPoly(C, length(p1.vars), 1:30, 2:5, true, C(100), info2, true).first
            # mult = _multiplySMP(p1, p2)

            # q, r = dividePolynomials(mult, p2)
            # #_onFail(@test true) do # TODO: rem
            # _onFail(@test isapprox(q, p1) && isapprox(r, zero(p1))) do
            #     _logInfo(mult, p2, info2[1], q, p1, info1[1], r, C, eps(C))
            # end

        end
    end

    @testset "dividePolynomials() Int" begin
        for C in [Int32, Int64, BigInt]
            info = GeneratedPolyInfo[]
            f = buildRandomPoly(C, 2:5, 2:30, 2:5, true, C(100), info, true).first
            g = buildRandomPoly(C, info[1].nvars, 2:30, 2:5, true, C(100), info, true).first

            m = _multiplySMP(f,g)
            q, r = dividePolynomials(m, g)

            #_onFail(@test true) do
            _onFail(@test q == f && iszero(r)) do
                _logInfo(m, f, g, info[1], q, r)
            end

            f = SparsePolynomial(Rational{C}, f);
            g = SparsePolynomial(Rational{C}, g);

            m = _multiplySMP(f,g)
            m = addPolynomials(m, f)
            q,r = dividePolynomials(m, g)

            qg = _multiplySMP(q,g);
            qgr = addPolynomials(qg, r);

            _onFail(@test m == qgr) do
                _logInfo(m, f, g, info[1], q, r)
            end

        end
    end
end


#=

@testset verbose = true "Unit tests for SparsePolyDiv.jl" begin
    # unit tests for divPolyByScalar()
    @testset "divPolyByScalar() tests" begin
        # test divPolyByScalar() functionality
        @testset "Test divPolyByScalar() functionality" begin
            testCases = []
            for C in [Float64] #floatCoefTypes
                info = GeneratedPolyInfo[]
                p = buildRandomPoly(C, 1:5, 1:30, 2:5, true, C(100), info, true).first
                s = _genRandCoef(C)
                mult = s * p
                q = divPolyByScalar(mult, s)

                #_onFail(@test true) do
                _onFail(@test isapprox(q, p)) do
                    _logInfo(mult, s, q, p, info[1], C, eps(C))
                end

                # 2 consts
                c = _genRandCoef(C)
                q = divPolyByScalar(SparsePolynomial([Term{C,2}(c)]), s)
                _onFail(@test isapprox(q, SparsePolynomial([Term{C,2}(c / s)]))) do
                    _logInfo(c, s, q, (c / s), C, eps(C))
                end
            end
        end

        # test divPolyByScalar() error trapping
        @testset "Test divPolyByScalar() error trapping" begin
            # divide by 0 - throw err
            p = giveMeAPolynomial(Float64)
            _onFail(@test_throws DivideError divPolyByScalar(p, Float64(0))) do
                println("div by 0 did not throw DivideError")
            end

            # p has undef terms - AssertionError thrown
            p = allocPoly(Float64, [Symbol("x"),Symbol("y"),Symbol("z")], 10)
            _onFail(@test_throws AssertionError divPolyByScalar(p, _genRandCoef(Float64))) do
                println("poly w undef terms passed to divPolyByScalar() doesn't throw AssertionError")
            end
        end
    end

    # unit tests for divideBySingleTerm()
    @testset "divideBySingleTerm() tests" begin
        # test divideBySingleTerm() functionality
        @testset "Test divideBySingleTerm() functionality" begin
            testCases = []
            for C in floatCoefTypes
                vars = [Symbol("x"),Symbol("y")]
                p = parsePoly(C, "4x^2 - 8x^1*y + 4x + 16y^2 - 8y")
                t = Term{C, 2}(C(2), makeExpVec([1,0]))

                testCases = [[p, # basic use case
                    t,
                    (parsePoly(C, "2x - 4y + 2"),
                    SparsePolynomial(parsePoly(C, "x + 16y^2 - 8y").terms[2:end], vars))],
                [parsePoly(C, "47x^3*y^2 + 12x^3 + 7x^2*y^4 - 47x^2*y - 71x^1*y^6 + 70x^1*y^3 - 15x^1*y - 82y^6 + 30y^3 - 25"), # larger basic use case
                    Term{C, 2}(C(1), makeExpVec([1,1])),
                    (parsePoly(C, "47x^2*y + 7x^1*y^3 - 47x - 71y^5 + 70y^2 - 15"),
                    parsePoly(C, "12x^3 - 82y^6 + 30y^3 - 25"))],
                [SparsePolynomial([Term{C, 2}(C(20))], vars), # 2 const terms
                    Term{C, 2}(C(2)),
                    (divPolyByScalar(SparsePolynomial([Term{C, 2}(C(20))], vars), C(2)),
                    zero(p))],
                [p, # divisor = const term
                    Term{C, 2}(C(2)),
                    (divPolyByScalar(p, C(2)),
                    zero(p))],
                [SparsePolynomial([Term{C, 2}(C(20))], vars), # dividend = const term
                    t,
                    (zero(p),
                    SparsePolynomial([Term{C, 2}(C(20))], vars))],
                [zero(p), # div 0 poly
                    t,
                    (zero(p),
                    zero(p))],
                [p, # no div necessary
                    Term{C, 2}(C(4), makeExpVec([3,1])),
                    (zero(p),
                    p)],
                [SparsePolynomial([t], vars), # a = b
                    t,
                    (SparsePolynomial([Term{C, 2}(C(1))], vars),
                    zero(p))]]

                # run test cases
                for x = 1:0 # TODO: rem
                #for x = 1:length(testCases)
                    q, r = divideBySingleTerm(testCases[x][1], testCases[x][2])

                    _onFail(@test q == testCases[x][3][1] && r == testCases[x][3][2]) do
                        println("Test case $x:")
                        _logInfo(testCases[x][1], testCases[x][2], q, testCases[x][3][1], r, testCases[x][3][2])
                    end
                end

                # div by coef
                c = _genRandCoef(C)
                q, r = divideBySingleTerm(p, Term{C,2}(c))
                div = divPolyByScalar(p, c)
                _onFail(@test isapprox(q, div) && r == zero(p)) do
                    _logInfo(p, c, q, div, r, C, eps(C))
                end

                # run rand test case using mult
                info = GeneratedPolyInfo[]
                p = buildRandomPoly(C, 1:5, 1:30, 2:5, true, C(100), info, true).first
                tpoly = buildRandomPoly(C, length(p.vars), 1:30, 2:5, true, C(100), nothing, true)
                t = tpoly.first.terms[rand(1:length(tpoly.first.terms))]
                mult = _multiplySMP(p, SparsePolynomial([t]))
                q, r = divideBySingleTerm(mult, t)

                #_onFail(@test true) do # TODO: rem
                _onFail(@test isapprox(q, p) && r == zero(p)) do
                    _logInfo(mult, t, tpoly.second, info[1], q, p, r, C, eps(C))
                end
            end
        end

        # test divideBySingleTerm() error trapping
        @testset "Test divideBySingleTerm() error trapping" begin
            # p has undef terms - AssertionError thrown
            p = allocPoly(Float64, [Symbol("x"),Symbol("y"),Symbol("z")], 10)
            _onFail(@test_throws AssertionError divideBySingleTerm(p, buildRandomPoly(Float64, 3, 1).terms[1])) do
                println("poly w undef terms passed to divideBySingleTerm() doesn't throw AssertionError")
            end
        end
    end

    # unit tests for dividePolynomials()
    @testset "dividePolynomials() tests" begin
        # test dividePolynomials() functionality
        @testset "Test dividePolynomials() functionality" begin
            for C in floatCoefTypes
                vars = [Symbol("x"),Symbol("y")]
                p = parsePoly(C, "4x^2 - 8x^1*y + 4x + 16y^2 - 8y")

                testCases = [[parsePoly(C, "4 + 4x^3 + 4y^2 + 4x^2 + 4y"), # basic use case
                    parsePoly(C, "x + 2y"),
                    (parsePoly(C, "4x^2 - 8x^1*y + 4x + 16y^2 - 8y"),
                    SparsePolynomial(parsePoly(C, "x - 32y^3 + 20y^2 + 4y + 4").terms[2:end], vars))],
                [parsePoly(C, "47x^3*y^2 + 12x^3 + 7x^2*y^4 - 47x^2*y - 71x^1*y^6 + 70x^1*y^3 - 15x^1*y - 82y^6 + 30y^3 - 25y"), # larger basic use case
                    parsePoly(C, "x^3 + 52x^2*y^4 + x^2*y^2 - 86x^1*y^5 - 30x^1*y^3 + x + 90y^4 + 96y^2 + 37"),
                    (SparsePolynomial(parsePoly(C, "x + 47y^2 + 12").terms[2:end], vars),
                    parsePoly(C, "-2444x^2*y^6 - 664x^2*y^4 - 12x^2*y^2 - 47x^2*y + 4042x^1*y^7 - 71x^1*y^6 + 2442x^1*y^5 + 430x^1*y^3 - 47x^1*y^2 - 15x^1*y - 12x - 4312y^6 - 5592y^4 + 30y^3 - 2891y^2 - 25y - 444"))],
                [zero(p), # div 0 poly
                    p,
                    (zero(p),
                    zero(p))],
                [p, # b = single term
                    SparsePolynomial([Term{C, 2}(C(2), makeExpVec([1,0]))], vars),
                    divideBySingleTerm(p, Term{C, 2}(C(2), makeExpVec([1,0])))],
                [p, # divisor = const term
                    SparsePolynomial([Term{C, 2}(C(2))], vars),
                    (divPolyByScalar(p, C(2)),
                    zero(p))],
                [SparsePolynomial([Term{C, 2}(C(20))], vars), # dividend = const term
                    p,
                    (zero(p),
                    SparsePolynomial([Term{C, 2}(C(20))], vars))],
                [parsePoly(C, "-46x^2*y^3 + 45x^2 - 51x^1*y^4 + 51x^1*y + 91y^5 - 20y^2 + 11"), # no div necessary
                    parsePoly(C, "61x^3 + 52x^2*y^4 - x^2*y^2 - 86x^1*y^5 - 30x^1*y^3 - x + 90y^4 + 96y^2 + 37"),
                    (zero(p),
                    parsePoly(C, "-46x^2*y^3 + 45x^2 - 51x^1*y^4 + 51x^1*y + 91y^5 - 20y^2 + 11"))],
                [p, # a = b
                    p,
                    (SparsePolynomial([Term{C, 2}(C(1))], vars),
                    zero(p))]]

                # run test cases
                for x = 1:0 # TODO: rem
                #for x = 1:length(testCases)
                    q, r = dividePolynomials(testCases[x][1], testCases[x][2])

                    _onFail(@test q == testCases[x][3][1] && r == testCases[x][3][2]) do
                        println("Test case $x:")
                        _logInfo(testCases[x][1], testCases[x][2], q, testCases[x][3][1], r, testCases[x][3][2])
                    end
                end

                # div by coef
                c = _genRandCoef(C)
                q, r = dividePolynomials(p, SparsePolynomial([Term{C,2}(c)]))
                div = divPolyByScalar(p, c)
                _onFail(@test isapprox(q, div) && r == zero(p)) do
                    _logInfo(p, c, q, div, r, C, eps(C))
                end

                # run rand test case using mult
                info1 = GeneratedPolyInfo[]
                info2 = GeneratedPolyInfo[]
                p1 = buildRandomPoly(C, 1:5, 1:30, 2:5, true, C(100), info1, true).first
                p2 = buildRandomPoly(C, length(p1.vars), 1:30, 2:5, true, C(100), info2, true).first
                mult = _multiplySMP(p1, p2)

                q, r = dividePolynomials(mult, p2)
                #_onFail(@test true) do # TODO: rem
                _onFail(@test isapprox(q, p1) && isapprox(r, zero(p1))) do
                    _logInfo(mult, p2, info2[1], q, p1, info1[1], r, C, eps(C))
                end

                q, r = dividePolynomials(mult, p1)
                #_onFail(@test true) do # TODO: rem
                _onFail(@test isapprox(q, p2) && r == zero(p1)) do
                    _logInfo(mult, p1, info1[1], q, p2, info2[1], r, C, eps(C))
                end
            end
        end

        # test dividePolynomials() error trapping
        @testset "Test dividePolynomials() error trapping" begin
            # divide by 0 - throw err
            p = giveMeAPolynomial(Float64)
            _onFail(@test_throws DivideError dividePolynomials(p, zero(p))) do
                println("div by 0 did not throw DivideError")
            end

            # p1 +/or p2 has undef terms - AssertionError thrown
            p = buildRandomPoly(Float64, 3, 10)
            p1 = allocPoly(Float64, [Symbol("x"),Symbol("y"),Symbol("z")], 10)
            p2 = allocPoly(Float64, [Symbol("x"),Symbol("y"),Symbol("z")], 10)

            _onFail(@test_throws AssertionError dividePolynomials(p1, p2)) do
                println("polys w undef terms passed to dividePolynomials() doesn't throw AssertionError")
            end

            _onFail(@test_throws AssertionError dividePolynomials(p1, p)) do
                println("polys w undef terms passed to dividePolynomials() doesn't throw AssertionError")
            end

            _onFail(@test_throws AssertionError dividePolynomials(p, p2)) do
                println("polys w undef terms passed to dividePolynomials() doesn't throw AssertionError")
            end
        end
    end
end

=#