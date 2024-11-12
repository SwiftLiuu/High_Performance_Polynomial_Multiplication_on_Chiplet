
"""

Super simple "sanity" tests for SparsePolynomial{C,N}.
Essentially enough code coverege to ensure successful compilation.

- Julia version:
- Authors: Alex Brandt
- Date: 2024-06-11

```jldoctest
julia>
```
"""


using Test

using .SparsePolynomials
using .SparsePolynomials.ExpVecs

@testset verbose = true "Sanity tests for SparsePolynomial" begin

    CoefTypes = [Base.BitInteger_types..., BigInt, Float32, Float64]
    # CoefTypes = [BigInt]

    @testset "Term" begin
        for C in CoefTypes
            for N in 1:10
            # let N = 3
                @test Term{C,N}() isa Term{C,N}
                @test Term{C,N}(one(C)) isa Term{C,N}
                @test Term{C,N}(one(C), ExpVec{UInt64,N}()) isa Term{C,N}

                t1 = Term{C,N}(C(1))
                t2 = Term{C,N}(C(2))
                @test Base.isless(t1, t2)
                @test t1 < t2
                @test Base.isless(t2, t1) == false
                @test (t2 < t1) == false
                @test Base.isequal(t1,t1)
                @test Base.isequal(t1,t2) == false
                @test t1 == t1
                @test (t1 == t2) == false
            end
        end
    end
``
    @testset "Constructors" begin
        for C in CoefTypes
            for N in 1:10

                @test SparsePolynomial(3, [Symbol('x'), Symbol('y'), Symbol('z')], [Term{C,3}(one(C))]) isa SparsePolynomial{C,3}

                terms = [Term(C(1), makeExpVec([1 for _ in 1:N])), Term(C(1), makeExpVec([0 for _ in 1:N]))]
                @test SparsePolynomial(terms) isa SparsePolynomial{C,N}
                @test SparsePolynomial(Term{C,N}[]) isa SparsePolynomial{C,N}

                coefs = C[C(1), C(2)]
                exps = stack([[1 for _ in 1:N], [0 for _ in 1:N]], dims=1)
                @test SparsePolynomial(coefs, exps) isa SparsePolynomial{C,N}
                @test SparsePolynomial(coefs, exps, ["x"*string(i) for i in 1:N]) isa SparsePolynomial{C,N}
                @test SparsePolynomial(coefs, exps, [Symbol("x"*string(i)) for i in 1:N]) isa SparsePolynomial{C,N}

            end
        end
    end

    @testset "Zero / One" begin
        for C in CoefTypes
            for N in 1:10
                z = zero(SMP{C,N})
                o = one(SMP{C,N})

                @test z isa SMP{C,N}
                @test zero(z) isa SMP{C,N}
                @test z == zero(z)
                @test iszero(z)
                @test iszero(o) == false

                @test o isa SMP{C,N}
                @test one(o) isa SMP{C,N}
                @test o == one(o)
                @test isone(o)
                @test isone(z) == false

            end
        end
    end

    @testset "Length / Size" begin
        for C in CoefTypes
            for N in 1:10
                z = zero(SMP{C,N})
                o = one(z)
                @test length(z) == 1
                @test size(z) == 1
                @test length(o) == 1
                @test size(o) == 1
            end
        end
    end

    @testset "Push" begin
        for C in CoefTypes
            for N in 1:10
                t = Term(C(1), makeExpVec([1 for _ in 1:N]))
                t2 = Term(C(1), makeExpVec([0 for _ in 1:N]))
                p = SparsePolynomial([t])
                oldSize = size(p)
                push!(p, t2)

                @test size(p) == oldSize + 1
            end
        end
    end

    @testset "Partial Degrees" begin
        for C in CoefTypes
            for N in 1:10
                t = Term(C(1), makeExpVec([1 for _ in 1:N]))
                p = SparsePolynomial([t])

                ans = ones(Int, N)
                @test partialDegrees(p) == ans
            end
        end
    end

    @testset "Parse" begin
        p1 = parsePoly(BigInt, "1 + a + b + c")
        p2 = parsePoly(Int64, "1 + a + b + c")
        p3 = parse(SparsePolynomial{BigInt, 3}, "1 + x + y + z")

        @test p1 isa SparsePolynomial{BigInt, 3}
        @test length(p1) == 4
        @test partialDegrees(p1) == [1,1,1]

        @test p2 isa SparsePolynomial{Int64, 3}
        @test length(p2) == 4
        @test partialDegrees(p2) == [1,1,1]

        @test p3 isa SparsePolynomial{BigInt,3}
        @test length(p3) == 4
        @test partialDegrees(p3) == [1,1,1]
    end


    @testset "Multiplication" begin
        p = parsePoly(Int64, "1 + x + y + z")
        q = parsePoly(Int64, "3 + x^2 + y^2 + z^2")

        prod = SparsePolynomials._multiplySMP(p,q)
        @test prod isa SparsePolynomial{Int64, 3}
        @test partialDegrees(prod) == [3,3,3]

        prod2 = SparsePolynomials._multiplySMP_SDMP(p, q, 4)
        @test prod2 isa SparsePolynomial{Int64, 3}
        @test prod == prod2
    end

end

