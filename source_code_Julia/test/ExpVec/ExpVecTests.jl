
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

using .SparsePolynomials.ExpVecs

@testset verbose = true "Unit tests for ExpVec" begin

    Encodings = [UInt8, UInt16, UInt32, UInt64, UInt128]
    MaxNvar = [4, 8, 16, 32, 64]

    @testset "Constructors" begin
        for (T,maxN) in zip(Encodings, MaxNvar)
            for N in 1:maxN
                v = ExpVec{T,N}()
                @test v isa ExpVec{T,N}
                @test iszero(v)

                v2 = makeExpVec([T(1) for _ in 1:N])
                @test v2 isa ExpVec{T,N}
                @test !iszero(v2)

                v3 = ExpVec(T, [1 for _ in 1:N])
                @test v3 isa ExpVec{T,N}
                @test !iszero(v3)
                @test v2 == v3
                @test v3 == makeExpVec(T, [1 for _ in 1:N])


                v4 = ExpVec(typeof(v3), [1 for _ in 1:N])
                @test v4 isa ExpVec{T,N}
                @test !iszero(v4)
                @test v4 == v3
            end
        end

    end

    @testset "Zero" begin
        for (T,maxN) in zip(Encodings, MaxNvar)
            for N in 1:maxN
                z = zero(ExpVec{T,N})
                z2 = zero(z)
                @test z isa ExpVec{T,N}
                @test iszero(z)
                @test z2 isa ExpVec{T,N}
                @test iszero(z2)
            end
        end
    end


    @testset "Comparisons" begin
        #TODO!

    end



end

