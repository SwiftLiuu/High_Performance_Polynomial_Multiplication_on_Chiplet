

module SparsePolynomials

export Term, ExpVec, SparsePolynomial, SMP, printPoly, parsePoly, partialDegrees, allocPoly

include("ExpVec/ExpVec.jl")
export ExpVec, getExp, makeExpVec, unpackExpVec, printMonomial, isDivisible

include("Common/Utils/Utils.jl")
export R_NVARS, R_NTERMS, R_SPARSITY, BIG_INT_MAX

using .ExpVecs, .Utils
const DEFAULT_EXP_T = UInt64

include("Term.jl")
include("SMP.jl")

include("Arithmetic/SparsePolyMult.jl")
export multPolyByScalar, negate, negateCoef

include("Arithmetic/SparsePolyAdd.jl")
export addPolynomials, subPolynomials

include("Arithmetic/SparsePolyDiv.jl")
export dividePolynomials, divideBySingleTerm, divPolyByScalar, isCoefZero

include("Common/TestGenerator/TestGenerator.jl")
export buildRandomPoly, generateTests, giveMeAPolynomial, giveMeAUniPolynomial, giveMeADensePolynomial, giveMeAMultiPolynomial, GeneratedPolyInfo
using .TestGenerator
end #end module

