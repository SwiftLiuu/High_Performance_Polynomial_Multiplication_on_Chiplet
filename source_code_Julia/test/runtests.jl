#init the modules
include("../src/SparsePolynomials.jl")

#call test groups
include("ExpVec/ExpVecTests.jl")
include("SparsePolynomialTests.jl")
