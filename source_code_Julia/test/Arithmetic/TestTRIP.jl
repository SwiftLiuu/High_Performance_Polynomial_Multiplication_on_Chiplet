include("../src/SparsePolynomial.jl")
include("../src/TRIP.jl")
using .SparsePoly
using .TRIP;

p = parsePoly(BigInt, "1 + x + y + z + u + w")
println(p)
println("*")
println(p)
println("=")

prod = TRIP. _multiplyTRIP(p,p)

stats = @timed TRIP._multiplyTRIP(prod, prod)
prod = stats.value
time = stats.time
println(prod)
println("Time: ", time)

