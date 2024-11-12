
using XInts
using TimerOutputs
const to = TimerOutput()

include("../../src/SparsePolynomials.jl")

using .SparsePolynomials;



println("Threads: ", Threads.nthreads())
p = parsePoly(XInt, "1 + x + y^2 + z^3 + u^4 + w^5")
# p = parsePoly(BigInt, "1 + x + y + z")
# p = parsePoly(BigInt, "12x^1y + 24x^2y^2")
# q = parsePoly(BigInt, "24x*y + 7x^2")
# io = IOBuffer();
# print("\n\n");
# print(io, p);
# print(io, "\n")
# print(io, q);
# print(String(take!(io)))
println(p)
println("*")
println(p)
println("=")
# println(q)

# min = zero(p.terms[1].exp)
# max = ExpVec([1,0,0,0])
# prod = SparsePolynomials._multiplySMP_range(p, p, min, max)
# SparsePolynomials._multiplySMP_Chiplet(p,p,[2,2])
# SparsePolynomials._multiplySMP_Chiplet(p,p,[1,1])

p2 = SparsePolynomials._multiplySMP(p,p);
p4 = SparsePolynomials._multiplySMP(p2,p2);
p8 = SparsePolynomials._multiplySMP(p4,p4);
p10 = SparsePolynomials._multiplySMP(p2, p8);
p12 = SparsePolynomials._multiplySMP(p10,p2);
p16 = SparsePolynomials._multiplySMP(p8,p8);

@timeit to "SDMP" SparsePolynomials._multiplySMP_SDMP(p2, p2, 2); #for compilation

reset_timer!(to)
println("Serial: ")
@timeit to "Serial" SparsePolynomials._multiplySMP(p12, p12)
show(to)

reset_timer!(to)
println("SDMP: ")
@timeit to "SDMP" SparsePolynomials._multiplySMP_SDMP(p12, p12, Threads.nthreads())
show(to)


# println(length(p16))
# prod = SparsePolynomials._multiplySMP_SDMPEdge(p, p)
# prod3 = SparsePolynomials._multiplySMP_Chiplet(p16, p8, [1 for _ in 1:8])
# @time prod4 = SparsePolynomials._multiplySMP_Chiplet(p16, p8, [2 for _ in 1:4])
# @time prod5 = SparsePolynomials._multiplySMP_Chiplet(p16, p8, [4 for _ in 1:2])

# @time prod2 = SparsePolynomials._multiplySMP(p16,p8)

# println("")
# println(prod)
# println(prod3 ==prod2)

# print(length.(prod))


#=

# prod = p;
inp = SparsePolynomials._multiplySMP(p,p)
inp = SparsePolynomials._multiplySMP(inp,inp)
inp = SparsePolynomials._multiplySMP(inp,inp)
inp = SparsePolynomials._multiplySMP(inp,inp)
# inp = SparsePolynomials._multiplySMP(inp,inp)
p2 = SparsePolynomials._multiplySMP_SDMP(p,p,1)
# println(inp)
# inp = SparsePolynomials._multiplySMP(inp, inp);
# inp = SparsePolynomials._multiplySMP(inp, inp);
# inp = SparsePolynomials._multiplySMP(inp, inp);
# inp = SparsePolynomials._multiplySMP(inp, inp);
# prod = SparsePolynomials._multiplySMP_SDMP(prod, prod,2);
# prod = SparsePolynomials._multiplySMP(prod, prod);
stats_serial = @timed SparsePolynomials._multiplySMP(inp, inp);
#prod = SparsePolynomials._multiplySMP(prod, prod);
# stats = @timed SparsePolynomials._multiplySMP(prod, prod)
println("start")
stats_para = @timed SparsePolynomials._multiplySMP_SDMP(inp, inp, Threads.nthreads())

println("Serial:   ", stats_serial.time)
println("Serial:   ", length(stats_serial.value))

println("Parallel: ", stats_para.time)
println("Parallel: ", length(stats_para.value))

if stats_serial.value != stats_para.value
	println(stats_serial.value, "\n")
	println(stats_para.value)
	error("TEST FAILED!!!")
else
	println("Test passed!")
end

# prod = stats.value
# time = stats.time
# prod = SparsePolynomials._multiplySMP(prod, prod)
# fp = open("test4.txt", "w")
# println(fp, prod)
# close(fp)
# println(length(prod))

# println("Time: ", time)

#println("All stats:\n", stats)


=#


