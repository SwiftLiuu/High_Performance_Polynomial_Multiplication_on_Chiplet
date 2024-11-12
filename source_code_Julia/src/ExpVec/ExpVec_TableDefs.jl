
include("ExpVec_8Bit.jl")
include("ExpVec_16Bit.jl")
include("ExpVec_32Bit.jl")
include("ExpVec_64Bit.jl")
include("ExpVec_128Bit.jl")

# For each mask/offset array, there is a fixed index
# at which the masks/offsets start for a particular number of variables. 
# Cache that index here to make table lookups easier.
const _nvarToTableIndex::Array{UInt16,1} = Array{UInt16,1}([
	1,
	3,
	7,
	13,
	21,
	31,
	43,
	57,
	73,
	91,
	111,
	133,
	157,
	183,
	211,
	241,
	273,
	307,
	343,
	381,
	421,
	463,
	507,
	553,
	601,
	651,
	703,
	757,
	813,
	871,
	931,
	993,
	1057,
	1123,
	1191,
	1261,
	1333,
	1407,
	1483,
	1561,
	1641,
	1723,
	1807,
	1893,
	1981,
	2071,
	2163,
	2257,
	2353,
	2451,
	2551,
	2653,
	2757,
	2863,
	2971,
	3081,
	3193,
	3307,
	3423,
	3541,
	3661,
	3783,
	3907,
	4033
]);











