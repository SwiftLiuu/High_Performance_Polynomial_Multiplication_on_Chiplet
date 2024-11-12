########################
# 8-bit exp vectors
########################

# Interleave bit masks and shift amounts to extract 
# partial degrees from a packed monomial (i.e. exponent vector)
const _uint8MasksOffs::Array{UInt8,1} = Array{UInt8,1}([
	#NVAR 1
	0xff,
	0,
	#NVAR 2
	0xf0,
	4,
	0x0f,
	0,
	#NVAR 3
	0xc0,
	6,
	0x38,
	3, 
	0x07,
	0,
	#NVAR 4	
	0xc0,
	6, 
	0x30,
	4, 
	0x0c,
	2, 
	0x03,
	0
]);

