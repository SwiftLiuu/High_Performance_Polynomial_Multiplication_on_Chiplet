########################
# 16-bit exp vectors
########################

# Interleave bit masks and shift amounts to extract 
# partial degrees from a packed monomial (i.e. exponent vector)
const _uint16MasksOffs::Array{UInt16,1} = Array{UInt16, 1}([
	#NVAR  1
	0xffff,
	0,
	#NVAR  2
	0xff00,
	8,
	0xff,
	0,
	#NVAR  3
	0xf800,
	11,
	0x7c0,
	6,
	0x3f,
	0,
	#NVAR  4
	0xf000,
	12,
	0xf00,
	8,
	0xf0,
	4,
	0xf,
	0,
	#NVAR  5
	0xe000,
	13,
	0x1c00,
	10,
	0x380,
	7,
	0x70,
	4,
	0xf,
	0,
	#NVAR  6
	0xc000,
	14,
	0x3000,
	12,
	0xe00,
	9,
	0x1c0,
	6,
	0x38,
	3,
	0x7,
	0,
	#NVAR  7
	0xc000,
	14,
	0x3000,
	12,
	0xc00,
	10,
	0x300,
	8,
	0xc0,
	6,
	0x38,
	3,
	0x7,
	0,
	#NVAR  8
	0xc000,
	14,
	0x3000,
	12,
	0xc00,
	10,
	0x300,
	8,
	0xc0,
	6,
	0x30,
	4,
	0xc,
	2,
	0x3,
	0
]);

