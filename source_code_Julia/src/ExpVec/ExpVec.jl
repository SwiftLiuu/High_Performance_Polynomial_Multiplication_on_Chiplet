


module ExpVecs

export ExpVec, getExp, makeExpVec, unpackExpVec, printMonomial, isDivisible
# export getExp, getExp_uint8, getExp_uint16, getExp_uint32, getExp_uint64, getExp_uint128
# export getMaxPartialDegree, setPartialDegree!
# export makeExpVec, unpackExpVec, fillExpVec!
# export printMonomial



#=
An encoding of packed exponent vectors.

Goals:
 - expand data representations as needed on overflow
 - use parametric types to ensure combined exponent vectors encode same nvar


=#

# Import the lookup tables for masks and offsets
include("ExpVec_TableDefs.jl")


# N is an integer that indicates the number of variables stored in the ExpVec
# It makes the number of variables part of the type itself.
mutable struct ExpVec{T<:Unsigned, N}
	vec::T
end
ExpVec{T,N}() where {T,N} = ExpVec{T,N}(0);


##############
# CONSTRUCTORS
##############

ExpVec(e::ExpVec{T,N}) where {T,N} = ExpVec{T,N}(e.vec);
ExpVec(exps::Vector{T}) where {T<:Integer} = makeExpVec(exps)
ExpVec(::Type{T1}, exps::Vector{T2}) where {T1<:Unsigned,T2<:Integer} = fillExpVec!(ExpVec{T1,length(exps)}(), exps);
ExpVec(::Type{ExpVec{T,N}}, exps::Vector{T2}) where {T,N,T2<:Integer} = fillExpVec!(ExpVec{T,N}(), exps);

makeExpVec(exps::Vector{T}) where {T<:Integer} = fillExpVec!(ExpVec{unsigned(T),length(exps)}(), exps);
#captured by above rule, but keeping to show use of Union in where
# makeExpVec(exps::Vector{T}) where {T<:Union{UInt128,Int128}} = fillExpVec!(ExpVec{UInt128,length(exps)}(), exps);
makeExpVec(::Type{T1}, exps::Vector{T2}) where {T1<:Unsigned,T2<:Integer} = fillExpVec!(ExpVec{T1,length(exps)}(), exps);
makeExpVec(::Type{ExpVec{T,N}}, exps::Vector{T2}) where {T,N,T2<:Integer} = fillExpVec!(ExpVec{T,N}(), exps);

import Base.convert, Base.+, Base.-;
# convert(::Type{ExpVec{T1,N}}, e::ExpVec{T2,N}) = ...
convert(::Type{ExpVec{T1,N}}, e::ExpVec{T2,N}) where {T1,T2,N} = makeExpVec(convert(Vector{T1},unpackExpVec(e)));
# convert(::Type{ExpVec{UInt128,N}}, e::ExpVec{T2,N}) where {N,T2} = makeExpVec(unpackExpVec)

ExpVec{T1,N}(other::ExpVec{T2,N}) where {T1<:Unsigned,T2<:Unsigned,N} = makeExpVec(convert(Vector{T1}, unpackExpVec(other)));


function Base.isequal(e1::ExpVec{T,N}, e2::ExpVec{T,N}) where {T<:Unsigned, N}
	return e1.vec == e2.vec;
end

function Base.isless(e1::ExpVec{T,N}, e2::ExpVec{T,N}) where {T<:Unsigned, N}
	return e1.vec < e2.vec;
end

function Base.iszero(e::ExpVec{T,N}) where {T<:Unsigned, N}
	return iszero(e.vec)
end

function Base.zero(e::ExpVec{T,N}) where {T<:Unsigned, N}
	return ExpVec{T,N}();
end

function Base.zero(e::Type{ExpVec{T,N}}) where {T<:Unsigned, N}
	return ExpVec{T,N}();
end


getNumVariables(e::ExpVec{T,N}) where {T<:Unsigned, N} = N

import Base.==, Base.<;
==(e1::ExpVec{T,N}, e2::ExpVec{T,N}) where {T<:Unsigned, N} = Base.isequal(e1,e2);
<(e1::ExpVec{T,N}, e2::ExpVec{T,N}) where {T<:Unsigned, N} = Base.isless(e1,e2);

"""
Function which determines whether a term is divisible by another
term by comparing their exponents.

@param e1 - Exponent of the dividend.
@param e2 - Exponent of the divisor.
@return Whether the dividend is divisible by the divisor.
"""
function isDivisible(e1::ExpVec{T,N}, e2::ExpVec{T,N})::Bool where {T,N}
	ret = true;
	i = 1
	while (ret && i <=N)
		ret = ret && (getExp(e1, i) >= getExp(e2, i))
		i+=1
	end
	return ret
end

function _getExp(m::ExpVec{T,N}, i::Integer, maskTable, idxTable) where {T,N}
	#to compute the correct index within a 1-based indexing scheme, convert i to 0-based.
	idx = idxTable[N] + (2*(i-1));
	mask = maskTable[idx];
	off = maskTable[idx+1];
	# print(idx, " ", mask, " ", off, "\n");
	return (m.vec & mask) >> off;
end

function _getMaxPartialDeg(m::ExpVec{T,N}, i::Integer, maskTable, idxTable) where {T,N}
	#to compute the correct index within a 1-based indexing scheme, convert i to 0-based.
	idx = idxTable[N] + (2*(i-1));
	mask = maskTable[idx]
	while mask & 1 == 0
		mask = mask >> 1
	end
	return mask
end

function _setExp!(m::ExpVec{T,N}, i::Integer, v::Integer, maskTable, idxTable) where {T,N}
	#to compute the correct index within a 1-based indexing scheme, convert i to 0-based.
	idx = idxTable[N] + (2*(i-1));
	mask = maskTable[idx];
	off = maskTable[idx+1];
	m.vec = m.vec & (~mask); #remove any previous value with inverted mask
	m.vec = m.vec | (T(v) << off); #note T(v) is crucial for shift
	# print(idx, " ", mask, " ", off, "\n");
	return m
end

function _setAllExp!(m::ExpVec{T,N}, v::Vector{Integer}, maskTable, idxTable) where {T,N}
	packedVec = T(0);
	for (index, val) in enumerate(v)
		#enumerate uses 1-based indexing, so compute maskTable idx as index-1.
		idx = idxTable[N] + (2*(index-1));
		off = maskTable[idx+1]; #+1 to get offset; mask is at +0.
		packedVec = packedVec | (T(val) << off);
	end
	m.vec = packedVec;
end


getExp_uint8(m::ExpVec{UInt8,N}, i::Integer)     where {N} = _getExp(m, i, _uint8MasksOffs, _nvarToTableIndex);
getExp_uint16(m::ExpVec{UInt16,N}, i::Integer)   where {N} = _getExp(m, i, _uint16MasksOffs, _nvarToTableIndex);
getExp_uint32(m::ExpVec{UInt32,N}, i::Integer)   where {N} = _getExp(m, i, _uint32MasksOffs, _nvarToTableIndex);
getExp_uint64(m::ExpVec{UInt64,N}, i::Integer)   where {N} = _getExp(m, i, _uint64MasksOffs, _nvarToTableIndex);
getExp_uint128(m::ExpVec{UInt128,N}, i::Integer) where {N} = _getExp(m, i, _uint128MasksOffs, _nvarToTableIndex);

##extract nvar from ExpVec type and call underlying method
getExp( m::ExpVec{UInt8,N}, i::Integer)   where {N} = getExp_uint8(m, i)
getExp( m::ExpVec{UInt16,N}, i::Integer)  where {N} = getExp_uint16(m, i)
getExp( m::ExpVec{UInt32,N}, i::Integer)  where {N} = getExp_uint32(m, i)
getExp( m::ExpVec{UInt64,N}, i::Integer)  where {N} = getExp_uint64(m, i)
getExp( m::ExpVec{UInt128,N}, i::Integer) where {N} = getExp_uint128(m, i)


getMaxPartialDegree_uint8(m::ExpVec{UInt8,N}, i::Integer)     where {N} = _getMaxPartialDeg(m, i, _uint8MasksOffs, _nvarToTableIndex);
getMaxPartialDegree_uint16(m::ExpVec{UInt16,N}, i::Integer)   where {N} = _getMaxPartialDeg(m, i, _uint16MasksOffs, _nvarToTableIndex);
getMaxPartialDegree_uint32(m::ExpVec{UInt32,N}, i::Integer)   where {N} = _getMaxPartialDeg(m, i, _uint32MasksOffs, _nvarToTableIndex);
getMaxPartialDegree_uint64(m::ExpVec{UInt64,N}, i::Integer)   where {N} = _getMaxPartialDeg(m, i, _uint64MasksOffs, _nvarToTableIndex);
getMaxPartialDegree_uint128(m::ExpVec{UInt128,N}, i::Integer) where {N} = _getMaxPartialDeg(m, i, _uint128MasksOffs, _nvarToTableIndex);

getMaxPartialDegree(m::ExpVec{UInt8,N}, i::Integer)   where {N} = getMaxPartialDegree_uint8(m, i);
getMaxPartialDegree(m::ExpVec{UInt16,N}, i::Integer)   where {N} = getMaxPartialDegree_uint16(m, i);
getMaxPartialDegree(m::ExpVec{UInt32,N}, i::Integer)   where {N} = getMaxPartialDegree_uint32(m, i);
getMaxPartialDegree(m::ExpVec{UInt64,N}, i::Integer)   where {N} = getMaxPartialDegree_uint64(m, i);
getMaxPartialDegree(m::ExpVec{UInt128,N}, i::Integer)   where {N} = getMaxPartialDegree_uint128(m, i);

setPartialDegree_uint8!(m::ExpVec{UInt8,N}, i::Integer, val::Integer)     where {N} = _setExp!(m, i, val, _uint8MasksOffs, _nvarToTableIndex);
setPartialDegree_uint16!(m::ExpVec{UInt16,N}, i::Integer, val::Integer)   where {N} = _setExp!(m, i, val, _uint16MasksOffs, _nvarToTableIndex);
setPartialDegree_uint32!(m::ExpVec{UInt32,N}, i::Integer, val::Integer)   where {N} = _setExp!(m, i, val, _uint32MasksOffs, _nvarToTableIndex);
setPartialDegree_uint64!(m::ExpVec{UInt64,N}, i::Integer, val::Integer)   where {N} = _setExp!(m, i, val, _uint64MasksOffs, _nvarToTableIndex);
setPartialDegree_uint128!(m::ExpVec{UInt128,N}, i::Integer, val::Integer) where {N} = _setExp!(m, i, val, _uint128MasksOffs, _nvarToTableIndex);

_setPartialDegree!(m::ExpVec{UInt8,N}, i::Integer, val::Integer)   where {N} = setPartialDegree_uint8!(m, i, val);
_setPartialDegree!(m::ExpVec{UInt16,N}, i::Integer, val::Integer)   where {N} = setPartialDegree_uint16!(m, i, val);
_setPartialDegree!(m::ExpVec{UInt32,N}, i::Integer, val::Integer)   where {N} = setPartialDegree_uint32!(m, i, val);
_setPartialDegree!(m::ExpVec{UInt64,N}, i::Integer, val::Integer)   where {N} = setPartialDegree_uint64!(m, i, val);
_setPartialDegree!(m::ExpVec{UInt128,N}, i::Integer, val::Integer)   where {N} = setPartialDegree_uint128!(m, i, val);


function setPartialDegree!(m::ExpVec{T,N}, idx::Integer, val::Integer) where {T,N}
	maxDeg = getMaxPartialDegree(m, idx);
	if (maxDeg < val)
		throw(OverflowError("The partial degree " * string(val) *
			" is too large for the exponent vector of type: " *
			string(typeof(m))));
	end

	return _setPartialDegree!(m, idx, val)
end


function fillExpVec!(m::ExpVec{T,N}, exps::Vector{I}) where {T,N,I<:Integer}
	m.vec = T(0); #reset the vector
	for idx in 1:min(N,length(exps))
		setPartialDegree!(m, idx, exps[idx]);
	end
	return m
end

function unpackExpVec(m::ExpVec{T,N}) where {T,N}
	vals = zeros(T,N);
	for i in 1:N
		vals[i] = getExp(m,i);
	end
	return vals;
end

function addExpVecs(m1::ExpVec{T,N}, m2::ExpVec{T,N}) where {T,N}
	return ExpVec{T,N}(m1.vec + m2.vec)
end

addExpVecs(m1::ExpVec{T1,N}, m2::ExpVec{T2,N}) where {T1, T2, N} = addExpVecs(promote(m1,m2)...)
+(m1::ExpVec{T1,N}, m2::ExpVec{T2,N}) where {T1,T2,N} = addExpVecs(m1,m2);

function subExpVecs(m1::ExpVec{T,N}, m2::ExpVec{T,N}) where {T,N}
	return ExpVec{T,N}(m1.vec - m2.vec)
end
subExpVecs(m1::ExpVec{T1,N}, m2::ExpVec{T2,N}) where {T1, T2, N} = subExpVecs(promote(m1,m2)...)
-(m1::ExpVec{T1,N}, m2::ExpVec{T2,N}) where {T1,T2,N} = subExpVecs(m1,m2);


## Promotion rules for encoding type ##

import Base.promote_rule;
promote_rule(::Type{ExpVec{T1,N}}, ::Type{ExpVec{T2,N}}) where {T1<:Unsigned,T2<:Unsigned,N} = ExpVec{promote_type(T1,T2),N}
# promote_rule(::Type{ExpVec{T1,N1}}, ::Type{ExpVec{T2,N2}}) where {T1<:Unsigned,T2<:Unsigned,N1,N2} = ExpVec{promote_type(T1,T2),max(N1,N2)}
## TODO Promotion rules for number of variables??? ##



## OUTPUT functions

function Base.show(io::IO, exp::ExpVec{T,N}) where {T,N}
	compact = get(io, :compact, true)
	print_object(io, exp, compact)
end

function Base.show(io::IO, mime::MIME"text/plain", exp::ExpVec{T,N}) where {T,N}
	compact = get(io, :compact, false)
	print_object(io, exp, compact)
end

function print_object(io::IO, exp::ExpVec{T,N}, compact::Bool) where {T,N}
	if compact
		print(io, convert(Vector{Int}, unpackExpVec(exp)))
	else
		print(io, typeof(exp));
		print(io, ":\n  ")
		print(io, convert(Vector{Int}, unpackExpVec(exp)))
	end
end

# function printMonomial(io::IO, exp::ExpVec{T,N}, vars::Vector{String}, withZeros::Bool = false) where {T,N}
# 	exps = unpackExpVec(exp)
# 	for (v,e) in zip(vars,exps)
# 		if withZeros || e > 0
# 			print(io, v * "^" * string(e))
# 		end
# 	end
# end

function printMonomial(exp::ExpVec{T,N}, vars::Vector{S}, withZerosAndOnes::Bool = false) where
	{T,N,S<:Union{Symbol, AbstractChar, AbstractString}}

	printMonomial(stdout, exp, map(string,vars));
end

function printMonomial(io::IO, exp::ExpVec{T,N}, vars::Vector{S}, withZerosAndOnes::Bool = false) where
	{T,N,S<:Union{Symbol, AbstractChar, AbstractString}}

	exps = unpackExpVec(exp)
	needsMult = false;
	for (v,e) in zip(vars,exps)
		if withZerosAndOnes || e > 1
			if needsMult
				print(io, "*")
				needsMult = false;
			end
			print(io, string(v) * "^" * string(e))
		elseif e == 1
			if needsMult
				print(io, "*")
				needsMult = false;
			end
			print(io, string(v))
			needsMult = true;
		end
	end
	# printMonomial(io, exp, map(string,vars));
end




## Metaprogramming
# for i in 1:32
# 	symbolname = Symbol("get_exp_$(i)")
# 	@eval function $(symbolname)(m::ExpVec, i::Integer, n::Integer)
# end

end;