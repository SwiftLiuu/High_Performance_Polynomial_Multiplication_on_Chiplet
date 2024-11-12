

"""
# module CircularBuffers

- Julia version:
- Author: Alex Brandt
- Date: 2024-06-17

# Examples

```jldoctest
julia>
```

"""
module CircularBuffers

import Base: lock, unlock

export CircularBuffer, BufferNode, CBuff, CBuffNode, write!, read!, capacity, is_full, close!, is_closed

using ..SparsePolynomials.ExpVecs

struct BufferNode{C,T}
    coefficient::C
    monomial::ExpVec{T,N} where {N}
end

const DEFAULT_BUFFER_SIZE::Int = 4096

mutable struct CircularBuffer{C,T}
	r::Int
    buffer::Vector{Union{BufferNode{C,T}, Nothing}}
    N::Int
    closed::Bool
    isFull::Bool
    lock::ReentrantLock
    w::Int


    CircularBuffer{C,T}() where {C,T} = new(1, fill(nothing, DEFAULT_BUFFER_SIZE), DEFAULT_BUFFER_SIZE, false, false, ReentrantLock(), 1)
end
CBuff{C} = CircularBuffer{C,UInt64}
CBuffNode{C} = BufferNode{C,UInt64}



function write!(b::CircularBuffer{C,T}, coefficient::C, monomial::ExpVec{T,N}) where {C,T,N}
    # lock(b.lock) do
        if is_full(b)
            throw("Buffer is full")
        end
        b.buffer[b.w] = BufferNode(coefficient, monomial)
        b.w = mod1(b.w + 1, b.N)
        b.isFull = (b.w == b.r)
    # end
end

function Base.read!(b::CircularBuffer{C,T})::BufferNode{C,T} where {C,T}
    # value = nothing
    # lock(b.lock) do
        if isempty(b)
            throw("Buffer is empty")
        end
        b.isFull = false
        value = b.buffer[b.r]

        b.r = mod1(b.r + 1, b.N)
    # end
    return value
end


function Base.isempty(b::CircularBuffer{C,T})::Bool where {C,T}
    # ret = false;
    # lock(b.lock) do
        ret = !b.isFull && (b.r == b.w)
    # end
    return ret
end

function is_full(b::CircularBuffer{C,T})::Bool where {C,T}
	# ret = false;
    # lock(b.lock) do
    	ret = b.isFull
    # end
    return ret;
end

@inline function capacity(b::CircularBuffer{C,T})::Int where {C,T}
	return b.N
end

function Base.size(b::CircularBuffer{C,T})::Int where {C,T}
    ret = b.N
    # println(Threads.threadid(), " N:", b.N, " w:", b.w, " r:", b.r, " f:", b.isFull)
    # lock(b.lock) do
    	if !b.isFull
	    	ret = (b.w - b.r) + (b.N * Int(b.w < b.r))
	    end
    # end
    return ret
end

function close!(b::CircularBuffer{C,T}) where {C,T}
    # lock(b.lock) do
        b.closed = true
    # end
end

function is_closed(b::CircularBuffer{C,T})::Bool where {C,T}
    # ret = false;
    # lock(b.lock) do
        ret = b.closed
    # end
    return ret
end


end #module