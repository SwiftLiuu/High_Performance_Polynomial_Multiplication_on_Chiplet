


"""
# module CircularBuffers

An implementation of circular buffers (i.e. ring buffers) using Julia `Channel`s.

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

using ..ExpVecs

struct BufferNode{T}
    coefficient::BigInt
    monomial::ExpVec{T,N} where {N}
end

const DEFAULT_BUFFER_SIZE::Int = 4096

mutable struct CircularBuffer{T}
    buffer::Channel{BufferNode{T}}
    # N::Int
    # r::Int
    # w::Int
    closed::Bool
    # isFull::Bool
    # lock::ReentrantLock

    CircularBuffer{T}() where {T} = new(Channel{BufferNode{T}}(DEFAULT_BUFFER_SIZE), false)
end
CBuff = CircularBuffer{UInt64}
CBuffNode = BufferNode{UInt64}



function write!(b::CircularBuffer{T}, coefficient::BigInt, monomial::ExpVec{T,N}) where {T,N}
    put!(b.buffer, BufferNode{T}(coefficient, monomial))
end

function Base.read!(b::CircularBuffer{T})::BufferNode{T} where {T}
    if isready(b.buffer)
        return take!(b.buffer)
    else
        return nothing
    end
end


function Base.isempty(b::CircularBuffer{T})::Bool where {T}
    return length(b.buffer.data) == 0
end

function is_full(b::CircularBuffer{T})::Bool where {T}
	return length(b.buffer.data) >= b.buffer.sz_max
end

@inline function capacity(b::CircularBuffer{T})::Int where {T}
	return b.buffer.sz_max
end

function Base.size(b::CircularBuffer{T})::Int where {T}
    length(b.buffer.data)
end

function close!(b::CircularBuffer{T}) where {T}
    close(b.buffer)
    b.closed = true;
end

function is_closed(b::CircularBuffer{T})::Bool where {T}
    return b.closed
end


end #module