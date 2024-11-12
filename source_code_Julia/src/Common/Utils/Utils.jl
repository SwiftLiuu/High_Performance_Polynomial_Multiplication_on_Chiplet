"""
# Utils file

- Julia version:
- Authors: Alex Brandt, Mars Semenova
- Date: 2024-05-30

# Examples

```jldoctest
julia>
```

NEXT STEPS (TODO)
- calculate a good value for the upper bound of R_NTERMS (inclusive) based on the
upper bounds (inclusive) of R_NVARS and R_SPARSITY
- validate preset bounds
"""

module Utils

export ObjectStream
export R_NVARS, R_NTERMS, R_SPARSITY, BIG_INT_MAX

# consts (TODO: verify vals)
R_NVARS::Pair = Pair(1, 10) # TODO: -10 good range
R_NTERMS::Pair = Pair(0, 200) # TODO: based on R_NVARS, R_SPARSITY
R_SPARSITY::Pair = Pair(2, 100) # TODO: -100 good range
BIG_INT_MAX = big"2"^256 # TODO: good max

"""
Helper function to validate a param passed to buildRandomPoly().

@param param - Parameter to validate.
@param defRange - Associated default range for the parameter.
@return Whether param is valid.
"""
_isValidParam(param, range::Pair)::Bool = param isa Int && param >= range.first && param <= range.second

"""
Helper function which makes a 2D vector a proper matrix by filling in values with
the type's zero element.

@param listOfLists - 2D Vector to expand.
@return Altered listOfLists.
"""
function _expandRaggedArray!(listOfLists::Vector{Vector{T}}) where {T}
    nCols = maximum(length.(listOfLists))

    for l in listOfLists
        if length(l) < nCols
            append!(l, [zero(T) for _ in length(l)+1:nCols]);
        end
    end

    return listOfLists
end

"""
Helper method to trim a Vector.

@param vec - Vector to trim.
@param nsize - Size of trimmed array.
@return Trimmed vector.
"""
function _truncate!(vec, nsize::Int)
	if nsize >= 0 && nsize < size(vec, 1)
		resize!(vec, nsize)
    end

	return vec
end

using DataStructures
using Base

export ObjectStream

"""
ObjectStream{T} provides a simple thread-safe read-only queue (or "stream") of objects.
Suitable for having multiple consumers simultaneously retrieve (and process) objects from a central list.

@see getTask
"""
struct ObjectStream{T}
    queue::Deque{T}
    lock::ReentrantLock

    function ObjectStream(objs::Vector{T}) where {T}
        queue = Deque{T}()
        for obj in objs
            push!(queue, obj)
        end
        tq = new{T}(queue, ReentrantLock())
        return tq
    end
end

"""
Retrieve the next object from the ObjectStream or `nothing` if the stream is empty.
"""
function Base.pop!(stream::ObjectStream{T}) where {T}
    if (isempty(stream.queue))
        return nothing
    end
    lock(stream.lock)
    obj = nothing
    if !isempty(stream.queue)
        obj = popfirst!(stream.queue)
    end
    unlock(stream.lock)
    return obj
end

end # end module
