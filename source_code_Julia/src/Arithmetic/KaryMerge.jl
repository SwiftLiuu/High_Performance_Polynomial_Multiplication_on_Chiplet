


"""
# module KaryMarge

Implements a k-way merge of sorted lists using a heap.
Each sorted list is a represented by a linked list: a HeapNode head of the list
and an implementation of getnext(::HeapNode).

- Julia version:
- Author: Alex Brandt
- Date: 2024-06-17

# Examples

```jldoctest
julia>
```
"""
module KaryMerge

export HeapNode, getnext, KaryMerger, removeAndTryInsert!
using DataStructures


abstract type HeapNode end



function Base.isless(a::HeapNode, b::HeapNode)
	return a < b;
end

function getnext(a::HeapNode)::Union{Nothing,HeapNode}
	return nothing
end



"""
	KaryMerger facilitates k-ary merge of sorted lists.
	This implementation is specialized to support
	up to and including a 7-way merge
	(i.e. a complete binary tree with 2 levels).
"""
mutable struct KaryMerger{T<:HeapNode}
	entries::Vector{Union{Nothing,T}}
	size::Int

	KaryMerger{T}() where {T} = new(fill(nothing, 7), 0)
end

function _swim!(heap::KaryMerger{T}, idx::Int) where {T<:HeapNode}
	parentIdx = idx >> 1;
	if parentIdx > 0 && isless(heap.entries[parentIdx], heap.entries[idx])
		temp = heap.entries[parentIdx]
		heap.entries[parentIdx] = heap.entries[idx]
		heap.entries[idx] = temp
		idx = parentIdx
		parentIdx = parentIdx >> 1
	end

	if parentIdx > 0 && isless(heap.entries[parentIdx], heap.entries[idx])
		temp = heap.entries[parentIdx]
		heap.entries[parentIdx] = heap.entries[idx]
		heap.entries[idx] = temp
	end
end

function _sinkRoot!(heap::KaryMerger{T}) where {T<:HeapNode}
	if heap.size == 1
		return
	elseif heap.size == 2
		if isless(heap.entries[1], heap.entries[2])
			temp = heap.entries[1]
			heap.entries[1] = heap.entries[2]
			heap.entries[2] = temp;
		end
		return
	end

	curIdx = 1;
	nextIdx = 2;
	if isless(heap.entries[2], heap.entries[3])
		nextIdx = 3;
	end

	if isless(heap.entries[curIdx], heap.entries[nextIdx])
		temp = heap.entries[curIdx]
		heap.entries[curIdx] = heap.entries[nextIdx]
		heap.entries[nextIdx] = temp;
		curIdx = nextIdx
		nextIdx = nextIdx << 1
	else
		return
	end

	if heap.size >= nextIdx + 1 && isless(heap.entries[nextIdx], heap.entries[nextIdx+1])
		nextIdx = nextIdx + 1
	end

	if heap.size >= nextIdx && isless(heap.entries[curIdx], heap.entries[nextIdx])
		temp = heap.entries[curIdx];
		heap.entries[curIdx] = heap.entries[nextIdx]
		heap.entries[nextIdx] = temp;
	end
end


@inline function Base.insert!(heap::KaryMerger{T}, elem::T) where {T<:HeapNode}
	if heap.size == 7
		error("KaryMerger is full and no elements can be inserted.");
	end

	heap.size += 1
	heap.entries[heap.size] = elem
	_swim!(heap, heap.size)
end

@inline Base.isempty(heap::KaryMerger{T}) where {T} = heap.size == 0

@inline Base.length(heap::KaryMerger{T}) where {T} = heap.size

function peek(heap::KaryMerger{T}) where {T<:HeapNode}
	if heap.size == 0
		return nothing;
	end

	return heap.entries[1]
end

function removeAndTryInsert!(heap::KaryMerger{T}) where {T<:HeapNode}
	if heap.size == 0
		return (nothing, false)
	end

	# println("before")
	# println(heap.entries[1].exp, " ", heap.entries[2].exp, " ", heap.entries[3].exp, " ", heap.entries[4].exp)

	ret = heap.entries[1]
	nextNode = getnext(heap.entries[1])
	success = !isnothing(nextNode)
	if !success
		# @assert false
		nextNode = heap.entries[heap.size]
		heap.entries[heap.size] = nothing
		heap.size -= 1
	end
	# println("next: ", nextNode.exp)

	heap.entries[1] = nextNode
	if heap.size > 1
		_sinkRoot!(heap)
	end
	# println("after")
	# println(heap.entries[1].exp, " ", heap.entries[2].exp, " ", heap.entries[3].exp, " ", heap.entries[4].exp)
	# println(heap.entries[1].exp, " ", heap.entries[2].exp)
	# println("")

	return (ret, success)
end


###########
# Interface implementation for BinaryMaxHeap
###########


@inline function Base.insert!(heap::BinaryMaxHeap{T}, elem::T) where {T<:HeapNode}
	push!(heap, elem)
end

peek(heap::BinaryMaxHeap{T}) where {T<:HeapNode} = first(heap)

function removeAndTryInsert!(heap::BinaryMaxHeap{T}) where {T<:HeapNode}
	if length(heap) == 0
		return (nothing, false)
	end

	ret = pop!(heap)
	nextNode = getnext(ret)
	success = !isnothing(nextNode)
	if success
		push!(heap, nextNode)
	end

	return (ret,success)
end









end #module

# Quick and dirty debugging
# using .KaryMerge

# struct MyNode <: HeapNode
# 	a::Int
# end

# Base.isless(a::MyNode, b::MyNode) = a.a < b.a

# heap = KaryMerger{MyNode}()
# insert!(heap, MyNode(5))
# insert!(heap, MyNode(3))
# insert!(heap, MyNode(2))
# insert!(heap, MyNode(5))
