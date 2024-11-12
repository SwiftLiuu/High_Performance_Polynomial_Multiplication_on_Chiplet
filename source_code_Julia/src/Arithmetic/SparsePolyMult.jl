#######################
# Functions related to multiplication of SparsePoylnomials.
# This internal module should not be exported, but rather
# individual functions should be exported from the SparsePolynomials module.
#
# Author: Alex Brandt
# Date: 2024-06-17
#
#######################

# include("CircularBuffer-Channel.jl")
include("CircularBuffer.jl")
using .CircularBuffers

include("KaryMerge.jl")
import .KaryMerge: KaryMerger, getnext, peek, removeAndTryInsert!

using Base.Threads

import Base: trylock, unlock

mutable struct HeapChainNode
	aidx::Int
	bidx::Int
	next::Union{HeapChainNode,Nothing}
end

mutable struct HeapNode{T,N}
	# exp::ExpVec{UInt64,5}
	exp::ExpVec{T,N}
	chain::HeapChainNode
end

mutable struct ProductHeap{T,N}
	entries::Vector{Union{HeapNode{T,N},Nothing}}
	heapSize::Int
	heapAlloc::Int
	lastB::Int #only needed for division
end

Base.size(heap::ProductHeap{T,N}) where {T,N} = heap.heapSize;

function _prodHeapInit(a::SMP{C,N}, b::SMP{C,N}) where {C,N}

	T = _getExpVecT(a)
	entries = Vector{Union{HeapNode{T,N},Nothing}}(nothing, size(a))
	heapSize = 0
	heapAlloc = size(a)
	return ProductHeap(entries, heapSize, heapAlloc, -1)
end

function _prodHeapPrint(heap::ProductHeap{T,N}) where {T,N}
	for i in 1:heap.heapSize
		print("Heap[$i]: ")
		node = heap.entries[i]
		print(node.exp, " ")
		chain = node.chain
		print("($(chain.aidx), $(chain.bidx))")
		chain = chain.next
		while !isnothing(chain)
			print(" -> ($(chain.aidx), $(chain.bidx))")
			chain = chain.next
		end
		println()
	end
end


###
# Increase the product heap's allocaton size.
###
function _prodHeapGrow!(heap::ProductHeap{T,N}, newSize::Int) where {T,N}
    if newSize > heap.heapAlloc
        resize!(heap.entries, newSize)
        heap.heapAlloc = newSize
    end
    return heap #for convenience
end

function _prodHeapInsert!(heap::ProductHeap{T,N}, a::SMP{C,N}, i::Integer, b::SMP{C,N}, j::Integer) where {C,N,T}
	_prodHeapInsert!(heap, HeapChainNode(i, j, nothing), a.terms[i].exp + b.terms[j].exp)
end

function _prodHeapInsert!(heap::ProductHeap{T,N}, chain::HeapChainNode, exp::ExpVec{T,N}) where {T,N}
	if heap.heapSize == 0
		heap.entries[1] = HeapNode(exp, chain);
		heap.heapSize = 1
		return
	end

	if heap.entries[1].exp == exp
		chain.next = heap.entries[1].chain
		heap.entries[1].chain = chain
		return
	end


	s = heap.heapSize;
	entries = heap.entries;

	#=
	 In a classic implementation of a binary heap insert,
	 you probably would insert the new entry at the bottom and swim it up.
	 Here, with chaining, we want to look for an existing node with matching exp first.
	 In the following while loop, we simulate a "swim" and keep track of the path as we go.
	 If we ever find a match, insert by adding on to a chain and return.
	 Otherwise, we continue until the "swim" finishes.
	 Once the swim finished, insert the new node where the swim stopped,
	 and "push" the entire path we followed down to maintain heap order.

	 i represents the index of the parent node and j represents the
	 index of the node being swam (i.e. the current insertion point)

	 *NOTE* This is based on some bit hacks for a 0-based indexing;
	 so i, j, path should always be +1 when actually accessing heap.entries.
 	=#
	i = (s-1) >> 1;
	j = s;
	path = UInt64(1);

	while j > 0
		if entries[i+1].exp == exp
			chain.next = entries[i+1].chain
			entries[i+1].chain = chain;
			return
		elseif entries[i+1].exp < exp
			path = path << 1;
			if !Bool(j & 1)
				#set the trailing bit to 1 distinguish left/right child
				path += 1
			end
			j = i;
			i = (i-1) >> 1
		else # entries[i+1].exp > exp
			break;
		end
	end

	#j is now the insertion point
	newNode = HeapNode(exp, chain);
	while (j <= s)
		j1 = j+1 #see above comment block
		temp = entries[j1]
		entries[j1] = newNode;
		newNode = temp;
		j = (j << 1) + 1 + (path & 1);
		path = path >> 1;
	end


	heap.heapSize = s+1;

end


function _prodHeapRemove!(heap::ProductHeap{T,N}) where {T,N}

	heap.heapSize -= 1;
	s = heap.heapSize

	entries = heap.entries;
	maxChain = entries[1].chain

	# bit-hacks for heap indices use 0-based indexing, use +1 for actual indexing into entries
	i = 0
	j = 1
	while (j < s)
		if j+1 < s && entries[j+1].exp < entries[j+2].exp
			j += 1
		end
		entries[i+1] = entries[j+1]
		i = j
		j = (j << 1) + 1;
	end

	j = (i-1) >> 1;
	while i > 0
		if entries[s+1].exp < entries[j+1].exp
			break
		end
		entries[i+1] = entries[j+1]
		i = j
		j = (j-1) >> 1;
	end

	#final swap outside of the loop
	#yes, s+1. This is technically out of bounds based on heap.heapSize, but the data is still there.
	entries[i+1] = entries[s+1]

	return maxChain
end

function _prodHeapPeek(heap::ProductHeap{T,N}) where {T,N}
	if heap.heapSize == 0
		return nothing
	else
		return heap.entries[1].exp
	end
end


function _willProductOverflow(a::SMP{C,N}, b::SMP{C,N}) where {C,N}
	pdegsA = partialDegrees(a)
	pdegsB = partialDegrees(b)
	pdegsA = pdegsA .+ pdegsB

	#TODO
	return false;
end

#TODO write a wrapper to ensure a and b have same var ordering
function _multiplySMP(a::SMP{C,N}, b::SMP{C,N}) where {C,N}

	if (iszero(a))
		return zero(b)
	elseif iszero(b)
		return zero(a)
	end

	# #TODO check for overflow in product
	# if _willProductOverflow(a,b)
	# 	#TODO replace with correct value of T
	# 	throw(OverflowError("The product of these two polynomials will overflow the exponent vector of type T"))
	# end

	#swap so a is the smaller element
	if size(a) > size(b)
		temp = a
		a = b
		b = temp
	end

	if length(b) == 1 && iszero(b.terms[1].exp) # mult by coef
		if length(a) == 1 && iszero(a.terms[1].exp) # mult 2 coefs
			return SparsePolynomial([Term{C,N}(a.terms[1].coef * b.terms[1].coef)])
		else
			return multPolyByScalar(a, b.terms[1].coef)
		end
	end

	heap = _prodHeapInit(a,b)
	chain = HeapChainNode(1, 1, nothing)
	exp = a.terms[1].exp + b.terms[1].exp
	# entries[1] = HeapNode(exp,chain)
	_prodHeapInsert!(heap, chain, exp)



	allocC = size(a) * size(b);
	c = allocPoly(C, a.vars, allocC)
	cTerms = c.terms;
	aTerms = a.terms;
	bTerms = b.terms;

	k = 1;
	lastA = size(a)
	lastB = size(b)
	firstB = 1 #first index of b terms

	while size(heap) > 0
		nextDegs = _prodHeapPeek(heap)
		cTerms[k] = Term{C,N}(zero(C),nextDegs)

		maxElem = nothing
		while nextDegs == cTerms[k].exp
			oldMax = maxElem;
			maxElem = _prodHeapRemove!(heap)
			while !isnothing(maxElem)
				cTerms[k].coef += aTerms[maxElem.aidx].coef * bTerms[maxElem.bidx].coef

				#if we extract the term corresponding to a_i*b_1, we need to insert a_{i+1}*b_1
				if maxElem.bidx == firstB && maxElem.aidx != lastA
					oldMax = HeapChainNode(maxElem.aidx+1, firstB, oldMax)
				end

				#add the head of maxElem to oldMax list
				nextMax = maxElem.next
				if maxElem.bidx != lastB
					maxElem.bidx += 1
					maxElem.next = oldMax
					oldMax = maxElem
				end

				maxElem = nextMax

			end #while maxElem

			maxElem = oldMax #reset maxElem to be head of linked list
			nextDegs = _prodHeapPeek(heap)

		end #while nextDegs

		#commit condensed like term to product
		if !iszero(cTerms[k].coef)
			k += 1
		end

		#insert all successors back to heap
		while !isnothing(maxElem)
			nextMax = maxElem.next
			maxElem.next = nothing
			_prodHeapInsert!(heap, maxElem, aTerms[maxElem.aidx].exp + bTerms[maxElem.bidx].exp)
			maxElem = nextMax
		end

		if k >= allocC
			#TODO resize C
		end

	end #while heap nonempty

	#TODO call the proper resize
	resize!(cTerms, k-1) #k-1 since 1-based indexing.

	return c

end



struct SDMP_GlobalNode{C<:Number} <: KaryMerge.HeapNode
	bufferRef::Ref{CBuff{C}}
	coef::C
	exp::ExpVec{UInt64,N} where {N}
end

@inline Base.isless(a::SDMP_GlobalNode, b::SDMP_GlobalNode) = a.exp < b.exp

function getnext(node::SDMP_GlobalNode{C}) where {C}
	buff = node.bufferRef[]
	if isempty(buff)
		return nothing
	end

	buffNode = read!(buff)
	if isnothing(buffNode)
		return nothing
	end

	return SDMP_GlobalNode(node.bufferRef, buffNode.coefficient, buffNode.monomial)
end

mutable struct SDMP_GlobalHeap{C<:Number,N}
	lock::ReentrantLock
	closedBuffers::Vector{Ref{CBuff{C}}}
	# allBuffers::Vector{Ref{CBuff}} #not needed right now
	emptyBuffers::Vector{Ref{CBuff{C}}} #buffers that are open but not in the merger yet
	merger #::KaryMerger{SDMP_GlobalNode{C}} #leave type undefined to swap out different mergers
	poly::SMP{C,N}
	polyAlloc::Int
	polySize::Int
end

trylock(heap::SDMP_GlobalHeap) = trylock(heap.lock)
unlock(heap::SDMP_GlobalHeap) = unlock(heap.lock)

function _multiplySMP_SDMP(a::SMP{C,N}, b::SMP{C,N}, nThreads::Int) where {C<:Number,N}

	# nThreads = 1
	bufferSet = [CBuff{C}() for _ in 1:nThreads]
	buffRefs = [Ref(bufferSet[i]) for i in 1:nThreads]

	allocC = size(a) * size(b);
	allocC = allocC > 100000000 ? 100000000 : allocC;
	prod = SMP(a.nvar, a.vars, Vector{Term{C,N}}(undef, allocC))
	globalHeap = SDMP_GlobalHeap{C,N}(ReentrantLock(), Vector{Ref{CBuff{C}}}[], buffRefs, KaryMerger{SDMP_GlobalNode{C}}(), prod, allocC, 0)

	#TODO
	#NO NO NO NO NO cannot use @threads for these non-independent tasks
	# @threads :static for i in 1:nThreads
	if nThreads == 1
		local_heap_merge(a,b,bufferSet,1,nThreads,globalHeap)
	else
		@sync begin
			for i in 1:nThreads
				errormonitor(@spawn local_heap_merge(a, b, bufferSet, i, nThreads, globalHeap))
			end
		end
	end

	global_heap_merge(globalHeap, length(a) * length(b))
	resize!(globalHeap.poly.terms, globalHeap.polySize)

	return globalHeap.poly

end

function local_heap_merge(a::SMP{C,N}, b::SMP{C,N}, buffers::Vector{CBuff{C}}, bufferIdx::Integer, nThreads::Int, globalHeap::SDMP_GlobalHeap{C,N}) where {C,N}

	buffer = buffers[bufferIdx]
	localHeap = _prodHeapInit(a,b)
	bufferCapacity = CircularBuffers.capacity(buffer)
	k = div(bufferCapacity, nThreads) #not to be confused with k in the serial algo
	lastA = size(a)
	lastB = size(b)
	firstB = 1 #first index of b terms

	startIdx = bufferIdx
	_prodHeapInsert!(localHeap, a, startIdx, b, 1);
	while size(localHeap) > 0
		nextDegs = _prodHeapPeek(localHeap)

		coef = zero(C)
		degs = nextDegs;

		maxElem = nothing
		while nextDegs == degs
			oldMax = maxElem;
			maxElem = _prodHeapRemove!(localHeap)

			while !isnothing(maxElem)
				coef += a.terms[maxElem.aidx].coef * b.terms[maxElem.bidx].coef
				#if we extract the term corresponding to a_i*b_1, we need to insert a_{i+nThreads}*b_1
				if maxElem.bidx == firstB && maxElem.aidx + nThreads <= lastA
					oldMax = HeapChainNode(maxElem.aidx+nThreads, firstB, oldMax)
				end

				#add the head of maxElem to oldMax list
				nextMax = maxElem.next
				if maxElem.bidx != lastB
					maxElem.bidx += 1
					maxElem.next = oldMax
					oldMax = maxElem
				end

				maxElem = nextMax

			end #while chain has elements

			maxElem = oldMax #reset maxElem to be head of linked list
			nextDegs = _prodHeapPeek(localHeap)
		end #while like terms

		#insert all successors back to heap
		while !isnothing(maxElem)
			nextMax = maxElem.next
			maxElem.next = nothing
			_prodHeapInsert!(localHeap, maxElem, a.terms[maxElem.aidx].exp + b.terms[maxElem.bidx].exp)
			maxElem = nextMax
		end

		#commit condensed like term to product
		if !iszero(coef)
			CircularBuffers.write!(buffer, coef, degs)
		end

		k -= 1

		while k == 0
			k = size(buffer)
			if trylock(globalHeap)
				madeProgress = global_heap_merge(globalHeap, min(k, div(bufferCapacity, nThreads)))
				unlock(globalHeap)
				if !madeProgress && k == bufferCapacity
					sleep(0.001)
				end
			elseif k == bufferCapacity
				sleep(0.001) #if buffer is full, try lock again
			end

			k = min(bufferCapacity - size(buffer), div(bufferCapacity,nThreads))
		end
	end #while heap nonempty

	close!(buffer)
end

function global_heap_merge(heap::SDMP_GlobalHeap{C,N}, maxIterations::Int) where {C,N}

	i = 0
	while i < maxIterations
		i += 1

		#ensure heap contains one element from each open buffer
		while !isempty(heap.emptyBuffers)
			buffRef = popfirst!(heap.emptyBuffers)
			if isempty(buffRef[])
				if is_closed(buffRef[])
					push!(heap.closedBuffers, buffRef)
				else
					push!(heap.emptyBuffers, buffRef)
					return (i!= 1)
				end
			else
				node = CircularBuffers.read!(buffRef[])
				insert!(heap.merger, SDMP_GlobalNode(buffRef, node.coefficient, node.monomial))

			end

		end

		if isempty(heap.merger)
			return (i != 1)
		end

		coef = C(0)
		exp = peek(heap.merger).exp
		while !isempty(heap.merger) && exp == peek(heap.merger).exp
			node, success = removeAndTryInsert!(heap.merger)
			coef += node.coef

			buff = node.bufferRef[]
			if !success
				#can NOT assert this. Buffer may have more in it after insert fails
				# @assert isempty(buff)
				if is_closed(buff)
					push!(heap.closedBuffers, node.bufferRef)
				else
					push!(heap.emptyBuffers, node.bufferRef)
				end
			end
		end

		if !iszero(coef)
			heap.polySize += 1
			if heap.polySize > heap.polyAlloc
				heap.polyAlloc += 10000000
				resize!(heap.poly.terms, heap.polyAlloc)
			end
			heap.poly.terms[heap.polySize] = Term(coef, exp)
			# push!(heap.poly, Term(coef, exp))
			# push is quite slow, so manually prealloc
			# https://stackoverflow.com/questions/34751225/how-efficient-are-push-and-append-methods-in-julia
		end
	end #while iterations
	return true;
end

"""
Helper function to allow addition functionality to be extended for division
by setting the operation between a polynomial and a scalar via an arg.

@param p - Polynomial.
@param s - Scalar. Must be the same type as p's coefficients.
@param op - Operation to apply to the coefficients of the passed polynomials
	when their terms' exps are the same.
@return Resulting SparsePolynomial{C,N}.
"""
function _iterWithScalarOp(p::SparsePolynomial{C,N}, s::C, op::Function) where {C,N}
	_assertCanonical(p)

	if s == C(0)
		if op(2, 2) == 4
			return zero(p)
		else
			throw(DivideError())
		end
	end

	if iszero(p)
		return zero(p)
	end

	if length(p) == 1 && iszero(p.terms[1].exp) # div 2 coefs
		np = SparsePolynomial([Term{C,N}(op(p.terms[1].coef, s))], p.vars)
		return np
	end

	np = SparsePolynomial(p)
	for x = 1:length(np)
		np.terms[x].coef = op(np.terms[x].coef, s)
	end

	return np
end

"""
Function which multiplies a polynomial with a scalar.

@param p - Polynomial.
@param s - Scalar. Must be the same type as p's coefficients.
@return Resulting SparsePolynomial{C,N}.
"""
function multPolyByScalar!(p::SparsePolynomial{C,N}, s::C) where {C,N}
    # TODO
    # TODO: replace implementation of multPolyByScalar w copy of inplace
end

function multPolyByScalar(p::SparsePolynomial{C,N}, s::C) where {C,N}
     return _iterWithScalarOp(p, s, ((c,x) -> c*x))
end

# Base defs
import Base.*
*(s::C, p::SparsePolynomial{C,N}) where {C,N} = multPolyByScalar(p, s)
*(p::SparsePolynomial{C,N}, s::C) where {C,N} = multPolyByScalar(p, s)

"""
Function which negated a polynomial.

@param p - Polynomial. Can't be of an Unsigned type.
@return Resulting SparsePolynomial{C,N}.
"""
function negate!(p::SparsePolynomial{C,N}) where {C,N}
    # TODO
    # TODO: replace implementation of negate w copy of inplace
end

function negate(p::SparsePolynomial{C,N}) where {C,N}
    return C(-1) * p
end

function negate(p::SparsePolynomial{C,N}) where {C<:Unsigned,N}
    _assertCanonical(p)

	if iszero(p)
		return zero(p)
	end

	np = SparsePolynomial(p)
	for x = 1:length(np)
		np.terms[x].coef = typemax(C) - np.terms[x].coef + 1
	end

	return np
end

"""
Function which negates a coefficient.

@param c - Coefficient of type C to negate.
@return Negated coefficient.
"""
function negateCoef(c::C) where {C} # TODO: inplace impl?
	return C(-1) * c
end

function negateCoef(c::C) where {C<:Unsigned}
	return C(0) - c
end

include("SparsePolyMult-Chiplet.jl")

