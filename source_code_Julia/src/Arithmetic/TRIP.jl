

module TRIP

using ..ExpVecs
using ..SparsePoly
using DataStructures
using Base.Threads

mutable struct HeapChainNode
	aidx::Int
	bidx::Int
	next::Union{HeapChainNode,Nothing}
end

mutable struct HeapNode
	exp::ExpVec{T,N} where {T,N}
	chain::HeapChainNode
end

mutable struct ProductHeap
	entries::Vector{Union{HeapNode,Nothing}}
	heapSize::Int
	heapAlloc::Int
	lastB::Int #only needed for division
end

struct Task
    startIndex::ExpVec{T,N} where {T,N}
    endIndex::ExpVec{T,N} where {T,N}
    resultIndex::Int
end

Base.size(heap::ProductHeap) = heap.heapSize;


function _willProductOverflow(a::SMP{C}, b::SMP{C}) where {C}
	pdegsA = partialDegrees(a)
	pdegsB = partialDegrees(b)
	pdegsA = pdegsA .+ pdegsB

	#TODO
	return false;
end

function computeSetS(l::Int, a::SMP{C}, b::SMP{C}) where {C}
    na = length(a)
    nb = length(b)
    na_l = div(na,l)
    nb_l = div(nb,l)
    aTerms = a.terms;
	bTerms = b.terms;

    T = typeof(aTerms[1].exp).parameters[1]
    N = typeof(aTerms[1].exp).parameters[2]
    SStar = ExpVec{T, N}[]


    for i = 1:na_l:na
        for j = calculateJ0i(i, na, nb, l)+1:nb_l:nb
            push!(SStar, aTerms[i].exp + bTerms[j].exp)
        end
    end

    for j = 1:nb_l:nb
        push!(SStar, aTerms[end].exp + bTerms[j].exp)
    end

    for i = 1:na_l:na
        push!(SStar, aTerms[i].exp + bTerms[end].exp)
    end

    sort!(SStar, rev=true)

    SetS = remove_duplicates(SStar)

    largestExponent = aTerms[1].exp + bTerms[1].exp
    smallestExponent = aTerms[end].exp + bTerms[end].exp

    if SetS[1] != largestExponent
        prepend!(SetS, [largestExponent])
    end
    if SetS[end] != smallestExponent
        push!(SetS, smallestExponent)
    end

    return SetS
end

#Function that serves for computeSetS
function calculateJ0i(i, na, nb, l)
    return 1 + ((i ÷ (na ÷ l)) % 2) * div(nb, 2*l)
end

function remove_duplicates(array::Vector{ExpVec{UInt64, N}}) where N
    unique_array = Vector{ExpVec{UInt64, N}}()

    for elem in array
        if !(elem in unique_array)
            push!(unique_array, elem)
        end
    end

    return unique_array
end


function findEdge(a::SMP{C}, b::SMP{C}, startIndex::ExpVec{T,N}, endIndex::ExpVec{T,N})where {C,T,N}
    na = length(a)
    nb = length(b)
    aTerms = a.terms;
	bTerms = b.terms;

    Lmin = fill(-1, na)
    Lmax = fill(-1, na)

    for i in na:-1:1

        if i == na || Lmin[i + 1] == -1
            start = 1  # Julia arrays are 1-indexed
        else
            start = Lmin[i + 1]
        end

        for j in start:nb
            sumExp = aTerms[i].exp + bTerms[j].exp

            if sumExp < endIndex
                break
            end
            if sumExp <= startIndex
                if Lmin[i] == -1
                    Lmin[i] = j
                end
                Lmax[i] = max(Lmax[i], j)
            end
        end
    end
    return [Lmin, Lmax]
end

Base.isless(a::HeapNode, b::HeapNode) = a.exp < b.exp

function heapBasedMultiplication(a::SMP{C}, b::SMP{C}, edge)where {C}
    heap = BinaryMaxHeap{HeapNode}()  # Initialize the max heap for HeapNode
    aTerms = a.terms;
	bTerms = b.terms;


    N = typeof(aTerms[1]).parameters[2]
    c = SMP{C,N}(a.nvar, a.vars, Vector{Term{C, N}}())
	cTerms = c.terms;

    # Build initial heap from the edge array
    for i in 1:length(edge[1])
        if edge[1][i] != -1
            indexA = i
            indexB = edge[1][i]
            exp = aTerms[indexA].exp + bTerms[indexB].exp
            chain = HeapChainNode(indexA, indexB, nothing)
            push!(heap, HeapNode(exp,chain))
        end
    end

    # Process the heap and build the result
    while !isempty(heap)
        maxTerm = pop!(heap)
        Coef = aTerms[maxTerm.chain.aidx].coef * bTerms[maxTerm.chain.bidx].coef
        M = maxTerm.exp

        insertNext!(a, b, heap, maxTerm, edge)  # Insertion right after extraction

        # Aggregate terms with the same monomial
        while !isempty(heap) && top(heap).exp == M
            maxTerm = pop!(heap)
            Coef += aTerms[maxTerm.chain.aidx].coef * bTerms[maxTerm.chain.bidx].coef
            insertNext!(a, b, heap, maxTerm, edge)
        end

        # Add the aggregated coefficient and monomial to the result
        if Coef != 0
            term = Term(Coef, M)
            push!(cTerms, term)
        end
    end

    return c
end


function insertNext!(a, b, heap, node, edge)
    aTerms = a.terms
	bTerms = b.terms
    if node.chain.bidx + 1 <= edge[2][node.chain.aidx]
        indexA = node.chain.aidx
        indexB = node.chain.bidx + 1
        exp = aTerms[indexA].exp + bTerms[indexB].exp
        chain = HeapChainNode(indexA, indexB, nothing)
        nextNode = HeapNode(exp, chain)
        push!(heap, nextNode)
    end
end

#This is the 'main' function where multiple dispatch takes place
function _multiplyTRIP(a::SMP{C}, b::SMP{C}) where {C}
    if (iszero(a))
		return zero(b)
	elseif iszero(b)
		return zero(a)
	end

	#TODO check for overflow in product
	if _willProductOverflow(a,b)
		#TODO replace with correct value of T
		throw(OverflowError("The product of these two polynomials will overflow the exponent vector of type T"))
	end

	#swap so a is the smaller element
	if size(a) > size(b)
		temp = a
		a = b
		b = temp
	end


    l = 2 # l determines the number of intervals (tasks), and should be tailored for different cases
    numberOfThreads = 2 # numberOfThreads should be tailored for different cases as well
    tasks = Queue{Task}()

    SetS = computeSetS(l, a, b)
    N = a.nvar
    expOne = zero(a.terms[1].exp)
    setPartialDegree!(expOne, N, 1)

    for i in 1:length(SetS) - 1
        enqueue!(tasks, Task(SetS[i], SetS[i + 1]+expOne, i))
    end

    containers = Array{Any}(undef, length(SetS) - 1)

    # Multithreading task execution
    @sync begin  # Ensure all spawned tasks complete before proceeding
        for i in 1:numberOfThreads
            @spawn begin
                while !isempty(tasks)
                    task = try
                        dequeue!(tasks)
                    catch e
                        break  # Exit if the queue is empty
                    end

                    # Task execution
                    edge = findEdge(a, b, task.startIndex, task.endIndex)
                    result = heapBasedMultiplication(a, b, edge)
                    containers[task.resultIndex] = result

                end
            end
        end
    end

    prodTerms = Vector{Term{C,N}}()
    for i in 1:length(containers)
        append!(prodTerms, containers[i].terms)
    end

    result = zero(a)
    result.terms = prodTerms
    return result

end


end
#module end