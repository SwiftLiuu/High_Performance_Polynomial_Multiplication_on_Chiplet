
#######################
# Functions related to multiplication of SparsePoylnomials specialized to chiplets.
# This internal module should not be exported, but rather
# individual functions should be exported from the SparsePolynomials module.
#
# This file expects to be included from within the SparsePolyMult module.
#
# Author: Alex Brandt
# Date: 2024-06-17
#
#######################

using .ExpVecs
using .CircularBuffers

using Base.Threads

using DataStructures

struct EdgeInfo
    minIdxList::Vector{Int}
    maxIdxList::Vector{Int}
end

# function _expInRange(exp::ExpVec{T,N}, start::ExpVec{T,N})

function _findEdge(a::SMP{C,N}, b::SMP{C,N}, minExp::ExpVec{T,N}, maxExp::ExpVec{T,N}) where {C,T,N}

    nA = size(a);
    nB = size(b);
    minIdxList = [-1 for _ in 1:nA];
    maxIdxList = [-1 for _ in 1:nA];

    #be inclusive if minExp is the constant term
    minCmp = >
    if (iszero(minExp))
        minCmp = >=
    end

    aTerms = a.terms;
    bTerms = b.terms;

    #find exp such that startExp >= exp > endExp
    # recall that a,b are sorted in decreasing order so startExp is actually the larger
    startIdx = nB;
    foundAnything = true;
    for i in 1:nA
        for j in startIdx:-1:1
            prodExp = aTerms[i].exp + bTerms[j].exp
            # println("i: ", i, " j: ", j, " prod: ", prodExp)
            if maxIdxList[i] == -1 && minCmp(prodExp, minExp)
                maxIdxList[i] = j;
                startIdx = j;
            end

            if prodExp <= maxExp
                minIdxList[i] = j;
            else
                #if prodExp is > startExp, we can stop processing
                #this row since A[i]*b[j-1] must also be bigger than prodExp.
                break
            end
        end #j

        #if, in the previous j loop, every prodExp was < endExp, then we can stop
        #processing as A[i+1]*b[j] must also be smaller than A[i]*b[j]
        if maxIdxList[i] == -1
            break
        end
    end #i


    return EdgeInfo(minIdxList, maxIdxList)



end

function calculateJ0i(i, na, nb, l)
    return 1 + ((div(i, div(na, l)) % 2) * div(nb, 2*l))
end

function remove_duplicates(vec::Vector{T}) where {T}
    unique = Vector{T}()
    for elem in vec
        if !(elem in unique)
            push!(unique, elem)
        end
    end
    return unique
end

function _gridOfExps(l::Int, a::SMP{C,N}, b::SMP{C,N}) where {C,N}
    na = length(a)
    nb = length(b)
    na_l = div(na,l)
    nb_l = div(nb,l)
    aTerms = a.terms;
	bTerms = b.terms;

    T = typeof(aTerms[1].exp).parameters[1]
    SStar = ExpVec{T, N}[]

    for i = 1:na_l:na
        for j = calculateJ0i(i, na, nb, l):nb_l:nb
            # println("(", i, ", ", j, ")")
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
    smallestExponent = zero(largestExponent)

    if SetS[1] != largestExponent
        prepend!(SetS, [largestExponent])
    end
    if SetS[end] != smallestExponent
        push!(SetS, smallestExponent)
    end

    # println(SetS)
    return SetS
end


@inline function _isInsideEdge(edge::EdgeInfo, row::Int, idx::Int)
    return edge.minIdxList[row] > 0 &&
            edge.minIdxList[row] <= idx && idx <= edge.maxIdxList[row]
end


function _prodHeapInit(a::SMP{C,N}, b::SMP{C,N}, edge::EdgeInfo) where {C,N}
    heap = SparsePolyMult._prodHeapInit(a,b)

    for i in eachindex(edge.minIdxList)
        if _isInsideEdge(edge, i, edge.minIdxList[i])
            SparsePolyMult._prodHeapInsert!(heap, a, i, b, edge.minIdxList[i])
        end
    end

    return heap
end

###
#Initialize the product heap using edgeinfo to find the first term of each "stream" (a_i * b).
#But, only use stream i if (i mod nThreads) + 1  == threadIdx. +1 since indices are 1-based in julia.
#Use nThreads as the stepsize between streams to be handled by this product heap: a_i and a_{i+nThreads}
###
function _prodHeapInit(a::SMP{C,N}, b::SMP{C,N}, edge::EdgeInfo, threadIdx::Int, nThreads::Int) where {C,N}
    heap = SparsePolyMult._prodHeapInit(a,b)
    for i in eachindex(edge.minIdxList)
        if (i % nThreads) + 1 == threadIdx
            if _isInsideEdge(edge, i, edge.minIdxList[i])
                SparsePolyMult._prodHeapInsert!(heap, a, i, b, edge.minIdxList[i])
            end
        end
    end

    return heap
end

function _multiplySMP_TRIP(a::SMP{C,N}, b::SMP{C,N}) where {C,N}

    l = 3;
    S = _gridOfExps(l, a, b)
    prods = SMP{C,N}[];
    for i in 1:length(S)-1
        p = _multiplySMP_range(a,b,S[i+1],S[i])
        push!(prods, p)
    end

    return SparsePolynomials._concatenatePolys(prods);

end

function _parseThreadConfig(threadConfig)
    nThreads = 0
    if threadConfig isa Vector{Int}
        nThreads = sum(threadConfig)
    elseif threadConfig isa Integer && threadConfig > 0
        nThreads = threadConfig
        threadConfig = [nThreads] #TODO use hardware topology to find clusters
    else
        nThreads = length(Sys.cpu_info())
        threadConfig = [nThreads]
    end

    return (nThreads, threadConfig)

end

function _getGridStepSize(a::SMP{C,N}, b::SMP{C,N}, threadConfig::Vector{Int}) where {C,N}
    return 3; #TODO!!!!
end


struct ChipletTask
    taskIdx::Int
    minExp::ExpVec
    maxExp::ExpVec
end


"""

threadConfig: If threadConfig is a `Vector` of integers, then this is interpretted as
              requesting length(threadConfig) number of top-level threads t_i
              and threadConfig[i] number of helper threads associated with t_i (including counting t_i itself).
              E.g. [4,4,8] means 16 threads in total will be used: two groups of 4 and one group of 8.
              If threadConfig is a positive integer, this method will use a total of threadConfig number of threads,
              potentially organized into groups based on the topology of the executing processor(s).
              If thread config is a non-positive integer or anything else, use as many threads as possible.

"""
function _multiplySMP_Chiplet(a::SMP{C,N}, b::SMP{C,N}, threadConfig) where {C,N}

    nThreads, threadConfig = _parseThreadConfig(threadConfig)
    println("Thread config: ", nThreads, " ", threadConfig)
    l = _getGridStepSize(a,b,threadConfig)
    S = _gridOfExps(l, a, b)
    results = Vector{SMP{C,N}}(undef, length(S)-1)

    count = 0;
    leaders = Int[]
    for i in 1:length(threadConfig)
        push!(leaders, count+1)
        count += threadConfig[i]
    end

    #Each task is the interval S[i] -> S[i+1] (where S[i] > S[i+1])
    tasks = ObjectStream([ChipletTask(i, S[i+1], S[i]) for i in 1:length(S)-1])

    #:static here is crucial. We cannot have task migration of the tasks being executed the leader threads
    asyncTasks = Task[]
    for i in 2:length(threadConfig)
        taskfunc() = _multiply_leaderThreadChiplet(a,b,tasks,results,leaders[i],threadConfig[i])
        t = Task(taskfunc)
        ccall(:jl_set_task_tid, Cvoid, (Any, Cint), t, leaders[i]-1)
        t = errormonitor(t)
        schedule(t)
        push!(asyncTasks, t)
    end
    _multiply_leaderThreadChiplet(a,b,tasks,results,leaders[1],threadConfig[1])
    for task in asyncTasks
        wait(task)
    end

    # @threads :static for i in 1:length(threadConfig)
    #     println("static thread ", i)
    #     _multiply_leaderThreadChiplet(a,b,tasks,results,leaders[i],threadConfig[i])
    # end
    return SparsePolynomials._concatenatePolys(results);

end

function _multiply_leaderThreadChiplet(a::SMP{C,N}, b::SMP{C,N}, tasks::ObjectStream{ChipletTask},
                                       results::Vector{SMP{C,N}}, leaderIdx::Int, nThreads::Int) where {C,N}

    #pin this thread to a particular CCX and get list of core IDs for the CCX
    #TODO pinthread(First core of CCX with index leaderIdx)

    coreIDs = [leaderIdx + i for i in 0:nThreads-1]
    task = pop!(tasks)
    while !isnothing(task)
        results[task.taskIdx] = _multiplySMP_SDMP_range(a,b,task.minExp,task.maxExp,coreIDs)
        task = pop!(tasks)
    end
end

## Broken, Was used for testing. use multiple worker threads to execute each TRIP task, but only one thread executing tasks.
# function _multiplySMP_SDMPEdge(a::SMP{C,N}, b::SMP{C,N}) where {C,N}

#     nThreads = 2;
#     l = 3;
#     S = _gridOfExps(l, a, b)
#     prods = SMP{C,N}[];
#     for i in 1:length(S)-1
#         p = _multiplySMP_SDMP_range(a, b, S[i+1], S[i], nThreads)
#         push!(prods, p)
#     end

#     return SparsePolynomials._concatenatePolys(prods);

# end

#multiply the subset of terms from a and b which produce product terms
# in the range [max, min).
function _multiplySMP_range(a::SMP{C,N}, b::SMP{C,N}, min::ExpVec{T,N}, max::ExpVec{T,N}) where {C,N,T}
    # println("multiple range: ", min, " ", max)
    edge = _findEdge(a,b,min,max)
    _multiplySMP_edgeTask(a, b, edge);
end

function _multiplySMP_edgeTask(a::SMP{C,N}, b::SMP{C,N}, edge::EdgeInfo) where {C,N}
    #assume size(a) < size(b) and they are non-zero.

    heap = _prodHeapInit(a,b,edge)
	allocC = size(a) * size(b);
	allocC = allocC > 60000000 ? 60000000 : allocC;
	c = SMP(a.nvar, a.vars, Vector{Term{C,N}}(undef, allocC))
	cTerms = c.terms;
	aTerms = a.terms;
	bTerms = b.terms;

	k = 1; #current index into c

	while size(heap) > 0
		nextDegs = SparsePolyMult._prodHeapPeek(heap)
		cTerms[k] = Term{C,N}(zero(C),nextDegs)

		maxElem = nothing
		while nextDegs == cTerms[k].exp
			oldMax = maxElem;
			maxElem = SparsePolyMult._prodHeapRemove!(heap)
			while !isnothing(maxElem)
				cTerms[k].coef += aTerms[maxElem.aidx].coef * bTerms[maxElem.bidx].coef

				#add the head of maxElem to oldMax list
				nextMax = maxElem.next
                if _isInsideEdge(edge, maxElem.aidx, maxElem.bidx + 1)
					maxElem.bidx += 1
					maxElem.next = oldMax
					oldMax = maxElem
				end

				maxElem = nextMax

			end #while maxElem

			maxElem = oldMax #reset maxElem to be head of linked list
			nextDegs = SparsePolyMult._prodHeapPeek(heap)

		end #while nextDegs

		#commit condensed like term to product
		if !iszero(cTerms[k].coef)
			k += 1
		end

		#insert all successors back to heap
		while !isnothing(maxElem)
			nextMax = maxElem.next
			maxElem.next = nothing
			SparsePolyMult._prodHeapInsert!(heap, maxElem, aTerms[maxElem.aidx].exp + bTerms[maxElem.bidx].exp)
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

function _multiplySMP_range(a::SMP{C,N}, b::SMP{C,N}, min::ExpVec{T,N}, max::ExpVec{T,N}) where {C,N,T}
    # println("multiple range: ", min, " ", max)
    edge = _findEdge(a,b,min,max)
    return _multiplySMP_edgeTask(a, b, edge);
end

struct SDMPLocalInfo{C<:Number,N,T}
    buffer::CBuff{C}
    bufferIdx::Int
    min::ExpVec{T,N}
    max::ExpVec{T,N}
    nThreads::Int
    globalHeap::SDMP_GlobalHeap{C,N}
    coreID::Int
end


function _multiplySMP_SDMP_range(a::SMP{C,N}, b::SMP{C,N}, min::ExpVec{T,N}, max::ExpVec{T,N}, coreIDs::Vector{Int}) where {C<:Number,N,T}
    edge = _findEdge(a,b,min,max)

    # coreIDs = [coreIDs[1]]
    nThreads = length(coreIDs)

    if nThreads == 1
       return _multiplySMP_range(a,b,min,max)
    end

	# nThreads = 1
	bufferSet = [CBuff{C}() for _ in 1:nThreads]
	buffRefs = [Ref(bufferSet[i]) for i in 1:nThreads]

	allocC = size(a) * size(b);
	allocC = allocC > 100000000 ? 100000000 : allocC;
	prod = SMP(a.nvar, a.vars, Vector{Term{C,N}}(undef, allocC))
	globalHeap = SDMP_GlobalHeap{C,N}(ReentrantLock(), Vector{Ref{CBuff{C}}}[], buffRefs, BinaryMaxHeap{SDMP_GlobalNode{C}}(), prod, allocC, 0)

    tasks = Task[];
    infos = Vector{SDMPLocalInfo}(undef, nThreads)
    for i in 2:nThreads
        infos[i] = SDMPLocalInfo(bufferSet[i], i, min, max, nThreads, globalHeap, coreIDs[i])
        taskfunc() = local_heap_merge_range(a,b,infos[i])
        t = Task(taskfunc)
        # t = @task local_heap_merge_range(a,b,info)
        ccall(:jl_set_task_tid, Cvoid, (Any, Cint), t, coreIDs[i]-1)
        t = errormonitor(t)
        push!(tasks, t)
        schedule(t)

    end
    infos[1] = SDMPLocalInfo(bufferSet[1], 1, min, max, nThreads, globalHeap, coreIDs[1])
    local_heap_merge_range(a,b,infos[1])
    for task in tasks
        Base._wait(task)
    end

    # @sync begin
    # 	for i in 1:nThreads
    #         info = SDMPLocalInfo(bufferSet[i], i, min, max, nThreads, globalHeap, i)
    # 		errormonitor(@spawn local_heap_merge_range(a, b, info))
    # 	end
    # end


    done = global_heap_merge(globalHeap, length(a) * length(b))
	resize!(globalHeap.poly.terms, globalHeap.polySize)

    # println("global heap poly: ", globalHeap.poly)
	return globalHeap.poly

end


function local_heap_merge_range(a::SMP{C,N}, b::SMP{C,N}, info::SDMPLocalInfo{C,N,T}) where {C,N,T}

    # coreID = info.coreID
    #TODO pin this thread onto the coreID

    #TODO don't call find edge here! Only once per leader.
    edge = _findEdge(a, b, info.min, info.max)
	buffer = info.buffer
    bufferIdx = info.bufferIdx
    nThreads = info.nThreads
    globalHeap = info.globalHeap

    localHeap = _prodHeapInit(a,b, edge, bufferIdx, nThreads)

	bufferCapacity = CircularBuffers.capacity(buffer)
	k = div(bufferCapacity, nThreads) #not to be confused with k in the serial algo; see 'k' in (Monagan & Pearce, 2009)

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

                #add the head of maxElem to oldMax list
				nextMax = maxElem.next
                if _isInsideEdge(edge, maxElem.aidx, maxElem.bidx + 1)
					maxElem.bidx += 1
					maxElem.next = oldMax
					oldMax = maxElem
				end

				maxElem = nextMax

			end #while maxElem

			maxElem = oldMax #reset maxElem to be head of linked list
			nextDegs = SparsePolyMult._prodHeapPeek(localHeap)

		end #while nextDegs

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

