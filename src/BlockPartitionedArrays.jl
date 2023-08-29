
"""
"""
struct BlockPRange{A} <: AbstractUnitRange{Int}
  ranges::Vector{PRange{A}}
  function BlockPRange(ranges::Vector{<:PRange{A}}) where A
    new{A}(ranges)
  end
end

Base.first(a::BlockPRange) = 1
Base.last(a::BlockPRange) = sum(map(last,a.ranges))

BlockArrays.blocklength(a::BlockPRange) = length(a.ranges)
BlockArrays.blocksize(a::BlockPRange) = (blocklength(a),)
BlockArrays.blockaxes(a::BlockPRange) = (Block.(Base.OneTo(blocklength(a))),)
BlockArrays.blocks(a::BlockPRange) = a.ranges

"""
"""
struct BlockPArray{T,N,A,B} <: BlockArrays.AbstractBlockArray{T,N}
  blocks::Array{A,N}
  axes::NTuple{N,B}

  function BlockPArray(blocks::Array{<:AbstractArray{T,N},N},
                       axes  ::NTuple{N,<:BlockPRange}) where {T,N}
    @check all(map(d->size(blocks,d)==blocklength(axes[d]),1:N))
    A = eltype(blocks)
    B = typeof(first(axes))
    new{T,N,A,B}(blocks,axes)
  end
end

const BlockPVector{T,A,B} = BlockPArray{T,1,A,B}
const BlockPMatrix{T,A,B} = BlockPArray{T,2,A,B}

@inline function BlockPVector(blocks::Vector{<:PVector},rows::BlockPRange)
  BlockPArray(blocks,(rows,))
end

@inline function BlockPVector(blocks::Vector{<:PVector},rows::Vector{<:PRange})
  BlockPVector(blocks,BlockPRange(rows))
end

@inline function BlockPMatrix(blocks::Matrix{<:PSparseMatrix},rows::BlockPRange,cols::BlockPRange)
  BlockPArray(blocks,(rows,cols))
end

@inline function BlockPMatrix(blocks::Matrix{<:PSparseMatrix},rows::Vector{<:PRange},cols::Vector{<:PRange})
  BlockPMatrix(blocks,BlockPRange(rows),BlockPRange(cols))
end

# AbstractArray API

Base.axes(a::BlockPArray) = a.axes
Base.size(a::BlockPArray) = Tuple(map(length,a.axes))

Base.IndexStyle(::Type{<:BlockPVector}) = IndexLinear()
Base.IndexStyle(::Type{<:BlockPMatrix}) = IndexCartesian()

function Base.similar(a::BlockPVector,::Type{T},inds::Tuple{<:BlockPRange}) where T
  vals = map(blocks(a),blocks(inds[1])) do ai,i
    similar(ai,T,i)
  end
  return BlockPArray(vals,inds)
end

function Base.similar(::Type{<:BlockPVector{T,A}},inds::Tuple{<:BlockPRange}) where {T,A}
  rows   = blocks(inds[1])
  values = map(rows) do r
    return similar(A,(r,))
  end
  return BlockPArray(values,inds)
end

function Base.similar(a::BlockPMatrix,::Type{T},inds::Tuple{<:BlockPRange,<:BlockPRange}) where T
  vals = map(CartesianIndices(blocksize(a))) do I
    rows = inds[1].ranges[I[1]]
    cols = inds[2].ranges[I[2]]
    similar(a.blocks[I],T,(rows,cols))
  end
  return BlockPArray(vals,inds)
end

function Base.similar(::Type{<:BlockPMatrix{T,A}},inds::Tuple{<:BlockPRange,<:BlockPRange}) where {T,A}
  rows = blocks(inds[1])
  cols = blocks(inds[2])
  values = map(CartesianIndices((length(rows),length(cols)))) do I
    i,j = I[1],I[2]
    return similar(A,(rows[i],cols[j]))
  end
  return BlockPArray(values,inds)
end

function Base.getindex(a::BlockPArray{T,N},inds::Vararg{Int,N}) where {T,N}
  @error "Scalar indexing not supported"
end
function Base.setindex(a::BlockPArray{T,N},v,inds::Vararg{Int,N}) where {T,N}
  @error "Scalar indexing not supported"
end

function Base.show(io::IO,k::MIME"text/plain",data::BlockPArray{T,N}) where {T,N}
  v = first(blocks(data))
  s = prod(map(si->"$(si)x",blocksize(data)))[1:end-1]
  map_main(partition(v)) do values
      println(io,"$s-block BlockPArray{$T,$N}")
  end
end

function Base.zero(v::BlockPArray)
  return mortar(map(zero,blocks(v)))
end

function Base.copyto!(y::BlockPVector,x::BlockPVector)
  @check blocklength(x) == blocklength(y)
  for i in blockaxes(x,1)
    copyto!(y[i],x[i])
  end
  return y
end

function Base.copyto!(y::BlockPMatrix,x::BlockPMatrix)
  @check blocksize(x) == blocksize(y)
  for i in blockaxes(x,1)
    for j in blockaxes(x,2)
      copyto!(y[i,j],x[i,j])
    end
  end
  return y
end

function Base.fill!(a::BlockPVector,v)
  map(blocks(a)) do a
    fill!(a,v)
  end
  return a
end

# AbstractBlockArray API

BlockArrays.blocks(a::BlockPArray) = a.blocks

function Base.getindex(a::BlockPArray,inds::Block{1})
  a.blocks[inds.n...]
end
function Base.getindex(a::BlockPArray{T,N},inds::Block{N}) where {T,N}
  a.blocks[inds.n...]
end
function Base.getindex(a::BlockPArray{T,N},inds::Vararg{Block{1},N}) where {T,N}
  a.blocks[map(i->i.n[1],inds)...]
end

function BlockArrays.mortar(blocks::Vector{<:PVector})
  rows = map(b->axes(b,1),blocks)
  BlockPVector(blocks,rows)
end

function BlockArrays.mortar(blocks::Matrix{<:PSparseMatrix})
  rows = map(b->axes(b,1),blocks[:,1])
  cols = map(b->axes(b,2),blocks[1,:])

  function check_axes(a,r,c)
    A = PartitionedArrays.matching_local_indices(axes(a,1),r)
    B = PartitionedArrays.matching_local_indices(axes(a,2),c)
    return A & B
  end
  @check all(map(I -> check_axes(blocks[I],rows[I[1]],cols[I[2]]),CartesianIndices(size(blocks))))

  return BlockPMatrix(blocks,rows,cols)
end

# PartitionedArrays API

Base.wait(t::Array)  = map(wait,t)
Base.fetch(t::Array) = map(fetch,t)

function PartitionedArrays.assemble!(a::BlockPArray)
  map(assemble!,blocks(a))
end

function PartitionedArrays.consistent!(a::BlockPArray)
  map(consistent!,blocks(a))
end

function PartitionedArrays.partition(a::BlockPArray)
  vals = map(partition,blocks(a)) |> to_parray_of_arrays
  return map(mortar,vals)
end

function PartitionedArrays.to_trivial_partition(a::BlockPArray)
  vals = map(to_trivial_partition,blocks(a))
  return mortar(vals)
end

# LinearAlgebra API

function LinearAlgebra.mul!(y::BlockPVector,A::BlockPMatrix,x::BlockPVector)
  o = one(eltype(A))
  for i in blockaxes(A,2)
    fill!(y[i],0.0)
    for j in blockaxes(A,2)
      mul!(y[i],A[i,j],x[j],o,o)
    end
  end
end

function LinearAlgebra.dot(x::BlockPVector,y::BlockPVector)
  return sum(map(dot,blocks(x),blocks(y)))
end

function LinearAlgebra.norm(v::BlockPVector)
  block_norms = map(norm,blocks(v))
  return sqrt(sum(block_norms.^2))
end

function LinearAlgebra.fillstored!(a::BlockPMatrix,v)
  map(blocks(a)) do a
    fillstored!(a,v)
  end
  return a
end

# Broadcasting

struct BlockPBroadcasted{A,B}
  blocks :: A
  axes   :: B
end

BlockArrays.blocks(b::BlockPBroadcasted) = b.blocks
BlockArrays.blockaxes(b::BlockPBroadcasted) = b.axes

function Base.broadcasted(f, args::Union{BlockPVector,BlockPBroadcasted}...)
  a1 = first(args)
  @boundscheck @assert all(ai -> blockaxes(ai) == blockaxes(a1),args)
  
  blocks_in = map(blocks,args)
  blocks_out = map((largs...)->Base.broadcasted(f,largs...),blocks_in...)
  
  return BlockPBroadcasted(blocks_out,blockaxes(a1))
end

function Base.broadcasted(f, a::Number, b::Union{BlockPVector,BlockPBroadcasted})
  blocks_out = map(b->Base.broadcasted(f,a,b),blocks(b))
  return BlockPBroadcasted(blocks_out,blockaxes(b))
end

function Base.broadcasted(f, a::Union{BlockPVector,BlockPBroadcasted}, b::Number)
  blocks_out = map(a->Base.broadcasted(f,a,b),blocks(a))
  return BlockPBroadcasted(blocks_out,blockaxes(a))
end

function Base.broadcasted(f,
                        a::Union{BlockPVector,BlockPBroadcasted},
                        b::Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{0}})
  Base.broadcasted(f,a,Base.materialize(b))
end

function Base.broadcasted(
  f,
  a::Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{0}},
  b::Union{BlockPVector,BlockPBroadcasted})
  Base.broadcasted(f,Base.materialize(a),b)
end

function Base.materialize(b::BlockPBroadcasted)
  blocks_out = map(Base.materialize,blocks(b))
  return mortar(blocks_out)
end

function Base.materialize!(a::BlockPVector,b::BlockPBroadcasted)
  map(Base.materialize!,blocks(a),blocks(b))
  return a
end

