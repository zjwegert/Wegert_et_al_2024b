abstract type SchurMethod end

struct AdditiveSchurMethod <: SchurMethod end
struct ExactSchurMethod <: SchurMethod end

const MultiFieldSpaceTypes = Union{<:MultiField.MultiFieldFESpace,<:GridapDistributed.DistributedMultiFieldFESpace}
const TriangulationTypes = Union{<:Triangulation,<:GridapDistributed.DistributedTriangulation}

"""
    struct BlockSchurPreconditioner{A,B,C,D,E,F,G} <: Gridap.Algebra.LinearSolver

Create a block diagonal/triangular preconditioner solver for saddle-point
problems of the form
```
|  A    B | * | x | = | b |
| -Bᵀ   C |   | y |   | c |
```
For example, the block triangular preconditioner is
```
P = |  A    0 |
    | -Bᵀ   S |
```
where `S = C + BᵀA⁻¹B` is the Schur complement.

The Schur complement can be built either exactly or via an approximation using the
`ExactSchurMethod` or `AdditiveSchurMethod`, respectively. The first is dense
and should only be used for very small problem. The latter is built
using the additive Schur complement approximation discussed by Cao and Neytcheva, 2021
(10.1002/nla.2362). In particular, we locally compute `Sₑ = Cₑ + BₑᵀAₑ⁻¹Bₑ` for each
element and then assemble the sparse approximate Schur complement from `Sₑ`. Note
that `Aₑ` is usually singular and it is common to perturb the diagonal by `h²`
where `h` is the side length of an element. We further note that this needs to be
adjusted for higher order finite elements as the lengths between DoFs decreases as
the order increases.
"""
struct BlockSchurPreconditioner{A,B,C,D,E,F,G} <: Gridap.Algebra.LinearSolver
  schur_method    :: A
  schur_assembler :: B
  a               :: C
  spaces          :: D
  trian           :: E
  solver          :: F
  diag_perturb    :: G

  function BlockSchurPreconditioner(
      solver          :: Union{GridapSolvers.BlockTriangularSolver,GridapSolvers.BlockDiagonalSolver},
      a               :: Function,
      UP              :: MultiFieldSpaceTypes, # Currently assume this is ordering of spaces
      VQ              :: MultiFieldSpaceTypes,
      trian           :: TriangulationTypes,
      diag_perturb    :: Number;
      schur_method    :: SchurMethod = AdditiveSchurMethod(),
      schur_assembler :: Assembler = SparseMatrixAssembler(UP[2],VQ[2])
    )
    @check length(UP) == length(VQ) == 2
    spaces = (UP,VQ)
    # _diag_peturb = UniformScaling(diag_perturb);
    U_dof_length = _get_dof_length(first(VQ))
    _diag_peturb = copyto!(zeros(U_dof_length,U_dof_length),UniformScaling(diag_perturb))

    A = typeof(schur_method)
    B = typeof(schur_assembler)
    C = typeof(a)
    D = typeof(spaces)
    E = typeof(trian)
    F = typeof(solver)
    G = typeof(_diag_peturb)

    new{A,B,C,D,E,F,G}(schur_method,schur_assembler,a,spaces,trian,solver,_diag_peturb)
  end
end

function _get_dof_length(V)
  length(first(get_cell_dof_ids(V)))
end

function _get_dof_length(V::GridapDistributed.DistributedFESpace)
  U_dof_length = map(_get_dof_length,local_views(V))
  return PartitionedArrays.getany(U_dof_length)
end

struct BlockSchurPreconditionerSS <: Gridap.Algebra.SymbolicSetup
  solver
end

function Gridap.Algebra.symbolic_setup(solver::BlockSchurPreconditioner, K::AbstractBlockMatrix)
  @check blocksize(K) == (2,2) "A 2x2 block structure is expected"
  return BlockSchurPreconditionerSS(solver)
end

mutable struct BlockSchurPreconditionerNS <: Gridap.Algebra.NumericalSetup
  solver
  K
  P_ns
  S
end

## Lazy additive Schur complement
struct SchurMap <: Map
  diag_perturb
end

function Gridap.Arrays.return_cache(map::SchurMap,Ki)
  return copy(Ki[2,2])
end

function Gridap.Arrays.evaluate!(cache,map::SchurMap,Ki)
  Si = cache
  axpy!(1,map.diag_perturb,Ki[1,1])
  ldiv!(lu!(Ki[1,1]),Ki[1,2])
  mul!(Ki[2,2],Ki[2,1],Ki[1,2],-1,1)
  copyto!(Si,Ki[2,2])
  return Si
end

function schur_complement_dom_contribs(
    solver::BlockSchurPreconditioner{<:AdditiveSchurMethod},
    K::BlockMatrix
  )
  a = solver.a
  UP, VQ = solver.spaces
  Ω = solver.trian
  diag_perturb = solver.diag_perturb

  du = get_trial_fe_basis(UP);
  dv = get_fe_basis(VQ);
  dom_cont = a(du,dv)

  return lazy_schur_complement(dom_cont,diag_perturb,Ω)
end

function schur_complement_dom_contribs(
    solver::BlockSchurPreconditioner{<:AdditiveSchurMethod},
    K::GridapDistributed.BlockPMatrix
  )
  a = solver.a
  UP, VQ = solver.spaces
  Ω = solver.trian
  diag_perturb = solver.diag_perturb

  du = get_trial_fe_basis(UP);
  dv = get_fe_basis(VQ);
  dom_cont = a(du,dv)

  cache = map(local_views(dom_cont),local_views(Ω)) do dom_cont,Ω
    lazy_schur_complement(dom_cont,diag_perturb,Ω)
  end;

  return GridapDistributed.DistributedDomainContribution(cache)
end

function lazy_schur_complement(dom_cont,diag_perturb,trian)
  Ke = get_array(dom_cont);
  schur_map = SchurMap(diag_perturb)
  dom_S = DomainContribution()
  dom_S.dict[trian] = lazy_map(schur_map,Ke);
  return dom_S
end

function assemble_schur_complement(solver,K)
  UP, VQ = solver.spaces
  P, Q = UP[2], VQ[2]
  schur_assem = solver.schur_assembler
  S_dom_contrib = schur_complement_dom_contribs(solver,K)
  matdata = collect_cell_matrix(P,Q,S_dom_contrib)
  return assemble_matrix(schur_assem,matdata);
end

function assemble_schur_complement!(S,solver,K)
  UP, VQ = solver.spaces
  P, Q = UP[2], VQ[2]
  schur_assem = solver.schur_assembler
  S_dom_contrib = schur_complement_dom_contribs(solver,K)
  matdata = collect_cell_matrix(P,Q,S_dom_contrib)
  return assemble_matrix!(S,schur_assem,matdata);
end

function assemble_schur_complement(::BlockSchurPreconditioner{<:ExactSchurMethod},K::BlockMatrix)
  A = blocks(K)[1,1];
  mBᵀ = blocks(K)[2,1];
  B = blocks(K)[1,2];
  C = blocks(K)[2,2];

  err = "The exact Schur complement is dense and thus
  should not be used for large problems. Instead use an
  AdditiveSchurMethod, for example."
  @assert any(size(A) .< 1000) err

  S = C - mBᵀ*(Matrix(A)\B)

  return S
end

function assemble_schur_complement!(S,solver::BlockSchurPreconditioner{<:ExactSchurMethod},K::BlockMatrix)
  _S = assemble_schur_complement(solver,K)
  copyto!(S,_S)
end

function Gridap.Algebra.numerical_setup(ss::BlockSchurPreconditionerSS,K::AbstractBlockMatrix)
  solver = ss.solver

  # Build Schur complement
  S = assemble_schur_complement(solver,K)

  # Setup block solver
  block_solver = solver.solver
  A = blocks(K)[1,1];
  B = blocks(K)[1,2];
  mBᵀ = blocks(K)[2,1];
  @check typeof(blocks(K)[2,1]) == typeof(S)

  P_mat = mortar((A,B),(mBᵀ,S));
  P_ns = numerical_setup(symbolic_setup(block_solver,P_mat),P_mat)
  return BlockSchurPreconditionerNS(solver,K,P_ns,S)
end

function Gridap.Algebra.numerical_setup!(ns::BlockSchurPreconditionerNS,K::AbstractBlockMatrix)
  solver = ns.solver
  S = ns.S
  P_ns = ns.P_ns

  # Build additive Schur complement
  assemble_schur_complement!(S,solver,K)

  # Setup block solver
  A = blocks(K)[1,1];
  B = blocks(K)[1,2];
  mBᵀ = blocks(K)[2,1];
  P_mat = mortar((A,B),(mBᵀ,S));
  numerical_setup!(P_ns,P_mat)
  ns.K = K
  return ns
end

function Gridap.Algebra.numerical_setup!(ns::BlockSchurPreconditionerNS,K::AbstractBlockMatrix,X::AbstractBlockVector)
  solver = ns.solver
  S = ns.S
  P_ns = ns.P_ns

  # Build additive Schur complement
  assemble_schur_complement!(S,solver,K)

  # Setup block solver
  A = blocks(K)[1,1];
  B = blocks(K)[1,2];
  mBᵀ = blocks(K)[2,1];
  P_mat = mortar((A,B),(mBᵀ,S));
  numerical_setup!(P_ns,P_mat,X)
  ns.K = K
  return ns
end

function Gridap.Algebra.solve!(X::AbstractBlockVector,ns::BlockSchurPreconditionerNS,B::AbstractBlockVector)
  solve!(X,ns.P_ns,B)
end

## Bug fix in GridapSolvers (tmp)
function Gridap.Algebra.numerical_setup!(ns::GridapSolvers.BlockSolvers.BlockTriangularSolverNS,mat::AbstractBlockMatrix)
  solver       = ns.solver
  mat_blocks   = blocks(mat)
  block_caches = map(GridapSolvers.BlockSolvers.update_block_cache!,ns.block_caches,solver.blocks,mat_blocks)
  map(diag(solver.blocks),ns.block_ns,diag(block_caches)) do bi, nsi, ci
    numerical_setup!(nsi,ci)
  end
  return ns
end