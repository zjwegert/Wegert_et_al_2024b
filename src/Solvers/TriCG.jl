using Gridap, Gridap.Algebra, Gridap.Helpers
using GridapSolvers
using GridapSolvers: ConvergenceLog, SolverTolerances, init!, update!, finalize!
using GridapDistributed
using LinearAlgebra
using PartitionedArrays
using AbstractTrees
using BlockArrays

macro kswap(x, y)
  quote
    local tmp = $(esc(x))
    $(esc(x)) = $(esc(y))
    $(esc(y)) = tmp
  end
end

"""
    struct TriCGSolver <: LinearSolver

Solver representing a TriCG block solver based on Montoison and Orban. See
the manuscript at https://doi.org/10.1137/20M1363030.

This solves the saddle-point problem

```
|  M    A | * | x | = | b |
|  Aᵀ   P |   | y |   | c |
```

where b and c must both be nonzero, M is positive definite, and
P is negative definite. Note that parts of this implementation are based on the
one available in Krylov.jl. Below we take N = -P to match the
description of TriCG published by Montoison and Orban (2021).

# Properties:
- `M_ls :: Gridap.Algebra.LinearSolver`: the linear solver for solving v↦M\\v.
- `N_ls :: Gridap.Algebra.LinearSolver`: the linear solver for solving u↦P\\u.
- `log  :: ConvergenceLog{Float64}`: solver log
"""
struct TriCGSolver <: Gridap.Algebra.LinearSolver
  M_ls :: Gridap.Algebra.LinearSolver
  N_ls :: Gridap.Algebra.LinearSolver
  log  :: ConvergenceLog{Float64}
end

AbstractTrees.children(s::TriCGSolver) = [s.M_ls,s.N_ls]

function TriCGSolver(M_ls,N_ls;maxiter=1000,atol=1e-12,rtol=1.e-6,verbose=false,name="TriCG")
  tols = SolverTolerances{Float64}(maxiter=maxiter,atol=atol,rtol=rtol)
  log  = ConvergenceLog(name,tols,verbose=verbose)
  return TriCGSolver(M_ls,N_ls,log)
end

struct TriCGSymbolicSetup <: Gridap.Algebra.SymbolicSetup
  solver
end

function Gridap.Algebra.symbolic_setup(solver::TriCGSolver, K::AbstractBlockMatrix)
  @check blocksize(K) == (2,2) "A 2x2 block structure is expected"
  return TriCGSymbolicSetup(solver)
end

mutable struct TriCGNumericalSetup <: Gridap.Algebra.NumericalSetup
  solver
  K
  M_ns
  N_ns
  caches
end

function get_solver_caches(solver::TriCGSolver,K::AbstractBlockMatrix)
  M = blocks(K)[1,1]; N = blocks(K)[2,2]

  Mvₖ₋₁  = allocate_in_domain(M)
  Mvₖ    = allocate_in_domain(M)
  q      = allocate_in_domain(M)
  gx₂ₖ₋₁ = allocate_in_domain(M)
  gx₂ₖ   = allocate_in_domain(M)
  vₖ     = allocate_in_domain(M)

  Nuₖ₋₁  = allocate_in_domain(N)
  Nuₖ    = allocate_in_domain(N)
  p      = allocate_in_domain(N)
  gy₂ₖ₋₁ = allocate_in_domain(N)
  gy₂ₖ   = allocate_in_domain(N)
  uₖ     = allocate_in_domain(N)

  return (Mvₖ₋₁,Mvₖ,q,gx₂ₖ₋₁,gx₂ₖ,vₖ),(Nuₖ₋₁,Nuₖ,p,gy₂ₖ₋₁,gy₂ₖ,uₖ)
end

function Gridap.Algebra.numerical_setup(ss::TriCGSymbolicSetup,K::AbstractBlockMatrix)
  solver = ss.solver
  M = blocks(K)[1,1]; N = blocks(K)[2,2];
  M_ns = numerical_setup(symbolic_setup(solver.M_ls,M),M)
  N_ns = numerical_setup(symbolic_setup(solver.N_ls,N),N)
  caches = get_solver_caches(solver,K)
  return TriCGNumericalSetup(solver,K,M_ns,N_ns,caches)
end

function Gridap.Algebra.numerical_setup!(ns::TriCGNumericalSetup,K::AbstractBlockMatrix)
  M = blocks(K)[1,1]; N = blocks(K)[2,2];
  numerical_setup!(ns.M_ns,M)
  numerical_setup!(ns.N_ns,N)
  ns.K = K
  return ns
end

function Gridap.Algebra.numerical_setup!(ns::TriCGNumericalSetup,K::AbstractBlockMatrix,X::AbstractBlockVector)
  M = blocks(K)[1,1]; N = blocks(K)[2,2];
  x = blocks(X)[1]; y = blocks(X)[2]
  numerical_setup!(ns.M_ns,M,x)
  numerical_setup!(ns.N_ns,N,y)
  ns.K = K
  return ns
end

# Parts of the below are based on Krylov.jl TriCG implementation
function Gridap.Algebra.solve!(X::AbstractBlockVector,ns::TriCGNumericalSetup,B::AbstractBlockVector)
  @check blocklength(X) == blocklength(B)
  T = eltype(X)
  ## Unpack
  solver, K, M_ns, N_ns, caches = ns.solver, ns.K, ns.M_ns, ns.N_ns, ns.caches
  log = solver.log
  M_cache, N_cache = caches
  Mvₖ₋₁, Mvₖ, q, gx₂ₖ₋₁, gx₂ₖ, vₖ = M_cache
  Nuₖ₋₁, Nuₖ, p, gy₂ₖ₋₁, gy₂ₖ, uₖ = N_cache
  vₖ₊₁ = vₖ; uₖ₊₁ = uₖ

  A = blocks(K)[1,2]; Aᵀ = blocks(K)[2,1]
  xₖ = blocks(X)[1]; yₖ = blocks(X)[2]
  b = blocks(B)[1]; c = blocks(B)[2]

  ## Initialise
  iter = 0
  fill!(xₖ,zero(T)); fill!(yₖ,zero(T));
  fill!(Mvₖ₋₁,zero(T)); fill!(Mvₖ,zero(T)); fill!(q,zero(T))
  fill!(gx₂ₖ₋₁,zero(T)); fill!(gx₂ₖ,zero(T)); fill!(vₖ,zero(T))
  fill!(Nuₖ₋₁,zero(T)); fill!(Nuₖ,zero(T)); fill!(p,zero(T))
  fill!(gy₂ₖ₋₁,zero(T)); fill!(gy₂ₖ,zero(T)); fill!(uₖ,zero(T))
  d₂ₖ₋₃ = d₂ₖ₋₂ = zero(T)
  π₂ₖ₋₃ = π₂ₖ₋₂ = zero(T)
  δₖ₋₁ = zero(T)

  # Tolerance for breakdown detection <-
  btol = eps()^(3/4)

  ## Initialise triorthogonalisation
  # Compute β₁Mv₁=b
  copyto!(Mvₖ, b)
  solve!(vₖ, M_ns, Mvₖ)
  βₖ = sqrt(dot(vₖ, Mvₖ)) # β₁ = ‖v₁‖_M
  if βₖ ≠ 0
    Mvₖ ./= βₖ
    vₖ  ./= βₖ
  else
    error("b must be nonzero")
  end
  # Compute γ₁Nu₁=c
  copyto!(Nuₖ, c)
  solve!(uₖ, N_ns, Nuₖ)
  uₖ .*= -1 # N = -P from (2,2) matrix block
  γₖ = sqrt(dot(uₖ, Nuₖ)) # γ₁ = ‖u₁‖_N
  if γₖ ≠ 0
    Nuₖ ./= γₖ
    uₖ  ./= γₖ
  else
    error("c must be nonzero")
  end
  # Compute ‖r₀‖² = (γ₁)² + (β₁)²
  res = sqrt(γₖ^2 + βₖ^2)
  done = init!(log, res)

  while !done
    iter = iter + 1

    ## Continue triorthogonalisation
    # Compute q = Auₖ - γₖMvₖ₋₁ and p = Aᵀvₖ - βₖNuₖ₋₁
    mul!(q, A , uₖ); mul!(p, Aᵀ, vₖ)
    if iter ≥ 2
      q .-= γₖ * Mvₖ₋₁; p .-= βₖ * Nuₖ₋₁
    end
    αₖ = dot(vₖ, q)
    # Compute q = q - αₖMvₖ and p = p - ̄αₖNuₖ
    q .-= αₖ * Mvₖ; p .-= conj(αₖ) * Nuₖ
    # Update Mvₖ₋₁ and Nuₖ₋₁
    Mvₖ₋₁ .= Mvₖ
    Nuₖ₋₁ .= Nuₖ

    ## Compute d,δ,σ,η,λ
    if iter == 1
      d₂ₖ₋₁ = 1
      δₖ    = conj(αₖ) / d₂ₖ₋₁
      d₂ₖ   = -1 - abs2(δₖ) * d₂ₖ₋₁
    else
      σₖ    = βₖ / d₂ₖ₋₂
      ηₖ    = γₖ / d₂ₖ₋₃
      λₖ    = -(ηₖ * conj(δₖ₋₁) * d₂ₖ₋₃) / d₂ₖ₋₂
      d₂ₖ₋₁ = 1 - abs2(σₖ) * d₂ₖ₋₂
      δₖ    = (conj(αₖ) - λₖ * conj(σₖ) * d₂ₖ₋₂) / d₂ₖ₋₁
      d₂ₖ   = -1 - abs2(ηₖ) * d₂ₖ₋₃ - abs2(λₖ) * d₂ₖ₋₂ - abs2(δₖ) * d₂ₖ₋₁
    end
    # Solve LₖDₖpₖ = (β₁e₁ + γ₁e₂)
    if iter == 1
      π₂ₖ₋₁ = βₖ / d₂ₖ₋₁
      π₂ₖ   = (γₖ - δₖ * βₖ) / d₂ₖ
    else
      π₂ₖ₋₁ = -(σₖ * d₂ₖ₋₂ * π₂ₖ₋₂) / d₂ₖ₋₁
      π₂ₖ   = -(δₖ * d₂ₖ₋₁ * π₂ₖ₋₁ + λₖ * d₂ₖ₋₂ * π₂ₖ₋₂ + ηₖ * d₂ₖ₋₃ * π₂ₖ₋₃) / d₂ₖ
    end

    ## Update gₖ
    if iter == 1
      @. gx₂ₖ₋₁ = vₖ
      @. gx₂ₖ   = - conj(δₖ) * gx₂ₖ₋₁
      @. gy₂ₖ   = uₖ
    else
      @. gx₂ₖ₋₁ = conj(ηₖ) * gx₂ₖ₋₁ + conj(λₖ) * gx₂ₖ
      @. gy₂ₖ₋₁ = conj(ηₖ) * gy₂ₖ₋₁ + conj(λₖ) * gy₂ₖ
      @. gx₂ₖ = vₖ - conj(σₖ) * gx₂ₖ
      @. gy₂ₖ =    - conj(σₖ) * gy₂ₖ
      @. gx₂ₖ₋₁ =    - gx₂ₖ₋₁ - conj(δₖ) * gx₂ₖ
      @. gy₂ₖ₋₁ = uₖ - gy₂ₖ₋₁ - conj(δₖ) * gy₂ₖ
      # g₂ₖ₋₃ == g₂ₖ and g₂ₖ₋₂ == g₂ₖ₋₁
      @kswap(gx₂ₖ₋₁, gx₂ₖ)
      @kswap(gy₂ₖ₋₁, gy₂ₖ)
    end

    ## Compute xₖ and yₖ
    xₖ .+= π₂ₖ₋₁ * gx₂ₖ₋₁ + π₂ₖ * gx₂ₖ
    yₖ .+= π₂ₖ₋₁ * gy₂ₖ₋₁ + π₂ₖ * gy₂ₖ

    ## Compute vₖ₊₁, uₖ₊₁, βₖ₊₁, γₖ₊₁
    solve!(vₖ₊₁, M_ns, q) # βₖ₊₁vₖ₊₁ = MAuₖ  - γₖvₖ₋₁ - αₖvₖ
    solve!(uₖ₊₁, N_ns, p) # γₖ₊₁uₖ₊₁ = NAᴴvₖ - βₖuₖ₋₁ - ᾱₖuₖ
    uₖ₊₁ .*= -1 # N = -P from (2,2) matrix block
    βₖ₊₁ = sqrt(dot(vₖ₊₁, q))  # βₖ₊₁ = ‖vₖ₊₁‖_E
    γₖ₊₁ = sqrt(dot(uₖ₊₁, p))  # γₖ₊₁ = ‖uₖ₊₁‖_F
    # βₖ₊₁ ≠ 0
    if βₖ₊₁ > btol
      q    ./= βₖ₊₁
      vₖ₊₁ ./= βₖ₊₁
    end
    # γₖ₊₁ ≠ 0
    if γₖ₊₁ > btol
      p    ./= γₖ₊₁
      uₖ₊₁ ./= γₖ₊₁
    end

    ## Update Mvₖ and Nuₖ
    Mvₖ .= q
    Nuₖ .= p

    # Compute ‖rₖ‖² = |γₖ₊₁ζ₂ₖ₋₁|² + |βₖ₊₁ζ₂ₖ|²
    ζ₂ₖ₋₁ = π₂ₖ₋₁ - conj(δₖ) * π₂ₖ
    ζ₂ₖ   = π₂ₖ
    res = sqrt(abs2(γₖ₊₁ * ζ₂ₖ₋₁) + abs2(βₖ₊₁ * ζ₂ₖ))
    done = update!(log,res)

    # Update βₖ, γₖ, π₂ₖ₋₃, π₂ₖ₋₂, d₂ₖ₋₃, d₂ₖ₋₂, δₖ₋₁, vₖ, uₖ.
    βₖ    = βₖ₊₁
    γₖ    = γₖ₊₁
    π₂ₖ₋₃ = π₂ₖ₋₁
    π₂ₖ₋₂ = π₂ₖ
    d₂ₖ₋₃ = d₂ₖ₋₁
    d₂ₖ₋₂ = d₂ₖ
    δₖ₋₁  = δₖ
  end

  finalize!(log,res)
  return X
end