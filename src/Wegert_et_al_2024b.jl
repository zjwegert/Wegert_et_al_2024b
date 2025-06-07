module Wegert_et_al_2024b

using Gridap
using GridapDistributed
using GridapPETSc
using GridapSolvers
using PartitionedArrays
using GridapTopOpt
using SparseMatricesCSR

using Gridap.TensorValues
using Gridap.Helpers
using Gridap.MultiField
using Gridap.CellData
using Gridap.FESpaces
using Gridap.Arrays

using GridapSolvers: ConvergenceLog, SolverTolerances, init!, update!, finalize!
using LinearAlgebra
using AbstractTrees
using BlockArrays

using Gridap.FESpaces: interpolate!
using GridapTopOpt: get_deriv_space, get_aux_space, get_state, update_descent_direction!, finished, print_msg
import Base.iterate
import Gridap.Algebra.numerical_setup!

__init__() = begin
  function numerical_setup!(ns::GridapSolvers.BlockSolvers.BlockTriangularSolverNS,mat::AbstractBlockMatrix)
    solver       = ns.solver
    mat_blocks   = blocks(mat)
    block_caches = map(GridapSolvers.BlockSolvers.update_block_cache!,ns.block_caches,solver.blocks,mat_blocks)
    map(diag(solver.blocks),ns.block_ns,diag(block_caches)) do bi, nsi, ci
      numerical_setup!(nsi,ci)
    end
    return ns
  end

  function iterate(m::GridapTopOpt.HilbertianProjection,state)
    it, J, C, θ, dJ, dC, uh, φh, vel, φ_tmp, γ, os_it = state
    history, params = m.history, m.params

    ## Periodicially call GC
    iszero(it % 50) && GC.gc();

    ## Check stopping criteria
    if finished(m)
      return nothing
    end

    ## Oscillation detection
    if (γ > 0.001) && m.has_oscillations(m,os_it)
      os_it = it + 1
      γ    *= params.os_γ_mult
      print_msg(m.history,"   Oscillations detected, reducing γ to $(γ)\n",color=:yellow)
    end

    ## Line search
    U_reg = get_deriv_space(m.problem.state_map)
    V_φ   = get_aux_space(m.problem.state_map)
    interpolate!(FEFunction(U_reg,θ),vel,V_φ)

    ls_enabled = params.ls_enabled
    reinit_mod = params.reinit_mod
    ls_max_iters,δ_inc,δ_dec = params.ls_max_iters,params.ls_δ_inc,params.ls_δ_dec
    ξ, ξ_reduce, ξ_reduce_tol = params.ls_ξ, params.ls_ξ_reduce_coef, params.ls_ξ_reduce_abs_tol
    γ_min, γ_max = params.ls_γ_min,params.ls_γ_max

    ls_it = 0; done = false
    φ = get_free_dof_values(φh); copy!(φ_tmp,φ)
    while !done && (ls_it <= ls_max_iters)
      # Advect  & Reinitialise
      evolve!(m.ls_evolver,φ,vel,γ)
      iszero(it % reinit_mod) && reinit!(m.ls_evolver,φ,params.γ_reinit)

      if ~ls_enabled
        J, C, dJ, dC = Gridap.evaluate!(m.problem,φh)
        break
      end

      # Calcuate new objective and constraints
      J_interm, C_interm, dJ_interm, dC_interm = Gridap.evaluate!(m.problem,φh)

      # Reduce line search parameter if constraints close to saturation
      _ξ = all(Ci -> abs(Ci) < ξ_reduce_tol, C_interm) ? ξ*ξ_reduce : ξ

      # Accept/reject
      if (J_interm < J + _ξ*abs(J)) || (γ <= γ_min)
        γ = min(δ_inc*γ, γ_max)
        done = true
        print_msg(history,"  Accepted iteration with γ = $(γ) \n";color=:yellow)
        J, C, dJ, dC = J_interm, C_interm, dJ_interm, dC_interm
      else
        γ = max(δ_dec*γ, γ_min)
        copy!(φ,φ_tmp)
        print_msg(history,"  Reject iteration with γ = $(γ) \n";color=:red)
      end
    end

    ## Calculate objective, constraints, and shape derivatives after line search
    uh = get_state(m.problem)

    ## Hilbertian extension-regularisation
    project!(m.vel_ext,dJ)
    project!(m.vel_ext,dC)
    θ = update_descent_direction!(m.projector,dJ,C,dC,m.vel_ext.K)

    ## Update history and build state
    push!(history,(J,C...,γ))
    state = (;it=it+1, J, C, θ, dJ, dC, uh, φh, vel, φ_tmp, γ, os_it)
    vars  = params.debug ? (it,uh,φh,state) : (it,uh,φh)
    return vars, state
  end
end

include("Helpers.jl")
export PZT5A
export cgamg_ksp_setup

include("Functionals.jl")
export piezo_hom_weak_form
export homogenised_mat_tensor_functionals

include("Solvers/Solvers.jl")

end