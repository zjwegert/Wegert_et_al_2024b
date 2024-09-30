## Non-dimensionalisation - see Functionals.jl for discussion
function non_dimensionalise(C,e,κ;ld=10^-5,τ_φ=10^4,c1=10^10,e1=1,κ1=10^-10)
  C = C/c1; e = e/e1; κ = κ/κ1
  A₁₁ = c1*ld^2
  A₁₂ = e1*ld*τ_φ
  A₂₂ = κ1*τ_φ^2

  return (;C,e,κ,A₁₁,A₁₂,A₂₂,c1,e1,κ1)
end

## PZT-5A
function PZT5A(D::Int;nondim = (C,e,κ) -> non_dimensionalise(C,e,κ))
  ε_0 = 8.854e-12;
  if D == 3
    C_Voigt = [12.0400e10  7.52000e10  7.51000e10  0.0        0.0        0.0
          7.52000e10  12.0400e10  7.51000e10  0.0        0.0        0.0
          7.51000e10  7.51000e10  11.0900e10  0.0        0.0        0.0
          0.0         0.0         0.0          2.1000e10  0.0        0.0
          0.0         0.0         0.0          0.0        2.1000e10  0.0
          0.0         0.0         0.0          0.0        0.0       2.30e10]
    e_Voigt = [0.0       0.0       0.0        0.0       12.30000  0.0
                0.0       0.0       0.0        12.30000  0.0       0.0
                -5.40000  -5.40000   15.80000   0.0       0.0       0.0]
    K_Voigt = [540*ε_0     0        0
                  0      540*ε_0    0
                  0        0    830*ε_0]
  elseif D == 2
    C_Voigt = [12.0400e10  7.51000e10     0.0
                7.51000e10  11.0900e10     0.0
                  0.0         0.0      2.1000e10]
    e_Voigt = [0.0       0.0       12.30000
                -5.40000   15.80000   0.0]
    K_Voigt = [540*ε_0        0
                  0    830*ε_0]
  else
    @notimplemented
  end
  C = voigt2tensor4(C_Voigt)
  e = voigt2tensor3(e_Voigt)
  κ = voigt2tensor2(K_Voigt)
  return nondim(C,e,κ)
end

"""
  Given a material constant given in Voigt notation,
  return a SymFourthOrderTensorValue using ordering from Gridap
"""
function voigt2tensor4(A::Array{M,2}) where M
  if isequal(size(A),(3,3))
    return SymFourthOrderTensorValue(A[1,1], A[3,1], A[2,1],
                                     A[1,3], A[3,3], A[2,3],
                                     A[1,2], A[3,2], A[2,2])
  elseif isequal(size(A),(6,6))
    return SymFourthOrderTensorValue(A[1,1], A[6,1], A[5,1], A[2,1], A[4,1], A[3,1],
                                     A[1,6], A[6,6], A[5,6], A[2,6], A[4,6], A[3,6],
                                     A[1,5], A[6,5], A[5,5], A[2,5], A[4,5], A[3,5],
                                     A[1,2], A[6,2], A[5,2], A[2,2], A[4,2], A[3,2],
                                     A[1,4], A[6,4], A[5,4], A[2,4], A[4,4], A[3,4],
                                     A[1,3], A[6,3], A[5,3], A[2,3], A[4,3], A[3,3])
  else
      @notimplemented
  end
end

"""
  Given a material constant given in Voigt notation,
  return a ThirdOrderTensorValue using ordering from Gridap
"""
function voigt2tensor3(A::Array{M,2}) where M
  if isequal(size(A),(2,3))
    return ThirdOrderTensorValue(A[1,1], A[2,1], A[1,3], A[2,3], A[1,3], A[2,3], A[1,2], A[2,2])
  elseif isequal(size(A),(3,6))
    return ThirdOrderTensorValue(
      A[1,1], A[2,1], A[3,1], A[1,6], A[2,6], A[3,6], A[1,5], A[2,5], A[3,5],
      A[1,6], A[2,6], A[3,6], A[1,2], A[2,2], A[3,2], A[1,4], A[2,4], A[3,4],
      A[1,5], A[2,5], A[3,5], A[1,4], A[2,4], A[3,4], A[1,3], A[2,3], A[3,3])
  else
    @notimplemented
  end
end

"""
  Given a material constant given in Voigt notation,
  return a SymTensorValue using ordering from Gridap
"""
function voigt2tensor2(A::Array{M,2}) where M
  return TensorValue(A)
end

## Low-level PETSc CG+AMG solver for Hilbertian extension-regularisation or
#   conductivity part of piezo problem
function cgamg_ksp_setup(;rtol=10^-8,maxits=100)

  function ksp_setup(ksp)
    pc = Ref{GridapPETSc.PETSC.PC}()

    rtol = PetscScalar(rtol)
    atol = GridapPETSc.PETSC.PETSC_DEFAULT
    dtol = GridapPETSc.PETSC.PETSC_DEFAULT
    maxits = PetscInt(maxits)

    @check_error_code GridapPETSc.PETSC.KSPSetType(ksp[],GridapPETSc.PETSC.KSPCG)
    @check_error_code GridapPETSc.PETSC.KSPSetTolerances(ksp[], rtol, atol, dtol, maxits)
    @check_error_code GridapPETSc.PETSC.KSPGetPC(ksp[],pc)
    @check_error_code GridapPETSc.PETSC.PCSetType(pc[],GridapPETSc.PETSC.PCGAMG)
    @check_error_code GridapPETSc.PETSC.KSPView(ksp[],C_NULL)
  end

  return ksp_setup
end

using Gridap.FESpaces: interpolate!
using GridapTopOpt: get_deriv_space, get_aux_space, get_state, update_descent_direction!, finished, print_msg
import Base.iterate

# The below iterate method for the Hilbertian projection optimiser is faster when NOT using AD.
#   The version in GridapTopOpt is faster for AD problems as the adjoint isn't required at every
#   LS iteration.
function Base.iterate(m::GridapTopOpt.HilbertianProjection,state)
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