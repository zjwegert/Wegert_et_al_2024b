using Wegert_et_al_2024b,GridapTopOpt,Gridap,Gridap.TensorValues,JSON3,Comonicon,Printf
# Distributed Dependencies
using GridapDistributed, GridapPETSc, GridapSolvers, PartitionedArrays,
  SparseMatricesCSR, GridapSolvers.BlockSolvers, Gridap.MultiField

using GridapTopOpt_Piezo: bulk_modulus_3d, hydrostatic_coupling_3d
using GridapTopOpt: forward_solve!,benchmark#,benchmark_forward_problem

function benchmark_forward_problem(m::AffineFEStateMap, φh, ranks; nreps = 10)
  function f(m, φh)
    forward_solve!(m,φh)
  end
  function reset!(m,φh)
    x = m.fwd_caches[4]
    fill!(x,zero(eltype(x)))
  end
  return benchmark(f, (m,φh), ranks; nreps, reset!)
end

function main(ranks,params)
  mesh_partition,n,order,lsf_fn,path,file_name,
    solver_method,nreps,elast_solver_rtol,cond_solver_rtol = params

  el_size = (n,n,n)
  xmax,ymax,zmax=(1.0,1.0,1.0)
  dom = (0,xmax,0,ymax,0,zmax)
  max_steps = floor(Int,order*minimum(el_size)/10)
  tol = 1/(5order^2)/minimum(el_size)
  η_coeff = 2
  i_am_main(ranks) && mkpath(path)

  sym = isequal(solver_method,:Schur) ? false : true

  ## FE Setup
  model = CartesianDiscreteModel(ranks,mesh_partition,dom,el_size,isperiodic=(true,true,true))
  el_Δ = get_el_Δ(model)
  f_Γ_D(x) = iszero(x)
  update_labels!(1,model,f_Γ_D,"origin")

  ## Triangulations and measures
  Ω = Triangulation(model)
  dΩ = Measure(Ω,2*order)

  ## Spaces
  reffe = ReferenceFE(lagrangian,VectorValue{3,Float64},order)
  reffe_scalar = ReferenceFE(lagrangian,Float64,order)
  V = TestFESpace(model,reffe;conformity=:H1,dirichlet_tags=["origin"])
  U = TrialFESpace(V,VectorValue(0.0,0.0,0.0))
  Q = TestFESpace(model,reffe_scalar;conformity=:H1,dirichlet_tags=["origin"])
  P = TrialFESpace(Q,0)
  mfs = BlockMultiFieldStyle()
  UP = MultiFieldFESpace([U,P];style=mfs)
  VQ = MultiFieldFESpace([V,Q];style=mfs)

  V_φ = TestFESpace(model,reffe_scalar)
  V_reg = TestFESpace(model,reffe_scalar)
  U_reg = TrialFESpace(V_reg)

  ## Create FE functions
  φh = interpolate(lsf_fn,V_φ)

  ## Weak form
  mat = PZT5A(3)
  interp = SmoothErsatzMaterialInterpolation(η = η_coeff*maximum(el_Δ),ϵ=10^-6)
  a,l,_,_ = piezo_hom_weak_form(mat,interp;sym)

  ## Finite difference solver and level set function
  stencil = HamiltonJacobiEvolution(FirstOrderStencil(3,Float64),model,V_φ,tol,max_steps)
  reinit!(stencil,get_free_dof_values(φh),0.5)

  ## Setup solver and FE operators
  Tm = SparseMatrixCSR{0,PetscScalar,PetscInt}
  Tv = Vector{PetscScalar}
  solver_u = ElasticitySolver(V;rtol=elast_solver_rtol)
  solver_ϕ = PETScLinearSolver(cgamg_ksp_setup(rtol=cond_solver_rtol))

  if isequal(solver_method,:Schur)
    P_solver = GridapSolvers.BlockTriangularSolver([solver_u,solver_ϕ],half=:lower);
    schur_preconditioner = BlockSchurPreconditioner(P_solver,(u,v)->a(u,v,φh,dΩ),UP,VQ,Ω,10^-10;
      schur_assembler = SparseMatrixAssembler(Tm,Tv,P,Q))
    solver = GridapSolvers.LinearSolvers.GMRESSolver(100;Pr=schur_preconditioner,rtol=1.e-8,
      verbose=i_am_main(ranks))
  elseif isequal(solver_method,:TriCG)
    solver = TriCGSolver(solver_u,solver_ϕ;rtol=1.e-8,verbose=i_am_main(ranks))
  elseif isequal(solver_method,:BlockGMRES)
    P = GridapSolvers.BlockDiagonalSolver([solver_u,solver_ϕ])
    solver = GridapSolvers.LinearSolvers.GMRESSolver(100;Pr=P,rtol=1.e-8,verbose=i_am_main(ranks))
  else
    error("Solver ($solver_method) not defined")
  end

  state_map = AffineFEStateMap(a,l[1],UP,VQ,V_φ,U_reg,φh,dΩ;
    assem_U = SparseMatrixAssembler(Tm,Tv,UP,VQ),
    assem_adjoint = SparseMatrixAssembler(Tm,Tv,VQ,UP),
    assem_deriv = SparseMatrixAssembler(Tm,Tv,U_reg,U_reg),
    ls = solver,adjoint_ls = solver
  )
  bfwd = benchmark_forward_problem(state_map, φh, ranks; nreps)

  ns = state_map.fwd_caches[1]
  if isequal(solver_method,:Schur)
    u_ns = ns.Pr_ns.P_ns.block_ns[1]
    ϕ_ns = ns.Pr_ns.P_ns.block_ns[2]
  elseif isequal(solver_method,:TriCG)
    u_ns = ns.M_ns
    ϕ_ns = ns.N_ns
  elseif isequal(solver_method,:BlockGMRES)
    u_ns = ns.Pr_ns.block_ns[1]
    ϕ_ns = ns.Pr_ns.block_ns[2]
  else
    error("Solver ($solver_method) not defined")
  end

  u_niters = Ref{PetscInt}()
  phi_niters = Ref{PetscInt}()
  @check_error_code GridapPETSc.PETSC.KSPGetIterationNumber(u_ns.ksp[],u_niters)
  @check_error_code GridapPETSc.PETSC.KSPGetIterationNumber(ϕ_ns.ksp[],phi_niters)
  outer_its = solver.log.num_iters
  outer_res = solver.log.residuals[solver.log.num_iters+1]/solver.log.residuals[1]

  if i_am_main(ranks)
    open(path*file_name*".txt","w") do f
      bcontent = "bfwd,u_its,phi_its,outer_its,outer_res\n"
      for i = 1:nreps
        bcontent *= "$(bfwd[i]),$(u_niters[]),$(phi_niters[]),$outer_its,$outer_res\n"
      end
      write(f,bcontent)
    end
  end
end

"""
Run the above solver benchmark. Note that the default value of
write_dir requires that the environment variable PROJECT_DIR is defined.

# Arguments

# Options

- `--px <arg>`: partition along x-axis.
- `--py <arg>`: partition along y-axis.
- `--pz <arg>`: partition along z-axis.
- `--n <arg>` : mesh partition along each axis.
- `--order <arg>`: mesh partition in each axis
- `--lsf_func <arg>`: initial level set function. This gets parsed and evaluated so it should be something like
  initial_lsf(2,0.2) or an anonymous function x->-sqrt((x[1]-0.5)^2+(x[2]-0.5)^2)+0.55
- `--write_dir <arg>`: directory where results are written
- `--file_name <arg>`: Output file name
- `--solver_method <arg>`: solver to be used (Schur,TriCG,BlockGMRES)
- `--nreps <arg>`: Number of benchmark repetitions
- `--petsc_opts <arg>`: options for PETScLinearSolver
- `--elast_solver_rtol <arg>`: rtol for elasticity solver
- `--cond_solver_rtol <arg>`: rtol for conductivity solver

# Flags

"""
@main function run(;
    px::Int,
    py::Int,
    pz::Int,
    n::Int=100,
    order::Int=1,
    lsf_func="initial_lsf(2,0.4)",
    solver_method::String="Schur",
    nreps::Int=10,
    write_dir::String=ENV["PROJECT_DIR"]*"/results/solver_benchmarks/",
    petsc_opts::String="-ksp_converged_reason",
    elast_solver_rtol::Float64=10^-4,
    cond_solver_rtol::Float64=10^-4,
    file_name = "3d_ncpus=$(px*py*pz)_n=$(n)_o=$(order)_solver=$(solver_method)_Urtol=$(@sprintf "%.0e" elast_solver_rtol)_Phirtol=$(@sprintf "%.0e" cond_solver_rtol)")
  mesh_partition = (px,py,pz)
  with_mpi() do distribute
    ranks = distribute(LinearIndices((prod(mesh_partition),)))
    solver_method = Symbol(solver_method)
    if i_am_main(ranks)
      println("--------------------------------------------------")
      println("Parsed args for $solver_method solver benchmark:")
      println("       (px,py,pz) => ", mesh_partition)
      println("                n => ", n)
      println("            order => ", order)
      println("         lsf_func => ", lsf_func)
      println("        write_dir => ", write_dir)
      println("    solver_method => ", solver_method)
      println("            nreps => ", nreps)
      println("       petsc_opts => ", petsc_opts)
      println("elast_solver_rtol => ", elast_solver_rtol)
      println(" cond_solver_rtol => ", cond_solver_rtol)
      println("--------------------------------------------------")
    end
    lsf_fn(x) = invokelatest(eval(Meta.parse(lsf_func)),x)
    params = (;mesh_partition,n,order,lsf_fn,path=write_dir,file_name,
      solver_method,nreps,elast_solver_rtol,cond_solver_rtol)
    GridapPETSc.with(args=split(petsc_opts)) do
      main(ranks,params)
    end
  end
end