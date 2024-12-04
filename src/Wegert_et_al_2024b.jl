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

include("Helpers.jl")
export PZT5A
export cgamg_ksp_setup

include("Functionals.jl")
export piezo_hom_weak_form
export homogenised_mat_tensor_functionals

include("Solvers/Solvers.jl")

end