module Floq_ED_v0

include("ED_Floq_Direct.jl")
include("ED_FloqPauliString.jl")
include("Operator_Expt_SVD.jl")

#######
#'To add any additional file I have to include("file.jl") and then I have to export all the functions (that will be used in the main file) in 
#that file here.'
######
export xyz_matrix,
    magnetic_field_matrix,
    Direct_Floq_Sys,
    Eig_Floq_Ham_Square,
    Ext_Floq_Sym_Asym_AuxQubit,
    #ADAPT_Aux_matrices
    PauliString_to_Matrix,
    generate_pauli_string,
    Aux_Diag_Str_Coeff,
    #Str_Coeff_Floq_Diag,
    Str_Coeff_Floq_OffDiag,
    Ham_Aux_Sym_Part,
    All_Aux_Ham,
    #Full Floquet Hamiltonian
    Floq_xyz_model,
    Floq_Magnetic_Field_model,
    ED_Floq_Ham_XYZ_B,
    ## State analysis
    Split_Ext_Vector,
    Op_Vec_Expectation,
    Op_Vec_Expectation_Low,
    One_Spin_Sum_Operator,
    Two_Spin_Sum_Operator,
    SVD_Entaglement_Entropy
end

