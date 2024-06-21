"""
This file deals with- 
1. Physical Floquet state from the extended Floquet states
2. Expectation value of an operator as a function of time.  
3.SVD, Entanglement entropy between auxialiary and physical space
"""

using LinearAlgebra 
using Kronecker
"""
Splitting a vector in extended state into an array of vectors.
"""
function Split_Ext_Vector(x_Vec,L,na)
  NFZ = 2^(na-1) - 1
  n = [2^L for j in 1:2^na] #Sorted in 
  ϕ_res = Vector{Vector{eltype(x_Vec)}}()
  start = firstindex(x_Vec)
  for len in n
    push!(ϕ_res, x_Vec[start:(start + len - 1)])
    start += len
  end
  return ϕ_res
end
#Split_Ext_Vector(x_Vec,L,na)
# Different Operators
# function Single_Spin_Sum_Operator(L,site)
# pdictOp0= Dict("X" => s1* , "Y" => , "Z" =>  )
# return pdictOp0
# end
function One_Spin_Sum_Operator(L)
  pdictOp1 = Dict("SumX" => magnetic_field_matrix(L,1.0,0.0,0.0),"SumY" =>  magnetic_field_matrix(L,0.0,1.0,0.0), "SumZ" =>  magnetic_field_matrix(L,0.0,0.0,1.0))
  return pdictOp1
end
function Two_Spin_Sum_Operator(L,PBCs)
  pdictOp2 = Dict("SumXX" => xyz_matrix(L,1.0,0.0,0.0,PBCs),"SumYY" => xyz_matrix(L,0.0,1.0,0.0,PBCs), "SumZZ" => xyz_matrix(L,0.0,0.0,1.0,PBCs))
  return pdictOp2
end

## Operator Expectation Vlaue:: Extnd_Floq_Vec
function Op_Vec_Expectation(Op,Vect,na,Omg,t)
  NFZ = 2^(na-1) - 1
  ϕ_lowest = sum(Vect[n+NFZ+1]*exp(im*n*Omg*t) for n in -NFZ:NFZ)  
  Op_expt = Adjoint(ϕ_lowest)*Op*ϕ_lowest
  return Op_expt  #ϕ_lowest
end
####
#function Op_ExtdVec_Expt(Op,Extd_Floq_Vec,L,na,t)
#   NFZ = 2^(na-1) - 1
#    Split_Ext_Vector(Extd_Floq_Vec,L,na)
#   ϕ_lowest = sum(Extd_Floq_Vec[n+NFZ+1]*exp(im*n*XΩ*t) for n in -NFZ:NFZ)  
#   Op_expt = Adjoint(ϕ_lowest)*Op*ϕ_lowest
#   return Op_expt  #ϕ_lowest
#end
####
function Op_Vec_Expectation_Low(Op,Vect,na,Omg,t)
  NFZ = 2^(na-1) - 1
  ϕ_lowest = sum(Vect[n]*exp(im*(n-2*NFZ-1)*Omg*t) for n in 1:(2*NFZ+1))  
  Op_expt = Adjoint(ϕ_lowest)*Op*ϕ_lowest
  return Op_expt  #ϕ_lowest
end

# Function for operator expectation values and entanglement
# In this cell we define a single and two spin operators. We also define a function that split  a large extended vector into array of small vectors.
# One_Spin_Sum_Operator(L)["SumX"]
# Note: There might be some asymmetric consideration  in -na to na :: Look at the structure of phi_n 

"""
SVD_Entaglement_Entropy(x_Vec,L,na)
SVD and Entanglement entropy between auxialiary and physical space. L= phyiscal sites, na=auxiliary sites.
"""
function SVD_Entaglement_Entropy(x_Vec,L,na)
  NFZ = 2^(na-1) - 1
  ϕntilde = Split_Ext_Vector(x_Vec, L, NFZ)
  CMat = zeros(ComplexF64, 2*NFZ+1, 2^L)
  for j in 1:(2*NFZ+1)
    CMat[j, :] = ϕntilde[j]
  end
  SVDϕ = svd(CMat)
  Svals0 = SVDϕ.S
  Svals = Svals0[Svals0 .> 0.0]
  Ent_bi = -sum(Svals[j]^2 * log(Svals[j]^2) for j in 1:length(Svals))
  return Svals0, Ent_bi
end
# #######
# #######
# input_folder="/Users/abhishek22/Documents/Floq_Exact_Diag/Floq_ED_v0/Examples/Output_Data_Plots/Data"
# Output_folder="/Users/abhishek22/Documents/Floq_Exact_Diag/Floq_ED_v0/Examples/Output_Data_Plots/Plots"
# ####
# input_file1="ED_UF_L_Nt_mesh_Omg_JxJyJz_dJxdJydJz_BxByBz_dBxdBydBz(3, 5000, 6.0, [3.7, 2.7, 1.4], [1.0, 1.4, 1.0], 2.3, 3.0, 4.0, 0.0, 0.0, 0.0, 1)_results.dat"
# input_file=joinpath(input_folder, input_file1)
# result1 = JSON.parsefile(input_file)
# #########
# input_file3="ED_HF_ADAPT_L_na_Omg_JxJyJz_dJxdJydJz_BxByBz_dBxdBydBz(3, 4, 6.0, [3.7, 2.7, 1.4], [1.0, 1.4, 1.0], 2.3, 3.0, 4.0, 0.0, 0.0, 0.0)_results.dat"
# input_file=joinpath(input_folder, input_file3)
# result3 = JSON.parsefile(input_file)
# #########
# na=4
# Ω0=6.0
# L=3
# Final_states=ComplexF64
# chosen_ops = results["paulis"][2:end]
# op_pars = results["opt_pars"][2:end]
# energies = results["energy"]
# max_grads = results["max_grad"]
# num_iters = length(results["opt_numevals"])-1
# Final_states = result3["Eigvecs_HFSqr"];
# display(Final_states[1])
#xArr=[(Final_states[j]["re"]+im*Final_states[j]["im"]) for j in 1:2^(L)]
##

#ϕntilde= Split_Ext_Vector(xArr, sites, na)
#ϕntilde= Split_Array(xArr, n)
#display(ϕntilde[1])
####
#Expectation value of an operator as a function of time.
####
# ϕntilde
# function Op_Expectation(Op, na, t)
#     ϕ_lowest = sum(ϕntilde[j+2^(na-1)]*exp(-im*j*Ω0*t) for j in (-2^(na-1)+1):(2^(na-1)-1))  
#     Op_expt = transpose(conj(ϕ_lowest))*Op*ϕ_lowest
#     return Op_expt  #ϕ_lowest
# end

# OperatorX3=One_Spin_Sum_Operator(3)["SumX"]
# OperatorY3=One_Spin_Sum_Operator(3)["SumY"]
# OperatorZ3=One_Spin_Sum_Operator(3)["SumZ"]
# OperatorXX3=Two_Spin_Sum_Operator(3,false)["SumXX"]
# OperatorYY3=Two_Spin_Sum_Operator(3,false)["SumYY"]
# OperatorZZ3=Two_Spin_Sum_Operator(3,false)["SumZZ"]

#Expectation value of an operator as a function of time.
# function Op_Expectation(Op, NFZ, t)
# ϕ_lowest = sum(ϕntilde[n+NFZ+1] * exp(im * n * XΩ * t) for n in -NFZ:NFZ)  
# Op_expt = Adjoint(ϕ_lowest) * Op * ϕ_lowest
# return Op_expt  #ϕ_lowest
# end
#There might be some asymmetric consideration  in -na to na :: Look at the structure of phi_n 
# Ω0=XΩ
# OperatorZ4=(kron(I2,I2,I2,s3)+kron(I2,I2,s3,I2)+kron(I2,s3,I2,I2)+kron(s3,I2,I2,I2))
# OperatorX4=(kron(I2,I2,I2,s1)+kron(I2,I2,s1,I2)+kron(I2,s1,I2,I2)+kron(s1,I2,I2,I2))
# OperatorY4=(kron(I2,I2,I2,s2)+kron(I2,I2,s2,I2)+kron(I2,s2,I2,I2)+kron(s2,I2,I2,I2))

# OperatorZZ4=(kron(I2,I2,s3,s3)+kron(I2,s3,I2,s3)+kron(s3,I2,I2,s3))
# OperatorXX4=(kron(I2,I2,s1,s1)+kron(I2,s1,I2,s1)+kron(s1,I2,I2,s1))
# OperatorYY4=(kron(I2,I2,s2,s2)+kron(I2,s2,I2,s2)+kron(s2,I2,I2,s2))

# to_plot3=[Op_Expectation(OperatorZ4, na, (2π/Ω0)*t) for t in 0.0:0.01:2.0]
# to_plot2=[Op_Expectation(OperatorY4, na, (2π/Ω0)*t) for t in 0.0:0.01:2.0]
# to_plot1=[Op_Expectation(OperatorX4, na, (2π/Ω0)*t) for t in 0.0:0.01:2.0]
# ##
# Mesht = [t for t in 0:0.01:2.0]
# scalefontsizes(1.1)
# labels = ["Expectation_Value"]
# #plot(Mesht, to_plot1, xlabel="Time", ylabel="Expectation_Value", xticks = Mesht[1]:0.1:Mesht[end], label = labels)
# #savefig("Expectation_Value_vs_(t/T)"*string(sites)*".pdf")

# to_plot3zz=[Op_Expectation(OperatorZZ4, na, (2π/Ω0)*t) for t in 0.0:0.01:2.0];
# to_plot2yy=[Op_Expectation(OperatorYY4, na, (2π/Ω0)*t) for t in 0.0:0.01:2.0];
# to_plot1xx=[Op_Expectation(OperatorZZ4, na, (2π/Ω0)*t) for t in 0.0:0.01:2.0];

# OperatorZ4=kron(sx,sx)
# Op_Expectation(Op,NFZ,t)

