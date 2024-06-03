# using Random
# using LinearAlgebra
import JSON
using DelimitedFiles
using Plots

#Basic parameters
# Ω0=10.0
# na=2
# Mesht= (2π/Ω0)*[t for t in -1:0.01:1]

# read in ADAPT results
input_folder="/Users/abhishek22/Documents/Floq_Exact_Diag/Floq_ED_v0/Examples/Output_Data_Plots/Data"
Output_folder="/Users/abhishek22/Documents/Floq_Exact_Diag/Floq_ED_v0/Examples/Output_Data_Plots/Plots"
##
input_file1="ED_UF_L_Nt_mesh_Omg_JxJyJz_dJxdJydJz_BxByBz_dBxdBydBz(3, 5000, 6.0, [3.7, 2.7, 1.4], [1.0, 1.4, 1.0], 2.3, 3.0, 4.0, 0.0, 0.0, 0.0, 1)_results.dat"
input_file=joinpath(input_folder, input_file1)
result1 = JSON.parsefile(input_file)

input_file2="ED_Ext_HFSqr_L_NFZ_Omg_JxJyJz_dJxdJydJz_BxByBz_dBxdBydBz(3, 7, 6.0, [3.7, 2.7, 1.4], [1.0, 1.4, 1.0], 2.3, 3.0, 4.0, 0.0, 0.0, 0.0)_results.dat"
input_file=joinpath(input_folder, input_file2)
result2 = JSON.parsefile(input_file)

input_file3="ED_HF_ADAPT_L_na_Omg_JxJyJz_dJxdJydJz_BxByBz_dBxdBydBz(3, 4, 6.0, [3.7, 2.7, 1.4], [1.0, 1.4, 1.0], 2.3, 3.0, 4.0, 0.0, 0.0, 0.0)_results.dat"
input_file=joinpath(input_folder, input_file3)
result3 = JSON.parsefile(input_file)

XL=3
energies_ED_UF = [result1["Quasienergies"][j]["re"] for j in 1:length(result1["Quasienergies"])]
# energies_ED_HF = [result2["Eigvals_HFSqr"][j]["re"] for j in 1:length(result1["Quasienergies"])]
# energies_ED_HF_ADAPT = [result3["Quasienergies"][j]["re"] for j in 1:length(result1["Quasienergies"])]
# display(energies)
energies_Sqr= sort((energies_ED_UF).^2)

iters = 1:length(energies_Sqr)
to_plot = [result2["Eigvals_HFSqr"][1:2^XL], result3["Eigvals_HFSqr"][1:2^XL], energies_Sqr]
#to_plot= energies
scalefontsizes(1.1)
labels = ["Ext_HFSqr" "ADAPT_Ext_HFSqr" "Direct_UF_Sqr"]
#labels = ["Quasienergy_UF" "Quasienergy(HF)" "Quasienergy(HFSqr_Expt)"]
plot(iters, to_plot, xlabel="Index", ylabel="energy", xticks = iters[1]:iters[end], lw = 3, label = labels,markershape=:circle, markersize=4)
Output_plot = joinpath(Output_folder, "Quasienergy_Square_comparison_0.pdf")
savefig(Output_plot)

# # Function for splitting array into subarray as a vector 
# Ω0=10.0
# Mesht= (2π/Ω0)*[t for t in -1:0.01:1]

# # Function for splitting array into subarray as a vector 
# x=rand(16)
# n=[2^(L-na) for j in 1:2^na];
# function Split_Array(x, n)
#          result = Vector{Vector{eltype(x)}}()
#          start = firstindex(x)
#          for len in n
#            push!(result, x[start:(start + len - 1)])
#            start += len
#          end
#          result
#        end
# ϕntilde= Split_Array(x, n)
# ####
# #Expectation value of an operator as a function of time.
# ####
# ϕntilde
# function Op_Expectation(Op, na, t)
#     ϕ_lowest = sum(ϕntilde[n+na+1]*exp(im*n*Ω0*t) for n in -na:(na-1))  
#     Op_expt = transpose(conj(ϕ_lowest))*Op*ϕ_lowest
#     return Op_expt  #ϕ_lowest
# end
# [Op_Expectation(s1,2, t) for t in -1:0.2:1]

# # create array of chosen operators
# #chosen_op_array = []
# #for pauli_str in chosen_ops
# #    push!(chosen_op_array,pauli_string_to_pauli(pauli_str))
# #end

# # current_state = deepcopy(init_state)
# # pauli_mult!(chosen_op_array[1], init_state, current_state)
# # println(current_state)
# # for i in 1:num_iters
# # end