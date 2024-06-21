using Random
using LinearAlgebra
import JSON
using DelimitedFiles
using Plots
using Floq_ED_v0
using LaTeXStrings

# read in ADAPT results
input_folder="/Users/abhishek22/Documents/Floq_Exact_Diag/Floq_ED_v0/Examples/Output_Data_Plots/Data"
Output_folder="/Users/abhishek22/Documents/Floq_Exact_Diag/Floq_ED_v0/Examples/Output_Data_Plots/Plots"
##
input_file1="ED_UF_L_Nt_mesh_Omg_JxJyJz_dJxdJydJz_BxByBz_dBxdBydBz(3, 5000, 6.0, [3.7, 2.7, 1.4], [1.0, 1.4, 1.0], 2.3, 3.0, 4.0, 0.0, 0.0, 0.0, 1)_results.dat"
input_file=joinpath(input_folder, input_file1)
result1 = JSON.parsefile(input_file)

#L_na_Omg_JxJyJz_dJxdJydJz_BxByBz_dBxdBydBz(3, 4, 6.0, 3.7, 2.7, 1.4, 1.0, 1.4, 1.0, 2.3, -3.0, 4.0, 0.0, 0.0, 0.0)
input_file2="ED_Ext_HFSqr_L_NFZ_Omg_JxJyJz_dJxdJydJz_BxByBz_dBxdBydBz(3, 7, 6.0, [3.7, 2.7, 1.4], [1.0, 1.4, 1.0], 2.3, 3.0, 4.0, 0.0, 0.0, 0.0)_results.dat"
input_file=joinpath(input_folder, input_file2)
result2 = JSON.parsefile(input_file)

input_file3="ED_HF_ADAPT_L_na_Omg_JxJyJz_dJxdJydJz_BxByBz_dBxdBydBz(3, 4, 6.0, [3.7, 2.7, 1.4], [1.0, 1.4, 1.0], 2.3, 3.0, 4.0, 0.0, 0.0, 0.0)_results.dat"
input_file=joinpath(input_folder, input_file3)
result3 = JSON.parsefile(input_file)

##################
############################################################################################################
XL=3
Xna=4
XL_tot=XL+Xna
Low_Cen_FloqZone=(2^(Xna-1)-1)*2^XL +1
High_Cen_FloqZone=(2^(Xna-1))*2^XL

  E_ED_UF = [result1["Quasienergies"][j]["re"] for j in 1:length(result1["Quasienergies"])]
  E_ED_dir = [result2["Eigvals_HF"][j] for j in Low_Cen_FloqZone:High_Cen_FloqZone]
  E_ED_ADAPT = [result3["Eigvals_HF"][j] for j in Low_Cen_FloqZone:High_Cen_FloqZone]
  #DeltaE_Eqr
  DeltaE_ED_Dir=E_ED_dir - E_ED_UF
  DeltaE_ED_ADAPT=E_ED_ADAPT - E_ED_UF
  println("DeltaE_Dir : ", DeltaE_ED_Dir)
  println("DeltaE_ADAPT : ", DeltaE_ED_ADAPT)
  iters = 1:2^XL
  to_plot = [DeltaE_ED_Dir DeltaE_ED_ADAPT]
  scalefontsizes(1.1)
  labels = ["DeltaE_ED_Dir" "DeltaE_ED_ADAPT"]
  plot(iters, to_plot, xlabel="Index", ylabel="Delta_Quasienergy", xticks = iters[1]:iters[end], lw = 3, label = labels,markershape=:circle, markersize=4)
  Output_plot = joinpath(Output_folder, "Delta_Quasienergy_comparison1.pdf")
  savefig(Output_plot)
#DeltaE_Eqr
  E_Sqr_UF = [result1["Eig_Values_Sqr"][j]["re"] for j in 1:2^XL]
  E_Sqr_ED_dir = [result2["Eigvals_HFSqr"][j] for j in 1:2^XL]
  E_Sqr_ED_ADAPT = [result3["Eigvals_HFSqr"][j] for j in 1:2^XL]

#DeltaE_Eqr
DeltaESqrED_Dir=E_Sqr_ED_dir - E_Sqr_UF
DeltaESqrED_ADAPT=E_Sqr_ED_ADAPT - E_Sqr_UF
println("DeltaESqr3 : ", DeltaESqrED_Dir)
println("DeltaESqr3 : ", DeltaESqrED_ADAPT)
iters = 1:2^XL
to_plot = [DeltaESqrED_Dir DeltaESqrED_ADAPT]
scalefontsizes(1.1)
labels = ["DeltaESqrED_Dir" "DeltaESqrED_ADAPT"]
plot(iters, to_plot, xlabel="Index", ylabel="Delta_Quasienergy_Square", xticks = iters[1]:iters[end], lw = 3, label = labels,markershape=:circle, markersize=4)
Output_plot = joinpath(Output_folder, "Delta_Quasienergy_Square_comparison1.pdf")
savefig(Output_plot)


# iters = 1:length(energies_Sqr)
# to_plot = [result2["Eigvals_HFSqr"][1:2^XL], result3["Eigvals_HFSqr"][1:2^XL], energies_Sqr]
# to_plot1 = [ DeltaESqr]
# #to_plot= energies
# scalefontsizes(1.1)
# labels = ["Ext_HFSqr" "ADAPT_Ext_HFSqr" "Direct_UF_Sqr"]
# #labels = ["Quasienergy_UF" "Quasienergy(HF)" "Quasienergy(HFSqr_Expt)"]
# plot(iters, to_plot, xlabel="Index", ylabel="energy", xticks = iters[1]:iters[end], lw = 3, label = labels,markershape=:circle, markersize=4)
# Output_plot = joinpath(Output_folder, "Quasienergy_Square_comparison_1.pdf")
# savefig(Output_plot)
# labels = ["DeltaESqr(DirectUF-ADAPT)"]
# #labels = ["Quasienergy_UF" "Quasienergy(HF)" "Quasienergy(HFSqr_Expt)"]
# # plot(iters, to_plot, xlabel="Index", ylabel="energy", xticks = iters[1]:iters[end], lw = 3, label = labels,markershape=:circle, markersize=4)
# # Output_plot = joinpath(Output_folder, "Quasienergy_Square_comparison_0.pdf")
# # savefig(Output_plot)

# ###################                                        ############################
# ############## Floquet State and Expectation value of an operator as a function of time.
# #####################                                       ############################



Final_states = result3["Eigvecs_HFSqr"];
xArr=[(Final_states[1][j]["re"]+im*Final_states[1][j]["im"]) for j in 1:2^(XL_tot)]
#display(Final_states[1])
#display(xArr)

# function Split_Ext_Vector(x_Vec,L,na)
#   NFZ = 2^(na-1) - 1
#   n = [2^L for j in 1:2^na] #Sorted in 
#   ϕ_res = Vector{Vector{eltype(x_Vec)}}()
#   start = firstindex(x_Vec)
#   for len in n
#     push!(ϕ_res, x_Vec[start:(start + len - 1)])
#     start += len
#   end
#   return ϕ_res
# end
phi_n= Split_Ext_Vector(xArr,XL,Xna)
# display(phi_n[1])
# display(One_Spin_Sum_Operator(2)["SumX"])
# display(Two_Spin_Sum_Operator(3,false)["SumXX"])

SumX3=One_Spin_Sum_Operator(3)["SumX"]
SumY3=One_Spin_Sum_Operator(3)["SumY"]
SumZ3=One_Spin_Sum_Operator(3)["SumZ"]

#######
XOmg=6.0
Mesht= (2π/XOmg)*[t for t in -1:0.01:1]

data1=[Op_Vec_Expectation(SumX3,phi_n,Xna,XOmg,Mesht[j]) for j in 1:length(Mesht)]

display(data1[1:10])

iter1 = 1:length(Mesht)
to_plot1 = real(data1)
#to_plot= energies
scalefontsizes(1.1)
labels = L"\sum_{j}{X_j}"
#labels = ["Quasienergy_UF" "Quasienergy(HF)" "Quasienergy(HFSqr_Expt)"]
plot(Mesht, to_plot1, xlabel="time", ylabel="Expt_SumX3", xticks = -1.0:0.5:1.0, lw = 3, label = labels,markershape=:circle, markersize=4)
#Output_plot1 = joinpath(Output_folder, "Expt_SumX3.pdf")
#savefig(Output_plot1)

#########
Floq_state_arr = result1["Eig_Vectors_Sqr"]
Floq_state=[(Floq_state_arr[1][j]["re"]+im*Floq_state_arr[1][j]["im"]) for j in 1:2^XL]
println("Floquet State is: ", Floq_state)


# Floq_state_arr = result1["Eigvecs"]
# Floq_state=[(Floq_state_arr[1][j]["re"]+im*Floq_state_arr[1][j]["im"]) for j in 1:2^XL]
# #println("Floquet State is: ", Floq_state)

E_Sqr0 = result1["Eig_Values_Sqr"]
E_Sqr=[(E_Sqr0[j]["re"]+im*E_Sqr0[j]["im"]) for j in 1:2^XL]
#println(" Squared Quasienergies are: ", E_Sqr)

function Op_Vec_Ut_Expt(Op,Vec)
    Ut_Mesht_arr = result1["Ut_Mesht"]
    Op_Expt_arr = []
    for j in 1:length(Ut_Mesht_arr)
        ED_Ut_Mesht = Ut_Mesht_arr[j]
        rows = length(ED_Ut_Mesht)
        cols = length(ED_Ut_Mesht[1])
        #ED_Ut_Mesht_matrix = Array{Complex{Float64}}(undef, rows, cols)
        ED_Ut_Mesht_matrix= zeros(Complex{Float64}, rows, cols)
        for i in 1:rows
            for j in 1:cols
                re1 = ED_Ut_Mesht[i][j]["re"]
                im1 = ED_Ut_Mesht[i][j]["im"]
                ED_Ut_Mesht_matrix[i, j] = Complex(re1, im1)
            end
        end
        #display(adjoint(ED_Ut_Mesht_matrix)*ED_Ut_Mesht_matrix)
        Op_Expt0=adjoint(Vec)*adjoint(transpose(ED_Ut_Mesht_matrix))*Op*transpose(ED_Ut_Mesht_matrix)*Vec
        #Op_Expt0=adjoint(Vec)*transpose(ED_Ut_Mesht_matrix)*Op*adjoint(transpose(ED_Ut_Mesht_matrix))*Vec
        push!(Op_Expt_arr, Op_Expt0)  
    end
    return Op_Expt_arr
end
    dataX3= Op_Vec_Ut_Expt(SumX3,Floq_state)
    Nt=length(dataX3)
    iter0 = 1:Nt
    #iter3=((2π/(XOmg*Nt))*iter0)
    iter3=((1/Nt)*iter0)
    to_plot3 = real(dataX3)
    scalefontsizes(1.1)
    labels = "ED_UF"
    #labels = ["Quasienergy_UF" "Quasienergy(HF)" "Quasienergy(HFSqr_Expt)"]
    plot(iter3, to_plot3, xlabel="time(t/T)", ylabel= L"\langle \psi(t)|\sum_{j}{X_j} |\psi(t)\rangle", xticks = 0.0:0.2:1.0, lw = 2, label = labels,markershape=:circle, markersize=4)
    Output_plot1 = joinpath(Output_folder, "Expt_SumX3_Ut1.pdf")
    savefig(Output_plot1)

    dataZ3= Op_Vec_Ut_Expt(SumZ3,Floq_state)
    Nt=length(dataZ3)
    iter0 = 1:Nt
    iter3=((1/Nt)*iter0)#((2π/(XOmg*Nt))*iter0)
    to_plot3 = real(dataZ3)
    scalefontsizes(1.1)
    labels = L"\sum_{j}{Z_j}"
    #labels = ["Quasienergy_UF" "Quasienergy(HF)" "Quasienergy(HFSqr_Expt)"]
    plot(iter3, to_plot3, xlabel="time(t/T)", ylabel= L"\langle \psi(t)|\sum_{j}{Z_j} |\psi(t)\rangle", xticks = 0.0:0.2:1.0, lw = 2, label = labels,markershape=:circle, markersize=4)
    Output_plot1 = joinpath(Output_folder, "Expt_SumZ3_Ut1.pdf")
    savefig(Output_plot1)

"""
  Flouqet state and  Expectation value of an operator as a function of time. 
"""
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

# xyz_matrix

 #Split_Ext_Vector(x_Vec, L, na),
    # Op_Vec_Expectation(op, NFZ, t),
#     Op_Vec_Expectation_Low,
#     One_Spin_Sum_Operator(L),
#     Two_Spin_Sum_Operator(L,PBCs)

"""
Change in Quasienergies and Quasienergies squared between different Hamiltonians of different physical sites.
"""














"""
SVD and Entanglement entropy between auxialiary and physical space. L= phyiscal sites, na=auxiliary sites.
"""

