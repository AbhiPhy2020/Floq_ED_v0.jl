# Go into the directory where you want Run_test_Data.jl is and the write following script in mac terminal
#julia --project="~/Documents/Floq_Exact_Diag/Floq_ED_v0" Run_test_Data.jl

using Floq_ED_v0  # MyPackage
using JSON
using JSON3
using JLD
using Plots

# Define a function to run source code from the package
"""
run_ED_Floq_Sym(L,na,Nt_mesh,Omg,JbarVec,dJVec,Bx,By,Bz,dBx,dBy,dBz,No_Fourier_Ham,PBCs)
"""
function run_ED_Floq_Sym(L,na,Nt_mesh,Omg,JbarVec,dJVec,Bx,By,Bz,dBx,dBy,dBz,No_Fourier_Ham,PBCs)
    # myparams = Dict("PBCs" => PBCs) #Do I need dic for all the parameters?
    Xtol=1.0e-12
    NFZ= (2^(na-1) - 1) #No_Floq_Zones
    rJ=rB=No_Fourier_Ham
    (Jx,Jy,Jz)= (JbarVec[1],JbarVec[2],JbarVec[3]) 
    (dJx,dJy,dJz)= (dJVec[1],dJVec[2],dJVec[3])
    
    #Call the function from the package
    res_ED_UF0= @timed Direct_Floq_Sys(L,Nt_mesh,Omg,JbarVec,dJVec,Bx,By,Bz,dBx,dBy,dBz,PBCs)
    res_ED_UF = res_ED_UF0[1] #Abhi: This is the output of the adapt_vqe
    elapsedtime = res_ED_UF0[2]
    println("time for ED calculation using UF is ",elapsedtime,"s")
    
    res_ED_HFSqr0= @timed Eig_Floq_Ham_Square(L,NFZ,Omg,JbarVec,dJVec,Bx,By,Bz,dBx,dBy,dBz,No_Fourier_Ham,PBCs) 
    res_ED_HFSqr = res_ED_HFSqr0[1] #Abhi: This is the output of the adapt_vqe
    elapsedtime = res_ED_HFSqr0[2]
    println("time for ED calculation using HF and HFSqr in extended Floquet method is ",elapsedtime,"s")

    res_ED_ADAPT0= @timed ED_Floq_Ham_XYZ_B(L,na,rJ,rB,Omg,Jx,Jy,Jz,dJx,dJy,dJz,Bx,By,Bz,dBx,dBy,dBz,Xtol,PBCs)
    res_ED_ADAPT = res_ED_ADAPT0[1] 
    elapsedtime = res_ED_ADAPT0[2]
    println("time for ED calculation using HF and HFSqr (ADAPT_style) extended Floquet method is ",elapsedtime,"s")
    
    # Save the data in the following folder
    folder_path="/Users/abhishek22/Documents/Floq_Exact_Diag/Floq_ED_v0/Examples/Output_Data_Plots/Data"

    # To store data in a file
    res1_string = JSON.json(Dict("Quasienergies" => res_ED_UF[1], "Eigvecs" =>res_ED_UF[2],"Ut_Mesht" => res_ED_UF[3]))
    file_name = "ED_UF_L_Nt_mesh_Omg_JxJyJz_dJxdJydJz_BxByBz_dBxdBydBz"*string((L,Nt_mesh,Omg,JbarVec,dJVec,Bx,By,Bz,dBx,dBy,dBz,No_Fourier_Ham))*"_results.dat"
    file_path = joinpath(folder_path, file_name)
    io1 = open(file_path,"w")
    write(io1,res1_string*"\n")
    close(io1)

    res2_string = JSON.json(Dict("Eigvals_HF" => res_ED_HFSqr[1], "Eigvecs_HF" => res_ED_HFSqr[2],"Eigvals_HFSqr" => res_ED_HFSqr[3],
    "Eigvecs_HFSqr" => res_ED_HFSqr[4], "Eigvals_HFSqrExpt" => res_ED_HFSqr[5],"Eigvecs_HFSqrExpt" => res_ED_HFSqr[6],"HF" => res_ED_HFSqr[7]))
    file_name= "ED_Ext_HFSqr_L_NFZ_Omg_JxJyJz_dJxdJydJz_BxByBz_dBxdBydBz"*string((L,NFZ,Omg,JbarVec,dJVec,Bx,By,Bz,dBx,dBy,dBz))*"_results.dat"
    file_path = joinpath(folder_path, file_name)
    io2= open(file_path,"w")
    write(io2,res2_string*"\n")
    close(io2)

    res3_string = JSON.json(Dict("Eigvals_HF" => res_ED_ADAPT[1], "Eigvecs_HF" => res_ED_ADAPT[2],"Eigvals_HFSqr" => res_ED_ADAPT[3],
    "Eigvecs_HFSqr" => res_ED_ADAPT[4], "Eigvals_HFSqrExpt" => res_ED_ADAPT[5],"Eigvecs_HFSqrExpt" => res_ED_ADAPT[6],"HF" => res_ED_ADAPT[7]))
    file_name= "ED_HF_ADAPT_L_na_Omg_JxJyJz_dJxdJydJz_BxByBz_dBxdBydBz"*string((L,na,Omg,JbarVec,dJVec,Bx,By,Bz,dBx,dBy,dBz))*"_results.dat"
    file_path = joinpath(folder_path, file_name)
    io3= open(file_path,"w")
    write(io3,res3_string*"\n")
    close(io3)
end

"""
In this execution, we will run the code for different Hamiltonians and compare the results.
"""

println("Start of the execution...")

"""
XXZ Hamiltonian with and wthout magnetic field
"""

#(XL,XJbarVec,XdJVec,XBx,XBy,XBz,XdBx,XdBy,XdBz,XΩ,XNF_Ham, XNF_Zones)=(3,[3.7,2.7,1.4],[1.0,1.4,1.0],2.3,3.0,0.0,0.0,0.0,0.0,6.0,1,7)
(XL,Xna,XJbarVec,XdJVec,XBx,XBy,XBz,XdBx,XdBy,XdBz,XOmg,XNF_Ham)=(3,4,[3.7,2.7,1.4],[1.0,1.4,1.0],2.3,3.0,4.0,0.0,0.0,0.0,6.0,1)
XNt_mesh=5000
run_ED_Floq_Sym(XL,Xna,XNt_mesh,XOmg,XJbarVec,XdJVec,XBx,XBy,XBz,XdBx,XdBy,XdBz,XNF_Ham,false)
#run_ED_Floq_Sym(L,No_Floq_Zones,Nt_mesh,Omg,JbarVec,dJVec,Bx,By,Bz,dBx,dBy,dBz,No_Fourier_Ham,PBCs)
#run_ED_Floq_Sym(L,JbarVec,dJVec,Bx,By,Bz,dBx,dBy,dBz,Omg,No_Fourier_Ham, No_Floq_Zones,Nt_mesh,PBCs)
# """
# XYZ Hamiltonian with and wthout magnetic field
# """

# #(XL,XJbarVec,XdJVec,XBx,XBy,XBz,XdBx,XdBy,XdBz,XΩ,XNF_Ham, XNF_Zones)=(3,[3.7,2.7,1.4],[1.0,1.4,1.0],2.3,3.0,0.0,0.0,0.0,0.0,6.0,1,7)
# (XL,XJbarVec,XdJVec,XBx,XBy,XBz,XdBx,XdBy,XdBz,XΩ,XNF_Ham, XNF_Zones)=(3,[0.0,0.0,0.0],[0.0,0.0,0.0],2.3,3.0,4.0,0.0,0.0,0.0,6.0,1,7)
# run_ED_Floq_Sym(L,JbarVec,dJVec,Bx,By,Bz,dBx,dBy,dBz,Omg,No_Fourier_Ham, No_Floq_Zones,PBCs)

# """
# Heisenberg Hamiltonian with and wthout magnetic field
# """

# #(XL,XJbarVec,XdJVec,XBx,XBy,XBz,XdBx,XdBy,XdBz,XΩ,XNF_Ham, XNF_Zones)=(3,[3.7,2.7,1.4],[1.0,1.4,1.0],2.3,3.0,0.0,0.0,0.0,0.0,6.0,1,7)
# (XL,XJbarVec,XdJVec,XBx,XBy,XBz,XdBx,XdBy,XdBz,XΩ,XNF_Ham, XNF_Zones)=(3,[0.0,0.0,0.0],[0.0,0.0,0.0],2.3,3.0,4.0,0.0,0.0,0.0,6.0,1,7)
# run_ED_Floq_Sym(L,JbarVec,dJVec,Bx,By,Bz,dBx,dBy,dBz,Omg,No_Fourier_Ham, No_Floq_Zones,PBCs)

# """
# Ising Hamiltonian with and wthout magnetic field
# """

# #(XL,XJbarVec,XdJVec,XBx,XBy,XBz,XdBx,XdBy,XdBz,XΩ,XNF_Ham, XNF_Zones)=(3,[3.7,2.7,1.4],[1.0,1.4,1.0],2.3,3.0,0.0,0.0,0.0,0.0,6.0,1,7)
# (XL,XJbarVec,XdJVec,XBx,XBy,XBz,XdBx,XdBy,XdBz,XΩ,XNF_Ham, XNF_Zones)=(3,[0.0,0.0,0.0],[0.0,0.0,0.0],2.3,3.0,4.0,0.0,0.0,0.0,6.0,1,7)
# run_ED_Floq_Sym(L,JbarVec,dJVec,Bx,By,Bz,dBx,dBy,dBz,Omg,No_Fourier_Ham, No_Floq_Zones,PBCs)


##############
##############
#####States analysis


"""
Working examples
Comparison between direct and extended Floquet Hamiltonian  
Different level of checks
1. For large frequency, without taking any time periodic term. Start with simple Ising and then take arbitrary XYZ matrix
2. Add time periodic term and then compare the result in the two cases.Start with simple Ising and then take arbitrary XYZ matrix
3. 
"""
