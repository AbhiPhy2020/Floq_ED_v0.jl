"""
Go into the directory where you want Run_test_aux_matrix_Data.jl is and the write following script in mac terminal
julia --project="~/Documents/Floq_Exact_Diag/Floq_ED_v0" Run_test_Data.jl

We will test the EDADAPT and ED_Direct in this file
"""

using Floq_ED_v0  # Replace "MyPackage" with the name of the package you want to use
using JSON
using LinearAlgebra
using Random

function run_my_Julia_package2(L,na,JbarVec,dJVec,Bx,By,Bz,dBx,dBy,dBz,Omg,No_Fourier_Ham,PBCs,tol)
    L_tot=L+na
    rJ=rB=No_Fourier_Ham
    NFZ= (2^(na-1) - 1) #No_Floq_Zones
    (Jx,Jy,Jz,dJx,dJy,dJz)=(JbarVec[1],JbarVec[2],JbarVec[3],dJVec[1],dJVec[2],dJVec[3])
    result_ED_Ext_Ham_dir=Eig_Floq_Ham_Square(L,NFZ,Omg,JbarVec,dJVec,Bx,By,Bz,dBx,dBy,dBz,No_Fourier_Ham,PBCs)
    #result_ED_Floq_xyz_ADAPT=Floq_xyz_model(L,na,r,Omg,JbarVec[1],JbarVec[2],JbarVec[3],dJVec[1],dJVec[2],dJVec[3],true,true,tol,PBCs)
    result_ED_Floq_xyz_ADAPT=ED_Floq_Ham_XYZ_B(L,na,rJ,rB,Omg,Jx,Jy,Jz,dJx,dJy,dJz,Bx,By,Bz,dBx,dBy,dBz,tol,PBCs)

    Check_Ham_Floq_ADAPT= zeros(2^L_tot,2^L_tot)
    Check_Ham_Floq_ADAPT=result_ED_Floq_xyz_ADAPT[7] 

    new_matrix_ADAPT= Check_Ham_Floq_ADAPT[1:end-2^L, 1:end-2^L]
    new_matrix_dir=result_ED_Ext_Ham_dir[7]
    new_matrix_ADAPT_End= Check_Ham_Floq_ADAPT[(end+1-2^L):end, (end+1-2^L):end]
    #display(new_matrix_ADAPT1)

    # Compare the matrices using isapprox
    tolerance = 1e-12
    are_equal = isapprox(new_matrix_ADAPT, new_matrix_dir, atol=tolerance)
    println("The matrices are approximately equal: $are_equal")

    return new_matrix_dir,new_matrix_ADAPT ,new_matrix_ADAPT_End
end
##
println("Start of the execution...")

Xtol= 1.0e-12
(XL,Xna)=(2,3)
XNFZ= (2^(Xna-1) - 1) #No_Floq_Zones
Xr=1
XL_tot=Xna+XL
XOmg=6.0
(Jx0,Jy0,Jz0,dJx0,dJy0,dJz0)= (0.0,0.0,5.0,0.0,0.0,0.0)
arr_rand = 0.5*rand(6)
(Bx0,By0,Bz0,dBx0,dBy0,dBz0)= (arr_rand[1],arr_rand[2],arr_rand[3],arr_rand[4],arr_rand[5],arr_rand[6])
XNF_Ham=1
XLZ_tot= XNFZ*(2^XL)
# XJbarVec=[Jx0,Jy0,Jz0] #5.0 * rand(3)  #
# XdJVec= [dJx0,dJy0,dJz0] #5.0 * rand(3)  #
XJbarVec=5.0 * rand(3)  
XdJVec= 5.0 * rand(3)  
new_data=run_my_Julia_package2(XL,Xna,XJbarVec,XdJVec,Bx0,By0,Bz0,dBx0,dBy0,dBz0,XOmg,XNF_Ham,false,Xtol)
 #display(new_data[1])
 display(new_data[1]-new_data[2])
 #display(new_data[3])
 #array_floats = 5.0 * rand(6)
 #display(array_floats)
if all(new_data[1].==new_data[2])
    println("Check:Both matrices are equal")
else
    println("Check:Both matrices are not equal")
end 
Xtolerance = 1e-12
    are_equal = isapprox(new_data[1],new_data[2], atol=Xtolerance)
    println("The matrices are approximately equal: $are_equal")

   # display(eigvals(new_data[2]))    
"""
Testing all functions output 
"""
#Xtol= 1.0e-12
# display(PauliString_to_Matrix("XY"))
#display(Str_Coeff_Floq_Diag(4,1.0, Xtol))
# println("pauli matrices for 5 qubits",Str_Coeff_Floq_Diag(5,1.0, Xtol))
# println("pauli matrices for 6 qubits",Str_Coeff_Floq_Diag(6,1.0, Xtol))
#dd1= @timed Str_Coeff_Floq_Diag(8,1.0, Xtol)
# display(Str_Coeff_Floq_Diag(6,1.0, Xtol))
# display(Str_Coeff_Floq_Diag(7,1.0, Xtol))
# display(Str_Coeff_Floq_Diag(8,1.0, Xtol))
# println("time for Str_Coeff_Floq_Diag is ",dd1[2],"s")
# display(dd1[1])
# display(generate_pauli_string(3))
#display(Str_Coeff_Floq_OffDiag(Xna,Xr,Xtol))
#display(Ham_Aux_Sym_Part(Xna, Xr))
#display(All_Aux_Ham(Xna,Xr,XOmg,true, true,Xtol))
# dd2= @timed Aux_Diag_Str_Coeff(8,1.0)
# println("time for Aux_Diag_Str_Coeff is ",dd2[2],"s")
# display(dd2[1])