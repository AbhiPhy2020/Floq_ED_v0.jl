# Go into the directory where you want Run_test_Data.jl is and the write following script in mac terminal
#julia --project="~/Documents/Floq_Exact_Diag/Floq_ED_v0" Run_test_Data.jl

using Floq_ED_v0  # Replace "MyPackage" with the name of the package you want to use
using JSON

function run_my_Julia_package2(L,JbarVec,dJVec,Bx,By,Bz,dBx,dBy,dBz,Omg,No_Fourier_Ham, No_Floq_Zones,PBCs,na,r,tol)
    L_tot=L+na
    result_ED_Ext_Ham_dir= Eig_Floq_Ham_Square(L,JbarVec,dJVec,Bx,By,Bz,dBx,dBy,dBz,Omg,No_Fourier_Ham,No_Floq_Zones,PBCs)  # Replace "some_function" with the actual function name
    result_ED_Floq_xyz_ADAPT=Floq_xyz_model(L,na,r,Omg,JbarVec[1],JbarVec[2],JbarVec[3],dJVec[1],dJVec[2],dJVec[3],true,true,tol,PBCs)
    Check_Ham_Floq_ADAPT= zeros(2^L_tot,2^L_tot)
    for j in 1:length(result_ED_Floq_xyz_ADAPT[2])
        Check_Ham_Floq_ADAPT += result_ED_Floq_xyz_ADAPT[2][j]*PauliString_to_Matrix(result_ED_Floq_xyz_ADAPT[1][j]) 
    end
    new_matrix_ADAPT= Check_Ham_Floq_ADAPT[1:end-2^L, 1:end-2^L]
    new_matrix_dir=result_ED_Ext_Ham_dir[7]
    new_matrix_ADAPT1= Check_Ham_Floq_ADAPT[(end+1-2^L):end, (end+1-2^L):end]
    #display(new_matrix_ADAPT1)
    return new_matrix_dir,new_matrix_ADAPT ,new_matrix_ADAPT1
end
##
println("Start of the execution...")

Xtol= 1.0e-11
(XL,Xna)=(2,2)
Xr=1
XL_tot=Xna+XL
XOmg=6.0
(Jx0,Jy0,Jz0,dJx0,dJy0,dJz0)= (2.7,2.7,3.4,1.0,1.0,1.0)
(Bx0,By0,Bz0,dBx0,dBy0,dBz0)= (0.0,0.0,0.0,0.0,0.0,0.0)
XNFZ= 2^(Xna-1)-1
XLZ_tot= XNFZ*(2^XL)
XJbarVec=[Jx0,Jy0,Jz0]
XdJVec=[1.0,1.0,1.0]
# new_data=run_my_Julia_package2(XL,XJbarVec,XdJVec,Bx0,By0,Bz0,dBx0,dBy0,dBz0,XOmg,1,XNFZ,false,Xna,Xr,Xtol)
# display(new_data[1])
# display(new_data[1]-new_data[2])
# display(new_data[3])
# display(PauliString_to_Matrix("XY"))
# display(Str_Coeff_Floq_Diag(Xna,XOmg, Xtol))
# display(generate_pauli_string(3))
#display(Str_Coeff_Floq_OffDiag(Xna,Xr,Xtol))
#display(Ham_Aux_Sym_Part(Xna, Xr))
#display(All_Aux_Ham(Xna,Xr,XOmg,true, true,Xtol))
