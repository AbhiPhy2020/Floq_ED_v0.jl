""" 
In this file, we will implement the following functions:
1. Exact diagolization (ED) of time periodic Hamiltonian Direct calculation (using U_F); 
2. ED of Extended Floquet Hamiltonian and its square using extended Floquet Hamiltonian method.
3. ED of Extended Floquet Hamiltonian and its square using extended Floquet Hamiltonian method with an auxiliary qubit.
For all of them we will use a static XYZ Hamilton with static magnetic field which is written at the beginning.
"""

using BenchmarkTools
using BlockDiagonals
using Kronecker
using LinearAlgebra
using Plots

###########
I2=Id=[1.0 0.0 ; 0.0 1.0]
s1=sx=[0.0 1.0 ; 1.0 0.0]
s2=sy=[0.0 -im*1.0 ; im*1.0 0.0]
s3=sz=[1.0 0.0 ; 0.0 -1.0];
              
#################################
# Static Spin Hamiltonian  (A general two spin XYZ model)
# functions to construct matrix forms of models - for purposes of exact diagonalization
##################################

#1:Define a dictionary where the keys are strings and values are the Pauli matrices alongwith identity matrix.
pdict = Dict('I' => [1.0+0.0im 0.0+0.0im; 0.0+0.0im 1.0+0.0im], 'X' => [0.0+0.0im 1.0+0.0im; 1.0+0.0im 0.0+0.0im],
             'Y' => [0.0+0.0im 0.0-1.0im; 0.0+1.0im 0.0+0.0im], 'Z' => [1.0+0.0im 0.0+0.0im; 0.0+0.0im -1.0+0.0im])

#2:Converting strings into kronecker product of Pauli matrices alongwith identity matrix.
function paulis_matrix(pstr)
    res = copy(pdict[pstr[1]]) #copy 
    for ch in pstr[2:length(pstr)]
        res = kron(res,pdict[ch])
    end
    return res
end
#3:Creating a function with  the input parameters and agenerate a matrix. 
function xyz_matrix(L,Jx,Jy,Jz,PBCs)
    hamJ = zeros(ComplexF64,2^L,2^L)
    for site=1:L-1
        xx_str = "I"^(site-1)*"XX"*"I"^(L-(site+1))
        yy_str = "I"^(site-1)*"YY"*"I"^(L-(site+1))
        zz_str = "I"^(site-1)*"ZZ"*"I"^(L-(site+1))
        hamJ = hamJ + Jx*paulis_matrix(xx_str) + Jy*paulis_matrix(yy_str) + Jz*paulis_matrix(zz_str)
    end
    ##If we have inhomogeneous coupling then we can put Jx as a array with length L and multiply each time.
    if PBCs
        xx_str = "X"*"I"^(L-2)*"X"
        yy_str = "Y"*"I"^(L-2)*"Y"
        zz_str = "Z"*"I"^(L-2)*"Z"
        hamJ = hamJ + Jx*paulis_matrix(xx_str) + Jy*paulis_matrix(yy_str) + Jz*paulis_matrix(zz_str)
    end
    return hamJ
end
####
function xxz_matrix(L,Jxy,Jz,PBCs)
    return xyz_matrix(L,Jxy,Jxy,Jz,PBCs)
end

function heisenberg_matrix(L,J,PBCs)
    return xyz_matrix(L,J,J,J,PBCs)
end

##Hamiltonian related to the Magnectic field 
function magnetic_field_matrix(L,Bx,By,Bz)
    hamB = zeros(ComplexF64,2^L,2^L)
    if typeof(Bx) == Float64
        Bxvec = fill(Bx,L)
    else
        Bxvec = Bx
    end
    if typeof(By) == Float64
        Byvec = fill(By,L)
    else
        Byvec = By
    end
    if typeof(Bz) == Float64
        Bzvec = fill(Bz,L)
    else
        Bzvec = Bz
    end
    for site=1:L
        Bx_str = "I"^(site-1)*"X"*"I"^(L-site)
        By_str = "I"^(site-1)*"Y"*"I"^(L-site)
        Bz_str = "I"^(site-1)*"Z"*"I"^(L-site)
        hamB = hamB + Bxvec[site]*paulis_matrix(Bx_str) + Byvec[site]*paulis_matrix(By_str) + Bzvec[site]*paulis_matrix(Bz_str)
    end
    return hamB
end
      
"""
Direct_Floq_Sys(L,Nt_mesh,ω0,JbarVec,dJVec,Bx,By,Bz,dBx,dBy,dBz,PBCs)
    Direct method of calculation for a time periodic Hamiltonian. 
"""
function Direct_Floq_Sys(L,Nt_mesh,ω0,JbarVec,dJVec,Bx,By,Bz,dBx,dBy,dBz,PBCs)
    #Generation of the mesht and corresponding Bt
    dt=2π/(Nt_mesh*ω0)
    Mesht=[t for t in dt:dt:(2π/ω0)]
    Mesh_Jvec_Coeff=[(JbarVec+ dJVec*cos(ω0*t)) for t in Mesht] # Need to change this for differe periodic fucnction
    Mesh_Bvec_Coeff=[([Bx,By,Bz]+[dBx,dBy,dBz]*cos(ω0*t)) for t in Mesht] # Need to change this for different periodic fucnction
    len_T=length(Mesht)
    Ut0=collect(kronecker(Id,L))
    Ut_Mesh0=Array{Matrix{ComplexF64},1}(undef,0)
    Hj=Matrix{ComplexF64}
    for j in 1:len_T
        (Jxtj,Jytj,Jztj)=Mesh_Jvec_Coeff[j]
        (Bxtj,Bytj,Bztj)=Mesh_Bvec_Coeff[j]
        Hj=xyz_matrix(L,Jxtj,Jytj,Jztj,PBCs)+ magnetic_field_matrix(L,Bxtj,Bytj,Bztj)
        Ut0= exp(-im*dt*Hj)*Ut0 # deltaT # Does the change in ordering matter: Yes if there is no time reversal symmetry
        push!(Ut_Mesh0,Ut0)
    end
    HF= im*(ω0/(2π))*log(Ut_Mesh0[len_T])
    #display(length(Ut_Mesh0))
    #log(eigvals(Ut_Mesh0[len_T]))
    eigsys= eigen(HF)
    Eig_HF_Vals= eigsys.values
    Eig_HF_Vec= eigsys.vectors
    Eig_HF_Vec= Eig_HF_Vec[:,sortperm(real(Eig_HF_Vals))] # To sort eigenvectors according to the eigenvalues
    return Eig_HF_Vals, Eig_HF_Vec, Ut_Mesh0 
end 

"""
Eig_Floq_Ham_Square(L,JbarVec,dJVec,Bx,By,Bz,dBx,dBy,dBz,Ω,No_Fourier_Ham, No_Floq_Zones,PBCs)
    The function calculates the eigenvalues and eigenvectors of the Floquet Hamiltonian and its square for a time 
    periodic Hamiltonian.
"""
function Eig_Floq_Ham_Square(L,No_Floq_Zones,Ω,JbarVec,dJVec,Bx,By,Bz,dBx,dBy,dBz,No_Fourier_Ham,PBCs)
    (Jxbar,Jybar,Jzbar)=(JbarVec)
    (Jx1,Jy1,Jz1)=0.5*(dJVec)
    (Bxbar,Bybar,Bzbar)=(Bx,By,Bz)
    (Bx1,By1,Bz1)=(0.5*dBx,0.5*dBy,0.5*dBz)
    Ω0=Ω
    NFHam=No_Fourier_Ham
    NFZ=No_Floq_Zones
    #Construction of the Floquet part of the matrix
    Ham_Aux_Ω0=Matrix{Float64}
    Ham_Aux_Ω0=diagm(0=>[(j-NFZ-1)*Ω0 for j in  1:(2*NFZ+1)])
    Ham_Aux=Array{Matrix{Float64},1}(undef,NFHam) # NFHam shows the number of matrices in Ham_Aux 
    for j in 1:NFHam
        Ham_Aux[j]=diagm(j=>[1.0 for i in 1:(2*NFZ+1-j)]) #It could be improved further for exact number of modes.
    end
    # Ham_Aux[j]=diagm(j=>[1.0 for i in  1:(2*NFZ+j-1)])
    Ham0= xyz_matrix(L,Jxbar,Jybar,Jzbar,PBCs)+ magnetic_field_matrix(L,Bxbar,Bybar,Bzbar)
    Ham1= xyz_matrix(L,Jx1,Jy1,Jz1,PBCs)+ magnetic_field_matrix(L,Bx1,By1,Bz1)
    HamF_tot =(kron(Matrix(I,2*NFZ+1,2*NFZ+1),Ham0) + kron(Ham_Aux_Ω0,kronecker(Id,L)) + kron(Ham_Aux[1],Ham1)+ kron(transpose(Ham_Aux[1]),Adjoint(Ham1)))
    #
    eigsys= eigen(HamF_tot)
    Eig_Values= eigsys.values
    Eig_Vectors= eigsys.vectors
    Eig_Vectors= Eig_Vectors[:, sortperm(real(Eig_Values))] # This sort the eigenvectors in terms of 
    #Squared Floquet Hamiltonian method
    HFSqur= (HamF_tot)*(HamF_tot)
    eigsys_Squr= eigen(HFSqur)
    Eig_Values_Squr= eigsys_Squr.values
    Eig_Vectors_Squr= eigsys_Squr.vectors
    Eig_Vectors_Squr= Eig_Vectors_Squr[:, sortperm(real(Eig_Values_Squr))] # This sort the eigenvectors in terms of 
    ####
    EigVal_Expt=[Adjoint(Eig_Vectors_Squr[:,j])*HamF_tot*Eig_Vectors_Squr[:,j] for j in 1:2^L]
    Eig_Vectors_Squr_Expt=[Eig_Vectors_Squr[:,j] for j in 1:2^L]
    #Eig_Vectors_Squr_Expt= Eig_Vectors_Squr_Expt[:, sortperm(real(EigVal_Expt))]
    return Eig_Values,Eig_Vectors, Eig_Values_Squr,Eig_Vectors_Squr, EigVal_Expt,Eig_Vectors_Squr_Expt, HamF_tot 
end

"""
Ext_Floq_Sym_Asym_AuxQubit(L,JbarVec,dJVec,Bx,By,Bz,dBx,dBy,dBz,Ω,No_Fourier_Ham, No_Floq_Zones,PBCs)
    The function calculates the eigenvalues and eigenvectors of the extended Floquet Hamiltonian for a time 
    periodic Hamiltonian when the system has an auxiliary qubit (auxilary space has 2^na dimension).
"""
function Ext_Floq_Sym_Asym_AuxQubit(L,JbarVec,dJVec,Bx,By,Bz,dBx,dBy,dBz,Ω,No_Fourier_Ham, No_Floq_Zones,PBCs)
    (Jxbar,Jybar,Jzbar)=(JbarVec)
    (Jx1,Jy1,Jz1)=0.5*(dJVec)
    (Bx1,By1,Bz1)=(0.5*dBx,0.5*dBy,0.5*dBz)
    (Bxbar,Bybar,Bzbar)=(Bx,By,Bz)
     Ω0=Ω
     NFHam=No_Fourier_Ham
     NFZ=No_Floq_Zones
    #Construction of the Floquet part of the matrix
    Ham_Aux_Ω0=Matrix{Float64}
    Ham_Aux_Ω0Mid=diagm(0=>[(j-NFZ-1)*Ω0 for j in  3:(2*NFZ-1)])
    Ham_Aux_Ω0Low=diagm(0=>[(NFZ-2)*Ω0 for j in  1:2])
    Ham_Aux_Ω0=BlockDiagonal([Ham_Aux_Ω0Low,Ham_Aux_Ω0Mid,Ham_Aux_Ω0Low])
    # Ham_Aux_Ω0=diagm(0=>[(j-NFZ-1)*Ω0 for j in  1:(2*NFZ+1)])# Here we are implementing modulo in the Frequency indices
    # Ham_Aux_Ω0=diagm(0=>[(j-NFZ-1)*Ω0 for j in  1:3*(2*NFZ+1)])# Here we are implementing the constant 
    ###
    Ham_Aux=Array{Matrix{Float64},1}(undef,NFHam)
    for j in 1:NFHam
    Ham_Aux[j]=diagm(j=>[1.0 for i in  1:(2*NFZ+j-1)])
    end
    ###
    Ham0= xyz_matrix(L,Jxbar,Jybar,Jzbar,PBCs)+ magnetic_field_matrix(L,Bxbar,Bybar,Bzbar)
    Ham1= xyz_matrix(L,Jx1,Jy1,Jz1,PBCs)+ magnetic_field_matrix(L,Bx1,By1,Bz1)
    ###
    HamF_tot = (kron(Matrix(I,(2*NFZ+1),(2*NFZ+1)),Ham0) + kron(Ham_Aux_Ω0,kronecker(Id,L)) + kron(Ham_Aux[1],Ham1)+ kron(transpose(Ham_Aux[1]),Adjoint(Ham1)))
    ###
    eigsys= eigen(HamF_tot)
    Eig_Values= eigsys.values
    Eig_Vectors= eigsys.vectors
    Eig_Vectors= Eig_Vectors[:, sortperm(real(Eig_Values))] # This sort the eigenvectors in terms of quasienergy values.
    return Eig_Values, Eig_Vectors   
end

#=
################################### This needs to be deleted once we make sure that everything in the square Hamiltonian
is correct.
## Floquet exetended Hamiltonian of time periodic Hamiltonian 
###################################
function Ext_Floq_Sys(L,JbarVec,dJVec,Bx,By,Bz,dBx,dBy,dBz,Ω,No_Fourier_Ham, No_Floq_Zones,PBCs)
    (Jxbar,Jybar,Jzbar)=(JbarVec)
    (Jx1,Jy1,Jz1)=0.5*(dJVec)
    (Bx1,By1,Bz1)=(0.5*dBx,0.5*dBy,0.5*dBz)
    (Bxbar,Bybar,Bzbar)=(Bx,By,Bz)
     Ω0=Ω
     NFHam=No_Fourier_Ham
     NFZ=No_Floq_Zones
    ###Construction of the Floquet part of the matrix
     Ham_Aux_Ω0=Matrix{Float64}
     Ham_Aux_Ω0=diagm(0=>[(j-NFZ-1)*Ω0 for j in  1:(2*NFZ+1)])
    ###
    Ham_Aux=Array{Matrix{Float64},1}(undef,NFHam) # NFHam shows the number of matrices in Ham_Aux 
    for j in 1:NFHam
        #Ham_Aux[j]=diagm(j=>[1.0 for i in  1:(2*NFZ)])
        Ham_Aux[j]=diagm(j=>[1.0 for i in  1:(2*NFZ+j-1)])
    end
    ###
    Ham0= xyz_matrix(L,Jxbar,Jybar,Jzbar,PBCs)+ magnetic_field_matrix(L,Bxbar,Bybar,Bzbar)
    Ham1= xyz_matrix(L,Jx1,Jy1,Jz1,PBCs)+ magnetic_field_matrix(L,Bx1,By1,Bz1)
    ###
    HamF_tot = (kron(Matrix(I,2*NFZ+1,2*NFZ+1),Ham0) + kron(Ham_Aux_Ω0,kronecker(Id,L)) + kron(Ham_Aux[1],Ham1) + kron(transpose(Ham_Aux[1]),Adjoint(Ham1)))
    ###
    eigsys= eigen(HamF_tot)
    Eig_Values= eigsys.values
    Eig_Vectors= eigsys.vectors
    Eig_Vectors= Eig_Vectors[:, sortperm(real(Eig_Values))] # This sort the eigenvectors in terms of quasienergy values.
    return Eig_Values, Eig_Vectors   
end =#
