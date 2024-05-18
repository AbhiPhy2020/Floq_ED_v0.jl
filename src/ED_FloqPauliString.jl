#using LinearAlgebra
#import Pkg; Pkg.add("Kronecker")
#import Pkg; Pkg.add("Plots");
#import Pkg; Pkg.add("BenchmarkTools");
#import Pkg; Pkg.add("BlockDiagonals");


#import Pkg; Pkg.add("PlotlyJS");
using Plots;
using BenchmarkTools;
using BlockDiagonals;
using Kronecker
#using PlotlyJS

#####################
## Here we take input as Pauli strings and then create a matrix from that. PauliString_to_Matrix("XYZI")=
function PauliString_to_Matrix(str::String)
    na=length(str)
   # display(na)
    PmatArr=Array{undef,1}
    #PmatArr=[zeros(2,2) for i in 1:na]
    PmatArr=[Matrix{}(undef,2,2) for _ in 1:na]
    #str0= "IZZ"
    for i in 1:na
        if str[i]== 'I' 
            PmatArr[i]= [1.0 0.0; 0.0 1.0]
        elseif str[i]== 'Z'
            PmatArr[i]= [1.0 0.0; 0.0 -1.0]
        elseif str[i]== 'X'
            PmatArr[i]= [0.0 1.0; 1.0 0.0]
        elseif str[i]== 'Y'
            PmatArr[i]= [0.0 -1*im; 1.0*im 0.0]
        else
            println("Invalid inputs as the diagonal parts should contain Z and I only")
        end
      # display(PmatArr[i])
    end
    MatKron=PmatArr[1];
for i in 2:na
MatKron= kronecker(MatKron,PmatArr[i]) #Here we get the big matrix
end
 return  MatKron  # This return a array of 
#return AmatArray
end
######################
PauliString_to_Matrix("XZ")*PauliString_to_Matrix("XY") #v = [ Matrix{ComplexF64}(undef,3,3) for _ in 1:5 ]
# For diagonal elements, Symmetric and assymetric Paulirepresentation does not matter. 
########################
function Str_Coeff_Floq_Diag(na::Int64,Om::Float64, tol::Float64)
    Ω0= Om # Omega=Om= frequency of the the time periodic Hamiltonian
    Str_ZI_Arr=Array{String, 1}(undef,2^na ) # this will used to all {Z_iI_i} strings
    Floq_Str_ZI_Arr=String[] # this will used to collect strings with only with non-zero coefficients (weights)
    Floq_Coeff_ZI_Arr=Float64[] # this will used to collect non-zero coefficients (weights)
    #1. Construction of Pauli string
    for j in 1:2^na
    Floqstr="Z"^0
    BitStr=reverse(bitstring(j-1))
     #display(BitStr)
for i in 1:na
   if BitStr[i]== '0'
            Floqstr *= 'Z'
        else
             Floqstr *= 'I'
    end       
end
    # display(Floqstr) ##to  display the strings of na length
     Str_ZI_Arr[j]=Floqstr
    end
    #2. Construction of the Floquet part of the matrix
    NFZ= 2^(na-1)
    Ham_Aux_Ω0=Matrix{Float64}
    Ham_Aux_Ω0=diagm(0=>[(j-NFZ)*Ω0 for j in  1:(2^na)])
    ##display(Ham_Aux_Ω0)
    #####################################
    #3.Generation of Strings with nonzero coefficients
    for j in 1:2^na
       Str_Coeff=tr(Ham_Aux_Ω0*PauliString_to_Matrix(Str_ZI_Arr[j]))/(2^na)
       if (abs(Str_Coeff)>tol)
        push!(Floq_Coeff_ZI_Arr,Str_Coeff)
        push!(Floq_Str_ZI_Arr,Str_ZI_Arr[j])
       end
    end
  #4. Checking the 
CheckHam_Aux_Ω0= zeros(size(PauliString_to_Matrix(Floq_Str_ZI_Arr[2])))
for j in 1:length(Floq_Coeff_ZI_Arr)
    CheckHam_Aux_Ω0 += Floq_Coeff_ZI_Arr[j]*PauliString_to_Matrix(Floq_Str_ZI_Arr[j]) 
end
#display(CheckHam_Aux_Ω0)
if all(CheckHam_Aux_Ω0.==Ham_Aux_Ω0)
 println("Check:The resconstructed matrix is equal to the initial auxiliary matrix")
else
 println("Check:The resconstructed matrix does not match the initial auxiliary matrix")
end 
    return Floq_Str_ZI_Arr, Floq_Coeff_ZI_Arr#, CheckHam_Aux_Ω0
end

#############
##############
Xtol= 1.0e-11;
Xna=5;
XOm=6.0;
Str_Coeff_Floq_Diag(Xna,XOm,Xtol)
### Off Diagonal terms using trace
function generate_pauli_string(N::Int)
    pauli_operators = ["I", "X", "Y"]
    num_operators = length(pauli_operators)
    # Generate all possible combinations of Pauli strings
    pauli_strings = String[]
    for i in 1:num_operators^N
        pauli_string = ""
        for j in 1:N
            # Map the index to the corresponding Pauli operator
            index = div(i - 1, num_operators^(j - 1)) % num_operators + 1
           # display(index)
            pauli_string=pauli_string*pauli_operators[index]
        end
        #display(pauli_string)
        push!(pauli_strings, pauli_string)
    end
    return pauli_strings
end
#########
##########
# For diagonal elements, Symmetric and assymetric Paulirepresentation does not matter. 
function Str_Coeff_Floq_OffDiag(na::Int64, r::Int64, tol::Float64)
    r0=r # Fourier coefficient of time periodic Hamiltonian
    #if(r0>=2^na-1)
      #  println("the number of modes can't be greater $2^na")
    Str_XYI_Arr=String[]#Array{String, 1}(undef,3^na) # this will used to all {Z_iI_i} strings
    Floq_Str_XYI_Arr=String[] # this will used to collect strings with only with non-zero coefficients (weights)
    Floq_Coeff_XYI_Arr=ComplexF64[] # this will used to collect non-zero coefficients (weights)
   # Str_Coeff::ComplexF64
    #1. Construction of Pauli string
    pauli_operators = ["I", "X", "Y"]
    num_ops = length(pauli_operators)
    # Generate all possible combinations of Pauli strings
    for i in 1:num_ops^na
        pauli_string = ""
        for j in 1:na
            # Map the index to the corresponding Pauli operator
            index = div(i - 1, num_ops^(j - 1)) % num_ops + 1
           # display(index)
            pauli_string=pauli_string*pauli_operators[index]
        end
        #display(pauli_string)
        push!(Str_XYI_Arr, pauli_string)
    end
   # display(Str_XYI_Arr)
    #2. Construction of the Floquet part of the matrix
    NFZ= 2^(na-1)
   Ham_AuxTryAsymr=diagm(r0=>[1.0 for i in  1:(2*NFZ-r0)])+im*zeros(2^na,2^na) # Collection of all A^(r)
    #####################################
    Str_Coeff=0.0+im*0.0
    #3.Generation of Strings with nonzero coefficients
    for j in 1:3^na
       Str_Coeff=tr(Ham_AuxTryAsymr*PauliString_to_Matrix(Str_XYI_Arr[j]))/(2^na)
       if (abs(Str_Coeff)>tol)
        push!(Floq_Coeff_XYI_Arr,Str_Coeff)
        push!(Floq_Str_XYI_Arr,Str_XYI_Arr[j])
      # display(Str_Coeff)
       end
    end 
  #4. Checking the 
CheckHam_AuxTryAsymr= zeros(size(PauliString_to_Matrix(Floq_Str_XYI_Arr[2])))
for j in 1:length(Floq_Coeff_XYI_Arr)
    CheckHam_AuxTryAsymr += Floq_Coeff_XYI_Arr[j]*PauliString_to_Matrix(Floq_Str_XYI_Arr[j]) 
end
#display(CheckHam_AuxTryAsymr)
if all(CheckHam_AuxTryAsymr.==Ham_AuxTryAsymr)
 println("Check:The resconstructed matrix is equal to the initial auxiliary matrix")
else
 println("Check:The resconstructed matrix does not match the initial auxiliary matrix")
end 
    return Floq_Str_XYI_Arr, Floq_Coeff_XYI_Arr#, CheckHam_Aux_Ω0
end
Xtol= 1.0e-11;
Xna=6;
Xr=1;
Str_Coeff_Floq_OffDiag(Xna, Xr, Xtol)[2]
### Symmetric Floquet zones of auxiliary matrices
################
##########
##############
function  Ham_Aux_Sym_Part(na::Int64,r::Int64) 
    r0=r
    #(na::Int64,r::Int64)
     k= (na - Int64(floor(log2(r)))-1)
 ####
     if r == 0
         return [0]
     end
     binary_array = Int[]
     while r > 0
         pushfirst!(binary_array, r % 2)
         r ÷= 2
     end
     ####
     input_string = ""
     ######
     for bit in binary_array
         if bit == 1
             input_string *= "A"
         else
             input_string *= "B"
         end
     end
     input_string= (("B")^k)*input_string
 ######
 replacements = Dict(
         'A' => ["X", "Y"],
         'B' => ["I", "Z"]
     )
 ######
     output_strings = [""]
     for char in input_string
         if char in keys(replacements)
             new_strings = []
             for output_str in output_strings
                 for new_char in replacements[char]
                     push!(new_strings, output_str * new_char)
                 end
             end
             output_strings = new_strings
         else
             for i in 1:length(output_strings)
                 output_strings[i] *= string(char)
             end
         end
     end
 #######
 coeff_arr=ComplexF64[]
     counts = Dict(
         'X' => 0,
         'Y' => 0,
         'Z' => 0,
         'I' => 0
     )
     for j in 1:2^na
         str=output_strings[j]
         counts['X'],counts['Y'],counts['Z'],counts['I']=[0,0,0,0]
     for char in str
         if char in keys(counts)
             counts[char] += 1
         end
     end
     C = (0.5)^counts['X'] * (0.5*im)^counts['Y'] * (-0.5)^counts['Z']*(0.5)^counts['I']
     #display(C)
     #display([counts['X'],counts['Y'],counts['Z'],counts['I']])
     push!(coeff_arr, C)
    end
####
Ham_AuxSymPart= im*zeros(2^na,2^na)
#display(r0)
Ham_AuxSymPart[2^na-r0,2^na]=1.0+im*0.0
#display(Ham_AuxSymPart)
#4. Checking the 
CheckHam_AuxSymPart= zeros(size(PauliString_to_Matrix(output_strings[1])))
for j in 1:length(coeff_arr)
    CheckHam_AuxSymPart += coeff_arr[j]*PauliString_to_Matrix(output_strings[j]) 
end
#display(CheckHam_AuxSymPart)
if all(CheckHam_AuxSymPart.==Ham_AuxSymPart)
 println("Check:The resconstructed matrix is equal to the initial auxiliary matrix")
else
 println("Check:The resconstructed matrix does not match the initial auxiliary matrix")
end 
####
     return output_strings, coeff_arr
 end
 
 #push!(coeff_arr, C)
Ham_Aux_Sym_Part(4, 1)[1]
### All Pauli strings and coeffients (auxiliary hamiltonian) in one function
###In general the following functions contains both Diagonal and offDiagonal. But I will comment out the symmetric and assymmetric 
function All_Aux_Ham(na::Int64,r,Omg::Float64,Aux_Sym_Zones::Bool,Real_Fourier::Bool,tol::Float64)
    Floq_Aux_Str_arr=Vector{String}[] # Corresponding to each Fourier mode 
    Floq_Aux_Coeff_arr=Vector{ComplexF64}[] # this will used to collect strings with only with non-zero coefficients (weights)
   ## For diagonal terms
  # push!(Floq_Aux_Str_arr, Str_Coeff_Floq_Diag(na,Omg, tol)[1])
   #push!(Floq_Aux_Coeff_arr, Str_Coeff_Floq_Diag(na,Omg, tol)[2])
   #display(Floq_Aux_Coeff_arr)
    ## For Off-diagonal terms (r=[1,3,5])
if typeof(r)==Int64
        rVec = fill(r,1)
        else
        rVec = r
end
####### 
## 1A. Asymmetric Floquet zones
 if Aux_Sym_Zones
    for j in 1: length(rVec)
push!(Floq_Aux_Str_arr,vcat(Str_Coeff_Floq_OffDiag(na,rVec[j],tol)[1],Ham_Aux_Sym_Part(na,rVec[j])[1]))
push!(Floq_Aux_Coeff_arr,vcat(Str_Coeff_Floq_OffDiag(na,rVec[j],tol)[2],-Ham_Aux_Sym_Part(na,rVec[j])[2]))#negative sign to subtract one element from the assymmetric one
    end 
 else
    for j in 1: length(rVec)
    push!(Floq_Aux_Str_arr, Str_Coeff_Floq_OffDiag(na,rVec[j],tol)[1])
    push!(Floq_Aux_Coeff_arr, Str_Coeff_Floq_OffDiag(na,rVec[j],tol)[2])
    end
end
 #2. For a real Forier Hamiltonian
if Real_Fourier
    for j in 1:length(rVec)
        indices_to_delete = findall(x -> abs(real(x)) < tol, Floq_Aux_Coeff_arr[j])
        sort!(indices_to_delete, rev=true)  # Sort indices in descending order to delete elements without affecting subsequent indices
        for idx in indices_to_delete
            splice!(Floq_Aux_Str_arr[j], idx)
            splice!(Floq_Aux_Coeff_arr[j], idx)
        end
        Floq_Aux_Coeff_arr[j]= (2.0).*Floq_Aux_Coeff_arr[j]
    end
end 
    return Floq_Aux_Str_arr,Floq_Aux_Coeff_arr ## this will return [A^(0), A^(1),..,A^(r)]
end

############
#############

### Checking different off diagonal auxilliary matrices

############
##################
#######
function Floq_xyz_model(L::Int64,na::Int64,r,Omg::Float64,Jx,Jy,Jz,dJx,dJy,dJz,Aux_Sym_Zones::Bool,Real_Fourier::Bool,tol::Float64, PBCs)
    Floq_Aux_Str_arr=String[] # Corresponding to each Fourier mode 
    Floq_Aux_Coeff_arr=ComplexF64[] # this will used to collect strings with only with non-zero coefficients (weights)
#############    
#### Diagonal terms 
 ##For Diagonal term related to A0
 Str_Coeff_Floq_Aux= Str_Coeff_Floq_Diag(Xna,XOm,Xtol) ## I guess here we have 
    for i=1:length(Str_Coeff_Floq_Aux[1])
        push!(Floq_Aux_Str_arr,Str_Coeff_Floq_Aux[1][i]*("I"^L))
        push!(Floq_Aux_Coeff_arr,Str_Coeff_Floq_Aux[2][i])
    end
 ##For Diagonal term related to H0: I^na*H0 terms
    for site=1:L-1
        xx_str = "I"^(site-1+na)*"XX"*"I"^(L-(site+1))
        yy_str = "I"^(site-1+na)*"YY"*"I"^(L-(site+1))
        zz_str = "I"^(site-1+na)*"ZZ"*"I"^(L-(site+1))
        push!(Floq_Aux_Str_arr,xx_str)
        push!(Floq_Aux_Coeff_arr,Jx)
        push!(Floq_Aux_Str_arr,yy_str)
        push!(Floq_Aux_Coeff_arr,Jy)
        push!(Floq_Aux_Str_arr,zz_str)
        push!(Floq_Aux_Coeff_arr,Jz)
    end
    if PBCs
        xx_str = "I"^(na)*"X"*"I"^(L-2)*"X"
        yy_str =  "I"^(na)*"Y"*"I"^(L-2)*"Y"
        zz_str = "I"^(na)*"Z"*"I"^(L-2)*"Z"
        push!(Floq_Aux_Str_arr,xx_str)
        push!(Floq_Aux_Coeff_arr,Jx)
        push!(Floq_Aux_Str_arr,yy_str)
        push!(Floq_Aux_Coeff_arr,Jy)
        push!(Floq_Aux_Str_arr,zz_str)
        push!(Floq_Aux_Coeff_arr,Jz)
    end
#############
   ##### For Off diagonal terms:    
##########
  Str_Coeff_OffDiag=All_Aux_Ham(na,r,Omg,Aux_Sym_Zones,Real_Fourier,tol)
############
if Real_Fourier
    for j in 1:length(r)
        for FloqJ in 1:length(Str_Coeff_OffDiag[1][j])
            for site=1:(L-1)
            xx_str = Str_Coeff_OffDiag[1][j][FloqJ]*"I"^(site-1)*"XX"*"I"^(L-(site+1))
            yy_str = Str_Coeff_OffDiag[1][j][FloqJ]*"I"^(site-1)*"YY"*"I"^(L-(site+1))
            zz_str = Str_Coeff_OffDiag[1][j][FloqJ]*"I"^(site-1)*"ZZ"*"I"^(L-(site+1))
            push!(Floq_Aux_Str_arr,xx_str)
            push!(Floq_Aux_Coeff_arr, Str_Coeff_OffDiag[2][j][FloqJ]*dJx[j])
            push!(Floq_Aux_Str_arr,yy_str)
            push!(Floq_Aux_Coeff_arr, Str_Coeff_OffDiag[2][j][FloqJ]*dJy[j])
            push!(Floq_Aux_Str_arr,zz_str)
            push!(Floq_Aux_Coeff_arr, Str_Coeff_OffDiag[2][j][FloqJ]*dJz[j])
            end
            if PBCs
            xx_str = Str_Coeff_OffDiag[1][j][FloqJ]*"I"^(na)*"X"*"I"^(L-2)*"X"
            yy_str = Str_Coeff_OffDiag[1][j][FloqJ]* "I"^(na)*"Y"*"I"^(L-2)*"Y"
            zz_str = Str_Coeff_OffDiag[1][j][FloqJ]*"I"^(na)*"Z"*"I"^(L-2)*"Z"
            push!(Floq_Aux_Str_arr,xx_str)
            push!(Floq_Aux_Coeff_arr,Str_Coeff_OffDiag[2][j][FloqJ]*dJx[j])
            push!(Floq_Aux_Str_arr,yy_str)
            push!(Floq_Aux_Coeff_arr,Str_Coeff_OffDiag[2][j][FloqJ]*dJy[j])
            push!(Floq_Aux_Str_arr,zz_str)
            push!(Floq_Aux_Coeff_arr,Str_Coeff_OffDiag[2][j][FloqJ]*dJz[j])
            end
        end  
    end        
   ###### When 
  #=   # Need attention: I have not implemented the complex coefficients in the following 
else
    for j in 1:length(r)
        for FloqJ in 1:length(Str_Coeff_OffDiag[1][j])
            for site=1:(L-1)
            xx_str = Str_Coeff_OffDiag[1][j][FloqJ]*"I"^(site-1)*"XX"*"I"^(L-(site+1))
            yy_str = Str_Coeff_OffDiag[1][j][FloqJ]*"I"^(site-1)*"YY"*"I"^(L-(site+1))
            zz_str = Str_Coeff_OffDiag[1][j][FloqJ]*"I"^(site-1)*"ZZ"*"I"^(L-(site+1))
            push!(Floq_Aux_Str_arr,xx_str)
            push!(Floq_Aux_Coeff_arr, Str_Coeff_OffDiag[2][j][FloqJ]*dJx[j])
            push!(Floq_Aux_Str_arr,xx_str)
            push!(Floq_Aux_Coeff_arr, Str_Coeff_OffDiag[2][j][FloqJ]*conj(dJx[j]))
            push!(Floq_Aux_Str_arr,yy_str)
            push!(Floq_Aux_Coeff_arr, Str_Coeff_OffDiag[2][j][FloqJ]*dJy[j])
            push!(Floq_Aux_Str_arr,yy_str)
            push!(Floq_Aux_Coeff_arr, Str_Coeff_OffDiag[2][j][FloqJ]*conj(dJy[j]))
            push!(Floq_Aux_Str_arr,zz_str)
            push!(Floq_Aux_Coeff_arr, Str_Coeff_OffDiag[2][j][FloqJ]*dJz[j])
            push!(Floq_Aux_Str_arr,zz_str)
            push!(Floq_Aux_Coeff_arr, Str_Coeff_OffDiag[2][j][FloqJ]*conj(dJz[j]))
            end
            if PBCs
            xx_str = Str_Coeff_OffDiag[1][j][FloqJ]*"I"^(na)*"X"*"I"^(L-2)*"X"
            yy_str = Str_Coeff_OffDiag[1][j][FloqJ]* "I"^(na)*"Y"*"I"^(L-2)*"Y"
            zz_str = Str_Coeff_OffDiag[1][j][FloqJ]*"I"^(na)*"Z"*"I"^(L-2)*"Z"
            push!(Floq_Aux_Str_arr,xx_str)
            push!(Floq_Aux_Coeff_arr,Str_Coeff_OffDiag[2][j][FloqJ]*dJx[j])
            push!(Floq_Aux_Str_arr,xx_str)
            push!(Floq_Aux_Coeff_arr,Str_Coeff_OffDiag[2][j][FloqJ]*conj(dJx[j]))
            push!(Floq_Aux_Str_arr,yy_str)
            push!(Floq_Aux_Coeff_arr,Str_Coeff_OffDiag[2][j][FloqJ]*dJy[j])
            push!(Floq_Aux_Str_arr,yy_str)
            push!(Floq_Aux_Coeff_arr,Str_Coeff_OffDiag[2][j][FloqJ]*conj(dJy[j]))
            push!(Floq_Aux_Str_arr,zz_str)
            push!(Floq_Aux_Coeff_arr,Str_Coeff_OffDiag[2][j][FloqJ]*dJz[j])
            push!(Floq_Aux_Str_arr,zz_str)
            push!(Floq_Aux_Coeff_arr,Str_Coeff_OffDiag[2][j][FloqJ]*conj(dJz[j]))
            end
        end  
    end       =#
end 
        ###
        ############ removing strings with zero weight
    for j in 1:length(r)
        indices_to_delete = findall(x -> abs(x) < tol, Floq_Aux_Coeff_arr)
        sort!(indices_to_delete, rev=true)  # Sort indices in descending order to delete elements without affecting subsequent indices
        for idx in indices_to_delete
            splice!(Floq_Aux_Str_arr, idx)
            splice!(Floq_Aux_Coeff_arr, idx)
        end
    end
        ###### 
    return Floq_Aux_Str_arr, Floq_Aux_Coeff_arr
end
XL=6
Xna=4
Xr=[1]
XOmg=6.0
Xtol=1.0e-8
(Jx0,Jy0,Jz0,dJx0,dJy0,dJz0)= (0.0,0.0,2.0,[1.0],[0.0],[0.0])
str_coeff1=Floq_xyz_model(XL,Xna,Xr,XOmg,Jx0,Jy0,Jz0,dJx0,dJy0,dJz0,true,true,Xtol,false)
str_coeff1[1]
#=## Various type of two spins interaction Hamiltonian 
function xxz_model(L,Jxy,Jz,PBCs)
    return xyz_model(L,Jxy,Jxy,Jz,PBCs)
end
#####
function heisenberg_model(L,J,PBCs)
    return xyz_model(L,J,J,J,PBCs)
end=#
##### Floquet Hamiltonian corresponding to magnetic field
######## Here the driven part of the magnetic field has been considered to be uniform/ homogeneous. We can make it inhomogeneous.
function Floq_Magnetic_Field_model(L::Int64,na::Int64,r,Omg::Float64,Bx,By,Bz,dBx,dBy,dBz,Aux_Sym_Zones::Bool,Real_Fourier::Bool,tol::Float64)
    Floq_Mag_Aux_Str_arr=String[] # Corresponding to each Fourier mode 
    Floq_Mag_Aux_Coeff_arr=ComplexF64[] #this will used to collect strings with only with non-zero coefficients (weights)
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
#############    
   #### Diagonal terms 
############
###For Diagonal term related to B0: I^na* terms
for site=1:L
    Bx_str = "I"^(site-1+na)*"X"*"I"^(L-site)
    By_str = "I"^(site-1+na)*"Y"*"I"^(L-site)
    Bz_str = "I"^(site-1+na)*"Z"*"I"^(L-site)
    if abs(Bxvec[site]) > 1.0e-10
        push!(Floq_Mag_Aux_Str_arr,Bx_str)
        push!(Floq_Mag_Aux_Coeff_arr,Bxvec[site])
    end
    if abs(Byvec[site]) > 1.0e-10
        push!(Floq_Mag_Aux_Str_arr,By_str)
        push!(Floq_Mag_Aux_Coeff_arr,Byvec[site])
    end
    if abs(Bzvec[site]) > 1.0e-10
        push!(Floq_Mag_Aux_Str_arr,Bz_str)
        push!(Floq_Mag_Aux_Coeff_arr,Bzvec[site])
    end
end
#############
  ##### For Off diagonal terms:    
##########
 Str_Coeff_OffDiag=All_Aux_Ham(na,r,Omg,Aux_Sym_Zones,Real_Fourier,tol)
############
if Real_Fourier
    for j in 1:length(r)
        for FloqJ in 1:length(Str_Coeff_OffDiag[1][j])
            for site=1:L
                Bx_str = Str_Coeff_OffDiag[1][j][FloqJ]*"I"^(site-1)*"X"*"I"^(L-site)
                By_str = Str_Coeff_OffDiag[1][j][FloqJ]*"I"^(site-1)*"Y"*"I"^(L-site)
                Bz_str = Str_Coeff_OffDiag[1][j][FloqJ]*"I"^(site-1)*"Z"*"I"^(L-site)
                push!(Floq_Mag_Aux_Str_arr,Bx_str)
                push!(Floq_Mag_Aux_Coeff_arr, Str_Coeff_OffDiag[2][j][FloqJ]*dBx[j])
                push!(Floq_Mag_Aux_Str_arr,By_str)
                push!(Floq_Mag_Aux_Coeff_arr, Str_Coeff_OffDiag[2][j][FloqJ]*dBy[j])
                push!(Floq_Mag_Aux_Str_arr,Bz_str)
                push!(Floq_Mag_Aux_Coeff_arr, Str_Coeff_OffDiag[2][j][FloqJ]*dBz[j])
           end
       end  
   end        
 ###### Complex part
 #=######
else
    for j in 1:length(r)
        for FloqJ in 1:length(Str_Coeff_OffDiag[1][j])
            for site=1:L
                Bx_str = Str_Coeff_OffDiag[1][j][FloqJ]*"I"^(site-1)*"X"*"I"^(L-site)
                By_str = Str_Coeff_OffDiag[1][j][FloqJ]*"I"^(site-1)*"Y"*"I"^(L-site)
                Bz_str = Str_Coeff_OffDiag[1][j][FloqJ]*"I"^(site-1)*"Z"*"I"^(L-site)
                push!(Floq_Mag_Aux_Str_arr,Bx_str)
                push!(Floq_Mag_Aux_Coeff_arr, Str_Coeff_OffDiag[2][j][FloqJ]*dBx[j])
                push!(Floq_Mag_Aux_Str_arr,Bx_str)
                push!(Floq_Mag_Aux_Coeff_arr, Str_Coeff_OffDiag[2][j][FloqJ]*conj(dBx[j]))
                push!(Floq_Mag_Aux_Str_arr,By_str)
                push!(Floq_Mag_Aux_Coeff_arr, Str_Coeff_OffDiag[2][j][FloqJ]*dBy[j])
                push!(Floq_Mag_Aux_Str_arr,By_str)
                push!(Floq_Mag_Aux_Coeff_arr, Str_Coeff_OffDiag[2][j][FloqJ]*conj(dBy[j]))
                push!(Floq_Mag_Aux_Str_arr,Bz_str)
                push!(Floq_Mag_Aux_Coeff_arr, Str_Coeff_OffDiag[2][j][FloqJ]*dBz[j])
                push!(Floq_Mag_Aux_Str_arr,Bz_str)
                push!(Floq_Mag_Aux_Coeff_arr, Str_Coeff_OffDiag[2][j][FloqJ]*conj(dBz[j]))
            end
        end  
    end  =#  
end
    ############ removing strings with zero weight
    for j in 1:length(r)
        indices_to_delete = findall(x -> abs(x) < tol, Floq_Mag_Aux_Coeff_arr)
        sort!(indices_to_delete, rev=true)  # Sort indices in descending order to delete elements without affecting subsequent indices
        for idx in indices_to_delete
            splice!(Floq_Mag_Aux_Str_arr, idx)
            splice!(Floq_Mag_Aux_Coeff_arr, idx)
        end
    end
 return Floq_Mag_Aux_Str_arr, Floq_Mag_Aux_Coeff_arr
end
XL=6
Xna=4
Xr=[1]
XOmg=6.0
Xtol=1.0e-10
(Bx0,By0,Bz0,dBx0,dBy0,dBz0)= (2.0,0.0,0.0,[1.0],[0.0],[0.0])
#str_coeff1=Floq_xyz_model(XL,Xna,Xr,XOmg,Jx0,Jy0,Jz0,dJx0,dJy0,dJz0,true,true,Xtol,false)
Mag_str_coeff1=Floq_Magnetic_Field_model(XL,Xna,Xr,XOmg,Bx0,By0,Bz0,dBx0,dBy0,dBz0,true,true,Xtol)
Mag_str_coeff1[1]
###### Matrix construction and eigenvalues
Xr=[1]
XOmg=6.0
Xtol=1.0e-10
(Jx0,Jy0,Jz0,dJx0,dJy0,dJz0)= (2.0,3.0,0.0,[1.0],[0.0],[0.0])
(Bx0,By0,Bz0,dBx0,dBy0,dBz0)= (1.0,2.0,2.0,[0.3],[0.2],[0.3])
#############
function ED_Floq_Ham_XYZ_B(L::Int64,na::Int64,tol::Float64)
       Ext_Floq_Ham_mat=Matrix{ComplexF64}
       #CheckHam_AuxSymPart= zeros(size(PauliString_to_Matrix(output_strings[1])))
       Total_Floq_Aux_Str_arr=String[] # Corresponding to each Fourier mode 
       Total_Floq_Aux_Coeff_arr=ComplexF64[] #this will used to collect strings with only with non-zero coefficients (weights)
        ####
        Floq_XYZ_Str_Coeff=Floq_xyz_model(L,na,Xr,XOmg,Jx0,Jy0,Jz0,dJx0,dJy0,dJz0,true,true,tol,false)
        Floq_Mag_Str_Coeff=Floq_Magnetic_Field_model(L,na,Xr,XOmg,Bx0,By0,Bz0,dBx0,dBy0,dBz0,true,true,tol)
        #######
       Total_Floq_Aux_Str_arr= vcat(Floq_XYZ_Str_Coeff[1],Floq_Mag_Str_Coeff[1])
       Total_Floq_Aux_Coeff_arr= vcat(Floq_XYZ_Str_Coeff[2],Floq_Mag_Str_Coeff[2])
       #Ext_Floq_Ham_mat=im*zeros(2^(L+na), 2^(L+na))
       Ext_Floq_Ham_mat= zeros(size(PauliString_to_Matrix(Total_Floq_Aux_Str_arr[1])))
       for j in 1:length(Total_Floq_Aux_Coeff_arr)
       Ext_Floq_Ham_mat += Total_Floq_Aux_Coeff_arr[j].*PauliString_to_Matrix(Total_Floq_Aux_Str_arr[j]) #+ Floq_Aux_Coeff_arr[j]*Floq_Aux_Str_arr[j]
       end
       eigsys= eigen(Ext_Floq_Ham_mat)
       Eig_Values= eigsys.values
       Eig_Vectors= eigsys.vectors
       Eig_Vectors= Eig_Vectors[:, sortperm(real(Eig_Values))] # This sort the eigenvectors in terms of quasienergy values. ##
     #check hermiticity:
       Ext_Floq_Ham_mat_Sqr= zeros(size(PauliString_to_Matrix(Total_Floq_Aux_Str_arr[1])))
       Ext_Floq_Ham_mat_Sqr=Ext_Floq_Ham_mat*Ext_Floq_Ham_mat
       eigsys_Sqr= eigen(Ext_Floq_Ham_mat_Sqr)
       Eig_Values_Sqr= eigsys_Sqr.values
       Eig_Vectors_Sqr= eigsys_Sqr.vectors
       Eig_Vectors_Sqr= Eig_Vectors_Sqr[:, sortperm(real(Eig_Values_Sqr))] # This sort the eigenvectors in terms of quasienergy values. ##
###########
    #HFSquare= (HamF_tot)*(HamF_tot)
    #eigsys= eigen(HFSquare)
    #Eig_Values= eigsys.values
    #Eig_Vectors= eigsys.vectors
    #Eig_Vectors= Eig_Vectors[:, sortperm(real(Eig_Values))] # This sort the eigenvectors in terms of 
    #####
    EigVal_Expt=[Adjoint(Eig_Vectors_Sqr[:,j])*Ext_Floq_Ham_mat*Eig_Vectors_Sqr[:,j] for j in 1:2^L]
    Eig_Vectors_Sqr_Expt=[Eig_Vectors_Sqr[:,j] for j in 1:2^L]
 return Eig_Values,Eig_Vectors,Eig_Values_Sqr,Eig_Vectors_Sqr, EigVal_Expt,Eig_Vectors_Sqr_Expt#Ext_Floq_Ham_mat, Ext_Floq_Ham_mat_Squred   
end
