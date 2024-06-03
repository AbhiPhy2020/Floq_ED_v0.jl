"""
This file deals with- 
1. Physical Floquet state from the extended Floquet states
2. Expectation value of an operator as a function of time.  
3.SVD, Entanglement entropy between auxialiary and physical space
"""

using LinearAlgebra 


"""
Splitting a vector in extended state into an array of vectors.
"""
function Split_Ext_Vector(x_Vec, L, na)
  NFZ = 2^(na-1) - 1
  n = [2^L for j in 1:(2*NFZ+1)]
  ϕ_res = Vector{Vector{eltype(x_Vec)}}()
  start = firstindex(x_Vec)
  for len in n
    push!(ϕ_res, x_Vec[start:(start + len - 1)])
    start += len
  end
  return ϕ_res
end

#####    Operator Expectation Vlaue
function Op_Vec_Expectation(Op,Vect,NFZ,t)
    ϕ_lowest = sum(Vect[n+NFZ+1]*exp(im*n*XΩ*t) for n in -NFZ:NFZ)  
    Op_expt = Adjoint(ϕ_lowest)*Op*ϕ_lowest
    return Op_expt  #ϕ_lowest
end

function Op_Vec_Expectation_Low(Op,Vect,NFZ,t)
    ϕ_lowest = sum(Vect[n]*exp(im*(n-2*NFZ-1)*XΩ*t) for n in 1:(2*NFZ+1))  
    Op_expt = Adjoint(ϕ_lowest)*Op*ϕ_lowest
    return Op_expt  #ϕ_lowest
end

#### Function for operator expectation values and entanglement
## In this cell we define a single and two spin operators. We also define a function that split  a large extended vector into array of small vectors.
########
function One_Spin_Sum_Operator(L)
  #Op=zeros(ComplexF64,2^L,2^L)
 pdictOp = Dict("SumX" => magnetic_field_matrix(L,1.0,0.0,0.0),"SumY" =>  magnetic_field_matrix(L,0.0,1.0,0.0), "SumZ" =>  magnetic_field_matrix(L,0.0,0.0,1.0))
 return pdictOp
end
########
function Two_Spin_Sum_Operator(L,PBCs)
 pdictOp2 = Dict("SumXX" => xyz_matrix(L,1.0,0.0,0.0,PBCs),"SumYY" => xyz_matrix(L,0.0,1.0,0.0,PBCs), "SumZZ" => xyz_matrix(L,0.0,0.0,1.0,PBCs))
 return pdictOp2
end
#Note: There might be some asymmetric consideration  in -na to na :: Look at the structure of phi_n 

"""
    SVD_Entaglement_Entropy(x_Vec,L,NFZ)

TBW
"""
function SVD_Entaglement_Entropy(x_Vec, L, na)
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

"""
Different level of checks
1. For large frequency, without taking any time periodic term. Start with simple Ising and then take arbitrary XYZ matrix
2. Add time periodic term and then compare the result in the two cases.Start with simple Ising and then take arbitrary XYZ matrix
3. 
"""

#Expectation value of an operator as a function of time.
function Op_Expectation(Op,NFZ, t)
    ϕ_lowest = sum(ϕntilde[n+NFZ+1]*exp(im*n*XΩ*t) for n in -NFZ:NFZ)  
    Op_expt = Adjoint(ϕ_lowest)*Op*ϕ_lowest
    return Op_expt  #ϕ_lowest
end
#=
s3=[1.0 0.0; 0.0 -1.0]
OperatorZ4=kron(s3,s3,s3,s3)
Op_Expectation(OperatorZ4, na, 0.3)
#transpose(conj(ϕntilde[1]))*OperatorZ4*ϕntilde[1] =#
XL
###
#There might be some asymmetric consideration  in -na to na :: Look at the structure of phi_n 
###
Ω0=XΩ
OperatorZ4=(kron(I2,I2,I2,s3)+kron(I2,I2,s3,I2)+kron(I2,s3,I2,I2)+kron(s3,I2,I2,I2))
OperatorX4=(kron(I2,I2,I2,s1)+kron(I2,I2,s1,I2)+kron(I2,s1,I2,I2)+kron(s1,I2,I2,I2))
OperatorY4=(kron(I2,I2,I2,s2)+kron(I2,I2,s2,I2)+kron(I2,s2,I2,I2)+kron(s2,I2,I2,I2))


OperatorZZ4=(kron(I2,I2,s3,s3)+kron(I2,s3,I2,s3)+kron(s3,I2,I2,s3))
OperatorXX4=(kron(I2,I2,s1,s1)+kron(I2,s1,I2,s1)+kron(s1,I2,I2,s1))
OperatorYY4=(kron(I2,I2,s2,s2)+kron(I2,s2,I2,s2)+kron(s2,I2,I2,s2))

to_plot3=[Op_Expectation(OperatorZ4, na, (2π/Ω0)*t) for t in 0.0:0.01:2.0]
to_plot2=[Op_Expectation(OperatorY4, na, (2π/Ω0)*t) for t in 0.0:0.01:2.0]
to_plot1=[Op_Expectation(OperatorX4, na, (2π/Ω0)*t) for t in 0.0:0.01:2.0]
######
######
Mesht = [t for t in 0:0.01:2.0]
scalefontsizes(1.1)
labels = ["Expectation_Value"]
#plot(Mesht, to_plot1, xlabel="Time", ylabel="Expectation_Value", xticks = Mesht[1]:0.1:Mesht[end], label = labels)
#savefig("Expectation_Value_vs_(t/T)"*string(sites)*".pdf")

to_plot3zz=[Op_Expectation(OperatorZZ4, na, (2π/Ω0)*t) for t in 0.0:0.01:2.0];
to_plot2yy=[Op_Expectation(OperatorYY4, na, (2π/Ω0)*t) for t in 0.0:0.01:2.0];
to_plot1xx=[Op_Expectation(OperatorZZ4, na, (2π/Ω0)*t) for t in 0.0:0.01:2.0];

OperatorZ4=kron(sx,sx)
Op_Expectation(Op,NFZ, t)
#### Function for operator expectation values and entanglement
## In this cell we define a single and two spin operators. We also define a function that split  a large extended vector into array of small vectors.
########
function One_Spin_Operator(L)
    #Op=zeros(ComplexF64,2^L,2^L)
   pdictOp = Dict("SumX" => magnetic_field_matrix(L,1.0,0.0,0.0),"SumY" =>  magnetic_field_matrix(L,0.0,1.0,0.0), "SumZ" =>  magnetic_field_matrix(L,0.0,0.0,1.0))
   return pdictOp
end
########
function Two_Spin_Operator(L,PBCs)
   pdictOp2 = Dict("SumXX" => xyz_matrix(L,1.0,0.0,0.0,PBCs),"SumYY" => xyz_matrix(L,0.0,1.0,0.0,PBCs), "SumZZ" => xyz_matrix(L,0.0,0.0,1.0,PBCs))
   return pdictOp2
end


#Note: There might be some asymmetric consideration  in -na to na :: Look at the structure of phi_n 
#####
# Defining the spin orthonormal basis:

