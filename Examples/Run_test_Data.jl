using Floq_ED_v0  # Replace "MyPackage" with the name of the package you want to use
using JSON
#######
#### Values of parameters
(XL,XJbarVec,XdJVec,XBx,XBy,XBz,XdBx,XdBy,XdBz,XΩ,XNF_Ham, XNF_Zones)=(2,[2.7,2.7,3.4],[1.0,1.0,1.0],0.0,0.0,0.0,0.0,0.0,0.0,6.0,1,3)
######
XNt_mesh=5000
######
XB0Vec=[XBx,XBy,XBz]
XB1Vec=[XdBx,XdBy,XdBz]
XJ0Vec=XJbarVec
XJ1Vec=XdJVec;

#@time data_HamF_Sqr1= Eig_Floq_Ham_Square(XL,XJbarVec,XdJVec,XBx,XBy,XBz,XdBx,XdBy,XdBz,XΩ,XNF_Ham, XNF_Zones,false)

# Define a function to run from the package
function run_my_Julia_package1(L,JbarVec,dJVec,Bx,By,Bz,dBx,dBy,dBz,Ω,No_Fourier_Ham, No_Floq_Zones,PBCs)
    # Call the function from the package
    result = Floq_ED_v0.Eig_Floq_Ham_Square(L,JbarVec,dJVec,Bx,By,Bz,dBx,dBy,dBz,Ω,No_Fourier_Ham, No_Floq_Zones,PBCs)  # Replace "some_function" with the actual function name
    return result
end

# Call the function
data = run_my_Julia_package1(XL,XJbarVec,XdJVec,XBx,XBy,XBz,XdBx,XdBy,XdBz,XΩ,XNF_Ham, XNF_Zones,false)
#adapt_results = timed_res[1]
   # elapsedtime = timed_res[2]
#print(data)
# Define the file path to store the data
#file_path = "/Users/abhishek22/Documents/Python Codes/Basic_Julia_Example/My Julia Packages/data.txt" 
 # You can change the file name and path as needed

# Convert to JSON and write to file
open("data.json", "w") do io
    JSON.print(io, data)
end
 println("start")
# # Read from JSON file
data = open("data.json", "r") do io
    JSON.parse(IOBuffer(read(io, String)))
end
#########
##########
#  io2 = open("zzFloq_L"*string(L)*"_pool_"*myparams["adapt_pool"]*"_initstate"*myparams["initial_state"]*"_results.dat","w")
#  write(io2,resstring*"\n")
#  close(io2)

#  parstring = JSON.json(myparams)

#     resstring = JSON.json(Dict("energy"=>adapt_results.energy,
#                              "max_grad"=>adapt_results.max_grad,
#                              "max_grad_ind"=>adapt_results.max_grad_ind,
#                              "grads"=>adapt_results.grads,
#                              "opt_pars"=>adapt_results.opt_pars,
#                              "paulis"=>chosen_op_array,
#                              "opt_numevals"=>adapt_results.opt_numevals))

#     io = open("zzFloq_L"*string(L)*"_pool_"*myparams["adapt_pool"]*"_initstate"*myparams["initial_state"]*"_params.dat","w")
#     write(io,parstring*"\n")
#     # write(io,resstring*"\n")
#     close(io)

# # Open the file in write mode
# #open(file_path, "w") do file
#     # Write the data to the file
#  #   write(file, data) =

# #####
######
##println("start")

# read data
# mydata = readdlm(inputfilename,Any,skipstart=2)

# Lvals = mydata[1:end,1]
# energyvals = mydata[1:end,4]

# if targetL in Lvals
#     indL = indexin(targetL,Lvals)[1]
#     run_z_adapt(targetL,inputpars["Jz"],inputpars["Bx"],inputpars["PBCs"],energyvals[indL],inputfilename)
# else
#     println("the desired system size is not present in the input file")
# end

# Example data
#data = [rand(3, 3) for _ in 1:5]

# # Convert to JSON and write to file
# open("data.json", "w") do io
#     JSON.print(io, data)
# end

# # Read from JSON file
# data = open("data.json", "r") do io
#     JSON.parse(IOBuffer(read(io, String)))
# end