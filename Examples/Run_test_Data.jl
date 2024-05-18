using Floq_ED_v0  # Replace "MyPackage" with the name of the package you want to use
#######
############

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

print(data)
# Define the file path to store the data
#file_path = "/Users/abhishek22/Documents/Python Codes/Basic_Julia_Example/My Julia Packages/data.txt" 
 # You can change the file name and path as needed

# Open the file in write mode
#open(file_path, "w") do file
    # Write the data to the file
 #   write(file, data)

