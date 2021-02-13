using Polymake
using Nemo
using AbstractAlgebra

include("parser.jl")
include("other_functions.jl")
include("myeval.jl")

#################################################################################################################################################
### Function to calculate matrix intersection
### Input: parsed matrix
### Output: intersected matrix
#################################################################################################################################################

function intersection_calc(parsed_matrix)
    row, col = size(parsed_matrix)
    
    #we create the canonical matrix of col dimension
    canon_matrix = zeros(col,col)
    for i in 1:col
       canon_matrix[i, i] = 1
    end

    #we create a zeros matrix that has twice the rows of the input matrix to add negation matrix.
    input_matrix = zeros(row*2, col)
    
    #now we fill the matrix with the parsed matrix and its negation matrix
    for i in 1:row
        for j in 1:col
            input_matrix[i, j] = parsed_matrix[i, j]
            input_matrix[i + row, j] = - parsed_matrix[i, j]
        end
    end

    #compute the cones
    matrix_cone = polytope.Cone(INPUT_RAYS=input_matrix)
    orthant_cone = polytope.Cone(INPUT_RAYS=canon_matrix)

    intersect_cone = polytope.intersection(matrix_cone, orthant_cone)
    intersect_matrix = intersect_cone.RAYS

    #Now turn all rational elements into integers
    
    intersect_matrix = Array{Rational{Int}}(intersect_matrix)

    for row in eachrow(intersect_matrix)
        denominator = 1
        lcm_ = 1
        for i in row
            if typeof(i) == Rational{Int64}
                #For some reason if i call denominator(i) I get ERROR: MethodError: objects of type Int64 are not callable even though i is Rational(Int64)
                #so I convert it to string then find the denominator with split. I will look into this in future project.
                i = string(i)
                denominator = parse(Int,split(i, "//")[2])
                lcm_ = lcm(lcm_, denominator)
            end
        end
        row .= row * lcm_
    end
    intersect_matrix = Array{Int}(intersect_matrix)
    #print(intersect_matrix)
    #R = FlintIntegerRing()
    #M = MatrixSpace(R, size(parsed_matrix)[1], size(parsed_matrix)[2])
    #parsed_matrix = Array{Int64,2}(parsed_matrix)
    #parsed_m = M(parsed_matrix::Array{Int64,2})
    #trans = transpose(parsed_m)
    #inver = inv(trans * parsed_m)
    #B,d = pseudo_inv(parsed_m)
    #trans_matrix = intersect_matrix * inv_matrix
    #println(inv)

    if size(intersect_matrix)[1] != polytope.dim(matrix_cone)
        error("Suitable matrix not found")
    else
        return intersect_matrix
    end
    
end



#################################################################################################################################################
### Function to find change of basis matrix 
### Input: original macro-variables text (to be changed later)
### Output: change of basis matrix and its inverse.
#################################################################################################################################################

function cob_calc(old_txt)
    # establish original and new macro-variable matrices (tb changed)
    parsed_matrix = Array{Int}(parse_matrix(old_txt))
    intersect_matrix = Array{Int}(intersection_calc(parsed_matrix))
    intersect_matrix = merge_sort_aux(intersect_matrix)

    # build rational ring and matrix space
    S = MatrixSpace(Nemo.QQ, size(parsed_matrix)[1], size(parsed_matrix)[2])

    # construct fmpq_poly matrices
    parsed_matrix = S(parsed_matrix)
    intersect_matrix = S(intersect_matrix)

    #compute change of basis matrix and its inverse
    cob_matrix = solve_left(intersect_matrix, parsed_matrix)
    cob_matrix_inverse = inv(cob_matrix)

    #print(cob_matrix_inverse)
    open(split(old_txt, ".txt")[1]*"_out.txt", "w") do io
        for line in 1:size(intersect_matrix)[1]
            for i in 1:size(intersect_matrix)[2]
                print(io, intersect_matrix[line,i])
                print(io, " ")
            end
            print(io, "\n")
        end
    end
    return cob_matrix, cob_matrix_inverse

end



#################################################################################################################################################
### Function to find new polynomials
### Input: original macro-variables text, polynomial text (to be changed later)
### Output: new polynomial system
#################################################################################################################################################

function poly_calc(old_txt, old_poly_txt)

    # first, we get the transition matrices, the array of old polynomials (without 0), and the counter for number of 0s
    cob_matrix, cob_matrix_inverse = cob_calc(old_txt)
    old_poly_system, counter = parse_polynomial(old_poly_txt)


    
    #We initialize the Ring again, and construct the dictionary of "variable" => variable
    variables_str = Array{String}([])
    for index in 0:size(cob_matrix)[1] - 1
        push!(variables_str, "y$index")
    end
    
    R, y = PolynomialRing(Nemo.QQ, variables_str)
    S = MatrixSpace(Nemo.QQ, size(cob_matrix)[1], size(cob_matrix)[1])


    #We first compute A^-1 y
    A_inv_y = Array{fmpq_mpoly}([])
    for i in 1:size(cob_matrix_inverse)[1]
        f = sum([cob_matrix_inverse[i, j] * y[j] for j in 1:size(cob_matrix_inverse)[2]])
        push!(A_inv_y, f)
    end

    
    #Now, we evaluate A^-1 y in f(y1,...yn) and we add back the 0 terms in the end. We need to 
    #do this because evaluate() works strictly for Array{fmpq_mpoly}.
    f_y = Array{Any}(map(p -> evaluate(p, A_inv_y), old_poly_system))

    for i in 0:counter
        push!(f_y, fmpq(0))
    end

    #print(f_y)
    #Finally, we multiply from the left by A
    new_poly = []
    for i in 1:size(cob_matrix)[1]
        y_prime = sum([cob_matrix[i, j] * f_y[j] for j in 1:size(cob_matrix)[2]])
        print(f_y[i])
        print('\n')
        print("=====================================")
        print('\n')
        print(y_prime)
        print('\n')
        print("+++++++++++++++++++++++++++++++++++++")
        print('\n')
        push!(new_poly, y_prime)
    end
   

    # printer
    open(split(old_poly_txt, ".txt")[1]*"_new.txt", "w") do io
        for i in 1:size(new_poly)[1]
            print(io, new_poly[i])
            print(io, '\n')
        end
    end

    #cob_matrix = Array{Int}(cob_matrix)
    open(split(old_txt, ".txt")[1]*"_trans.txt", "w") do io
        for line in 1:size(cob_matrix)[1]
            for i in 1:size(cob_matrix)[2]
                print(io, cob_matrix[line,i])
                print(io, " ")
            end
            print(io, "\n")
        end
    end
end






#################################################################################################################################################
#Utility functions#
#################################################################################################################################################
function printer(intersect_matrix,txt)
    open(split(txt, ".txt")[1]*"_out.txt", "w") do io
        for line in eachrow(intersect_matrix)
            for i in line
                print(io, i)
                print(io, " ")
            end
            print(io, "\n")
        end
    end
    print("\n")
    print("File " * split(txt, ".txt")[1] * "_out.txt" * " created!")
end


#To run the script inside julia with > include("poly.jl")   >run(arg)
function run(txt)
    parsed_matrix = parse_matrix(txt)
    intersect_matrix = intersection_calc(parsed_matrix)
    printer(intersect_matrix,txt)
end

#To run the script in terminal with commandline >julia poly.jl arg1
function main()
    parsed_matrix = parse_matrix(ARGS[1])
    intersect_matrix = intersection_calc(parsed_matrix)
    printer(intersect_matrix,ARGS[1])
end



if !isdefined(Base, :active_repl)
    main()
end

#################################################################################################################################################
#################################################################################################################################################
