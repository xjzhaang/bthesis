using Polymake
using Nemo
using AbstractAlgebra

include("parser.jl")
include("other_functions.jl")
include("myeval.jl")

#################################################################################################################################################
### Function to calculate matrix intersection
# Gleb: types (Array{Int}??)
### Input: parsed matrix
### Output: intersected matrix
# Gleb: what if no matrix?
#################################################################################################################################################

function intersection_calc(parsed_matrix)
    row, col = size(parsed_matrix)
    
    #we create the canonical matrix of col dimension
    canon_matrix = zeros(col,col)
    for i in 1:col
       canon_matrix[i, i] = 1
    end

    #we create a zeros matrix that has twice the rows of the input matrix to add negation matrix.
    input_matrix = zeros(row * 2, col)
    
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
    intersect_matrix = rational_to_int(intersect_matrix)
    #We now modify the matrix into upper triangular form using merge sort
    intersect_matrix = merge_sort_aux(intersect_matrix)

    if size(intersect_matrix)[1] != polytope.dim(matrix_cone)
        intersect_matrix = find_best_basis(intersect_matrix)
    end
        #throw(DimensionMismatch("No suitable matrix found"))
    print(size(intersect_matrix)[1])
    print("\n")
    print(polytope.dim(matrix_cone))

    return intersect_matrix
end



#################################################################################################################################################
### Function to find change of basis matrix 
### Input: old macrovariable matrix A, new macrovariable matrix B
### Output: change of basis matrix M such that A = M * B (?) and M^(-1)
# Gleb: types
#################################################################################################################################################

function cob_calc(parsed_matrix, intersect_matrix)
    # build rational ring and matrix space
    # Gleb: use MatrixSpace(Nemo.QQ, size(parsed_matrix)...)
    S = MatrixSpace(Nemo.QQ, size(parsed_matrix)[1], size(parsed_matrix)[2])

    # construct fmpq_poly matrices
    parsed_matrix = S(parsed_matrix)
    intersect_matrix = S(intersect_matrix)

    #compute change of basis matrix and its inverse
    cob_matrix = solve_left(intersect_matrix, parsed_matrix)
    cob_matrix_inverse = inv(cob_matrix)

    return cob_matrix, cob_matrix_inverse
end



#################################################################################################################################################
### Function to find new polynomials
### Input: cob_matrix, cob_matrix_inverse, old_poly_system, counter
# Gleb: remove counter by having zero polynomials
### Output: new polynomial system
#
# Gleb: types
#################################################################################################################################################

function poly_calc(cob_matrix, cob_matrix_inverse, old_poly_system, counter)
    
    # We initialize the Ring
    # Gleb: variables_str = ["y$index" for index in ???]
    variables_str = Array{String}([])
    for index in 0:size(cob_matrix)[1] - 1
        push!(variables_str, "y$index")
    end
    
    R, y = PolynomialRing(Nemo.QQ, variables_str)
    # Gleb: again using ...
    S = MatrixSpace(Nemo.QQ, size(cob_matrix)[1], size(cob_matrix)[1])


    # Gleb: what is A?
    # We first compute A^-1 y
    A_inv_y = Array{fmpq_mpoly}([])
    for i in 1:size(cob_matrix_inverse)[1]
        f = sum([cob_matrix_inverse[i, j] * y[j] for j in 1:size(cob_matrix_inverse)[2]])
        push!(A_inv_y, f)
    end

    # Gleb: hopefully, not necessary
    # Now, we evaluate A^-1 y in f(y1,...yn) and we add back the 0 terms in the end. We need to 
    # do this because evaluate() works strictly for Array{fmpq_mpoly}.
    f_y = Array{Any}(map(p -> evaluate(p, A_inv_y), old_poly_system))

    for i in 0:counter
        push!(f_y, fmpq(0))
    end

    #Finally, we multiply from the left by A
    new_poly = []
    # Gleb: use ncols or nrows
    for i in 1:size(cob_matrix)[1]
        y_prime = sum([cob_matrix[i, j] * f_y[j] for j in 1:size(cob_matrix)[2]])
        push!(new_poly, y_prime)
    end

    return new_poly
end


#################################################################################################################################################
### Function to find a best basis from linearly dependent rows
### Input: matrix of rank < dimension
### Output: matrix of rank = dimension
#################################################################################################################################################
#matrix = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 2 1 1 0 1 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 1 1 2 0 0 1 1 0 0 1 0 1 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 1 1 1 1 0 0 0 2 0 1 1 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 1 1 1 0 1 1 1 2 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]
function find_best_basis(matrix)
    basis_matrix = Array{Int64}(undef, 0, size(matrix)[2])
    return_matrix = Array{Int64}(undef, 0, size(matrix)[2])

    #problem: I cant index fmpz matrix, so I only transform Array{Int} to fmpz when calculating rank
    #problem2: MatrixSpace changest after adding another row, so I need to redefine S every time
    for row_index in 1:size(matrix)[1]
        S = MatrixSpace(Nemo.ZZ, size(basis_matrix)[1], size(matrix)[2])
        old_rank = rank(S(basis_matrix))
        basis_matrix = vcat(basis_matrix, reshape(matrix[row_index, :], 1, size(matrix)[2]))
        S1 = MatrixSpace(Nemo.ZZ, size(basis_matrix)[1], size(matrix)[2])
        new_rank = rank(S1(basis_matrix))

        #print((old_rank, new_rank))
        #print("\n")

        if old_rank != new_rank
            return_matrix = vcat(return_matrix, reshape(matrix[row_index, :], 1, size(matrix)[2]))
        end
    end
    return return_matrix
end


#################################################################################################################################################
#Utility functions#
#################################################################################################################################################
function new_matrix_printer(intersect_matrix, txt)
    open(split(txt, ".txt")[1] * "_new.txt", "w") do io
        for row_index in 1:size(intersect_matrix)[1]
            for col_index in 1:size(intersect_matrix)[2]
                print(io, intersect_matrix[row_index, col_index])
                print(io, " ")
            end
            print(io, "\n")
        end
    end
    print("\n")
    print("File " * split(txt, ".txt")[1] * "_new.txt" * " created!")
end

function new_poly_printer(new_poly, txt)
    open(split(txt, ".txt")[1] * "_new.txt", "w") do io
        for row_index in 1:size(new_poly)[1]
            print(io, new_poly[row_index])
            print(io, "\n")
        end
    end
    print("\n")
    print("File " * split(txt, ".txt")[1] * "_new.txt" * " created!")
end


#To run the script inside julia with > include("poly.jl")   >run(arg)
function get_new_matrix(matrix_txt)
    parsed_matrix = parse_matrix(matrix_txt)
    intersect_matrix = intersection_calc(parsed_matrix)
    #print(intersect_matrix)
    new_matrix_printer(intersect_matrix, matrix_txt)
end


function get_new_poly(matrix_txt, poly_txt)
    parsed_matrix = parse_matrix(matrix_txt)
    intersect_matrix = intersection_calc(parsed_matrix)
    cob_matrix, cob_matrix_inverse = cob_calc(parsed_matrix, intersect_matrix)

    old_poly_system, counter = parse_polynomial(poly_txt)

    new_poly = poly_calc(cob_matrix, cob_matrix_inverse, old_poly_system, counter)
    new_poly_printer(new_poly, poly_txt)
end

function run_all_PP(n)
    for i in 2:n
        get_new_matrix("PP/$i" * "m.txt")
        get_new_poly("PP/$i" * "m.txt", "PP/$i" * "p.txt")
    end
end

function run_all_fceri(n)
    for i in 1:n
        get_new_matrix("fceri/$i" * "m.txt")
        get_new_poly("fceri/$i" * "m.txt", "fceri/$i" * "p.txt")
    end
end

#To run the script in terminal with commandline >julia intersect.jl arg1
function main()
    run_all(ARGS[1])
end



if !isdefined(Base, :active_repl)
    main()
end

#################################################################################################################################################
#################################################################################################################################################
