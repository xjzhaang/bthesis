using Polymake
using Nemo
using AbstractAlgebra

include("parser.jl")
include("other_functions.jl")
include("myeval.jl")

#################################################################################################################################################
### Function to calculate matrix intersection
### Input: parsed matrix of type Array{Int}
### Output: intersected matrix of type Array{Int}
###         if no suitable matrix, throw DimensionMismatch 
#################################################################################################################################################

function intersection_calc(parsed_matrix::Array{Int})
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


    if size(intersect_matrix)[1] < polytope.dim(matrix_cone)
        throw(DimensionMismatch("No suitable matrix found"))
    elseif size(intersect_matrix)[1] > polytope.dim(matrix_cone)
        intersect_matrix = find_best_basis(intersect_matrix)
        intersect_matrix = merge_sort_aux(intersect_matrix)
    end
    print(size(intersect_matrix)[1])
    print('\n')
    print(polytope.dim(matrix_cone))
    return intersect_matrix
end



#################################################################################################################################################
### Function to find change of basis matrix 
### Input: old macrovariable matrix A, new macrovariable matrix B of type Array{Int}
### Output: change of basis matrix M such that A = M * B and M^(-1) of type fmpq_mat
#################################################################################################################################################

function cob_calc(parsed_matrix::Array{Int}, intersect_matrix::Array{Int})
    # build rational ring and matrix space
    S = MatrixSpace(Nemo.QQ, size(parsed_matrix)...)

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
### Input: cob_matrix, cob_matrix_inverse, old_poly_system 
### Output: new polynomial system of type Array{fmpq_mpoly}
### We use the equation y' = M f( M^(-1) y) where M, M^(-1) are the change of basis matrices
#################################################################################################################################################

function poly_calc(cob_matrix::fmpq_mat, cob_matrix_inverse::fmpq_mat, old_poly_system::Array{fmpq_mpoly})
    
    # We initialize the Ring
    variables_str = ["y$index" for index in 0:size(cob_matrix)[1] - 1]
    
    R, y = PolynomialRing(Nemo.QQ, variables_str)
    S = MatrixSpace(Nemo.QQ, size(cob_matrix)...)

    # We first compute M^(-1) y
    A_inv_y = Array{fmpq_mpoly}([])
    for i in 1:size(cob_matrix_inverse)[1]
        f = sum([cob_matrix_inverse[i, j] * y[j] for j in 1:size(cob_matrix_inverse)[2]])
        push!(A_inv_y, f)
    end

    # Now, we evaluate M^(-1) y in f(y1,...yn)
    f_y = Array{Any}(map(p -> evaluate(p, A_inv_y), old_poly_system))

    #Finally, we multiply from the left by M
    new_poly = []
    for i in 1:nrows(cob_matrix)
        y_prime = sum([cob_matrix[i, j] * f_y[j] for j in 1:size(cob_matrix)[2]])
        push!(new_poly, y_prime)
    end

    return new_poly
end


#################################################################################################################################################
### Function to find the best basis from linearly dependent rows using Kruskal's greedy algorithm
### Input: matrix of rank < dimension
### Output: matrix of rank = dimension
#################################################################################################################################################
function find_best_basis(matrix::Array{Int})
    #first we construct a list of edges where each edge is a tuple (weight, row1)
    col = size(matrix)[2]
    edges_list = []
    for i in 1:size(matrix)[1]
        weight = find_nonzero(matrix[i, :])
        push!(edges_list, (weight, i))
    end

    #the edges are sorted in increasing order of weight
    sort!(edges_list, by = x -> x[1])

    #we pop the first edge and construct the return matrix to be [row1, row2]
    return_matrix = reshape(matrix[edges_list[1][2], :], 1, col)
    visited_row = []
    popfirst!(edges_list)

    #for every edge in list of edges, we add to the return_matrix and see if it is linearly independent, if not, we remove it
    #after everyloop, we check if the rank of return_matrix is the same as the input matrix
    for edge_index in 1:size(edges_list)[1]
        push!(visited_row, edges_list[edge_index][2])
        S = MatrixSpace(Nemo.ZZ, size(return_matrix)...)
        old_rank = rank(S(return_matrix))
        return_matrix = vcat(return_matrix, vcat(reshape(matrix[edges_list[edge_index][2], :], 1, col)))
        S = MatrixSpace(Nemo.ZZ, size(return_matrix)...)
        new_rank = rank(S(return_matrix))
        if old_rank != new_rank - 1
            return_matrix = return_matrix[1:size(return_matrix)[1] - 1, :]
        end
        S = MatrixSpace(Nemo.ZZ, size(return_matrix)...)
        S1 = MatrixSpace(Nemo.ZZ, size(matrix)...)
        if rank(S(return_matrix)) == rank(S1(matrix))
            break
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

function fceri_macro_printer(intersect_matrix, fceri_varnames, txt)
    variables_dict = parse_varnames(fceri_varnames)
    new_names = Dict()
    open(split(txt, ".txt")[1] * "_macrovariables.txt", "w") do io
        for row_index in 0:size(intersect_matrix)[1] - 1
            str = "y" * "$row_index = "
            for col_index in 0:size(intersect_matrix)[2] - 1
                if intersect_matrix[row_index + 1, col_index + 1] == 0
                    continue
                elseif intersect_matrix[row_index + 1, col_index + 1] == 1
                    str = str * variables_dict[col_index] * " + "
                else
                    str = str * string(intersect_matrix[row_index + 1, col_index + 1])
                    str = str * variables_dict[col_index] * " + "
                end
            
            end
            str = str[1: length(str) - 3]
            if !occursin(" + ", str)
                new_names["y" * "$row_index"] = split(str, "= ")[2] 
            end
            print(io, str)
            print(io, "\n")
        end
    end
    return new_names
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

###########################################################################################
#Functions to run within julia
###########################################################################################
function get_new_matrix(matrix_txt, is_fceri)
    parsed_matrix = parse_matrix(matrix_txt)
    intersect_matrix = intersection_calc(parsed_matrix)
    new_matrix_printer(intersect_matrix, matrix_txt)
    if is_fceri
        new_names = fceri_macro_printer(intersect_matrix, "fceri/fceri_varnames.txt", matrix_txt)
    end
end


function get_new_poly(matrix_txt, poly_txt)
    parsed_matrix = parse_matrix(matrix_txt)
    intersect_matrix = intersection_calc(parsed_matrix)
    cob_matrix, cob_matrix_inverse = cob_calc(parsed_matrix, intersect_matrix)

    old_poly_system= parse_polynomial(poly_txt)

    new_poly = poly_calc(cob_matrix, cob_matrix_inverse, old_poly_system)
    new_poly_printer(new_poly, poly_txt)
end

function run_all_PP(n)
    for i in 2:n
        get_new_matrix("PP/$i" * "m.txt", false)
        get_new_poly("PP/$i" * "m.txt", "PP/$i" * "p.txt")
    end
end

function run_all_fceri(n)
    for i in 2:n
        get_new_matrix("fceri/$i" * "m.txt", true)
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


#PP_2 = ['S0', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10', 'S11', 'S12', 'S13', 'S14', 'S15', 'S16', 'S17', 'kCatE', 'kCatF', 'kOffE', 'kOffF', 'kOnE', 'kOnF']

