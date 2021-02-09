using Polymake
using LinearAlgebra
using Nemo

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
    println("positive orthant dimension: \n", polytope.dim(orthant_cone))
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


#We parse the txt file into matrix
function parse_matrix(txt)
    f = open(txt)
    lines = readlines(f)
    #first we create a zeros matrix of the correct dimension.
    matrix = zeros(length(lines), length(split(lines[1])))
    matrix = Array{Rational{Int64}}(matrix)
    for line_index in 1:length(lines)
        vector = split(lines[line_index])
        for v_index in 1:length(vector)
            max_deminominator = 1
            # first we check if the element is a fraction, if so we find its numerator and denominator and create a rational to represent it
            if occursin("/", vector[v_index])
                numerator = parse(Int,split(vector[v_index], "/")[1])
                denominator = parse(Int,split(vector[v_index], "/")[2])
                if denominator > max_deminominator
                    max_deminominator = denominator
                end
                rational = numerator // denominator
                matrix[line_index, v_index] = rational
            #if the element is an int, we parse the string to int and add it to the matrix
            else
                matrix[line_index, v_index] = parse(Int64, vector[v_index])
            end
        end
    end
    return matrix
end

#We parse the original polynomials into a dictionary of Expr
function parse_polynomial(txt)
    f = open(txt)
    lines = readlines(f)
    old_poly = Dict()
    for line_index in 1:length(lines)
        old_poly[line_index] = Meta.parse(lines[line_index])
    end
    return old_poly
end

#A function to convert matrix to rref
function rref_matrix(matrix)
    S = MatrixSpace(ZZ,size(matrix)[1],size(matrix)[2])
    mat = Array{Int}(matrix)
    flint_mat = S(mat)
    return rref(flint_mat)
end





#A function to convert between fmpz_mat and Array{Int} formats
function Base.convert(::Type{Array{Int}}, x::Nemo.fmpz_mat)
    m,n = size(x)
    mat = Int[x[i,j] for i = 1:m, j = 1:n]
    return mat
end

Base.convert(::Type{Array}, x::Nemo.fmpz_mat) = convert(Array{Int}, x)

#Main function to find new polynomial 
function polynomial_calc(txt2)
    #old_poly = parse_polynomial(txt1)
    parsed_matrix = parse_matrix(txt2)
    intersect_matrix = intersection_calc(parsed_matrix)
    #We first make new matrix into rref form
    ref_intersect_matrix = merge_sort_aux(intersect_matrix)


    # Transtion matrix
    transition_matrix = Array{Int}(cob_matrix_c(parsed_matrix))
 
    open(split(txt2, ".txt")[1]*"_bruh1.txt", "w") do io
        for line in eachrow(a)
            for i in line
                print(io, i)
                print(io, " ")
            end
            print(io, "\n")
        end
    end


end

function cob_matrix_c(matrix)
    cob_matrix = Array{Rational{Int64}}(zeros(size(matrix)[1], size(matrix)[1]))
    for row_index in 1:size(matrix)[1]
        if matrix[row_index, :] == matrix[row_index, :]
            cob_matrix[row_index,row_index] = 1
        end
    end
    for col_index in 1:size(matrix)[1]
        cob_matrix[3, col_index] = matrix[col_index, 3]
    end
    return cob_matrix
end

#################################################
#useless function
function find_cob_matrix(matrix1, matrix2)
    cob_matrix = Array{Rational{Int64}}(zeros(size(matrix1)[1], size(matrix1)[1]))
    list = []
    count = []
    c = 0
    m = Array{Int64}(zeros(size(matrix1)[1], size(matrix1)[2]))

    #if two rows are the same, the change of basis matrix' row is the identity
    for row_index in 1:size(matrix1)[1]
        if matrix1[row_index, :] == matrix2[row_index, :]
            cob_matrix[row_index,row_index] = 1
        else
            push!(list,row_index)
        end
    end
    for l in list
        for col_index in 1:size(matrix1)[2]
            for row_index in 1:size(matrix1)[1]
                if row_index in count
                    break
                else
                    if matrix2[row_index, col_index] != 0 && row_index != c
                        m[row_index, col_index] = cob_matrix[l, col_index] * matrix2[row_index, col_index]
                        push!(count,row_index)
                    end
                end
            end
            if c in count
                break
            else
                c += 1
                total1 = sum([m[row_index, col_index] for row_index in 1:size(matrix1)[1]])
                print(matrix2[c, col_index])
                cob_matrix[l, last(count)] = (matrix1[l, col_index] - total1) // matrix2[c, col_index]
            end
        end
    end
    return cob_matrix
end
###################################################################################

#To make intersected matrix into ref 
function merge_sort_aux(matrix)
    if size(matrix)[1] == 1
        return matrix
    end
    matrix1 = matrix[1:Int(floor(size(matrix)[1] / 2)), :]
    matrix2 = matrix[size(matrix1)[1] + 1 : size(matrix)[1], :]

    matrix1 = merge_sort_aux(matrix1)
    matrix2 = merge_sort_aux(matrix2)
    
    return merge_sort(matrix1,matrix2)
end
function merge_sort(matrix1,matrix2)
    res_matrix = Array{Int64,2}(undef,0,size(matrix1)[2])
    #print(res_matrix)
    while (size(matrix1)[1] != 0 && size(matrix2)[1] != 0) 
        if (findfirst(x -> x != 0, matrix1)[2] < findfirst(x -> x != 0, matrix2)[2])
            if size(matrix1)[1] > 1
                res_matrix = [res_matrix ; reshape(matrix1[1, :],(1,size(matrix1)[2]))]
                matrix1 = matrix1[setdiff(1:end, 1), :]
                #print(matrix1)
            else
                res_matrix = [res_matrix ; matrix1]
                matrix1 = matrix1[setdiff(1:end, 1), :] 
            end
        else
            if size(matrix2)[1] > 1
                res_matrix = [res_matrix ; reshape(matrix2[1, :],(1,size(matrix2)[2]))]
                matrix2 = matrix2[setdiff(1:end, 1), :]
            else
                res_matrix = [res_matrix ; matrix2]
                matrix2 = matrix2[setdiff(1:end, 1), :]
            end
        end
    end
    while (size(matrix1)[1] != 0)
        if size(matrix1)[1] > 1
            res_matrix = [res_matrix ; reshape(matrix1[1, :],(1,size(matrix1)[2]))]
            matrix1 = matrix1[setdiff(1:end, 1), :]
        else
            res_matrix = [res_matrix ; matrix1]
            matrix1 = matrix1[setdiff(1:end, 1), :] 
        end
    end
    while (size(matrix2)[1] != 0) 
        if size(matrix2)[1] > 1
            res_matrix = [res_matrix ; reshape(matrix2[1, :],(1,size(matrix2)[2]))]
            matrix2 = matrix2[setdiff(1:end, 1), :]
        else
            res_matrix = [res_matrix ; matrix2]
            matrix2 = matrix2[setdiff(1:end, 1), :]
        end
    end
    return res_matrix
end

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


