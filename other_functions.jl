
#################################################################################################################################################
### Merge Sort algorithm to show new variables in reduced echelon form ###
#################################################################################################################################################
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

#################################################################################################################################################
#################################################################################################################################################









#################################################################################################################################################
#Unused/Not yet used functions
#################################################################################################################################################

#A function to convert matrix to rref
function rref_matrix(matrix)
    S = MatrixSpace(ZZ,size(matrix)[1],size(matrix)[2])
    mat = Array{Int}(matrix)
    flint_mat = S(mat)
    return rref(flint_mat)
end

macro assert(ex)
    return :( $ex ? nothing : throw(AssertionError($(string(ex)))) )
end

#A function to convert between fmpz_mat and Array{Int} formats
function Base.convert(::Type{Array{Int}}, x::Nemo.fmpz_poly)
    m,n = size(x)
    mat = Rational{Int}[x[i,j] for i = 1:m, j = 1:n]
    return mat
end

Base.convert(::Type{Array}, x::Nemo.fmpz_poly) = convert(Array{Rationall{Int}}, x)


function polynomial_calc(txt)
    #old_poly = parse_polynomial(txt1)
    parsed_matrix = parse_matrix(txt)
    intersect_matrix = intersection_calc(parsed_matrix)
    #We first make new matrix into rref form
    ref_intersect_matrix = merge_sort_aux(intersect_matrix)
    trans = solve_left(ref_intersect_matrix, parsed_matrix)
    print(trans)
    # Transtion matrix
    transition_matrix = find_cob_matrix(parsed_matrix, ref_intersect_matrix)
    inv_transition_matrix = inv(transition_matrix)
    inv_A_Y = Array{Int}(inv_transition_matrix * ref_intersect_matrix)

    open(split(txt, ".txt")[1]*"_AinvY.txt", "w") do io
        for line in eachrow(inv_A_Y)
            for i in line
                print(io, i)
                print(io, " ")
            end
            print(io, "\n")
        end
    end

    verifier_matrix = Array{Int}(transition_matrix*ref_intersect_matrix)
    if verifier_matrix == parsed_matrix
        print("Transition matrix good")
    end

    
    open(split(txt, ".txt")[1]*"_trans.txt", "w") do io
        for line in eachrow(ref_intersect_matrix)
            for i in line
                print(io, i)
                print(io, " ")
            end
            print(io, "\n")
        end
    end

end

#Calculate transition matrix
#Takes two macro-variables, produce transition matrix
function find_cob_matrix(matrix1, matrix2)
    cob_matrix = Array{Rational{Int64}}(zeros(size(matrix1)[1], size(matrix1)[1]))
    list = []
    count = [0]
    m = Array{Rational{Int64}}(zeros(size(matrix1)[1], size(matrix1)[2]))
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
            sum = 0
            if matrix1[last(count)+1, col_index] != 0
                push!(count,last(count)+1)
                for row_index in 1:size(matrix1)[1]
                    if matrix2[row_index, col_index] != 0 && row_index != last(count)
                        m[row_index, col_index] = cob_matrix[l, row_index] * matrix2[row_index, col_index]
                        sum += m[row_index, col_index]
                    end
                end
                cob_matrix[l, last(count)] = (matrix1[l, col_index] - sum) // matrix2[last(count), col_index]
            end
        end
    end
    return cob_matrix
end


#################################################################################################################################################
#################################################################################################################################################