using Nemo
include("myeval.jl")

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



function parse_polynomial(txt)
    # parse the text file
    f = open(txt)
    lines = readlines(f)
    

    variables_str = Array{String}([])
    for index in 0:length(lines)-1
        push!(variables_str, "y$index")
    end

    #Create the PolynomiaL ring
    QQ = FlintQQ
    R, v = PolynomialRing(QQ, variables_str)
    S = MatrixSpace(R, length(lines), length(lines))

    #Create a dictionary sending symbol to v[i]
    expr_dict = Dict(Symbol("y0") => v[1])
    for index in 2:length(v)
        num = index-1
        expr_dict[Symbol("y$num")] = v[index]
    end
    
    counter = 0
    poly_system = Array{fmpq_mpoly}([])
    for line_index in 1:length(lines)
        #We need to have the Array as type fmpq_mpoly, therefore we need to remove the 0 terms. We add a counter for the Number
        # of 0 terms, as we will add them back when multiplying f(A^-1y) by A
        if myeval(Meta.parse(lines[line_index]), expr_dict) != fmpq(0)
            push!(poly_system, myeval(Meta.parse(lines[line_index]), expr_dict))
        else
            counter += 1
        end
    end
    #print(typeof(poly_system))
    return poly_system, counter
end
