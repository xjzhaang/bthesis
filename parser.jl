using Nemo

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
    

    variables = Array{String}([])
    for line_index in 1:length(lines)
        push!(variables, "y$line_index")
    end
    variables_tuple = Tuple(variables)

    #Create the PolynomiaL ring
    QQ = FlintQQ
    R, variables_tuple = PolynomialRing(QQ, variables)

    for line_index in 1:length(lines)
        @eval $(Symbol(:f, line_index)) = eval(Meta.parse(lines[line_index]))
    end
    print(f1)
end


