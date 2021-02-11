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
    

    variables_str = Array{String}([])
    variables = []
    for index in 0:length(lines)-1
        push!(variables_str, "y$index")
    end

    for index in 1:length(variables_str)
        #string_as_varname_function(variables_str[index],nothing)
        push!(variables, Symbol(variables_str[index]))
    end
    variables_tuple = Tuple(variables)

    #Create the PolynomiaL ring
    QQ = FlintQQ
    R, variables_tuple = PolynomialRing(QQ, variables_str)
    S = MatrixSpace(R, length(lines), length(lines))
    for line_index in 1:length(lines)
        a = Meta.parse(lines[line_index])
        print(eval(a))
        #@eval $(Symbol(:f, line_index)) = a
    end
    #print(f1)
end

function string_as_varname_function(s::AbstractString, v::Any)
	s = Symbol(s)
	@eval (($s) = ($v))
end

