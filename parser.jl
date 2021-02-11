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
    for index in 0:length(lines)-1
        push!(variables_str, "y$index")
    end
    #var_tuple = Tuple(variables_str)
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
    
    poly_system = []
    for line_index in 1:length(lines)
        push!(poly_system, myeval(Meta.parse(lines[line_index]), expr_dict))
    end
    return poly_system
end

function myeval(e::Union{Expr,Symbol,Number}, map::Dict{Symbol,fmpq_mpoly})
    try
        return f(e, map)
    catch ex
        println("Can't parse \"$e\"")
        rethrow(ex)
    end
end 

function f(s::Symbol, map::Dict{Symbol,fmpq_mpoly})
    if haskey(map, s)
        return map[s]
    else
        throw(UndefVarError(s))
    end
end    

# Numbers are converted to type Float64.
function f(x::Number, map::Dict{Symbol,fmpq_mpoly})
    return Int(x)
end    

# To parse an expression, convert the head to a singleton
# type, so that Julia can dispatch on that type.
function f(e::Expr, map::Dict{Symbol,fmpq_mpoly})
    return f(Val(e.head), e.args, map)
end

# Call the function named in args[1]
function f(::Val{:call}, args, map::Dict{Symbol,fmpq_mpoly})
    return f(Val(args[1]), args[2:end], map)
end

# Addition
function f(::Val{:+}, args, map::Dict{Symbol,fmpq_mpoly})
    x = 0
    for arg ∈ args
        x += f(arg, map)
    end
    return x
end

# Subtraction and negation
function f(::Val{:-}, args, map::Dict{Symbol,fmpq_mpoly})
    len = length(args)
    if len == 1
        return -f(args[1], map)
    else
        return f(args[1], map) - f(args[2], map)
    end
end    

# Multiplication
function f(::Val{:*}, args, map::Dict{Symbol,fmpq_mpoly})
    x = 1
    for arg ∈ args
        x *= f(arg, map)
    end
    return x
end    

# Division
function f(::Val{:/}, args, map::Dict{Symbol,fmpq_mpoly})
    return f(args[1], map) / f(args[2], map)
end    