using Test
include("../intersect.jl")

test_matrix1 = [1 0 0 0; 0 -2 1 0; 0 0 2 1]
test_matrix1_res = [1 0 0 0; 0 4 0 1; 0 0 2 1]
S = MatrixSpace(Nemo.QQ, 3, 3)
test_matrix2 = [-2 0 1; 0 1 0]

R, (y0,y1,y2) = PolynomialRing(Nemo.QQ, ["y0", "y1", "y2"])

@testset "Parsers test" begin
    #simple test
    @test parse_matrix("test files/test_matrix1.txt") == [1 0 0 0; 0 -2 1 0; 0 0 2 1]
    #test if file with fractions become int
    @test parse_matrix("test files/test_fraction.txt") == [1 0 0 0; 0 3 0 7; 0 0 -16 21]

    #test if polynomial parser is good
    @test parse_polynomial("test files/test_poly1.txt") == Array{fmpq_mpoly}([y0*y1 - y2, -y0*y2 + y1, R(0)])
    @test parse_polynomial("test files/test_poly2.txt") == Array{fmpq_mpoly}([y0*y1 - y0, y2*y0*y1, y1*y1 - y2*y2])
    println("\n")
end

@testset "Mergesort matrix test" begin
    #test if mergesorting for turning matrix into triangular form works 
    @test merge_sort_aux([0 -2 1 0; 1 0 0 0; 0 0 2 1]) == test_matrix1
    @test merge_sort_aux([0 0 2 1; 0 4 0 1; 1 0 0 0]) == test_matrix1_res
    println("\n")
end
@testset "Intersection test" begin
    #test type
    @test typeof(intersection_calc(test_matrix1)) == Array{Int64,2}

    #test intersect result
    @test intersection_calc(test_matrix1) == Array{Int64,2}(test_matrix1_res)
    @test_throws DimensionMismatch("No suitable matrix found") intersection_calc(test_matrix2)
    println("\n")
end

@testset "Change of basis test" begin
    cob_matrix, cob_matrix_inverse = cob_calc(test_matrix1, test_matrix1_res)
    #test type
    @test typeof(cob_matrix) == fmpq_mat
    @test typeof(cob_matrix_inverse) == fmpq_mat

    #test cob results
    @test cob_matrix == S([1 0 0; 0 -1//2 1//2; 0 0 1])
    @test cob_matrix_inverse == S([1 0 0; 0 -2 1; 0 0 1])
    println("\n")
end

@testset "Polynomial calc test" begin
    new_poly1 = poly_calc(S([1 0 0; 0 -1//2 1//2; 0 0 1]),  S([1 0 0; 0 -2 1; 0 0 1]), Array{fmpq_mpoly}([y0*y1 - y2, -y0*y2 + y1, R(0)]))
    new_poly2 = poly_calc(S([1 0 0; 0 -1//2 1//2; 0 0 1]),  S([1 0 0; 0 -2 1; 0 0 1]), Array{fmpq_mpoly}([y0*y1 - y0, y2*y0*y1, y1*y1 - y2*y2]))

    #test types
    @test typeof(new_poly1) == Array{Any,1}
    @test typeof(new_poly1[3]) == fmpq_mpoly
    @test typeof(new_poly2) == Array{Any,1}
    @test typeof(new_poly2[3]) == fmpq_mpoly

    #test results
    @test new_poly1 == [-2*y0*y1 + y0*y2 - y2, 1//2*y0*y2 + y1 - 1//2*y2, 0]
    @test new_poly2 == [-2*y0*y1 + y0*y2 - y0, y0*y1*y2 - 1//2*y0*y2^2 + 2*y1^2 - 2*y1*y2, 4*y1^2 - 4*y1*y2]
    println("\n")
end