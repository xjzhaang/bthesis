
# Bachelor Thesis code

intersect.jl is a program that computes new macro-variables.
Requires polymake, Nemo, AbstractAlgebra.jl packages

## How to Use

Inside Julia terminal:
```sh 
include("intersect.jl")
``` 
To get the new macro-variables in a [n]m_new.txt file of matrices:
```sh 
get_new_matrix("[n]m.txt") 
``` 

To get the new ODE system of polynomials in a [n]p_new.txt file: 
```sh
get_new_poly("[n]p.txt") 
``` 

To get both outputs for 1~i reduced models of case studies:
```sh
run_all_[case_study_name](i)
```
Available case study names: PP, fceri, Barua



