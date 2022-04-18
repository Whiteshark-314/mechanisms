"""
    getcondvects_n_k(n::Int, k=1, val_vector=0:n-1)
    getcondvects_n_k(n::Vector{Int}, k::Int=length(n), val_vector::Union{Tuple, Vector}=missing)

Produces condition vectors. When n=2 and it produces a k-input truth table.


# Examples
julia> getcondvects_n_k(2, 3)
8×3 Matrix{Any}:
 0  0  0
 0  0  1
 0  1  0
 0  1  1
 1  0  0
 1  0  1
 1  1  0
 1  1  1
julia> getcondvects_n_k(3, 3, ['x', 'y', 'z'])
27×3 Matrix{Any}:
 'x'  'x'  'x'
 'x'  'x'  'y'
 'x'  'x'  'z'
 'x'  'y'  'x'
 'x'  'y'  'y'
 'x'  'y'  'z'
 'x'  'z'  'x'
 ⋮         
 'z'  'y'  'x'
 'z'  'y'  'y'
 'z'  'y'  'z'
 'z'  'z'  'x'
 'z'  'z'  'y'
 'z'  'z'  'z'
julia> getcondvects_n_k([1,2,3], 3, [["x"], ['a','b'], [1, 2, 3]])
6×3 Matrix{Any}:
 "x"  'a'  1
 "x"  'a'  2
 "x"  'a'  3
 "x"  'b'  1
 "x"  'b'  2
 "x"  'b'  3

"""

function getcondvects_n_k(n::Int, k=1, val_vector=0:n-1)
	rows = n^k
	condvects=Matrix{Any}(missing, rows, k);
    for i = 1:k
        @views condvects[:, end+1-i] = reshape(repeat(permutedims(val_vector),
		Int(rows/n^(k-i+1)), Int(rows/n^i)), rows, 1);
    end
	return condvects
end

function getcondvects_n_k(n::Vector{Int}, k::Int=length(n),
	val_vector::Union{Tuple, Vector}=missing)
	@assert k == length(n) "n should have k entries"
	rows = prod(n)
	condvects=Matrix{Any}(missing, rows, k);
	if ismissing(val_vector)
		for i = 1:k
			@views condvects[:, i] = reshape(repeat(permutedims(1:n[i]),
			Int(rows/prod(n[1:i])), Int(rows/prod(n[i:end]))), rows, 1);
		end
	else
		@assert prod(length.(val_vector) .== n) "Length of each entry in val_vector should equal "
		for i = 1:k
	        @views condvects[:, i] = reshape(repeat(permutedims(val_vector[i]),
			Int(rows/prod(n[1:i])), Int(rows/prod(n[i:end]))), rows, 1);
		end
	end
	return condvects
end