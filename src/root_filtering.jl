# stuff for singular roots, if these appear
using HomotopyContinuation
using LinearAlgebra

function filter_singular(system, raw)
    J_mat = jacobian(system)

    sings = singular.(getindex.(raw, 1))
    !isempty(flatten(sings)) ? @info("Singular solutions found, beware!") : nothing

    total = []

    for res in raw # for each parameter set
        trials = solutions(res[1], only_nonsingular=true)
        p_vals = res[2]
        this_J(x_vals) = evaluate(J_mat, vcat(system.variables, system.parameters) => vcat(x_vals, p_vals))

        for path in filter( x -> x.singular && x.return_code == :success, path_results(res[1])) # go over singular paths found here
            sing_idx = singular_variables(this_J, solution(path)) # find spurious index
            
            # try and set spurious variable to zero, see if still a solution
            trial_soln = solution(path)
            trial_soln[sing_idx] = 0.0 
            e = norm(evaluate(system.expressions, vcat(system.variables, system.parameters) => vcat(trial_soln, p_vals)))

            if abs(e) < 3*path.accuracy 
                append!(trials, [trial_soln])
            end
        end
        append!(total, [Vector{Vector{ComplexF64}}(unique_points(trials))])
        #display.(Vector{Vector{ComplexF64}}(unique_points(trials)))
    end
    Vector{Vector{Vector{ComplexF64}}}(total)
end


function singular_variables(J, x_vals)
    eigs = eigen(J(x_vals))
    spur = findall( x -> abs(x) < 1e-20, eigs.values)
    length(spur) > 1 ? error("Singular solutions in more than 1 variable") : nothing
    max_elem = argmax(abs.(eigs.vectors[:,spur]))[1]
    max_elem
end


function _is_equivalent_root(r1, r2; spur_idx=[])
    diff = (r1-r2)[1:end .∉ spur_idx]
end

