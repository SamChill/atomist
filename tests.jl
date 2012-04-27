#!/usr/bin/env julia
load("src/init.jl")

let
    atoms = read_xyz("examples/lj38.xyz")
    nl = neighbor_list(1.5)
    nl.update(atoms)
    @assert sum(nl.coordination) == 288
end
