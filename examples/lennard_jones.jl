#!/usr/bin/env julia
load("src/init.jl")

atoms = read_xyz("examples/lj38.xyz")
atoms.calculator = lennard_jones_calculator()

U = get_potential_energy(atoms)
println("energy: $U")

atoms.positions += 0.1 * rand(length(atoms), 3)
U = get_potential_energy(atoms)
println("energy: $U")
