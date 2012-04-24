#!/usr/bin/env julia
function read_xyz(filename)
    fh = open(filename, "r")

    natoms = int(readline(fh))
    positions = zeros(natoms, 3)
    elements = zeros(Int, natoms)

    comment = chomp(readline(fh))

    for (line,i) in enumerate(each_line(fh))
        fields = split(line)
        elements[i] = element_properties[fields[1]].number
        positions[i,:] = map(float, fields[2:])
    end
    
    Atoms(elements, positions)
end

#let
#    load("lennard-jones.jl")
#    atoms = read_xyz("auh.xyz")
#    println(atoms)
#    atoms.calculator = lennard_jones_calculator(1.0, 1.0)
#    U = potential_energy(atoms)
#    println(U)
#end
