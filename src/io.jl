#!/usr/bin/env julia

function write_xyz(filename, atoms::Atoms)
    fh = open(filename, "w")
    
    N = length(atoms)

    write(fh, "$N\n")
    write(fh, "$atoms\n")

    for i = 1:N
        symbol = element_properties[atoms.elements[i]].symbol
        pos = atoms.positions[i,:]
        write(fh, "$symbol $(pos[1]) $(pos[2]) $(pos[3])\n")
    end

    close(fh)
end

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

    close(fh)
    
    Atoms(elements, positions)
end

