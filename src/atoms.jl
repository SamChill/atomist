type Atoms
    elements::Array{Int,1}
    positions::Array{Float64, 2}
    box::Array{Float64}
    pbc::Bool
    calculator::Calculator
end

function ==(a::Atoms, b::Atoms)
    if length(a) != length(b)
        return false
    elseif any(a.elements != b.elements)
        return false
    elseif any(a.positions != b.positions)
        return false
    elseif any(a.box != b.box)
        return false
    end

    return true
end

function copy(atoms::Atoms)
    Atoms(copy(atoms.elements), copy(atoms.positions), copy(atoms.box),
          copy(atoms.pbc), atoms.calculator)
end

function show(atoms::Atoms)
    count = HashTable()
    for z in atoms.elements
        count[z] = get(count, z, 0) + 1
    end

    formula = ""
    for (z,n) in count
        symbol = element_properties[z].symbol
        formula = strcat(formula, "$(symbol)$(n)")
    end
    print("Atoms: $formula")
end

function no_calculator()
    function f()
        error("no calculator was defined")
    end
    Calculator(f,f)
end

function Atoms(elements, positions)
    Atoms(elements, positions, zeros(1,3), false, no_calculator())
end

function Atoms(positions)
    elements = ones(Int, size(positions)[1])
    Atoms(elements, positions)
end

function Atoms()
    Atoms(Array(Float64, 1, 3))
end

function length(atoms::Atoms)
    length(atoms.elements)
end

function get_potential_energy(atoms::Atoms)
    atoms.calculator.get_potential_energy(atoms)
end

function get_forces(atoms::Atoms)
    atoms.calculator.get_forces(atoms)
end

function get_distance(atoms::Atoms, i::Int, j::Int)
    r  = atoms.positions[i,1:] - atoms.positions[j,1:]
    if atoms.pbc
        r = pbc(r, atoms.box)
    end
    sqrt(sum(r.^2))
end

function pbc(r::Array{Float64}, box::Array{Float64})
    r - int(r ./ box) .* box
end
