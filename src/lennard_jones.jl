function lennard_jones_calculator()
    lennard_jones_calculator(1.0, 1.0)
end

function lennard_jones_calculator(epsilon, sigma)
    energy = 0.0
    forces = zeros(Float64, 1, 3)
    
    function calculate(atoms::Atoms)
        N = length(atoms)
        energy = 0.0
        forces = zeros(Float64, N, 3)

        for i = 1:N
            for j = i+1:N
                r = get_distance(atoms, i, j)
                energy += 4 * epsilon * ( (sigma/r)^12.0 - (sigma/r)^6.0 )
                f = -24 * epsilon * ( 2*(sigma^12.0/r^13.0) - (sigma^6.0/r^7.0) )

                diff = atoms.positions[i,:] - atoms.positions[j,:]
                forces[i,:] -= f*diff/r
                forces[j,:] += f*diff/r
            end
        end
    end

    function get_potential_energy(atoms::Atoms)
        calculate(atoms)
        return energy
    end

    function get_forces(atoms::Atoms)
        calculate(atoms)
        return forces
    end

    return Calculator(get_potential_energy, get_forces)
end
