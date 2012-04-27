#!/usr/bin/env julia

type NeighborList
    update::Function
    distances::Array
    neighbors::Array
end

function neighbor_list(cutoff::Real)
    neighbors = Array[]
    distances = Array[]
    atoms_old = Atoms()
    function update(atoms::Atoms)
        if atoms_old == atoms
            return
        end
        N = length(atoms)
        for i = 1:N
            push(neighbors, Array[])
            push(distances, Float64[])
            for j = 1:N
                if i == j
                    continue
                end
                r = get_distance(atoms, i, j)
                if r < cutoff
                    push(neighbors[i], [i,j])
                    push(distances[i], r)
                end
            end
        end
    end
    NeighborList(update, distances, neighbors)
end
