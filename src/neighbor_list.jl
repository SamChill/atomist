#!/usr/bin/env julia

type NeighborList
    update::Function
    distances::Array
    neighbors::Array
    coordination::Array{Int}
end

function neighbor_list(cutoff::Real)
    neighbors = Array[]
    distances = Array[]
    coordination = Int[]
    atoms_old = Atoms()

    #This is a brute force method O(N^2). It could be improved with 
    #the use of verlet and cell lists.
    function update(atoms::Atoms)
        if atoms_old == atoms
            return
        else
            atoms_old = copy(atoms)
        end

        N = length(atoms)
        for i = 1:N
            push(neighbors, Array[])
            push(distances, Float64[])
            push(coordination, 0)
            for j = 1:N
                if i == j
                    continue
                end
                r = get_distance(atoms, i, j)
                if r < cutoff
                    push(neighbors[i], [i,j])
                    push(distances[i], r)
                    coordination[i] += 1
                end
            end
        end
    end
    NeighborList(update, distances, neighbors, coordination)
end
