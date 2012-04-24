#!/usr/bin/env julia
liblammps = dlopen("./liblammps")
libc = dlopen("libc")

function lammps_open()
    lammps_handle = Array(Uint, 1)
    #args = "lammps -log none -echo log -screen none"
    args = "lammps -log none -screen none"
    #args = "lammps"
    args = split(args)
    arr = Array(ByteString, length(args))
    for i = 1:length(args)
        arr[i] = cstring(args[i])
    end
    argv = Array(Ptr{Uint8}, length(arr)+1)
    for i = 1:length(arr)
        argv[i] = arr[i].data
    end
    argv[length(args)+1] = C_NULL
    ccall(dlsym(liblammps, :lammps_open_no_mpi), Ptr{Void}, 
          (Int, Ptr{Ptr{Uint8}}, Ptr{Uint}), length(argv)-1, argv, 
          lammps_handle)

    return lammps_handle[1]
end

function lammps_command(lammps_handle, cmd)
    println(cmd)
    ccall(dlsym(liblammps, :lammps_command), Ptr{Void}, (Uint, Ptr{Uint8}), 
          lammps_handle, cstring(cmd))
end

function lammps_put_coords(lammps_handle, positions)
    ccall(dlsym(liblammps, :lammps_put_coords), Ptr{Void}, 
          (Uint, Ptr{Float64}), lammps_handle, positions[:])
end

function lammps_extract_variable(lammps_handle, variable)
    lammps_extract_variable(lammps_handle, variable, "", 1)
end

function lammps_extract_variable(lammps_handle, variable, group, n)
    var = ccall(dlsym(liblammps, :lammps_extract_variable), Ptr{Void},
                (Uint, Ptr{Uint8}, Ptr{Uint8}), lammps_handle, 
                variable, group)
    var_copy = Array(Float64, n)
    ccall(dlsym(libc, :memcpy), Ptr{Void}, (Ptr{Void}, Ptr{Void}, Int),
          var_copy, var, n*sizeof(Float64))
    ccall(dlsym(libc, :free), Ptr{Void}, (Ptr{Void},), var)

    if n == 1
        return var_copy[1]
    else
        return var_copy
    end
end

function create_lammps(atoms::Atoms)
    handle = lammps_open()
    cmd(s) = lammps_command(handle, s)

    ntypes = 0
    atom_types = HashTable()
    j = 1
    for i = 1:length(atoms)
        z = atoms.elements[i]
        if !has(atom_types, z)  
            atom_types[z] = j
            j += 1
        end
    end
    ntypes = length(atom_types)

    #Gives units in Angstoms and eV
    cmd("units metal")
    cmd("atom_style	atomic")

    #Preserves atomic index ordering
    cmd("atom_modify map array sort 0 0")

    #Define periodic cell
    if atoms.pbc
        lx = ly = lz = 0
        mx = atoms.box[1]
        my = atoms.box[2]
        mz = atoms.box[3]
    else
        lx = min(atoms.positions[:,1])
        ly = min(atoms.positions[:,2])
        lz = min(atoms.positions[:,3])
        mx = max(atoms.positions[:,1])
        my = max(atoms.positions[:,2])
        mz = max(atoms.positions[:,3])
        cmd("boundary s s s")
    end

    cmd("region box block $lx $mx $ly $my $lz $mz units box")
    cmd("create_box $ntypes box")

    #Initialize the atoms and their types
    for i = 1:length(atoms)
        atom_type = atom_types[atoms.elements[i]]
        x = atoms.positions[i,1]
        y = atoms.positions[i,2]
        z = atoms.positions[i,3]
        cmd("create_atoms $atom_type single $x $y $z units box")
    end

    #We don't care about mass but have to set it
    cmd("mass * 1.0")

    #Define variables for force and energy so they can be extracted
    cmd("variable fx atom fx")
    cmd("variable fy atom fy")
    cmd("variable fz atom fz")
    cmd("variable pe equal pe")

    #lammps_put_coords(handle, atoms.positions)
    #cmd("run 0")
    #U = lammps_extract_variable(handle, "pe")
    #println("Potential Energy: $U")
    #forces = Array(Float64, length(atoms), 3)
    #for (var,i) in enumerate([ "fx", "fy", "fz" ])
    #    forces[:,i] = lammps_extract_variable(handle, var, "all", 
    #                                                length(atoms))
    #end
    #println("Forces: $forces")

    return handle
end

function pt_morse_calculator()
    prev_atoms = Atoms()
    energy = 0.0
    forces = zeros(Float64, 1, 3)
    lammps_handle = 0

    function calculate(atoms::Atoms)
        #if (atoms == prev_atoms)
        #    println("IDENTICAL")
        #    return
        #end

        cmd(s) = lammps_command(lammps_handle, s)
        if length(atoms) != length(prev_atoms)
            forces = Array(Float64, length(atoms), 3)
            #println("CREATING LAMMPS")
            lammps_handle = create_lammps(atoms)

            cmd("pair_style morse 9.5")
            cmd("pair_coeff * * 0.7102 1.6047 2.897")
            cmd("pair_modify shift yes")
        end

        lammps_put_coords(lammps_handle, atoms.positions)
        cmd("run 0")
        energy = lammps_extract_variable(lammps_handle, "pe")
        for (var,i) in enumerate([ "fx", "fy", "fz" ])
            forces[:,i] = lammps_extract_variable(lammps_handle, var, "all", 
                                                  length(atoms))
        end

        prev_atoms = copy(atoms)
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
