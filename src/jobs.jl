
function enrico_2stage_job(dir::String, str::Structure, state;
                           ecutwfc = 60,
                           ecutrho = 480,
                           mixing_beta = 0.1,
                           kpoints = (4, 4, 4, 0, 0, 0),
                           server = "daint",
                           pseudos = "oncvpsp_sr",
                           exec = "pw-occupations-10nodes")
    e = get_exec(exec; server = server)
    scf1 = Calculation{QE}("scf_1", :calculation => "scf",
                           :mixing_beta => 1,
                           :mixing_fixed_ns => 1,
                           :electron_maxstep => 1;
                           exec = e)
    scf2 = Calculation{QE}("scf", :calculation => "scf",
                           :mixing_beta => mixing_beta,
                           :electron_maxstep => 1000; exec = e)
    set_kpoints!(scf1, kpoints)
    set_kpoints!(scf2, kpoints)

    starting_ns = generate_starting_ns_eigenvalue(state)
    j = Job(splitpath(dir)[end], str, [scf1, scf2],
            :ecutwfc => ecutwfc,
            :ecutrho => ecutrho,
            :starting_ns_eigenvalue => starting_ns,
            :conv_thr => 1e-9;
            dir = dir,
            server = server,
            environment = "normal_10nodes")
    set_pseudos!(j, pseudos; server = "localhost")
    return j
end

function starting_ns_job(dir::String, str::Structure, state;
                         ecutwfc = 60,
                         ecutrho = 480,
                         mixing_beta = 0.1,
                         kpoints = (4, 4, 4, 0, 0, 0),
                         server = "daint",
                         pseudos = "oncvpsp_sr",
                         exec = "pw-occupations-10nodes")
    e = get_exec(exec; server = server)
    scf2 = Calculation{QE}("scf", :calculation => "scf",
                           :mixing_beta => mixing_beta,
                           :Hubbard_maxstep => 100,
                           :Hubbard_mixing_beta => 0.4,
                           :Hubbard_conv_thr => 0.1,
                           :electron_maxstep => 1000; exec = e)
    set_kpoints!(scf2, kpoints)

    starting_ns = generate_starting_ns_eigenvalue(state)
    j = Job(splitpath(dir)[end], str, [scf1, scf2],
            :ecutwfc => ecutwfc,
            :ecutrho => ecutrho,
            :starting_ns_eigenvalue => starting_ns,
            :conv_thr => 1e-9;
            dir = dir,
            server = server,
            environment = "normal_10nodes")
    set_pseudos!(j, pseudos; server = "localhost")
    return j
end

function set_Hubbard!(j, occ;
                      conv_thr = 0.1,
                      maxstep = 100,
                      mix = 0.4,
                      type = "energy",
                      strength = 1.0)
    if ndims(occ) == 4
        j[:Hubbard_occupations] = occ
    else
        j[:starting_ns_eigenvalue] = occ
    end
    j[:Hubbard_conv_thr] = conv_thr
    j[:Hubbard_maxstep] = maxstep
    j[:Hubbard_mixing_beta] = mix
    j[:Hubbard_constraint_type] = type
    j[:Hubbard_strength] = strength
    return j
end

delete_Hubbard!(j::Job) = delete_hubbard!.(j.calculations)
function delete_Hubbard!(c)
    delete!(c, :Hubbard_mixing_beta)
    delete!(c, :Hubbard_conv_thr)
    delete!(c, :Hubbard_occupations)
    delete!(c, :Hubbard_maxstep)
    delete!(c, :Hubbard_constraint_type)
    delete!(c, :Hubbard_strength)
    return c
end


function local_save(job::Job, local_dir::String)
    jdir = job.dir
    env = job.environment
    s = job.server

    job.dir = local_dir
    job.server = local_server().name
    job.environment = "default"
    suppress() do
        # We do save different versions for saving potential outputs that were already there
        # while the job is resubmitted
        return save(job; fillexecs = false)
    end
    job.dir = jdir
    job.environment = env
    job.server = s
    return job
end

function local_load(job::Job)
    job.dir = abspath(job.dir)
    return load(local_server(), job)
end
gencalc(job, settings) = throw(MethodError(gencalc, (job, settings)))
