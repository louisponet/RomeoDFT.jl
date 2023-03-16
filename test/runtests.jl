using Test
using RomeoDFT
using Pkg.Artifacts
using RomeoDFT: RemoteHPC
using RomeoDFT.RemoteHPC: configure_local, julia_main

using RomeoDFT: ProgressMeter

# Sets up the local server
tconfdir = tempname()
if ispath(tconfdir)
    rm(tconfdir; recursive = true)
end
import RomeoDFT.RemoteHPC: config_path
config_path(p...) = joinpath(tconfdir, p...)

paths = ["jobs",
         "logs/jobs",
         "logs/runtimes",
         "storage/servers",
         "storage/execs",
         "storage/environments"]
         
for p in paths
    mkpath(config_path(p))
end

redirect_stdin(devnull) do
    redirect_stderr(devnull) do
        redirect_stdout(devnull) do
            s = configure_local(; interactive = false)
            s.jobdir = tempname()
            RemoteHPC.SLEEP_TIME[] = 0.01
            !ispath(s.jobdir) && mkdir(s.jobdir)
            s.max_concurrent_jobs = 10
            save(s)
            return tsk = @async julia_main(sleep_time=0.01)
        end
    end
end

while !isalive(local_server())
    sleep(0.1)
end

const TESTSERVER = local_server()

env = Environment(name="default")
save(TESTSERVER, env)
# Adding the fake quantum espresso binary
function create_mockdata(case::String)
    bindir = tempname()
    mkpath(bindir) 
    
    testdir = joinpath(splitdir(pathof(RomeoDFT))[1], "../test/assets/$case")
    if !ispath(testdir)
        testdir = @artifact_str case
    end
    
    info = RomeoDFT.DFControl.FileIO.qe_parse_calculation(joinpath(testdir, "5/scf.in"))
    str = info.structure
    pw_path = joinpath(bindir, "pw.x")
    open(pw_path, "w") do f
        write(f, """
        #!$(joinpath(Sys.BINDIR, "julia")) --startup-file=no
        lines = readlines(stdin)
        pref_id = findfirst(x -> occursin("prefix", x), lines)
        if pref_id === nothing
            println("error")
        end

        prefix = strip(split(lines[pref_id])[end], ''')
        entity_str = match(r"id_(\\d+)", prefix)[1]
        entity_id = parse(Int, entity_str)

        jobdir = joinpath("$testdir",entity_str)
        calcid = findfirst(x -> occursin("calculation", x), lines)
        if calcid === nothing
            fname = "projwfc"
            files = filter(x -> occursin("pdos", x), readdir(jobdir))
            for f in files
                cp(joinpath(jobdir, f), f; force = true)
            end
        else
            calculation = strip(split(lines[calcid])[end], ''')
            if calculation == "bands"
                fname = "bands"
            elseif calculation == "nscf"
                fname = "nscf"
            elseif calculation == "scf"
                fname = "scf"
            end
        end

        print(read(joinpath(jobdir, fname * ".out"), String))
        """)
    end
    chmod(pw_path, 0o777)
    projwfc_path = joinpath(bindir, "projwfc.x")
    symlink(pw_path, projwfc_path)
    pwexec = Exec(name="pw_$case", path=pw_path)
    save(TESTSERVER, pwexec)
    save(TESTSERVER, Exec(name="projwfc_$case", path=projwfc_path))
    
    pseudo_dir = joinpath(bindir, "pseudos")
    mkpath(pseudo_dir)
    for a in unique(str.atoms)
        touch(joinpath(pseudo_dir, "$(a.element.symbol).UPF"))
    end
    pseudo_set = configure_pseudoset(TESTSERVER, case, pseudo_dir)
    set_pseudos!(str, pseudo_set)
    calc = RomeoDFT.DFControl.Calculation(name="scf", flags=info.flags, data=info.data, exec=pwexec)

    return ServerInfo(TESTSERVER.name, "pw_$case", "default", 5), Template(str, calc), length(readdir(testdir))
end


test_cases = ("Ni2O2",)
for case in test_cases
    @testset "$case" begin
        server_info, template, n_entities = create_mockdata(case)
    
        # ledger_dir = tempname()
        ledger_dir = "/tmp/Romeotest/"
        ispath(ledger_dir) && rm(ledger_dir, recursive=true)
        
        @info "Running at $ledger_dir"
        l = Searcher(rootdir=ledger_dir, sleep_time = 0.01)
        Entity(l, server_info)
        base_e = Entity(l, BaseCase(),
                           deepcopy(template),
                           Generation(1))
        sim_entity = Entity(l, RandomSearcher(10),
                               deepcopy(template),
                               Unique(1e-2, true),
                               IntersectionSearcher(0.25, 100),
                               StopCondition(0.1, 3),
                               Generation(1))
                               
        RomeoDFT.set_mode!(l, :search)
        save(l)
        @test ispath(joinpath(ledger_dir, string(RomeoDFT.DATABASE_VERSION), "ledger.jld2"))
        l = load(l)
        @test l.mode == :search
        @test l.sleep_time == 0.01
        l.loop = Threads.@spawn RomeoDFT.loop(l, verbosity=3)
        p = RomeoDFT.ProgressMeter.ProgressUnknown("Simulating $case...", spinner=true)        
        while l.loop !== nothing && !istaskfailed(l.loop)
            RomeoDFT.ProgressMeter.next!(p, showvalues=[("Processed:", length(l[Done]))])
            sleep(0.5)
        end

        @test l.loop === nothing
        if l.loop !== nothing
            display(l.loop)
        end
        
        @test isempty(l[Error])

        for e in @entities_in(l, Trial)
            @test !ispath(joinpath(l, "$(Entity(e).id)"))
            @test e in l[Results] && !isempty(l[Results][e].state.occupations)
        end
    end
end
