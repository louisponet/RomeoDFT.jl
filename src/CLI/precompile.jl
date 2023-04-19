using RomeoDFT
try
    RomeoDFT.Orchestrator.status()
catch
    nothing
end
l = Searcher(; rootdir = tempname(), sleep_time = 5)
save(l)
load(l)
RomeoDFT.command_main()
