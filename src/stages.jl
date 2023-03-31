creation_stage() = Stage(:creation, [JobCreator(),
                                     Relaxor(),
                                     HPCreator(),
                                     BandsCreator(),
                                     NSCFCreator(),
                                     ProjwfcCreator()])

postprocessing_stage() = Stage(:postprocess, [ResultsProcessor(),
                                              RelaxProcessor(),
                                              HPProcessor(),
                                              UniqueExplorer()])
    
full_server_interaction_stage() = Stage(:server_interaction, [[JobSubmitter(), Cleaner(), OutputPuller()],
                                                              SimJobRemover(),
                                                              JobMonitor()])

nosubmit_server_interaction_stage() = Stage(:server_interaction, [[Cleaner(), OutputPuller()],
                                                                   SimJobRemover(),
                                                                   JobMonitor()])

finalize_stage() = Stage(:finalize, [ErrorCorrector(), Rerunner(), Stopper()])

core_stage() = Stage(:core, [creation_stage(),
                             full_server_interaction_stage(),
                             postprocessing_stage(),
                             finalize_stage()])

cleanup_stage() = Stage(:cleanup, [creation_stage(),
                                   nosubmit_server_interaction_stage(),
                                   postprocessing_stage(),
                                   finalize_stage()])

intersection_stage() = Stage(:intersection, [Intersector(), RandomTrialGenerator()])
    
# firefly_stage() = Stage(:firefly, [FireFly(), PostFireflyExplorer(), Archiver()])

search_stage() = Stage(:main, [intersection_stage(), core_stage()])

function set_searcher_stages!(l::AbstractLedger, s::Symbol)
    if s ∈ (:postprocess, :manual)
        Overseer.ledger(l).stages = [core_stage()]
    elseif s == :search
        Overseer.ledger(l).stages = [search_stage()]
    elseif s == :cleanup
        Overseer.ledger(l).stages = [cleanup_stage()]
    else
        error("Searcher stage $s not recognized...")
    end
    prepare(l)
end

