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

search_stages() = [Stage(:main, [intersection_stage(), core_stage()])]
