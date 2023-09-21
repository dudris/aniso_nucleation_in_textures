function results_structure = load_results(resultsfile,results_name)
load(resultsfile,results_name)
results_structure = eval(results_name);
end