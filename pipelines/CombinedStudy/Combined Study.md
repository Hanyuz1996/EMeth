# Combined Study

Details of experiment design and results can be found in the supplementary materials. This pipeline include four files. 

- GetMethylationMatrix.R: preprocess data, filtering probes.
- EstimatingMethylation.R:  estimate mean and variance for each probe and construct the reference data.
- Simulation.R: functions to generate simulation data from mixture model and run different algorithms.
- batchsim.R: script to generate experiments with repetitions under different settings.