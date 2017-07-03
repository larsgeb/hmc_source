This is the first implementation of Hamiltonian Monte Carlo for synthetic point source inversion. The code has two parts: one for the generation of artificial data (far-field P waveforms), and one for Monte Carlo inversion.


Synthetic Data:
--------------

Compilation: g++ aux.cpp make_data.cpp -o make_data
Run: ./make_data

This generates far-field P waveforms for the source-receiver setup in INPUT/setup.txt and the (prior mean) parameters in INPUT/parameters.txt.


Monte Carlo Inversion:
—--------------------

Compilation: g++ aux.cpp mc.cpp sampling.cpp -o sampling
Run Metropolis-Hastings: ./sampling metropolis [# samples] [verbose]
Run HMC: ./sampling hamilton [# samples] [# time steps] [time increment] [verbose]

The output is written to OUTPUT/samples.txt and can then be analysed with the Python scripts plot_samples.py and plot_trajectory.py.