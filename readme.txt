This is the readme for the model file associated with the paper:

Hoshino O, Kameno, R, Watanabe K (2019) Reducing variability in motor cortex activity at a resting state by extracellular GABA for reliable perceptual decision-making. J Comput Neurosci

This c-program was contributed by O Hoshino. It was originally built using Microsoft Visual C++ and also works in Microsoft Visual Studio 2012 (create a new project, add the c file to it, and build and run).

More usage instructions:

1.    Set input current to the sensory network by giving a proper value to "int_inp0_3". "onset_0" and "period_0" define its onset time and duration.
2.    Set times for output data by giving values "OUT" (starting time) and "PERIOD" (recording time-period).
3.    Run. The default (as provided) settings of the model is for figure 2A.
4.    Output data files (***.dat) give membrane potentials of pyramidal cells (and local GABA concentrations in N_S that are kept at the basal level). 