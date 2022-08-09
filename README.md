# In-Silico Pace-mapping

If using any of these scripts, please cite

Sofia Monaci, Marina Strocchi, Cristobal Rodero, Karli Gillette, John Whitaker, Ronak Rajani, Christopher A. Rinaldi, Mark O'Neill, Gernot Plank, Andrew King, Martin J. Bishop. *In-silico pace-mapping using a detailed whole torso model and implanted electronic device electrograms for more efficient ablation planning*, Computers in Biology and Medicine, Volume 125, 2020.

## ABOUT

This computational platform can be used to generate in-silico pace-maps from clinical/simulated ECGs/EGMs paces and simulated/clinical VT ECGs/EGMs to **localise critical ablation targets of post-infarct VTs (e.g. exits, entrances, isthmuses) non-invasively**.

Python scripts only require basic libraries such as numpy, pandas, os, sys, argparse and matlab.engine (visit (here) [https://uk.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html] for instructions on how to set up matlab engine for python). However, the overall pipeline requires tools such as meshtool and GlVTKConvert
