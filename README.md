# Chaotic Winds from a Dying World
This repository contains code used in Chaotic winds from a dying world: a one-dimensional map for evolving atmospheres by Bromley & Chiang 2023. <https://ui.adsabs.harvard.edu/abs/2023MNRAS.521.5746B/abstract>

This code was primarily used to iterate the map and generate the figures used in the paper. In general, functions called `radiativeTransfer` or similar correspond to Map A, functions called `pierrehunbert` corespond to Map B and functions called `guillot` correspond to Map C.
These names are derived from the references where the radiative transfer equations are derived from. Additionally, the data for figures 10 and 14 are contained in files `chaoticPointsMC2.txt` and `chaoticPointsGuillot2.txt` respectively. 
These data are generated in the file `chaosClassifier.py` and plotted in `volumePlotter.py`. The format for the data in each row is `p1,p2,p3,lyExp`

Any questions about this code should be directed to Joshua Bromley (`j.bromley AT mail.utoronto.ca`)
