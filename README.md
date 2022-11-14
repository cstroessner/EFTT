Code to reproduce the figures in "Functional low-rank approximation based on the extended tensor train format" by C. Stroessner, B. Sun and D. Kressner.
(Arxiv: https://arxiv.org/abs/TODO.)

Matlab Dependencies:
Chebfun (download from https://github.com/chebfun/chebfun)
TT-Toolbox (download from https://github.com/oseledets/TT-Toolbox)
MATLAB Statistics and Machine Learning Toolbox
MATLAB Partial Differential Equation Toolbox

Python dependencies and installations are described in the notebooks (we use https://github.com/goroda/Compressed-Continuous-Computation and https://fenicsproject.org/).
After running MATLAB setup.m, you need to copy the file testPoints.mat from the test folder into the notebook1 and notebook2 folders.

Matlab File Descriptions:
setup.m: Initalization and correct setting of paths (needs to be run first).
generateData.m: Runs all numerical experiments not involving the solution of PDEs.
displayTables.m: Displays all results not involving PDEs as tables.
plotFigureX.m: Generates the plot for Figure X in the paper. 
generatePDEData.: Runs and displays the results for all numerical experiments involving the solution of PDEs.

Python File Descriptions:
Notebook1: Generates and displays the data for the table in the appendix and for Figure 3.
Notebook2: Generates and displays the data for Table 1 which involves the solution of PDEs.
