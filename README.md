# Modeldefinition and code repository for publication "Resource competition in neural protein synthesis"

This repository includes all code to reproduce all data and subfigures for the publication. The instruction below lists all needed packages and indicates which figures have been assembled using inkscape. 
No figure content was changed using inkscape. Alle indicators in the figures are part of the presented python code. Mathematical descriptions of the 9 model variations are given as SBML files in the directory *SBML_model_descriptions/*. 
Provided simulation data can be found in the directory *simulation_results/*. Any newly simulated data will be automatically safed to that directory aswell.
All figures and subfigures will be automatically safed to the directory *figures/*. 

### Python Dependencies and Package version 
- Python 3.9
- Matplotlib 3.5.1
- Numpy 1.25.3
- Sklearn 1.2.2
- Cmcrameri 1.8
  
*Cmcrameri* is a collection of colourmaps and styles. Information on how to install it can be found here:
https://github.com/callumrollo/cmcrameri 

**Optional**: If you do not want to download *Cmcrameri* you can change the color map in the .py scripts to use any other continuous colour map: 

CTRL+F for *palette* in the .py files: 
*palette = ... #your custom continuous colourmap*. 

### Shorthands for model variations

The code uses following shorthands for the different resource conditions: 

- ul : unlimited
- l : limited
- c : competitive

This modifiers are used in function names and followed by: 

- R : ribosomes / initiation resource
- Y : aminoacids / elongation resource


### Producing the Figures
We recommend running the python scripts in this directory instead of the jupyter notebook. Both include the same code. The jupyter notebooks are still provided for an easy starting point if anyone wishes to further develop the model.
#### Figure 1
Run the notebook *figure1_mock_data.py*. This will provide you with the mock data which is used in Figure 1. The illustrations were added postscript using inkscape. 

#### Figure 2
This figure is composed of 9 subfigures. To produce the subfigures run notebook *figure2_parameter_exploration_models.py*. We recommend to only rerun the simulated data with sufficient compute time (up to 48hrs). Data to plot the results is provided. The produced subfigures where imported to inkscape and assembled in the grid seen in Figure 2. The subfigure names indicate the model variation used in the simulation. 


#### Figure 3-6
The remaining figures are all produced with the same script. It also saves the simulation data automatically in the directory *simulation_results/*. This way figures can be reproduced easily from the provided simulation data as well. Run the notebok *figure3_to_6_simulation_all_model_variations.py* to produce Figures 3-6. Figures 3-5 have not been modified postscript. Figure 6 received indications and coloured underlay to visualise the proposed regulation theory using inkscape. The simulation provides the underlying curve plot.


#### Appendix Figure
The appendix figure can be plotted by running *appendix_figures.py*. We recommend to rerun the simulation only with sufficient compute time. You can choose the number of initialisation before rerunning. The paper figure is run with 100 initialisations. This is the only script without corresponding jupyter notebook. The code for the appendix figure is at the end of the Figure 3-6 notebook. 
