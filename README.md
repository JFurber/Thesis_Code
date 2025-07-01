# Analysis of Badger Movement - Thesis
Files for the analysis of GPS badger movement in relation to the Thesis by Jessica Furber.

Files include:
<ul>
  <li> Initial modelling strategies (Chapter 3) </li>
  <li> Conversion of coordinates from Latitude/Longitude to Easting/Northing (and back) (Chapter 4) - <i>R coding</i> </li>
  <li> Code to estimate metastable states (Chapter 5), including: </li>
    <ol>
      <li> Interpolation of coordinates for correct format of data for extended dynamic mode decomposition - <i>Matlab</i> </li>
      <li> Generation of metastable states via extended dynamic mode decomposition (EDMD) - <i> Python </i> </li>
    </ol>
  <li> Code to calculate the model (Chapter 6), including: </li>
    <ol>
      <li> Estimation of diffusion and statistical analysis - <i> Python, Matlab, and R </i> </li>
      <li> Creation of KDE and GMM potential - <i> Python</i> </li>
      <li> Estimation of mountain passes - <i> Python</i> </li>
      <li> Comparison of potentials  - <i> Python</i> </li>
      <li> Simulation (with Attraction and Repuslion (AR) and without) - <i> Python</i> </li>
    </ol>
</ul>

You will also require to download the folder d3s from https://github.com/sklus/d3s in order to use the EDMD.py file.

The GPS tracking data of badger movements, collected from Cornwall, are lodged on Movebank (www.movebank.org; Movebank Project 158275131). For the GPS tracking data of badger movements, collected from Woodchester, please contact data access at APHA on badgerdata@apha.gov.uk. The GPS tracking data of badger movements, collected from Northern Ireland, are lodged on Movebank (www.movebank.org; Movebank Project 5019193945).
