========
Overview
========


Package to calculate thermal discomfort severity under several thermal definitions (e.g. traditional thermal comfort, sleep comfort, occupant health and safety limits, etc.).

* Free software: MIT license

Installation
============

::

    pip install compy   
    

Instructions
============
With this package, you can get the discomfort severity associated with a large set of comfort parameters. Please follow these steps:

**Step 1** - After installing the package, please download the following files: 



    - `main <https://drive.google.com/file/d/1Hg5VSDoSRkicoWOsJpVOdjX93ajakcro/view?usp=sharing>`_
    - `PMV model_template <https://drive.google.com/file/d/10ZniYVqR-SPKyVC1qz7ml7JPQcSh-xPY/view?usp=sharing>`_ ----> *if working with PMV comfort model*
    - `Adaptive model_template <https://drive.google.com/file/d/1qEHnlmfTOpabHgXz10lTxmak254B5ogf/view?usp=sharing>`_ ----> *if working with Adaptive comfort model*
    
**Step 2** - Open the downloaded .csv file relative to your project. For instance, if you are working with adaptive thermal model, open *Adaptive model_template.csv* file. Columns represent comfort parameters related to the selected comfort model (i.e., PMV or Adaptive) along with date and time. It's important that the order of columns remains unchanged. For your reference each template contains some example data, please replace them with your data. If you do not have data for a certain parameter, simply fill the missing data with a representative value. Every row must contain a value for each column.  

**Step 3** - Run the *main* file. A file named *discomfort_severities_PMV.csv* or *discomfort_severities_Adaptive.csv* (depending on the selected comfort model) will be automatically saved in the same directory as the downloaded files are located. The file contains all the discomfort severities. 
