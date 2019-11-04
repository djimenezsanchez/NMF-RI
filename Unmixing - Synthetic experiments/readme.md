# NMF-RI: Blind spectral unmixing of highly mixed multispectral flow and image cytometry data.

**Experiments and code described in: "NMF-RI: Blind spectral unmixing of highly mixed multispectral flow and image cytometry"
Authors: Daniel Jimenez-Sanchez, Mikel Ariz, Jose Mario Morgado, Ivan Cortes-Dominguez, Carlos Ortiz-de-Solorzano.
Contact: codesolorzano@unav.es

**Each folder describes one section of the paper:**
- Synthetic Experiments
- Experiments with fluorescent microbeads
- Flow Cytometry (AURORA, Cytek)
- Tissue Cytometry (Vectra Polaris, Perkin Elmer)

**In each folder, execute "main.m" to:**
- Load the data from the related folders, and create the Observation Matrix, Y.
- Pre-Process data using "dataPreprocessing.m" using the loaded data and the sparseness threshold parameter.
- Unmix data using using "NMF-RI.m".
- Display results as described in the paper. "displayResults.m"
