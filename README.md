# NMF-RI: Blind spectral unmixing of highly mixed multispectral flow and image cytometry data.

This repository contains the official MATLAB implementation of:

**NMF-RI: Blind spectral unmixing of highly mixed multispectral flow and image cytometry data.**  
Daniel Jiménez-Sánchez*, Mikel Ariz*, José Mário Morgado, Iván Cortés-Domínguez, Carlos Ortiz-de-Solórzano. https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz751/5583734  (*: equal contribution)

contact: codesolorzano@unav.es, danijimnzs@gmail.com

![NMF-RI](nmf-ri.png)

**Abstract:**

**Motivation**
Recent advances in multiplex immunostaining and multispectral cytometry have opened the door to simultaneously visualizing an unprecedented number of biomarkers both in liquid and solid samples. Properly unmixing fluorescent emissions is a challenging task, which normally requires the characterization of the individual fluorochromes from control samples. As the number of fluorochromes increases, the cost in time and use of reagents becomes prohibitively high. Here we present a fully-unsupervised blind spectral unmixing method for the separation of fluorescent emissions in highly mixed spectral data, without the need for control samples. To this end, we extend an existing method based on Non-negative Matrix Factorization, and introduce several critical improvements: initialization based on the theoretical spectra, automated selection of ‘sparse’ data and use of a re-initialized multi-layer optimizer.

**Results**
Our algorithm is exhaustively tested using synthetic data to study its robustness against different levels of colocalization, signal to noise ratio, spectral resolution, and the effect of errors in the initialization of the algorithm. Then we compare the performance of our method to that of traditional spectral unmixing algorithms using novel multispectral flow and image cytometry systems. In all cases, we show that our blind unmixing algorithm performs robust unmixing of highly spatially and spectrally mixed data with an unprecedently low computational cost. In summary, we present the first use of a blind unmixing method in multispectral flow and image cytometry, opening the door to the widespread use of our method to efficiently pre-process multiplex immunostaining samples without the need of experimental controls.

### Run experiments

From the project's root folder for each experiment, simply run
```
main.m
```

### Use of NMF-RI on new data

To use NMF-RI on new data (tissue or flow cytometry) we recommend to use the ´main.m´ code distributed and load new cytometry data to the existing folder hierarchy.

### Cite
If you make use of this code in your own work, please cite our paper:
```
@article{10.1093/bioinformatics/btz751,
    author = {Jiménez-Sánchez, Daniel and Ariz, Mikel and Morgado, José Mário and Cortés-Domínguez, Iván and Ortiz-de-Solórzano, Carlos},
    title = "{NMF-RI: Blind spectral unmixing of highly mixed multispectral flow and image cytometry data}",
    journal = {Bioinformatics},
    year = {2019},
    month = {10},
    doi = {10.1093/bioinformatics/btz751},
}
```
