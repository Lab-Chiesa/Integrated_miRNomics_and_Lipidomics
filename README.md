This repository has been created to host files (code and data files) related to the integrated miRNomics and lipidomics datasets produced by the [Laboratory of Pharmacology of Dyslipidemias and Atherosclerosis](http://www.disfeb.unimi.it/ecm/home/ricerca/laboratori-ricerca/laboratorio-di-farmacologia-delle-dislipidemie-e-dellaterosclerosi), of the [University of Milan](http://www.unimi.it/).

**Authors:**
[Stefano Manzini](mailto:stefano.manzini@gmail.com)<sup>1</sup>*, Marco Busnelli1*, Mika Hilvo2, Matteo Chiara3, Cinzia Parolini1, Federica Dellera1, Giulia Sara Ganzetti1, Nadia Vaira1, David Stephen Horner3, Reijo Laaksonen2, and [Giulia Chiesa](mailto:giulia.chiesa@unimi.it)1

 _* these authors equally contributed to this work_

1.	Department of Pharmacological and Biomolecular Sciences, Università degli Studi di Milano, via Balzaretti 9, 20133 Milano, Italy;
2.	Zora Biosciences, 02150 Espoo, Finland;
3.	Department of Biosciences, Università degli Studi di Milano, Via Giovanni Celoria 26 - 20133 Milano, Italy;



Please note that, as the article is currently being written, data is uploaded in either _anonymized tables_ or in _password-protected zipfiles_. Final data and results will be available without restriction upon publication.

- **Reproducing the results**

The repository root contains two bash scripts that can be launched to reproduce the analyses described in the paper. Please consider that each reconciler.py run takes a few hours to complete (depending on your hardware). reconciler.py supports a wide range of paramenters and settings, please refer to the paper and to reconciler.py’s help to adapt the parameters to suit the needs of your analysis.

- **Repository content**

**/data**

Contains the two datasets: digested data from smallRNAseq and mass-spec lipidomics of murine tissues, in tab-delimited tables.

**/lib**

Contains a Python library of common functions called by other programs.

**/miRNA_targets**

Contains TargetScan mRNA predictions of miRNAs discussed in the paper, keyword-based GeneCards gene lists, and results about miRNA to mRNA predictions.

**/reconciler**

Contails reconciler.py, the main program that does the correlations. Extensive description of its design is contained in the paper’s Supplementary Materials and Methods.

**/results**

Contains the pre-computed output of reconciler.py, discussed in the paper.

**/source**

Contains other Python source files of programs and scripts used in the analysis workflow. It also contains pickled Python objects used by these programs.
