This repository has been created to host files (code and data files) related to the integrated miRNomics and lipidomics datasets produced by the **Laboratory of Experimental Pharmacology and Integrative Biology of Atherosclerosis** ([EN](https://disfeb-unimi-it.translate.goog/it/ricerca/risorse-e-luoghi-della-ricerca/laboratori-di-ricerca/farmacologia-sperimentale-e-biologia-integrata-dellaterosclerosi?_x_tr_sl=auto&_x_tr_tl=en&_x_tr_hl=en-US&_x_tr_pto=wapp)|[IT](https://disfeb.unimi.it/it/ricerca/risorse-e-luoghi-della-ricerca/laboratori-di-ricerca/farmacologia-sperimentale-e-biologia-integrata-dellaterosclerosi)), of the [University of Milan](https://www.unimi.it/en).

**Authors:**
[Stefano Manzini](mailto:stefano.manzini@gmail.com)<sup>1</sup>*, Marco Busnelli<sup>1*</sup>, Mika Hilvo<sup>2</sup>, Matteo Chiara<sup>3</sup>, Cinzia Parolini<sup>1</sup>, Federica Dellera<sup>1</sup>, Giulia Sara Ganzetti<sup>1</sup>, Nadia Vaira<sup>1</sup>, David Stephen Horner<sup>3</sup>, Reijo Laaksonen<sup>2</sup>, and [Giulia Chiesa](mailto:giulia.chiesa@unimi.it)<sup>1</sup>

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
