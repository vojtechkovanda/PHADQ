This file describes the MATLAB code and other files used in the article *Audio dequantization using instantaneous frequency.*

The code has been developed in MATLAB version R2025a

In the root folder, there is a main m-file `main.m`, and it runs a dequantization experiment. Within this file, it is possible to change the input audio file (`audiofile`) and the bit depth of quantization (`param.delta`) for the quantized observation. It is also possible to modify other parameters such as the DGT frame settings, number of iterations or general options of the CP numerical algorithm. The default values are those used in the experiments in the paper for the sake of reproducibility.

The other files in the root folder provide supporting functions and tools for the dequantization algorithms, including signal quantization and projection onto the set of feasible values.
 
Folder  `phase-correction` contains a set of MATLAB codes written as a supplementary material of the following tutorial paper explaining phase-related property of spectrograms:

-   Kohei Yatabe, Yoshiki Masuyama, Tsubasa Kusano and Yasuhiro Oikawa, "**Representation of complex spectrogram via phase conversion**,"  _Acoustical Science and Technology_, vol.40, no.3, May 2019. (Open Access) ([PHAIN/phase_correction at main · TomoroTanaka/PHAIN](https://github.com/TomoroTanaka/PHAIN/tree/main/phase_correction))

The MATLAB codes were used for time-frequency signal processing and phase correction.
