# README: Internal Std Correction

Lexi Jones-Kellett

1. internal_standards_quant.ipynb
	- Wrote after following "PicoGreen Fluorometric DNA Quantification Protocol" to estimate internal standard
          concentration in each batch.
	- Used to derive DNA concentration averages from replicate estimates. 

2. normalize_by_IS_quant.ipynb
	- Tested code to calculate absolute ASV abundance from internal standard reads.
	- Includes some QC steps to check that everything worked as planned

3. calc_ASV_abundance_from_IS_v2.py
	- Code applied to entire corrected merged 16S+18S relative abundance table to calculate absolute abundance. 
	- Needs to be run for each internal standard (Blautia_producta, Deinococcus_radiodurans, Thermus_thermophilus)

4. We averaged the counts estimated from each of the internal standards for the rest of analysis. 


NOTE: Much of this protocol followed Lin et al. 2019 (https://doi.org/10.1128/AEM.02634-18); PNGs include some of the 
constants and equations.
