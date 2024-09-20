# README: Internal Std Correction

Lexi Jones-Kellett

The scripts use a "conf.json" file to point to the local directories as follows:
	{
		"data_dir": "/path/to/file/",
		"ASV_count_dir": "/path/to/file"
	}


internal_standards_quant.ipynb
	- Wrote after following "PicoGreen Fluorometric DNA Quantification Protocol" to estimate internal standard
          concentration in each batch.
	- Used to derive DNA concentration averages from replicate estimates. 

calc_ASV_abundance_from_IS_v2.py
	- Code applied to entire corrected merged 16S+18S relative abundance table to calculate absolute 
	  abundance. 
	- Needs to be run for each internal standard (Blautia_producta, Deinococcus_radiodurans, 
	  Thermus_thermophilus)
	- methods described in "ASV_abundance_methodology.png"


NOTE: Much of this protocol followed Lin et al. 2019 (https://doi.org/10.1128/AEM.02634-18).
