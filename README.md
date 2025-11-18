Analysis code to study Timelike Compton Scattering and J/Psi photoproduction with the CLAS12 detector


# Run JPsi analysis

## Setup the path to the data
change the path to the data here:
https://github.com/PChatagnon/TCS_Analysis/blob/aec59ef95ab8ef469854883ef780eb2e5eb9a391/CS_Extraction/bib_CS_extraction/Sample_Class.h#L297

## Run the BSA analysis

```
cd CS_Extraction
root -l BSA_Extraction.C Results_CS/Config_Standard.dat
```

## Run the CS analysis

```
cd CS_Extraction
root -l CS_Extraction.C Results_CS/Config_Standard.dat
```

Pdf containing all plots in Results_CS
