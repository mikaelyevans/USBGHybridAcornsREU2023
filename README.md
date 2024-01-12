## Project Description

This github project analyzes oak seedlings to detect hybridization using
DNA markers. Oak trees are critical ecosystem creators and are of very
high conservation value. Additionally, oaks are an exceptional species,
meaning they cannot be seedbanked using traditional methods and can only
be conserved in living collections like botanic gardens or arboreta.
Maintaining living collections of rare and endangered oaks is a
challenging task because they may frequently hybridize with other oak
species. In other words, seedlings collected from a maternal tree in a
botanic garden may be unintentional hybrid offspring that are not ideal
for conservation purposes. This limits conservationists ability to use
seeds produced in ex situ collections for restoration of threatened and
endangered plant species unless the parentage is verified and the
seedlings are not hybrids.

We performed a study analyzing the parentage of acorns produced by
maternal Quercus muehlenbergii individuals at the Morton Arboretum to
identify if (1) any offspring individuals were hybrids and (2) what
factors contribute to the parentage of offspring produced in living
collections. We collected a total of 385 seeds, grew them to seedling
stage, and sampled leaf tissue from the seedlings. Using 14
microsatellite loci, we performed parentage analysis using CERVUS
software.

## Folder Descriptions

### Archive:

The Archive folder includes all files that are **not** integral to the
analysis of this project. Inside it, there are three folders:
Data\_Files\_Archive, UHA\_Attempted\_md\_Analysis, and
UHA\_nomd\_PolyPatEx\_Analysis. The organization of the Archive folder
is as follows:

-   Data\_Files\_Archive → This folder contains two duplicate data files
    that are not used in the final analysis.
-   UHA\_Attempted\_md\_Analysis → This folder contains two folders:
    UHA\_Cervus\_Attempts and UHA\_PolyPatEx\_Attempts. The analysis
    included in these folders was done with files that had too much
    missing data so they are not relevant to the rest of our study. The
    purpose of the files was only for practicing the types of analysis
    with the data we had at the time.
-   UHA\_nomd\_PolyPatEx\_Analysis → This folder contains the fully
    cleaned data set analysis using PolyPatEx. For this project, using
    the CERVUS analysis gave us better, more accurate parentage
    assignment results, so we decided against using the PolyPatEx
    results. All the data and code files can be found here for PolyPatEx
    analysis.

### Data\_Files:

The Data\_Files folder includes two subfolers: Clean\_Data\_Files and
Data\_Cleaning\_Code. The purpose of this folder is to include all the
clean data files and to show how they were cleaned. These cleaned data
files were used in the paternity analysis to produce our final results.

-   Clean\_Data\_Files → This folder contains two subfolders, Clean\_csv
    and Clean\_xlsx. The files are the same, just in two different
    forms. These files are the results from the Data\_Cleaning\_Code
    folder.
-   Data\_Cleaning\_Code → The Data\_Cleaning\_Code folder contains an
    important script: UHA\_datacleaning.R, written by Emily Schumacher.
    This R script cleans the data in order to make it useable for the
    analysis. The files directly involved with the data cleaning R
    script are “UHA\_Final\_Scores.arp”, “UHA\_Final\_Scores.gen”, and
    “UHA\_Final\_Scores\_nomd\_genalex.csv”. The other csv files that
    were not directly involved with this script have minor adjustments
    made to them so they can be used in analysis as well. These files
    are “UHA\_Final\_Scores\_nomd\_df.csv” and
    “UHA\_Final\_Scores\_GenAlEx.csv”.

### Paternity\_Analysis:

The Paternity\_Analysis folder includes the cleaned data file used for
analysis and code for the CERVUS analysis. The organization of the
paternity analysis is as follows:

-   Results → This folder contains a R Markdown script written by
    Mikaely Evans that was used to create the final figures for this
    project. These figures were used in her symposium presentation, but
    are now outdated because this project is still in process. This
    folder also contains “UHA\_parentage\_results.csv” which is the csv
    file used to create these figures. It includes variables like
    Offspring ID, Mother ID, Candidate father ID, and Candidate father
    species.

-   UHA\_nomd\_CERVUS\_Analysis → This folder contains the final no
    missing data parentage analysis input and output files from this
    study. CERVUS is a software made specifically for parentage analysis
    and it uses the three specific data sets mentioned in the
    CERVUS\_Input\_Files folder.
