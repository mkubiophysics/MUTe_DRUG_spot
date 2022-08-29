# Mute_DRUG_spot
A freely available tool for plotting mutations of a protein along with the corresponding drug against which they are known to confer resistance. Hence, goes the name Mute_Drug_spot. 

![Capture_app.png](https://github.com/mkubiophysics/Mute_Drug_spot/blob/main/Capture_app.png)

## Installation
You only need to have an installation of RStudio in order to run this application (but it is recommended that the dependencies are already present on your system; see the next section). All it takes to launch Mute_Drug_spot is 2 steps, which are as follows:
 - To install `shiny`, open an RStudio session and paste the command `install.packages("shiny")`, followed by `library(shiny)` to import the package into your session.
 - Once `shiny` is imported into the session, simply paste the command `runGitHub("MUTe_DRUG_spot", "mkubiophysics")`. This command directly downloads and runs the application within the RStudio session.
 - _BONUS:_ If you wish to verify and/or modify the code first, you can clone the repository on to your local machine. This can be achieved either by clicking _Code -> Download ZIP_, or by typing the command `git clone https://github.com/mkubiophysics/MUTe_DRUG_spot.git` on a terminal. Once the repository is downloaded on your local machine, open the `app.R` file in your RStudio session, and once you are done verifying the code, simply press _Run App_ on the top right corner of the editor window, or press `Ctrl+A` and `Ctrl+Enter` to launch the app.

## Dependencies
The application has the following dependencies:
 - shiny
 - xml2
 - BiocManager
 - tidyr
 - dplyr
 - stringr
 - shinyWidgets
 - shinycustomloader
 - RColorBrewer
 - bslib
 - trackViewer (installed via BiocManager)

Although the application will automatically install these packages if they are not found on your system, installing them may take a while. Thus, it is recommended to install them beforehand so that launching the application for the first time does not take a long time. To install the above packages (except `trackViewer`), type the command `install.packages({packagename})`, where {packagename} is the name of the package you are installing. Installing `trackViewer` requires `BiocManager` to be already installed on your system. Once you have it installed, paste the command `BiocManager::install("trackViewer")` to install `trackViewer`.

# Usage
Once you launch Mute_Drug_spot through RStudio, all you need to do is provide a UniProt ID of the protein and a file containing a list of mutations and their corresponding drug resistances. This file, ideally in the CSV format, should look like this:

```
Mutation, Resistance
A1C, "Drug 1"
D2E, "Drug 2"
F3G, "Drug 3"
```

Make sure that the drug names are enclosed inside double quotes ("). Mute_Drug_spot supports visualization of different types of mutations, includins mis-sense mutations, non-sense mutations, and frameshift mutations. These are displayed via different symbols on the plot, corresponding to:

|Mutation Type      | Shape   | Example |
|-------------------|---------|---------|
|Mis-sense mutation | Circle  |  A1C    |
|Non-sense mutation | Square  |  A1*    |
|Frameshift indel   | Diamond |  A1fs   |
|In-frame insertion | Triangle| 1A>AAA  |

The syntax for writing the different kinds of mutations is also available in the file `example_mutaions.csv` in this repository.

Once you have this list of mutations ready, the main steps for the usage of Mute_Drug_spot are:
 - On the top left corner of the app, type the UniProt ID of your protein of interest and press Enter. Mute_Drug_spot uses the UniProt API to fetch information about the protein, including its name, length, and domain information.
 - If your protein of interest does not have a UniProt ID, or if Mute_Drug_spot gives missing/incorrect information, you can manually edit all these parameters in the input fields below. For example, the length of the sequence and number of domains can be modified just below the UniProt ID input field (note that Mute_Drug_spot will always display at least one domain; if your protein does not have any domains, set the length of this domain equal to the total length of the protein).
 - Upload the list of mutations in the input box provided, and set the file delimiter from the dropdown box below (Mute_Drug_spot supports both CSV and TSV file formats, but CSV files are preferred).
 - Ensure that the plot displayed on the center of the screen is updated and all the parameters are displayed correctly. If the mutations or domains are displayed incorrectly, check the list of mutations and the domain information input fields to make sure there are no errors.
 - In order to change the color selected for each drug, go to the Customize Colors panel below the plot and set the colors manually.
 - The labels of the x- and y-axes can be changed through their corresponding input fields in the Customize Labels panel below the plot. The Title field is automatically set to the name of the protein determined by Mute_Drug_spot if you provide a UniProt ID for the protein. Any of these fields can be set to empty values in order to remove the corresponding label(s) from the plot.
 - If the title of the plot is displayed too far from, or too close to, the plot, you can change its display height using the Title Position slider.
 - The Shrink Plot checkbox reduces the height of the mutation spots in the plot. This can be useful if your list of mutations has a large number of unique drug names, in which case the legend on top of the plot will take a large amount of space and potentially overlap with the plot title.
 - Once all the parameters have been set to the desired values, click on Download Plot in order to save the generated plot as a PDF file.
 
## FAQs
 - _Error below the UniProt ID field_: Mute_Drug_spot uses the UniProt API to fetch information corresponding to the UniProt ID supplied by the user. It will display an error just below this field if it is unable to find any information on UniProt using your given ID. This value is updated only when Enter is pressed in the field (to reduce the number of calls made to the UniProt API), so any time you wish to modify the ID of the protein, make sure to press Enter to trigger the API call. If you wish to clear any information obtained from UniProt and set all the parameters manually, simply remove any value in the UniProt ID field and press Enter.
 - _Warning about number of drugs_: If there are more than 12 unique drug names in your list of mutations, Mute_Drug_spot will display a warning message. This is because `RColorBrewer`, the package used by Mute_Drug_spot to color the mutation spots, only supports a maximum of 12 colors for color palette generation. Although Mute_Drug_spot can automatically interpolate in-between this color palette to create a larger palette, the resulting colors may not be easily distinguishable from each other. In such a case, you can manually set the colors for each drug, as mentioned in the Usage section.
 - _There are files with random names in the Mute_Drug_spot folder_: Mute_Drug_spot saves a temporary version of your plot every time you make an update to any parameter value. This not just ensures that any previous version of the plot is available for recovery, but is in fact necessary for Mute_Drug_spot to be able to save the plots at all (due to the way R handles file downloads). In an ideal case, Mute_Drug_spot automatically deletes any temporary files older than 24 hours (note that Mute_Drug_spot needs to be running in order for this housekeeping to work), but if you wish to remove these files manually, you can safely remove any PDF files in the folder.
 
 ## Issues/Feedback
 If you have any other query, issue, or feedback regarding Mute_Drug_spot, feel free to open a new Issue on this Github repository. But make sure to first read through the code documentation in order to see if your query has already been documented here.
 
 ## Cite:
 If you used Mute_Drug_spot in your research work, please cite the tool as:
 {WIP}
