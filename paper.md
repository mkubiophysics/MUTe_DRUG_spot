title: 'MUTE_DRUG_spot : A flexible R based shiny application to plot mutation-drug relations'
tags:
  - R
  - mutations
  - data visulization
  - drug targets
  - 
authors:
  - name: Aakriti Jain
    orcid: 0000-0003-2451-4551
    affiliation: 1
  - name: Aayush Srivastava
    affiliation: 1
  - name: Manish Kumar
    orcid: 0000-0002-7936-9892
    corresponding: true 
    affiliation: 1
affiliations:
 - name: Department of Biophysics, University of Delhi South Campus, New Delhi, India. 
   index: 1
   ror: 04gzb2213
date:  17 April 2025
bibliography: paper.bib


# Summary
Mutations in drug targets is one of the fundamental mechanisms in development of drug resistance and therefore considered as 1:1 biomarkers of drug efficacy. Sequencing technologies have helped in generating massive amounts of data on gene-mutation-drug relations but the analysis and visualization of this data sometimes might be difficult for biomedical researchers with limited bioinformatics’ expertise. Till now most of the drug-resistance causing mutations are analyzed and visualized using their genomic sequences but there are many instances where amino-acid mutations lead to drug-resistance. 
With an objective to understand and analyse the latter  we have  developed, MUTe_DRUG_spot, a freely available R based shiny application.  It is easy to use, download and comprehend with no requirement of advanced computational skills. Once downloaded, the interface of our application (Capture_app.PNG) would open and mutation-drug data can be input in the required fields to obtain customized plots for your own mutation dataset. The user has to just provide a list of mutations and drugs along with the UniProt ID of the drug target. The tool automatically fetches all other information of the drug target. 
The generated plots () are actually 3D lollipop plots with drug target (protein) length on the x-axis , the mutations on the y-axis and the drugs against which resistance is conferred as a color of the lollipop. Additionally the type of mutations can be indicated using different shapes provided in a catalogue. The diversity of plots and the information which can be analysed through it vouches for the utility of this application.
This application can be a useful tool for studying mutation-drug relations from data collected through published literature or generated through experiments. It can also be combined with other text mining tools to assess the suitability of any potential drug target. This tool will be useful not only for plotting mutation-drug target relations of pathogenic diseases but will also come handy for visualization of cancer genomics and proteomics data thereby, advancing precision medicine.  



# Statement of need
Sequencing technologies have advanced to a great extent leading to reduction in costs and time and increase in the amount of data. The application of this surplus amount of data is subject to knowledge generation which requires compatible tools for easy visualization and analysis. Alteration in gene or protein sequences of target molecules are important markers for investigation in 
We present a tool called MUTe_Drug_Spot which offers an easy way to plot a set of point mutations/Single Nucleotide polymorphisms (SNPs) and an additional parameter in one single plot. It can be a useful tool for experimentalists, clinicians and scientists analyzing the effect of mutations/SNPs on a drug target and drug simultaneously. The tool is freely accessible at https://github.com/mkubiophysics/MUTe_DRUG_spot
Studying mutation-drug relations is a relevant piece of information of any potential drug target. 
Therefore, this application can be an excellent tool for visualizing mutation-drug relations of your protein from the data collected through published literature or generated through experiments. It can be combined with other text mining tools to retrieve mutation-drug data of various proteins and visualize it easily using our application.    

   




# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.


# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.



# References
