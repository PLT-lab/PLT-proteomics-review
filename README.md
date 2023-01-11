# **Proteomics-wise, how similar are mouse and human platelets (and megakaryocytes)?**

Proteomics data analysis performed in the following systematic review:

## Database search of published platelet proteomics studies

We performed a search on Pubmed using the following keywords: platelet
AND proteomics AND (lysate OR releasate OR secretome). After manually
screening the abstract and method section, 30 studies of whole platelet
lysate (27 in human, 3 in mouse )<sup>30</sup>, and 5 studies of
platelet releasate (4 in human, 1 in mouse) were selected<sup>34</sup>.
The selection criteria took into account that the proteomics analysis
was done using mass spectrometry, and that the protein digestion was
performed in solution. Differences on the sample processing, mass
spectrometer and search engine used were noted (see
data/metadata_datasets.xlsx). In addition, these studies had to have a
publicly available raw protein dataset (*i.e.*, datasets only showing
differentially expressed proteins were not included), where healthy
control samples had to be either clearly identified or easily deduced.

In addition, FASTA files of the human and mouse reference proteomes were
downloaded from UniProt (November 2022). These included only reviewed,
canonical proteins; with 20385 and 17127 entries, respectively.

## Dataset cleaning, filtering, and analysis

All datasets were stored in a single Excel file (data/all_datasets.xlsx,
one dataset per sheet), and imported into RStudio. For each dataset, the
following cleaning and filtering was performed, when needed: protein
and/or gene identifications (*i.e.*, UniProt IDs, gene names) were
cleaned, so that only the leading identified, canonical protein/gene
remained; additionally, a filtering step was performed to ensure that
the working proteins were expressed in at least two thirds of the
controls. Moreover, and if it was specified in the dataset, contaminants
and/or detected decoys were removed, as well as proteins with
low-confidence detection. Lastly, to homogenize all datasets,
protein/gene identifications were annotated so that they had an
accompanying UniProt IDs, ENTREZID, SYMBOL and ENSEMBL ID.

In order to select a reliable set of proteins that composed the core of
the platelet proteome/releasate, all identified proteins across any of
the respective datasets (*i.e.*, 27 datasets for the human platelet
lysate) were filtered so that only reviewed proteins were used, and the
remaining were merged, and their occurrence counted. As a rule of thumb,
those proteins that were detected in at least half of the datasets
(*i.e.*, in 12 of the 27 datasets, for the human platelet lysate, so
that it also reached a total count above 2000) were considered as
reliably expressed and part of the core proteome. In addition, orthologs
were obtained in each case (both human to mouse, and mouse to human), as
well as the overlap with the reference proteomes, all of it represented
as Venn diagrams.

Those datasets that presented with clear, reliable relative
quantification were used to study the protein distribution (x11 and x39
for the platelet lysate proteome, and s5 and s4 for the releasate
proteome; for human and mouse, respectively). Thus, controls were
selected, features were filtered based on missing values and the median
expression of each protein was calculated. Based on this value, proteins
were ranked for plotting, and the distribution was divided into three
subsets, based on its quartiles (first quartile, inter-quartile, and
third quartile). Each subset was further subject to a gene ontology
enrichment analysis, plus an extra simplification step to remove
redundancy of the resulting enriched GO terms, to determine which known
biological functions were over-represented in each one of them.

Lastly, to study the overlap and correlation between platelet
transcriptomics and proteomics, both in mouse and human, datasets from
the study by Rowley *et al*. were used<sup>35</sup>. Transcripts were
filtered based on their RPKM expression (RPKM \> 0.3), and the resulting
datasets were overlapped with the respective proteomics data (core PLT
proteome). Those common features further underwent a Pearson correlation
analysis against the proteomics datasets (x11 and x39 for mouse and
human, respectively), after log2-transformation of RPKM values.

All the data manipulation and analysis were conducted using R (R Core
Team, version 4.0.3). Handling of the data was performed using the
‘dplyr’<sup>36</sup> and ‘stringr’<sup>37</sup> libraries, and plotting
with the ‘ggplot2’<sup>38</sup>, ‘ggpubr’<sup>39</sup> and
‘eulerr’<sup>40</sup> libraries. The Bioconductor packages
‘AnnotationDbi’<sup>41</sup>, ‘org.Hs.eg.db’<sup>42</sup> and
‘org.Mm.eg.db’<sup>43</sup> were used to annotate the data,
‘gprofiler2’<sup>44</sup> to extract the orthologs, and
‘clusterProfiler’<sup>45</sup> to perform the enrichment analysis.

# References

<span class="csl-left-margin">1. </span><span
class="csl-right-inline">Rodrigues RM, Valim V de S, Berger M, et al.
[The proteomic and particle composition of human platelet lysate for
cell therapy products](https://doi.org/10.1002/jcb.30310). *Journal of
Cellular Biochemistry*. 2022;123(9):1495–1505. </span>

<span class="csl-left-margin">2. </span><span
class="csl-right-inline">Linge CP, Jern A, Tydén H, et al. Enrichment of
complement, immunoglobulins, and autoantibody targets in the proteome of
platelets from patients with systemic lupus erythematosus. *Thromb.
Haemost.* 2022;122(9):1486–1501. </span>

<span class="csl-left-margin">3. </span><span
class="csl-right-inline">Trugilho MRO, Azevedo-Quintanilha IG, Gesto
JSM, et al. Platelet proteome reveals features of cell death, antiviral
response and viral replication in covid-19. *Cell Death Discov.*
2022;8(1):324. </span>

<span class="csl-left-margin">4. </span><span
class="csl-right-inline">Tassi Yunga S, Gower AJ, Melrose AR, et al.
Effects of ex vivo blood anticoagulation and preanalytical processing
time on the proteome content of platelets. *J. Thromb. Haemost.*
2022;20(6):1437–1450. </span>

<span class="csl-left-margin">5. </span><span
class="csl-right-inline">Almeida LGN de, Young D, Chow L, et al.
Proteomics and metabolomics profiling of platelets and plasma mediators
of thrombo-inflammation in gestational hypertension and preeclampsia.
*Cells*. 2022;11(8):1256. </span>

<span class="csl-left-margin">6. </span><span
class="csl-right-inline">Su M, Chen C, Li S, et al. Gasdermin
d-dependent platelet pyroptosis exacerbates NET formation and
inflammation in severe sepsis. *Nat. Cardiovasc. Res.*
2022;1(8):732–747. </span>

<span class="csl-left-margin">7. </span><span
class="csl-right-inline">Allan HE, Hayman MA, Marcone S, et al. Proteome
and functional decline as platelets age in the circulation. *J. Thromb.
Haemost.* 2021;19(12):3095–3112. </span>

<span class="csl-left-margin">8. </span><span
class="csl-right-inline">Nebie O, Carvalho K, Barro L, et al. Human
platelet lysate biotherapy for traumatic brain injury: Preclinical
assessment. *Brain*. 2021;144(10):3142–3158. </span>

<span class="csl-left-margin">9. </span><span
class="csl-right-inline">Mirlashari MR, Vetlesen A, Nissen-Meyer LSH, et
al. Proteomic study of apheresis platelets made HLA class I deficient
for transfusion of refractory patients. *Proteomics Clin. Appl.*
2021;15(6):e2100022. </span>

<span class="csl-left-margin">10. </span><span
class="csl-right-inline">Bergemalm D, Ramström S, Kardeby C, et al.
Platelet proteome and function in x-linked thrombocytopenia with
thalassemia and in silico comparisons with gray platelet syndrome.
*Haematologica*. 2021;106(11):2947–2959. </span>

<span class="csl-left-margin">11. </span><span
class="csl-right-inline">Van Bergen MGJM, Marneth AE, Hoogendijk AJ, et
al. Specific proteome changes in platelets from individuals with GATA1-,
GFI1B-, and <span class="nocase">RUNX1-linked</span> bleeding disorders.
*Blood*. 2021;138(1):86–90. </span>

<span class="csl-left-margin">12. </span><span
class="csl-right-inline">Bianchetti A, Chinello C, Guindani M, et al. A
blood bank standardized production of human platelet lysate for
mesenchymal stromal cell expansion: Proteomic characterization and
biological effects. *Front. Cell Dev. Biol.* 2021;9:650490. </span>

<span class="csl-left-margin">13. </span><span
class="csl-right-inline">Bhat A, Das S, Yadav G, et al. Hyperoxidized
albumin modulates platelets and promotes inflammation through CD36
receptor in severe alcoholic hepatitis. *Hepatol. Commun.*
2020;4(1):50–65. </span>

<span class="csl-left-margin">14. </span><span
class="csl-right-inline">Sereni L, Castiello MC, Marangoni F, et al.
Autonomous role of Wiskott-Aldrich syndrome platelet deficiency in
inducing autoimmunity and inflammation. *J. Allergy Clin. Immunol.*
2018;142(4):1272–1284. </span>

<span class="csl-left-margin">15. </span><span
class="csl-right-inline">Salunkhe V, De Cuyper IM, Papadopoulos P, et
al. A comprehensive proteomics study on platelet concentrates: Platelet
proteome, storage time and mirasol pathogen reduction technology.
*Platelets*. 2019;30(3):368–379. </span>

<span class="csl-left-margin">16. </span><span
class="csl-right-inline">Nassa G, Giurato G, Cimmino G, et al. Splicing
of platelet resident <span class="nocase">pre-mRNAs</span> upon
activation by physiological stimuli results in functionally relevant
proteome modifications. *Sci. Rep.* 2018;8(1):498. </span>

<span class="csl-left-margin">17. </span><span
class="csl-right-inline">Rijkers M, Eshof BL van den, Meer PF van der,
et al. Monitoring storage induced changes in the platelet proteome
employing label free quantitative mass spectrometry. *Sci. Rep.*
2017;7(1):11045. </span>

<span class="csl-left-margin">18. </span><span
class="csl-right-inline">Trugilho MR de O, Hottz ED, Brunoro GVF, et al.
Platelet proteome reveals novel pathways of platelet activation and
platelet-mediated immunoregulation in dengue. *PLoS Pathog.*
2017;13(5):e1006385. </span>

<span class="csl-left-margin">19. </span><span
class="csl-right-inline">Lee H, Chae S, Park J, et al. Comprehensive
proteome profiling of platelet identified a protein profile predictive
of responses to an antiplatelet agent sarpogrelate. *Mol. Cell.
Proteomics*. 2016;15(11):3461–3472. </span>

<span class="csl-left-margin">20. </span><span
class="csl-right-inline">Thiele T, Braune J, Dhople V, et al. Proteomic
profile of platelets during reconstitution of platelet counts after
apheresis. *Proteomics Clin. Appl.* 2016;10(8):831–838. </span>

<span class="csl-left-margin">21. </span><span
class="csl-right-inline">Cimmino G, Tarallo R, Nassa G, et al.
Activating stimuli induce platelet <span class="nocase">microRNA</span>
modulation and proteome reorganisation. *Thromb. Haemost.*
2015;114(1):96–108. </span>

<span class="csl-left-margin">22. </span><span
class="csl-right-inline">Zeiler M, Moser M, Mann M. Copy number analysis
of the murine platelet proteome spanning the complete abundance range.
*Mol. Cell. Proteomics*. 2014;13(12):3435–3445. </span>

<span class="csl-left-margin">23. </span><span
class="csl-right-inline">Burkhart JM, Vaudel M, Gambaryan S, et al. The
first comprehensive and quantitative analysis of human platelet protein
composition allows the comparative analysis of structural and functional
pathways. *Blood*. 2012;120(15):e73–82. </span>

<span class="csl-left-margin">24. </span><span
class="csl-right-inline">Wright B, Stanley RG, Kaiser WJ, Mills DJ,
Gibbins JM. Analysis of protein networks in resting and collagen
receptor (<span class="nocase">GPVI)-stimulated</span> platelet
sub-proteomes. *Proteomics*. 2011;11(23):4588–4592. </span>

<span class="csl-left-margin">25. </span><span
class="csl-right-inline">Martı́nez-Botı́a P, Meinders M, De Cuyper IM, et
al. Dissecting platelet proteomics to understand the pathophysiology of
immune thrombocytopenia: Studies in mouse models. *Blood Adv.*
2022;6(11):3529–3534. </span>

<span class="csl-left-margin">26. </span><span
class="csl-right-inline">Solari FA, Mattheij NJA, Burkhart JM, et al.
Combined quantification of the global proteome, phosphoproteome, and
proteolytic cleavage to characterize altered platelet functions in the
human scott syndrome. *Mol. Cell. Proteomics*. 2016;15(10):3154–3169.
</span>

<span class="csl-left-margin">27. </span><span
class="csl-right-inline">Loroch S, Trabold K, Gambaryan S, et al.
Alterations of the platelet proteome in type I glanzmann thrombasthenia
caused by different homozygous delG frameshift mutations in ITGA2B.
*Thromb. Haemost.* 2017;117(3):556–569. </span>

<span class="csl-left-margin">28. </span><span
class="csl-right-inline">Stokhuijzen E, Koornneef JM, Nota B, et al.
Differences between platelets derived from neonatal cord blood and adult
peripheral blood assessed by mass spectrometry. *J. Proteome Res.*
2017;16(10):3567–3575. </span>

<span class="csl-left-margin">29. </span><span
class="csl-right-inline">Shah P, Yang W, Sun S, et al. Platelet
glycoproteins associated with aspirin-treatment upon platelet
activation. *Proteomics*. 2017;17(6):1600199. </span>

<span class="csl-left-margin">30. </span><span
class="csl-right-inline">Holten TC van, Bleijerveld OB, Wijten P, et al.
Quantitative proteomics analysis reveals similar release profiles
following specific PAR-1 or PAR-4 stimulation of platelets. *Cardiovasc.
Res.* 2014;103(1):140–146. </span>

<span class="csl-left-margin">31. </span><span
class="csl-right-inline">Au AE-L, Sashindranath M, Borg RJ, et al.
Activated platelets rescue apoptotic cells via paracrine activation of
EGFR and <span class="nocase">DNA-dependent</span> protein kinase. *Cell
Death Dis.* 2014;5(9):e1410. </span>

<span class="csl-left-margin">32. </span><span
class="csl-right-inline">Servais L, Wéra O, Dibato Epoh J, et al.
Platelets contribute to the initiation of colitis‐associated cancer by
promoting immunosuppression. *J. Thromb. Haemost.* 2018;16(4):762–777.
</span>

<span class="csl-left-margin">33. </span><span
class="csl-right-inline">Parsons MEM, Szklanna PB, Guerrero JA, et al.
Platelet releasate proteome profiling reveals a core set of proteins
with low variance between healthy adults. *Proteomics*.
2018;18(15):1800219. </span>

<span class="csl-left-margin">34. </span><span
class="csl-right-inline">Szklanna PB, Parsons ME, Wynne K, et al. The
platelet releasate is altered in human pregnancy. *Proteomics Clin.
Appl.* 2019;13(3):e1800162. </span>

<span class="csl-left-margin">35. </span><span
class="csl-right-inline">Rowley JW, Oler AJ, Tolley ND, et al.
Genome-wide <span class="nocase">RNA-seq</span> analysis of human and
mouse platelet transcriptomes. *Blood*. 2011;118(14):e101–11. </span>

<span class="csl-left-margin">36. </span><span
class="csl-right-inline">Wickham H, François R, Henry L, Müller K.
[Dplyr: A grammar of data
manipulation](https://CRAN.R-project.org/package=dplyr). 2022. </span>

<span class="csl-left-margin">37. </span><span
class="csl-right-inline">Wickham H. [Stringr: Simple, consistent
wrappers for common string
operations](https://CRAN.R-project.org/package=stringr). 2019. </span>

<span class="csl-left-margin">38. </span><span
class="csl-right-inline">Wickham H. [ggplot2: Elegant graphics for data
analysis](https://ggplot2.tidyverse.org). Springer-Verlag New York;
2016. </span>

<span class="csl-left-margin">39. </span><span
class="csl-right-inline">Kassambara A. [Ggpubr: ’ggplot2’ based
publication ready plots](https://CRAN.R-project.org/package=ggpubr).
2020. </span>

<span class="csl-left-margin">40. </span><span
class="csl-right-inline">Larsson J. [<span class="nocase">eulerr</span>:
Area-proportional euler and venn diagrams with
ellipses](https://CRAN.R-project.org/package=eulerr). 2022. </span>

<span class="csl-left-margin">41. </span><span
class="csl-right-inline">Pagès H, Carlson M, Falcon S, Li N.
[AnnotationDbi: Manipulation of SQLite-based annotations in
bioconductor](https://bioconductor.org/packages/AnnotationDbi). 2020.
</span>

<span class="csl-left-margin">42. </span><span
class="csl-right-inline">Carlson M. Org.hs.eg.db: Genome wide annotation
for human. 2020. </span>

<span class="csl-left-margin">43. </span><span
class="csl-right-inline">Carlson M. Org.mm.eg.db: Genome wide annotation
for mouse. 2020. </span>

<span class="csl-left-margin">44. </span><span
class="csl-right-inline">Kolberg L, Raudvere U, Kuzmin I, Vilo J,
Peterson H. gprofiler2– an r package for gene list functional enrichment
analysis and namespace conversion toolset g:profiler. *F1000Research*.
2020;9 (ELIXIR)(709): </span>

<span class="csl-left-margin">45. </span><span
class="csl-right-inline">Yu G, Wang L-G, Han Y, He Q-Y.
[clusterProfiler: An r package for comparing biological themes among
gene clusters](https://doi.org/10.1089/omi.2011.0118). *OMICS: A Journal
of Integrative Biology*. 2012;16(5):284–287. </span>
