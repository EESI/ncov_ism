# Informative Subtype Marker (ISM) Analysis Report
Informative Subtype Marker (ISM) is an efficient framework for genetic subtyping of a pandemic virus and implement it for SARS-CoV-2, the novel coronavirus that causes COVID-19.        
Drexel University EESI Lab, 2020        
Maintainer: Zhengqiao Zhao, zz374 at drexel dot edu  
Owner: Gail Rosen, gailr at ece dot drexel dot edu  

**Report created on 2021/02/02**
<!--- dividing line --->

## Abstract
The novel coronavirus responsible for COVID-19, SARS-CoV-2, expanded to reportedly 8.7 million confirmed cases worldwide by June 21, 2020. The global SARS-CoV-2 pandemic highlights the importance of tracking viral transmission dynamics in real-time. Through June 2020, researchers have obtained genetic sequences of SARS-CoV-2 from over 50 thousand samples from infected individuals worldwide. Since the virus readily mutates, each sequence of an infected individual contains useful information linked to the individual's exposure location and sample date. But, there are over 30,000 bases in the full SARS-CoV-2 genome, so tracking genetic variants on a whole-sequence basis becomes unwieldy. *ncov_ism* is a method to instead efficiently identify and label genetic variants, or "subtypes" of SARS-CoV-2. This method defines a compact set of nucleotide sites that characterize the most variable (and thus most informative) positions in the viral genomes sequenced from different individuals, called an Informative Subtype Marker or *ISM*. This tool defines viral subtypes for each ISM, and analyze the regional distribution of subtypes to track the progress of the pandemic.

## Entropy time series analysis
The following figure shows how entropy values at different positions change over time.

<img src="results/0_Entropy_time_series_analysis.png" alt="ents" width="800"/>

A few covarying positions are identified

<!--- covarying table starts --->
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>Covarying group</th>
      <th>Coverage</th>
      <th>NT configurations</th>
      <th>Representative position</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>28881;28882;28883</td>
      <td>0.991934</td>
      <td>GGG;AAC</td>
      <td>28881</td>
    </tr>
    <tr>
      <td>445;6286;21255;22227;26801;28932;29645</td>
      <td>0.947184</td>
      <td>TCGCCCG;CTCTGTT</td>
      <td>26801</td>
    </tr>
    <tr>
      <td>204;21614;27944</td>
      <td>0.898983</td>
      <td>GCC;TTT</td>
      <td>27944</td>
    </tr>
    <tr>
      <td>913;3267;5388;5986;6954;14676;15279;16176;2306...</td>
      <td>0.901591</td>
      <td>CCCCTCCTACCCCTGCGAGATC;TTATCTTCTAATCGCTTGCTAT</td>
      <td>23604</td>
    </tr>
    <tr>
      <td>1163;7540;16647;18555;22992;23401</td>
      <td>0.938080</td>
      <td>ATGCGG;TCTTAA</td>
      <td>22992</td>
    </tr>
    <tr>
      <td>10097;23731</td>
      <td>0.996828</td>
      <td>GC;AT</td>
      <td>10097</td>
    </tr>
    <tr>
      <td>18877;27964;28854</td>
      <td>0.889275</td>
      <td>CCC;CTC</td>
      <td>28854</td>
    </tr>
    <tr>
      <td>8782;28144</td>
      <td>0.993198</td>
      <td>CT;TC</td>
      <td>8782</td>
    </tr>
    <tr>
      <td>241;3037;14408;23403</td>
      <td>0.973077</td>
      <td>TTTG;CCCA</td>
      <td>23403</td>
    </tr>
  </tbody>
</table>
<!--- covarying table ends --->
<!--- dividing line --->

## ISM positions
The following table shows the annotations of ISM sites using the reference viral genome.

<!--- annotation table starts --->
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>Ref position</th>
      <th>Entropy</th>
      <th>Gene</th>
      <th>Is silent</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>204</td>
      <td>0.600872</td>
      <td>Non-coding</td>
      <td>True</td>
    </tr>
    <tr>
      <td>241</td>
      <td>0.327798</td>
      <td>Non-coding</td>
      <td>True</td>
    </tr>
    <tr>
      <td>445</td>
      <td>0.800201</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>True</td>
    </tr>
    <tr>
      <td>913</td>
      <td>0.362935</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>True</td>
    </tr>
    <tr>
      <td>1059</td>
      <td>0.627709</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>False</td>
    </tr>
    <tr>
      <td>1163</td>
      <td>0.262509</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>False</td>
    </tr>
    <tr>
      <td>3037</td>
      <td>0.330777</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>True</td>
    </tr>
    <tr>
      <td>3267</td>
      <td>0.367061</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>False</td>
    </tr>
    <tr>
      <td>5388</td>
      <td>0.361585</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>False</td>
    </tr>
    <tr>
      <td>5986</td>
      <td>0.374293</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>True</td>
    </tr>
    <tr>
      <td>6286</td>
      <td>0.804308</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>True</td>
    </tr>
    <tr>
      <td>6954</td>
      <td>0.363575</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>False</td>
    </tr>
    <tr>
      <td>7540</td>
      <td>0.197486</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>True</td>
    </tr>
    <tr>
      <td>8782</td>
      <td>0.151998</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>True</td>
    </tr>
    <tr>
      <td>10097</td>
      <td>0.165350</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>False</td>
    </tr>
    <tr>
      <td>11083</td>
      <td>0.352021</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>False</td>
    </tr>
    <tr>
      <td>14408</td>
      <td>0.333818</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>False</td>
    </tr>
    <tr>
      <td>14676</td>
      <td>0.368721</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>True</td>
    </tr>
    <tr>
      <td>14805</td>
      <td>0.228778</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>True</td>
    </tr>
    <tr>
      <td>15279</td>
      <td>0.364136</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>True</td>
    </tr>
    <tr>
      <td>16176</td>
      <td>0.361970</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>True</td>
    </tr>
    <tr>
      <td>16647</td>
      <td>0.205341</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>True</td>
    </tr>
    <tr>
      <td>18555</td>
      <td>0.213862</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>True</td>
    </tr>
    <tr>
      <td>18877</td>
      <td>0.327941</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>True</td>
    </tr>
    <tr>
      <td>20268</td>
      <td>0.390456</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>True</td>
    </tr>
    <tr>
      <td>21255</td>
      <td>0.816114</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>True</td>
    </tr>
    <tr>
      <td>21614</td>
      <td>0.504229</td>
      <td>YP_009724390.1: surface glycoprotein</td>
      <td>False</td>
    </tr>
    <tr>
      <td>22227</td>
      <td>0.803554</td>
      <td>YP_009724390.1: surface glycoprotein</td>
      <td>False</td>
    </tr>
    <tr>
      <td>22992</td>
      <td>0.325385</td>
      <td>YP_009724390.1: surface glycoprotein</td>
      <td>False</td>
    </tr>
    <tr>
      <td>23063</td>
      <td>0.379714</td>
      <td>YP_009724390.1: surface glycoprotein</td>
      <td>False</td>
    </tr>
    <tr>
      <td>23271</td>
      <td>0.364818</td>
      <td>YP_009724390.1: surface glycoprotein</td>
      <td>False</td>
    </tr>
    <tr>
      <td>23401</td>
      <td>0.215601</td>
      <td>YP_009724390.1: surface glycoprotein</td>
      <td>True</td>
    </tr>
    <tr>
      <td>23403</td>
      <td>0.328119</td>
      <td>YP_009724390.1: surface glycoprotein</td>
      <td>False</td>
    </tr>
    <tr>
      <td>23604</td>
      <td>0.398903</td>
      <td>YP_009724390.1: surface glycoprotein</td>
      <td>False</td>
    </tr>
    <tr>
      <td>23709</td>
      <td>0.364862</td>
      <td>YP_009724390.1: surface glycoprotein</td>
      <td>False</td>
    </tr>
    <tr>
      <td>23731</td>
      <td>0.147387</td>
      <td>YP_009724390.1: surface glycoprotein</td>
      <td>True</td>
    </tr>
    <tr>
      <td>24334</td>
      <td>0.276680</td>
      <td>YP_009724390.1: surface glycoprotein</td>
      <td>True</td>
    </tr>
    <tr>
      <td>24506</td>
      <td>0.362065</td>
      <td>YP_009724390.1: surface glycoprotein</td>
      <td>False</td>
    </tr>
    <tr>
      <td>24914</td>
      <td>0.366751</td>
      <td>YP_009724390.1: surface glycoprotein</td>
      <td>False</td>
    </tr>
    <tr>
      <td>25563</td>
      <td>0.803961</td>
      <td>YP_009724391.1: ORF3a protein</td>
      <td>False</td>
    </tr>
    <tr>
      <td>26735</td>
      <td>0.296868</td>
      <td>YP_009724393.1: membrane glycoprotein</td>
      <td>True</td>
    </tr>
    <tr>
      <td>26801</td>
      <td>0.845113</td>
      <td>YP_009724393.1: membrane glycoprotein</td>
      <td>True</td>
    </tr>
    <tr>
      <td>27944</td>
      <td>0.654932</td>
      <td>YP_009724396.1: ORF8 protein</td>
      <td>True</td>
    </tr>
    <tr>
      <td>27964</td>
      <td>0.320243</td>
      <td>YP_009724396.1: ORF8 protein</td>
      <td>False</td>
    </tr>
    <tr>
      <td>27972</td>
      <td>0.363867</td>
      <td>YP_009724396.1: ORF8 protein</td>
      <td>False</td>
    </tr>
    <tr>
      <td>28048</td>
      <td>0.361900</td>
      <td>YP_009724396.1: ORF8 protein</td>
      <td>False</td>
    </tr>
    <tr>
      <td>28111</td>
      <td>0.361117</td>
      <td>YP_009724396.1: ORF8 protein</td>
      <td>False</td>
    </tr>
    <tr>
      <td>28144</td>
      <td>0.150106</td>
      <td>YP_009724396.1: ORF8 protein</td>
      <td>False</td>
    </tr>
    <tr>
      <td>28280</td>
      <td>0.368327</td>
      <td>YP_009724397.2: nucleocapsid phosphoprotein</td>
      <td>False</td>
    </tr>
    <tr>
      <td>28281</td>
      <td>0.360936</td>
      <td>YP_009724397.2: nucleocapsid phosphoprotein</td>
      <td>False</td>
    </tr>
    <tr>
      <td>28282</td>
      <td>0.362545</td>
      <td>YP_009724397.2: nucleocapsid phosphoprotein</td>
      <td>False</td>
    </tr>
    <tr>
      <td>28854</td>
      <td>0.348248</td>
      <td>YP_009724397.2: nucleocapsid phosphoprotein</td>
      <td>False</td>
    </tr>
    <tr>
      <td>28881</td>
      <td>0.940754</td>
      <td>YP_009724397.2: nucleocapsid phosphoprotein</td>
      <td>False</td>
    </tr>
    <tr>
      <td>28882</td>
      <td>0.936498</td>
      <td>YP_009724397.2: nucleocapsid phosphoprotein</td>
      <td>True</td>
    </tr>
    <tr>
      <td>28883</td>
      <td>0.936610</td>
      <td>YP_009724397.2: nucleocapsid phosphoprotein</td>
      <td>False</td>
    </tr>
    <tr>
      <td>28932</td>
      <td>0.800105</td>
      <td>YP_009724397.2: nucleocapsid phosphoprotein</td>
      <td>False</td>
    </tr>
    <tr>
      <td>28977</td>
      <td>0.373717</td>
      <td>YP_009724397.2: nucleocapsid phosphoprotein</td>
      <td>False</td>
    </tr>
    <tr>
      <td>29645</td>
      <td>0.802626</td>
      <td>YP_009725255.1: ORF10 protein</td>
      <td>False</td>
    </tr>
  </tbody>
</table>
<!--- annotation table ends --->
<!--- dividing line --->

## ISM distribution worldwide
The following figure shows the major ISMs in selective countries/regions (in the legend next to each
country/region, we show the date when a major ISM was first sequenced in that country/region). 
ISMs with less than 5% abundance are plotted as “OTHER”. 

<img src="results/1_regional_ISM.png" alt="regional" width="800"/>

<!--- ![Fig 1](results/1_regional_ISM.png "Subtype composition in different locations worldwide") --->

## ISM distribution in US
The following figure shows the ISM distribution in the United States in 25 states. ISMs with less than 5% abundance are plotted as “OTHER”. 

<img src="results/2_intra-US_ISM.png" alt="states" width="800"/>

<!--- ![Fig 2](results/2_intra-US_ISM.png "Subtype composition in different locations in US") --->

## ISM abundance PCA analysis
ISM abundance table is constructed for regions with more than 150 submissions. We then visualize the pattern of viral genetic variation by the first two principle components of the pairwise Bray-Curtis dissimilarity matrix between regions. Regions with similar genetic variation pattern are grouped together in the PCA plot.

<img src="results/country_2d.png" alt="pca" width="800"/>


## The dynamic of ISM in different locations
1. The relative abundance (%) of ISMs in DNA sequences from the United Kingdom as sampled over time.

<img src="results/3_ISM_growth_United Kingdom.png" alt="UK" width="800"/>

2. The relative abundance (%) of ISMs in DNA sequences from France as sampled over time.

<img src="results/3_ISM_growth_France.png" alt="France" width="800"/>

3. The relative abundance (%) of ISMs in DNA sequences from the USA as sampled over time.

<img src="results/3_ISM_growth_USA.png" alt="USA" width="800"/>
<!--- ![Fig 3](results/3_ISM_growth_USA.png "the dynamic subtype composition in US over time") --->

## Reference
Research article in PLOS COMPUTATIONAL BIOLOGY [Genetic Grouping of SARS-CoV-2 Coronavirus Sequences using Informative Subtype Markers for Pandemic Spread Visualization](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008269).
If you find our work helpful, please cite:
```
@article{10.1371/journal.pcbi.1008269,
    author = {Zhao, Zhengqiao AND Sokhansanj, Bahrad A. AND Malhotra, Charvi AND Zheng, Kitty AND Rosen, Gail L.},
    journal = {PLOS Computational Biology},
    publisher = {Public Library of Science},
    title = {Genetic grouping of SARS-CoV-2 coronavirus sequences using informative subtype markers for pandemic spread visualization},
    year = {2020},
    month = {09},
    volume = {16},
    url = {https://doi.org/10.1371/journal.pcbi.1008269},
    pages = {1-32},
    number = {9},
    doi = {10.1371/journal.pcbi.1008269}
}
```
## Acknowledgement
We would like to thank [GISAID](www.gisaid.org) for sharing the sequence data and metadata. We also gratefully acknowledge the authors, originating and submitting laboratories of the sequences from GISAID’s EpiFlu Database on which this research is based. The list is detailed in [here](results/acknowledgement_table.txt). All submitters of data may be contacted directly via the [GISAID](www.gisaid.org) website.
