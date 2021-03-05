# Informative Subtype Marker (ISM) Analysis Report
Informative Subtype Marker (ISM) is an efficient framework for genetic subtyping of a pandemic virus and implement it for SARS-CoV-2, the novel coronavirus that causes COVID-19.        
Drexel University EESI Lab, 2020        
Maintainer: Zhengqiao Zhao, zz374 at drexel dot edu  
Owner: Gail Rosen, gailr at ece dot drexel dot edu  

**Report created on 2021/03/04**
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
      <th>NT configurations</th>
      <th>Coverage</th>
      <th>Representative position</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>28881;28882;28883</td>
      <td>GGG;AAC</td>
      <td>0.992289</td>
      <td>28881</td>
    </tr>
    <tr>
      <td>445;6286;21255;22227;26801;28932;29645</td>
      <td>TCGCCCG;CTCTGTT</td>
      <td>0.960428</td>
      <td>26801</td>
    </tr>
    <tr>
      <td>204;21614;27944</td>
      <td>GCC;TTT</td>
      <td>0.901312</td>
      <td>27944</td>
    </tr>
    <tr>
      <td>913;3267;5388;5986;6954;14676;15279;16176;23063;23271;23604;23709;24506;24914;27972;28048;28111;28280;28281;28282;28977</td>
      <td>CCCCTCCTACCCTGCGAGATC;TTATCTTCTAATGCTTGCTAT</td>
      <td>0.933089</td>
      <td>23604</td>
    </tr>
    <tr>
      <td>18877;27964;28854</td>
      <td>CCC;CTC</td>
      <td>0.887817</td>
      <td>27964</td>
    </tr>
    <tr>
      <td>17615;28095</td>
      <td>AA;GA</td>
      <td>0.958332</td>
      <td>17615</td>
    </tr>
    <tr>
      <td>10319;18424;21304;24334;25907;26735;28472;28869;28975</td>
      <td>CACCGCCCG;TGTCTCTTG</td>
      <td>0.850017</td>
      <td>28869</td>
    </tr>
    <tr>
      <td>10097;23731</td>
      <td>GC;AT</td>
      <td>0.996450</td>
      <td>10097</td>
    </tr>
    <tr>
      <td>1163;7540;16647;18555;22992;23401</td>
      <td>ATGCGG;ATGCAG</td>
      <td>0.955938</td>
      <td>22992</td>
    </tr>
    <tr>
      <td>8782;28144</td>
      <td>CT;TC</td>
      <td>0.993449</td>
      <td>8782</td>
    </tr>
    <tr>
      <td>241;3037;14408;23403</td>
      <td>TTTG;CCCA</td>
      <td>0.973927</td>
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
      <th>AA position</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>204</td>
      <td>0.527503</td>
      <td>Non-coding</td>
      <td>True</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>241</td>
      <td>0.252706</td>
      <td>Non-coding</td>
      <td>True</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>445</td>
      <td>0.729446</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>True</td>
      <td>V60V</td>
    </tr>
    <tr>
      <td>913</td>
      <td>0.694239</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>True</td>
      <td>S216S</td>
    </tr>
    <tr>
      <td>1059</td>
      <td>0.676355</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>False</td>
      <td>T265I</td>
    </tr>
    <tr>
      <td>1163</td>
      <td>0.193912</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>False</td>
      <td>I300F</td>
    </tr>
    <tr>
      <td>3037</td>
      <td>0.253332</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>True</td>
      <td>F924F</td>
    </tr>
    <tr>
      <td>3267</td>
      <td>0.699064</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>False</td>
      <td>T1001I</td>
    </tr>
    <tr>
      <td>5388</td>
      <td>0.698970</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>False</td>
      <td>A1708D</td>
    </tr>
    <tr>
      <td>5986</td>
      <td>0.702848</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>True</td>
      <td>F1907F</td>
    </tr>
    <tr>
      <td>6286</td>
      <td>0.734245</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>True</td>
      <td>T2007T</td>
    </tr>
    <tr>
      <td>6954</td>
      <td>0.697813</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>False</td>
      <td>I2230T</td>
    </tr>
    <tr>
      <td>7540</td>
      <td>0.142399</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>True</td>
      <td>T2425T</td>
    </tr>
    <tr>
      <td>8782</td>
      <td>0.118147</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>True</td>
      <td>S2839S</td>
    </tr>
    <tr>
      <td>10097</td>
      <td>0.136113</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>False</td>
      <td>G3278S</td>
    </tr>
    <tr>
      <td>10319</td>
      <td>0.377448</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>False</td>
      <td>L3352F</td>
    </tr>
    <tr>
      <td>11083</td>
      <td>0.293754</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>False</td>
      <td>L3606F</td>
    </tr>
    <tr>
      <td>14408</td>
      <td>0.258560</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>False</td>
      <td>P314L</td>
    </tr>
    <tr>
      <td>14676</td>
      <td>0.699718</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>True</td>
      <td>P403P</td>
    </tr>
    <tr>
      <td>14805</td>
      <td>0.264886</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>True</td>
      <td>Y446Y</td>
    </tr>
    <tr>
      <td>15279</td>
      <td>0.697850</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>True</td>
      <td>H604H</td>
    </tr>
    <tr>
      <td>16176</td>
      <td>0.696308</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>True</td>
      <td>T903T</td>
    </tr>
    <tr>
      <td>16647</td>
      <td>0.151641</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>True</td>
      <td>T1060T</td>
    </tr>
    <tr>
      <td>17615</td>
      <td>0.335089</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>False</td>
      <td>K1383R</td>
    </tr>
    <tr>
      <td>18424</td>
      <td>0.350202</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>False</td>
      <td>N1653D</td>
    </tr>
    <tr>
      <td>18555</td>
      <td>0.164434</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>True</td>
      <td>D1696D</td>
    </tr>
    <tr>
      <td>18877</td>
      <td>0.321750</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>True</td>
      <td>L1804L</td>
    </tr>
    <tr>
      <td>20268</td>
      <td>0.372304</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>True</td>
      <td>L2267L</td>
    </tr>
    <tr>
      <td>21255</td>
      <td>0.747558</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>True</td>
      <td>A2596A</td>
    </tr>
    <tr>
      <td>21304</td>
      <td>0.372334</td>
      <td>YP_009724389.1: ORF1ab polyprotein</td>
      <td>False</td>
      <td>R2613C</td>
    </tr>
    <tr>
      <td>21614</td>
      <td>0.457142</td>
      <td>YP_009724390.1: surface glycoprotein</td>
      <td>False</td>
      <td>L18F</td>
    </tr>
    <tr>
      <td>22227</td>
      <td>0.733969</td>
      <td>YP_009724390.1: surface glycoprotein</td>
      <td>False</td>
      <td>A222V</td>
    </tr>
    <tr>
      <td>22992</td>
      <td>0.296281</td>
      <td>YP_009724390.1: surface glycoprotein</td>
      <td>False</td>
      <td>S477N</td>
    </tr>
    <tr>
      <td>23063</td>
      <td>0.714322</td>
      <td>YP_009724390.1: surface glycoprotein</td>
      <td>False</td>
      <td>N501Y</td>
    </tr>
    <tr>
      <td>23271</td>
      <td>0.699519</td>
      <td>YP_009724390.1: surface glycoprotein</td>
      <td>False</td>
      <td>A570D</td>
    </tr>
    <tr>
      <td>23401</td>
      <td>0.165068</td>
      <td>YP_009724390.1: surface glycoprotein</td>
      <td>True</td>
      <td>Q613Q</td>
    </tr>
    <tr>
      <td>23403</td>
      <td>0.248856</td>
      <td>YP_009724390.1: surface glycoprotein</td>
      <td>False</td>
      <td>D614G</td>
    </tr>
    <tr>
      <td>23604</td>
      <td>0.743056</td>
      <td>YP_009724390.1: surface glycoprotein</td>
      <td>False</td>
      <td>P681H</td>
    </tr>
    <tr>
      <td>23709</td>
      <td>0.699827</td>
      <td>YP_009724390.1: surface glycoprotein</td>
      <td>False</td>
      <td>T716I</td>
    </tr>
    <tr>
      <td>23731</td>
      <td>0.116574</td>
      <td>YP_009724390.1: surface glycoprotein</td>
      <td>True</td>
      <td>T723T</td>
    </tr>
    <tr>
      <td>24334</td>
      <td>0.225993</td>
      <td>YP_009724390.1: surface glycoprotein</td>
      <td>True</td>
      <td>A924A</td>
    </tr>
    <tr>
      <td>24506</td>
      <td>0.696112</td>
      <td>YP_009724390.1: surface glycoprotein</td>
      <td>False</td>
      <td>S982A</td>
    </tr>
    <tr>
      <td>24914</td>
      <td>0.701545</td>
      <td>YP_009724390.1: surface glycoprotein</td>
      <td>False</td>
      <td>D1118H</td>
    </tr>
    <tr>
      <td>25563</td>
      <td>0.825306</td>
      <td>YP_009724391.1: ORF3a protein</td>
      <td>False</td>
      <td>Q57H</td>
    </tr>
    <tr>
      <td>25907</td>
      <td>0.349064</td>
      <td>YP_009724391.1: ORF3a protein</td>
      <td>False</td>
      <td>G172V</td>
    </tr>
    <tr>
      <td>26735</td>
      <td>0.294816</td>
      <td>YP_009724393.1: membrane glycoprotein</td>
      <td>True</td>
      <td>Y71Y</td>
    </tr>
    <tr>
      <td>26801</td>
      <td>0.778550</td>
      <td>YP_009724393.1: membrane glycoprotein</td>
      <td>True</td>
      <td>L93L</td>
    </tr>
    <tr>
      <td>27944</td>
      <td>0.579214</td>
      <td>YP_009724396.1: ORF8 protein</td>
      <td>True</td>
      <td>H17H</td>
    </tr>
    <tr>
      <td>27964</td>
      <td>0.417859</td>
      <td>YP_009724396.1: ORF8 protein</td>
      <td>False</td>
      <td>S24L</td>
    </tr>
    <tr>
      <td>27972</td>
      <td>0.699151</td>
      <td>YP_009724396.1: ORF8 protein</td>
      <td>False</td>
      <td>Q27_</td>
    </tr>
    <tr>
      <td>28048</td>
      <td>0.698559</td>
      <td>YP_009724396.1: ORF8 protein</td>
      <td>False</td>
      <td>R52I</td>
    </tr>
    <tr>
      <td>28095</td>
      <td>0.225796</td>
      <td>YP_009724396.1: ORF8 protein</td>
      <td>False</td>
      <td>K68_</td>
    </tr>
    <tr>
      <td>28111</td>
      <td>0.695749</td>
      <td>YP_009724396.1: ORF8 protein</td>
      <td>False</td>
      <td>Y73C</td>
    </tr>
    <tr>
      <td>28144</td>
      <td>0.114729</td>
      <td>YP_009724396.1: ORF8 protein</td>
      <td>False</td>
      <td>L84S</td>
    </tr>
    <tr>
      <td>28280</td>
      <td>0.705916</td>
      <td>YP_009724397.2: nucleocapsid phosphoprotein</td>
      <td>False</td>
      <td>D3H</td>
    </tr>
    <tr>
      <td>28281</td>
      <td>0.693491</td>
      <td>YP_009724397.2: nucleocapsid phosphoprotein</td>
      <td>False</td>
      <td>D3V</td>
    </tr>
    <tr>
      <td>28282</td>
      <td>0.693967</td>
      <td>YP_009724397.2: nucleocapsid phosphoprotein</td>
      <td>False</td>
      <td>D3E</td>
    </tr>
    <tr>
      <td>28472</td>
      <td>0.344767</td>
      <td>YP_009724397.2: nucleocapsid phosphoprotein</td>
      <td>False</td>
      <td>P67S</td>
    </tr>
    <tr>
      <td>28854</td>
      <td>0.341757</td>
      <td>YP_009724397.2: nucleocapsid phosphoprotein</td>
      <td>False</td>
      <td>S194L</td>
    </tr>
    <tr>
      <td>28869</td>
      <td>0.405139</td>
      <td>YP_009724397.2: nucleocapsid phosphoprotein</td>
      <td>False</td>
      <td>P199L</td>
    </tr>
    <tr>
      <td>28881</td>
      <td>0.981990</td>
      <td>YP_009724397.2: nucleocapsid phosphoprotein</td>
      <td>False</td>
      <td>R203K</td>
    </tr>
    <tr>
      <td>28882</td>
      <td>0.977984</td>
      <td>YP_009724397.2: nucleocapsid phosphoprotein</td>
      <td>True</td>
      <td>R203R</td>
    </tr>
    <tr>
      <td>28883</td>
      <td>0.979321</td>
      <td>YP_009724397.2: nucleocapsid phosphoprotein</td>
      <td>False</td>
      <td>G204R</td>
    </tr>
    <tr>
      <td>28932</td>
      <td>0.729689</td>
      <td>YP_009724397.2: nucleocapsid phosphoprotein</td>
      <td>False</td>
      <td>A220V</td>
    </tr>
    <tr>
      <td>28975</td>
      <td>0.292876</td>
      <td>YP_009724397.2: nucleocapsid phosphoprotein</td>
      <td>False</td>
      <td>M234I</td>
    </tr>
    <tr>
      <td>28977</td>
      <td>0.702618</td>
      <td>YP_009724397.2: nucleocapsid phosphoprotein</td>
      <td>False</td>
      <td>S235F</td>
    </tr>
    <tr>
      <td>29645</td>
      <td>0.731812</td>
      <td>YP_009725255.1: ORF10 protein</td>
      <td>False</td>
      <td>V30L</td>
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
