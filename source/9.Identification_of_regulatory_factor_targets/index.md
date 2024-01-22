# 9. Identification of Transcriptional Factor Targets

## 9.1. BeReTa


We integrated BeReTa into the ETGEMs framework. This algorithm, by incorporating a transcriptional regulatory network, identifies regulatory factors that have a controlling effect on increasing the flux of the target reaction. This chapter provides a detailed guide on the application of BeReTa within the ETGEMs framework. 

Taking the production of succinic acid by Corynebacterium glutamicum as an example


```python
import pandas as pd
from autoETGEMs.io.json import *
tcmodel = read_json_model("./autoETGEMs/model/cgmodel_test_BeReTa.json")
```

Input transcriptome data and transcriptional regulatory networks.


```python
exp_comp_genes_file='./RETGEMs/exp_comp_genes.csv'
exp_comp_genes = pd.read_csv(exp_comp_genes_file,index_col=0)
exp_comp={}
exp_comp['genes']=exp_comp_genes.index
exp_comp['expression']=exp_comp_genes
reg_net_file='./RETGEMs/reg_net.csv'
reg_net = pd.read_csv(reg_net_file)
```

The BeReTa_v1 function supports four models in the ETGEMs framework, and additional constraints can be optionally added based on the parameters.


```python
from autoETGEMs.analysis.BeReTa import *
result = BeReTa_v1(tcmodel,"ETM", reg_net, exp_comp, "succ_c")
```

The results are in dictionary format. The "reglist" contains transcription factors that exhibit significant regulatory effects on the target reaction. The "metrics" corresponds to the associated parameters.   

(1)Beneficial Score, indicating the extent to which each transcription factor regulates the target reaction. Positive values indicate a promoting effect, while negative values indicate an inhibitory effect.   

(2)P-value, which measures the significance of the Beneficial Score, indicating whether it significantly differs from random conditions.

(3)Effective Genes, representing the count of genes that have an impact on the target reaction under the regulation of the transcription factor. These genes may be associated with metabolic reactions in the model.   

(4)Number of Effective Reactions, indicating the count of reactions that are influenced by the target reaction under the regulation of the transcription factor. These reactions may interact with other parts of the metabolic network.

(5)Fraction of Effective Genes, representing the proportion of genes that have an impact on the target reaction under the regulation of the transcription factor, relative to all genes associated with that transcription factor. This column reflects how many genes in this regulatory network are related to the regulation of the target reaction.


```python
result["reglist"]
```




    ['Cgl1919', 'Cgl1931', 'Cgl0291', 'Cgl0369', 'Cgl1541', 'Cgl1920', 'Cgl3089']




```python
result["metrics"]
```




    [array([4.82300374, 0.        , 8.        , 8.        , 0.33333333]),
     array([-4.0713046 ,  0.        ,  9.        , 10.        ,  0.17647059]),
     array([-3.21287011e+00,  1.67000000e-02,  1.70000000e+01,  1.50000000e+01,
             1.22302158e-01]),
     array([-1.85285556,  0.0132    ,  6.        ,  6.        ,  0.16216216]),
     array([-0.44300455,  0.043     ,  1.        ,  2.        ,  0.25      ]),
     array([0.38487191, 0.0372    , 7.        , 4.        , 0.5       ]),
     array([0.03473783, 0.0439    , 2.        , 1.        , 1.        ])]


