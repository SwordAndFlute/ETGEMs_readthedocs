# 6. Bottleneck Reactions

Using ETGEMs, the identification of thermodynamic bottleneck reactions can be easily accomplished. This chapter primarily outlines the workflow for discovering bottleneck reactions in ETGEMs, as well as exploring potential coupling strategies targeted at these bottleneck reactions.


```python
from autoETGEMs.io.json import *
from autoETGEMs.analysis.flux_balance_analysis import *
import pandas as pd
import numpy as np
import datetime
from pyomo.environ import *
from pyomo.opt import SolverFactory
import pyomo.environ as pyo
from concurrent.futures import ProcessPoolExecutor, as_completed
from autoETGEMs.analysis.max_driving_force import *
from autoETGEMs.analysis.metabolite_concentration_variability_analysis import *
from autoETGEMs.analysis.enzyme_cost_variability_analysis import *
from autoETGEMs.plotting.plotting import *
from autoETGEMs.core.reaction import *
from autoETGEMs.analysis.coupling_strategy import *
from autoETGEMs.analysis.continuous_reaction import *
norm_model = read_json_model("./autoETGEMs/model/iML1515_etcmodel_20231108.json")
```

Using the simulation of arginine production as an example, the computational objective of ETGEMs is first modified to maximize the production of Serine.


```python
norm_model.add_demand_reaction("ser__L_c")
norm_model.obj_name = "DM_ser__L_c"
```

## 6.1. Search for Bottleneck Reactions

### 6.1.1. Calculate MDF and determine the B value (inflection point):  

(1) Determine the maximum product yield under no thermodynamic constraints using enzyme-constrained models and metabolic network models;  

(2) Starting from a minimum product yield of 0.1, calculate the MDF at this point. Then, compute the maximum yield at this MDF, which represents the first inflection point. Proceed by adding a very small value (0.02) to this yield, and repeat the above process until it exceeds the maximum value obtained in (1);  

(3) Obtain the yield and MDF at each point, and plot a step graph.


```python
MDF_result_file = './analysis/MDF_clc_result.csv'
MDF_list_ETM = calculate_MDF_list(norm_model, "ETM", MDF_result_file)
```

    0:06:43.595517



```python
import pandas as pd
import numpy as np

df = pd.read_csv(MDF_result_file, header=None, names=['x', 'y'])
df = df.drop(0)

def find_inflection_points(data):
    inflection_points = []
    for i in range(1, len(data)):
        if data[i-1, 1] > data[i, 1]:
            inflection_points.append((data[i-1, 0], data[i-1, 1]))
    return inflection_points
data_array = df.to_numpy()
inflection_points = find_inflection_points(data_array)
print("Inflection Points:")
for point in inflection_points:
    print(point)
```

    Inflection Points:
    (6.407, '9.764')
    (8.693, '6.387')
    (8.77, '4.661')
    (8.976, '3.879')
    (15.629, '3.767')
    (16.797, '3.035')
    (19.553, '1.844')
    (19.688, '0.856')



```python
MDF_dict = {"ETM" : MDF_list_ETM}
ytick = [10,8,6,4,2,0,-2]
png_file = "./analysis/max_MDF_by_four_model.png"
xlab = "Product rate (mmol/gDW/h)"
Draw_MDF_By_Product_rate(MDF_dict,0,xlab,20,1,ytick,png_file)
```


    
![png](6.Bottleneck_reactions_files/6.Bottleneck_reactions_10_0.png)
    



Taking the inflection point corresponding to a growth rate between 16 and 17 as an example for analysis, we search for the corresponding thermodynamic bottleneck reactions.


```python
set_MDF_substrate = 16
B_value = MDF_Calculation(norm_model,set_MDF_substrate,"DM_ser__L_c")
opt_etm = calculate_ETM(norm_model,"maximize",B_value)
obj_enz_constraint = opt_etm.obj()
```

### 6.1.2. Identifying the thermodynamic bottleneck reactions associated with inflection points

Minimum flux sum calculation（pFBA）


```python
[min_V,Concretemodel]=Min_Flux_Sum_Calculation(norm_model,"ETM",opt_etm.obj(),norm_model.obj_name,B_value)
print("Min flux amount : " +str(min_V))
```

    Min flux amount : 673.1761518200325



```python
use_result = Get_Results_Thermodynamics(norm_model,Concretemodel)
use_result = use_result[use_result['flux'] > 1e-10] 
use_result = use_result.sort_values(by = 'flux',axis = 0,ascending = False)
use_result["reaction"] = use_result.apply(lambda row: norm_model.get_reaction_by_id(row.name).build_reaction_string(), axis = 1)
use_result["gpr"] = use_result.apply(lambda row: norm_model.get_reaction_by_id(row.name).gene_reaction_rule, axis = 1)
use_result.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>flux</th>
      <th>z</th>
      <th>f</th>
      <th>enz</th>
      <th>met_concentration</th>
      <th>reaction</th>
      <th>gpr</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>CYTBO3_4pp</th>
      <td>36.080377</td>
      <td>1.0</td>
      <td>102.619807</td>
      <td>2.684710e-05</td>
      <td>;h2o_c : 1.0;h_c : 1.0;h_p : 1.0;o2_c : 4.9999...</td>
      <td>4.0 h_c + 0.5 o2_c + q8h2_c &lt;=&gt; h2o_c + 4.0 h_...</td>
      <td>b0429 and b0432 and b0431 and b0430</td>
    </tr>
    <tr>
      <th>ATPS4rpp_num2</th>
      <td>29.519810</td>
      <td>1.0</td>
      <td>-9999.000000</td>
      <td>1.035685e-01</td>
      <td>;adp_c : 0.0020000000008562916;atp_c : 0.02000...</td>
      <td>adp_c + 4.0 h_p + pi_c &lt;=&gt; atp_c + h2o_c + 3.0...</td>
      <td>b3732 and b3735 and b3731 and b3733 and b3734 ...</td>
    </tr>
    <tr>
      <th>H2Otex_reverse_num5</th>
      <td>26.432151</td>
      <td>1.0</td>
      <td>-9999.000000</td>
      <td>2.424768e-06</td>
      <td>;h2o_e : 1.0;h2o_p : 1.0</td>
      <td>h2o_p &lt;=&gt; h2o_e</td>
      <td>s0001</td>
    </tr>
    <tr>
      <th>EX_h2o_e</th>
      <td>26.432151</td>
      <td>1.0</td>
      <td>-9999.000000</td>
      <td>7.342338e-11</td>
      <td>;h2o_e : 1.0</td>
      <td>h2o_e &lt;=&gt;</td>
      <td></td>
    </tr>
    <tr>
      <th>H2Otpp_reverse_num1</th>
      <td>26.432151</td>
      <td>1.0</td>
      <td>-9999.000000</td>
      <td>7.342338e-11</td>
      <td>;h2o_c : 1.0;h2o_p : 1.0</td>
      <td>h2o_c &lt;=&gt; h2o_p</td>
      <td>s0001</td>
    </tr>
  </tbody>
</table>
</div>



 List extraction of candidate bottleneck reactions

(1) Bottleneck reactions represent the thermodynamic limitations within a pathway, determining the optimal level achievable for the thermodynamic driving force of the pathway. Identifying and optimizing the thermodynamic bottleneck reactions within a pathway can effectively enhance the thermodynamic feasibility of the pathway.  

(2) Specifically, set the lower flux limit of biomass synthesis reactions, then determine the MDF level achievable for the pathway. Subsequently, use this MDF value as a constraint to solve for the maximum growth rate level. Then, with both the maximum biomass synthesis rate and the MDF value as constraints, calculate the upper and lower bounds of the thermodynamic driving force for each reaction in the model. Finally, identify reactions as bottleneck reactions if both the upper and lower bounds are equal to the MDF value.   

(3) Based on the fact that bottleneck reactions are certain to participate in the pathway and their driving force (Dfi) is not variable, a variability analysis can be performed only on the reactions output by the pFBA method within the pathway, where the thermodynamic driving force values equal the MDF. This significantly reduces the size of the reaction list for variability analysis, thereby shortening the time required to identify bottleneck reactions.


```python
use_result_tmp=use_result[use_result['f']>-norm_model.K_value]
use_result_select=use_result_tmp[abs(use_result_tmp['f']-B_value)<=1e-05]
use_result_select.head(10)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>flux</th>
      <th>z</th>
      <th>f</th>
      <th>enz</th>
      <th>met_concentration</th>
      <th>reaction</th>
      <th>gpr</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>PGCD</th>
      <td>13.868163</td>
      <td>1.0</td>
      <td>3.034508</td>
      <td>2.422459e-02</td>
      <td>;3pg_c : 0.015054490135661808;3php_c : 4.99999...</td>
      <td>3pg_c + nad_c &lt;=&gt; 3php_c + h_c + nadh_c</td>
      <td>b2913</td>
    </tr>
    <tr>
      <th>PGK_reverse</th>
      <td>10.000000</td>
      <td>1.0</td>
      <td>3.034508</td>
      <td>7.717342e-05</td>
      <td>;13dpg_c : 0.0002540594417664234;3pg_c : 0.015...</td>
      <td>13dpg_c + adp_c &lt;=&gt; 3pg_c + atp_c</td>
      <td>b2926</td>
    </tr>
    <tr>
      <th>GAPD</th>
      <td>10.000000</td>
      <td>1.0</td>
      <td>3.034513</td>
      <td>9.879880e-05</td>
      <td>;13dpg_c : 0.0002540594417664234;g3p_c : 0.020...</td>
      <td>g3p_c + nad_c + pi_c &lt;=&gt; 13dpg_c + h_c + nadh_c</td>
      <td>b1779</td>
    </tr>
    <tr>
      <th>ACONTb_num1</th>
      <td>7.736326</td>
      <td>1.0</td>
      <td>3.034508</td>
      <td>9.908535e-04</td>
      <td>;acon_C_c : 0.0001740933479233313;h2o_c : 1.0;...</td>
      <td>acon_C_c + h2o_c &lt;=&gt; icit_c</td>
      <td>b0118</td>
    </tr>
    <tr>
      <th>ICL</th>
      <td>7.736326</td>
      <td>1.0</td>
      <td>3.034508</td>
      <td>2.314813e-05</td>
      <td>;glx_c : 4.516250333787466e-05;icit_c : 9.2371...</td>
      <td>icit_c &lt;=&gt; glx_c + succ_c</td>
      <td>b4015</td>
    </tr>
    <tr>
      <th>ACONTa_num1</th>
      <td>7.736326</td>
      <td>1.0</td>
      <td>3.034508</td>
      <td>5.586844e-08</td>
      <td>;acon_C_c : 0.0001740933479233313;cit_c : 0.02...</td>
      <td>cit_c &lt;=&gt; acon_C_c + h2o_c</td>
      <td>b0118</td>
    </tr>
    <tr>
      <th>GND</th>
      <td>7.439561</td>
      <td>1.0</td>
      <td>3.034508</td>
      <td>3.874772e-06</td>
      <td>;6pgc_c : 0.0036589255431281505;co2_c : 9.9999...</td>
      <td>6pgc_c + nadp_c &lt;=&gt; co2_c + nadph_c + ru5p__D_c</td>
      <td>b2029</td>
    </tr>
    <tr>
      <th>EDA</th>
      <td>7.087788</td>
      <td>1.0</td>
      <td>3.034508</td>
      <td>6.080807e-04</td>
      <td>;2ddg6p_c : 0.02000000000856292;g3p_c : 0.0200...</td>
      <td>2ddg6p_c &lt;=&gt; g3p_c + pyr_c</td>
      <td>b1850</td>
    </tr>
    <tr>
      <th>PGI_reverse</th>
      <td>4.527349</td>
      <td>1.0</td>
      <td>3.034508</td>
      <td>3.549192e-07</td>
      <td>;f6p_c : 0.02000000000856292;g6p_c : 0.0162562...</td>
      <td>f6p_c &lt;=&gt; g6p_c</td>
      <td>b4025</td>
    </tr>
    <tr>
      <th>GLXCL</th>
      <td>3.868163</td>
      <td>1.0</td>
      <td>3.034508</td>
      <td>3.683965e-06</td>
      <td>;2h3oppan_c : 0.0021110700049181337;co2_c : 9....</td>
      <td>2.0 glx_c + h_c &lt;=&gt; 2h3oppan_c + co2_c</td>
      <td>b0507</td>
    </tr>
  </tbody>
</table>
</div>



Determination of bottleneck reaction


```python
path_reac_list=list(use_result_select.index)
reaction_Df_list_fixed=Get_Max_Min_Df_Complete(norm_model,"maximize",path_reac_list,B_value,opt_etm.obj(),norm_model.obj_name)
reaction_Df_list_fixed
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>maximize</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>PGCD</th>
      <td>3.034513</td>
    </tr>
    <tr>
      <th>PGK_reverse</th>
      <td>3.034513</td>
    </tr>
    <tr>
      <th>GAPD</th>
      <td>3.034513</td>
    </tr>
    <tr>
      <th>ACONTb_num1</th>
      <td>16.494213</td>
    </tr>
    <tr>
      <th>ICL</th>
      <td>41.977534</td>
    </tr>
    <tr>
      <th>ACONTa_num1</th>
      <td>16.494213</td>
    </tr>
    <tr>
      <th>GND</th>
      <td>37.189908</td>
    </tr>
    <tr>
      <th>EDA</th>
      <td>21.717833</td>
    </tr>
    <tr>
      <th>PGI_reverse</th>
      <td>28.455846</td>
    </tr>
    <tr>
      <th>GLXCL</th>
      <td>73.808182</td>
    </tr>
    <tr>
      <th>ACKr_num1</th>
      <td>18.667088</td>
    </tr>
    <tr>
      <th>PTAr_reverse_num2</th>
      <td>29.007038</td>
    </tr>
    <tr>
      <th>GHMT2r_reverse</th>
      <td>29.476390</td>
    </tr>
    <tr>
      <th>ASPTA_reverse_num2</th>
      <td>28.814757</td>
    </tr>
    <tr>
      <th>MTHFC_reverse</th>
      <td>22.328721</td>
    </tr>
    <tr>
      <th>ASAD_reverse</th>
      <td>3.034513</td>
    </tr>
    <tr>
      <th>ASPK_num3</th>
      <td>3.034513</td>
    </tr>
    <tr>
      <th>TALA_num1</th>
      <td>52.415552</td>
    </tr>
    <tr>
      <th>TKT1_num1</th>
      <td>28.828726</td>
    </tr>
    <tr>
      <th>GK1</th>
      <td>33.367088</td>
    </tr>
    <tr>
      <th>NDPK5_reverse_num2</th>
      <td>19.417480</td>
    </tr>
    <tr>
      <th>PPM_reverse</th>
      <td>13.528721</td>
    </tr>
    <tr>
      <th>PUNP3_reverse_num2</th>
      <td>26.268671</td>
    </tr>
    <tr>
      <th>PUNP4_num2</th>
      <td>29.588427</td>
    </tr>
    <tr>
      <th>DGK1_reverse</th>
      <td>19.417480</td>
    </tr>
  </tbody>
</table>
</div>




```python
Bottleneck_reaction=reaction_Df_list_fixed[(reaction_Df_list_fixed['maximize']-B_value)<=0.01]
Bottleneck_reaction_list = list(Bottleneck_reaction.index)
Bottleneck_reaction_list
```




    ['PGCD', 'PGK_reverse', 'GAPD', 'ASAD_reverse', 'ASPK_num3']



## 6.2. Discovery of Coupling Strategies

### 6.2.1 Homoenzymatic Coupled Reactions

Based on the identified bottleneck reactions, search for reactions in ETGEMs that share the same gene-protein-reaction (GPR) association and have a standard Gibbs free energy less than 0. These reactions could be potential examples of homoenzymatic coupled reactions.


```python
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
Fusion_Protein_List = Get_Fusion_Protein(Bottleneck_reaction_list,norm_model)
Fusion_Protein_List
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Bottleneck_reaction</th>
      <th>Bottleneck_reaction_g0</th>
      <th>Bottleneck_reaction_gpr</th>
      <th>Bottleneck_reaction_equation</th>
      <th>Coupling_reaction</th>
      <th>Coupling_reaction_g0</th>
      <th>Coupling_reaction_gpr</th>
      <th>Coupling_reaction_equation</th>
      <th>final_g0</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>PGCD</td>
      <td>29.5</td>
      <td>b2913</td>
      <td>3pg_c + nad_c &lt;=&gt; 3php_c + h_c + nadh_c</td>
      <td>AHGDx_reverse</td>
      <td>-27.1</td>
      <td>b2913</td>
      <td>akg_c + h_c + nadh_c &lt;=&gt; S2hglut_c + nad_c</td>
      <td>2.4</td>
    </tr>
    <tr>
      <th>1</th>
      <td>GAPD</td>
      <td>0.1</td>
      <td>b1779</td>
      <td>g3p_c + nad_c + pi_c &lt;=&gt; 13dpg_c + h_c + nadh_c</td>
      <td>E4PD_num2</td>
      <td>-53.9</td>
      <td>b1779</td>
      <td>e4p_c + h2o_c + nad_c &lt;=&gt; 4per_c + 2.0 h_c + n...</td>
      <td>-53.8</td>
    </tr>
  </tbody>
</table>
</div>



You can merge the two reactions and re-solve for the MDF, then compare whether the MDF has improved.


```python
norm_model.merge_reaction("PGCD","AHGDx_reverse")
MDF_list_ETM_ACONT = calculate_MDF_list(norm_model,"ETM","./MDF_list.csv")
MDF_dict["PGCD+AHGDx_reverse"] = MDF_list_ETM_ACONT
```

    0:07:02.379270



```python
ytick = [ -2,0,2,4,6,8,10,12,14,16,18,20]
png_file = "./analysis/max_MDF_by_four_model.png"
xlab = "Product rate (mmol/gDW/h)"
Draw_MDF_By_Product_rate(MDF_dict,0,xlab, 20,1,ytick,png_file)
```


    
![png](6.Bottleneck_reactions_files/6.Bottleneck_reactions_29_0.png)
    


To split the merged reactions


```python
norm_model.split_reaction("PGCD+AHGDx_reverse")
```

### 6.2.1 Successive Reactions

Successive reactions refer to a scenario where the product of one reaction is consumed by another reaction. Based on the identified bottleneck reactions, search for successive reactions with a standard Gibbs free energy less than 0.


```python
Continuous_Reaction_List = Get_Continuous_Reaction(Bottleneck_reaction_list,norm_model)
Continuous_Reaction_List
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Bottleneck_reaction</th>
      <th>Bottleneck_reaction_g0</th>
      <th>Bottleneck_reaction_gpr</th>
      <th>Bottleneck_reaction_equation</th>
      <th>Coupling_reaction</th>
      <th>Coupling_reaction_g0</th>
      <th>Coupling_reaction_gpr</th>
      <th>Coupling_reaction_equation</th>
      <th>Coupling_metabolite</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>GAPD</td>
      <td>0.1</td>
      <td>b1779</td>
      <td>g3p_c + nad_c + pi_c &lt;=&gt; 13dpg_c + h_c + nadh_c</td>
      <td>TKT1_num1</td>
      <td>-1.50</td>
      <td>b2935</td>
      <td>r5p_c + xu5p__D_c &lt;=&gt; g3p_c + s7p_c</td>
      <td>g3p_c</td>
    </tr>
    <tr>
      <th>1</th>
      <td>ASPK_num3</td>
      <td>21.8</td>
      <td>b4024</td>
      <td>asp__L_c + atp_c &lt;=&gt; 4pasp_c + adp_c</td>
      <td>ASNN_num1</td>
      <td>-17.30</td>
      <td>b1767</td>
      <td>asn__L_c + h2o_c &lt;=&gt; asp__L_c + nh4_c</td>
      <td>asp__L_c</td>
    </tr>
    <tr>
      <th>2</th>
      <td>PGCD</td>
      <td>29.5</td>
      <td>b2913</td>
      <td>3pg_c + nad_c &lt;=&gt; 3php_c + h_c + nadh_c</td>
      <td>PSERT</td>
      <td>-10.90</td>
      <td>b0907</td>
      <td>3php_c + glu__L_c &lt;=&gt; akg_c + pser__L_c</td>
      <td>3php_c</td>
    </tr>
    <tr>
      <th>3</th>
      <td>ASAD_reverse</td>
      <td>25.4</td>
      <td>b3433</td>
      <td>4pasp_c + h_c + nadph_c &lt;=&gt; aspsa_c + nadp_c +...</td>
      <td>DHDPS</td>
      <td>-129.00</td>
      <td>b2478</td>
      <td>aspsa_c + pyr_c &lt;=&gt; 23dhdp_c + 2.0 h2o_c + h_c</td>
      <td>aspsa_c</td>
    </tr>
    <tr>
      <th>4</th>
      <td>GAPD</td>
      <td>0.1</td>
      <td>b1779</td>
      <td>g3p_c + nad_c + pi_c &lt;=&gt; 13dpg_c + h_c + nadh_c</td>
      <td>TRPS1</td>
      <td>-31.40</td>
      <td>b1260 and b1261</td>
      <td>3ig3p_c + ser__L_c &lt;=&gt; g3p_c + h2o_c + trp__L_c</td>
      <td>g3p_c</td>
    </tr>
    <tr>
      <th>5</th>
      <td>GAPD</td>
      <td>0.1</td>
      <td>b1779</td>
      <td>g3p_c + nad_c + pi_c &lt;=&gt; 13dpg_c + h_c + nadh_c</td>
      <td>TKT2_num1</td>
      <td>-10.50</td>
      <td>b2935</td>
      <td>e4p_c + xu5p__D_c &lt;=&gt; f6p_c + g3p_c</td>
      <td>g3p_c</td>
    </tr>
    <tr>
      <th>6</th>
      <td>ASPK_num3</td>
      <td>21.8</td>
      <td>b4024</td>
      <td>asp__L_c + atp_c &lt;=&gt; 4pasp_c + adp_c</td>
      <td>ASPabcpp</td>
      <td>-30.61</td>
      <td>b0655 and b0654 and b0652 and b0653</td>
      <td>asp__L_p + atp_c + h2o_c &lt;=&gt; adp_c + asp__L_c ...</td>
      <td>asp__L_c</td>
    </tr>
    <tr>
      <th>7</th>
      <td>PGCD</td>
      <td>29.5</td>
      <td>b2913</td>
      <td>3pg_c + nad_c &lt;=&gt; 3php_c + h_c + nadh_c</td>
      <td>PGM_num1</td>
      <td>-4.30</td>
      <td>b0755</td>
      <td>2pg_c &lt;=&gt; 3pg_c</td>
      <td>3pg_c</td>
    </tr>
    <tr>
      <th>8</th>
      <td>PGCD</td>
      <td>29.5</td>
      <td>b2913</td>
      <td>3pg_c + nad_c &lt;=&gt; 3php_c + h_c + nadh_c</td>
      <td>GLYCK</td>
      <td>-15.60</td>
      <td>b0514</td>
      <td>atp_c + glyc__R_c &lt;=&gt; 3pg_c + adp_c + h_c</td>
      <td>3pg_c</td>
    </tr>
    <tr>
      <th>9</th>
      <td>PGCD</td>
      <td>29.5</td>
      <td>b2913</td>
      <td>3pg_c + nad_c &lt;=&gt; 3php_c + h_c + nadh_c</td>
      <td>HPYRP</td>
      <td>-17.30</td>
      <td>b1813</td>
      <td>3php_c + h2o_c &lt;=&gt; hpyr_c + pi_c</td>
      <td>3php_c</td>
    </tr>
    <tr>
      <th>10</th>
      <td>ASPK_num3</td>
      <td>21.8</td>
      <td>b4024</td>
      <td>asp__L_c + atp_c &lt;=&gt; 4pasp_c + adp_c</td>
      <td>ASPTA_reverse_num1</td>
      <td>-2.80</td>
      <td>b0928</td>
      <td>glu__L_c + oaa_c &lt;=&gt; akg_c + asp__L_c</td>
      <td>asp__L_c</td>
    </tr>
    <tr>
      <th>11</th>
      <td>PGCD</td>
      <td>29.5</td>
      <td>b2913</td>
      <td>3pg_c + nad_c &lt;=&gt; 3php_c + h_c + nadh_c</td>
      <td>PGK_reverse</td>
      <td>-19.50</td>
      <td>b2926</td>
      <td>13dpg_c + adp_c &lt;=&gt; 3pg_c + atp_c</td>
      <td>3pg_c</td>
    </tr>
    <tr>
      <th>12</th>
      <td>GAPD</td>
      <td>0.1</td>
      <td>b1779</td>
      <td>g3p_c + nad_c + pi_c &lt;=&gt; 13dpg_c + h_c + nadh_c</td>
      <td>PGK_reverse</td>
      <td>-19.50</td>
      <td>b2926</td>
      <td>13dpg_c + adp_c &lt;=&gt; 3pg_c + atp_c</td>
      <td>13dpg_c</td>
    </tr>
    <tr>
      <th>13</th>
      <td>ASAD_reverse</td>
      <td>25.4</td>
      <td>b3433</td>
      <td>4pasp_c + h_c + nadph_c &lt;=&gt; aspsa_c + nadp_c +...</td>
      <td>HSDy_reverse_num1</td>
      <td>-18.60</td>
      <td>b3940</td>
      <td>aspsa_c + h_c + nadph_c &lt;=&gt; hom__L_c + nadp_c</td>
      <td>aspsa_c</td>
    </tr>
    <tr>
      <th>14</th>
      <td>GAPD</td>
      <td>0.1</td>
      <td>b1779</td>
      <td>g3p_c + nad_c + pi_c &lt;=&gt; 13dpg_c + h_c + nadh_c</td>
      <td>TKT1_num2</td>
      <td>-1.50</td>
      <td>b2465</td>
      <td>r5p_c + xu5p__D_c &lt;=&gt; g3p_c + s7p_c</td>
      <td>g3p_c</td>
    </tr>
    <tr>
      <th>15</th>
      <td>ASPK_num3</td>
      <td>21.8</td>
      <td>b4024</td>
      <td>asp__L_c + atp_c &lt;=&gt; 4pasp_c + adp_c</td>
      <td>ASNN_num2</td>
      <td>-17.30</td>
      <td>b0828</td>
      <td>asn__L_c + h2o_c &lt;=&gt; asp__L_c + nh4_c</td>
      <td>asp__L_c</td>
    </tr>
    <tr>
      <th>16</th>
      <td>GAPD</td>
      <td>0.1</td>
      <td>b1779</td>
      <td>g3p_c + nad_c + pi_c &lt;=&gt; 13dpg_c + h_c + nadh_c</td>
      <td>TKT2_num2</td>
      <td>-10.50</td>
      <td>b2465</td>
      <td>e4p_c + xu5p__D_c &lt;=&gt; f6p_c + g3p_c</td>
      <td>g3p_c</td>
    </tr>
    <tr>
      <th>17</th>
      <td>PGCD</td>
      <td>29.5</td>
      <td>b2913</td>
      <td>3pg_c + nad_c &lt;=&gt; 3php_c + h_c + nadh_c</td>
      <td>PGM_num2</td>
      <td>-4.30</td>
      <td>b3612</td>
      <td>2pg_c &lt;=&gt; 3pg_c</td>
      <td>3pg_c</td>
    </tr>
    <tr>
      <th>18</th>
      <td>ASPK_num3</td>
      <td>21.8</td>
      <td>b4024</td>
      <td>asp__L_c + atp_c &lt;=&gt; 4pasp_c + adp_c</td>
      <td>ASPTA_reverse_num2</td>
      <td>-2.80</td>
      <td>b4054</td>
      <td>glu__L_c + oaa_c &lt;=&gt; akg_c + asp__L_c</td>
      <td>asp__L_c</td>
    </tr>
    <tr>
      <th>19</th>
      <td>ASAD_reverse</td>
      <td>25.4</td>
      <td>b3433</td>
      <td>4pasp_c + h_c + nadph_c &lt;=&gt; aspsa_c + nadp_c +...</td>
      <td>HSDy_reverse_num2</td>
      <td>-18.60</td>
      <td>b0002</td>
      <td>aspsa_c + h_c + nadph_c &lt;=&gt; hom__L_c + nadp_c</td>
      <td>aspsa_c</td>
    </tr>
  </tbody>
</table>
</div>


