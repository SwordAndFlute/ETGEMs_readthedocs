# 7. Limitation Metabolites

This chapter explains how to use ETGEMs for the identification of limiting metabolites.


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

Before conducting the analysis of limiting metabolites, it is necessary to first identify bottleneck reactions. The detailed analysis of thermodynamic bottleneck reactions can be found in Chapter 6.


```python
norm_model.add_demand_reaction("ser__L_c")
norm_model.obj_name = "DM_ser__L_c"
Bottleneck_reaction_list = ['PGCD', 'PGK_reverse', 'GAPD', 'ASAD_reverse', 'ASPK_num3']
B_value = 3.034507565524991
obj_enz_constraint = 16.78392460294908
```

(1) Limiting metabolites refer to those whose concentrations reach the set boundaries or play dual roles as substrates and products in different bottleneck reactions. This dual role leads to their concentrations being constrained, eliminating variability. Therefore, there exists an optimal concentration level.

(2) Since limiting metabolites are certain to participate in bottleneck reactions, a list of metabolites involved in bottleneck reactions can be compiled. Subsequently, the concentration variability range of these metabolites is calculated, and metabolites with non-variable concentrations are identified as limiting metabolites.


```python
Bottleneck_reaction_met=[]
for rea in norm_model.reactions.values():
    if rea.id in Bottleneck_reaction_list:
        #print(rea)
        for met in rea.metabolites.keys():
            if met !='h_c' and met !='h2o_c':
                Bottleneck_reaction_met.append(met)
                

Bottleneck_reaction_met=list(set(Bottleneck_reaction_met))
Bottleneck_reaction_met
```




    ['aspsa_c',
     'nadph_c',
     'asp__L_c',
     'adp_c',
     '13dpg_c',
     '3pg_c',
     'nad_c',
     'nadp_c',
     '4pasp_c',
     'nadh_c',
     'g3p_c',
     'pi_c',
     '3php_c',
     'atp_c']




```python
max_min_concentration_list_fixed = Get_Max_Min_Met_Concentration_for_list(norm_model,B_value,obj_enz_constraint,norm_model.obj_name,Bottleneck_reaction_list,Bottleneck_reaction_met,4)
max_min_concentration_list_fixed
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
      <th>max_concentration</th>
      <th>min_concentration</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>nadph_c</th>
      <td>-3.912023</td>
      <td>-12.206073</td>
    </tr>
    <tr>
      <th>adp_c</th>
      <td>-6.214608</td>
      <td>-14.508658</td>
    </tr>
    <tr>
      <th>asp__L_c</th>
      <td>-3.912023</td>
      <td>-3.912025</td>
    </tr>
    <tr>
      <th>aspsa_c</th>
      <td>-14.508656</td>
      <td>-14.508658</td>
    </tr>
    <tr>
      <th>13dpg_c</th>
      <td>-8.277940</td>
      <td>-8.277942</td>
    </tr>
    <tr>
      <th>nad_c</th>
      <td>-3.912023</td>
      <td>-12.206073</td>
    </tr>
    <tr>
      <th>3pg_c</th>
      <td>-4.196077</td>
      <td>-4.196079</td>
    </tr>
    <tr>
      <th>nadp_c</th>
      <td>-6.214608</td>
      <td>-14.508658</td>
    </tr>
    <tr>
      <th>4pasp_c</th>
      <td>-11.238948</td>
      <td>-11.238950</td>
    </tr>
    <tr>
      <th>nadh_c</th>
      <td>-6.214608</td>
      <td>-14.508658</td>
    </tr>
    <tr>
      <th>g3p_c</th>
      <td>-3.912023</td>
      <td>-3.912025</td>
    </tr>
    <tr>
      <th>pi_c</th>
      <td>-5.453106</td>
      <td>-5.453108</td>
    </tr>
    <tr>
      <th>atp_c</th>
      <td>-3.912023</td>
      <td>-12.206073</td>
    </tr>
    <tr>
      <th>3php_c</th>
      <td>-14.508656</td>
      <td>-14.508658</td>
    </tr>
  </tbody>
</table>
</div>




```python
Limiting_metabolite = max_min_concentration_list_fixed[(max_min_concentration_list_fixed['max_concentration'] - max_min_concentration_list_fixed['min_concentration']) <= 0.001]
Limiting_metabolite
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
      <th>max_concentration</th>
      <th>min_concentration</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>asp__L_c</th>
      <td>-3.912023</td>
      <td>-3.912025</td>
    </tr>
    <tr>
      <th>aspsa_c</th>
      <td>-14.508656</td>
      <td>-14.508658</td>
    </tr>
    <tr>
      <th>13dpg_c</th>
      <td>-8.277940</td>
      <td>-8.277942</td>
    </tr>
    <tr>
      <th>3pg_c</th>
      <td>-4.196077</td>
      <td>-4.196079</td>
    </tr>
    <tr>
      <th>4pasp_c</th>
      <td>-11.238948</td>
      <td>-11.238950</td>
    </tr>
    <tr>
      <th>g3p_c</th>
      <td>-3.912023</td>
      <td>-3.912025</td>
    </tr>
    <tr>
      <th>pi_c</th>
      <td>-5.453106</td>
      <td>-5.453108</td>
    </tr>
    <tr>
      <th>3php_c</th>
      <td>-14.508656</td>
      <td>-14.508658</td>
    </tr>
  </tbody>
</table>
</div>


