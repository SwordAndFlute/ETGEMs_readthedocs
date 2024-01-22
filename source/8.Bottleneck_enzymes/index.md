# 8. Bottleneck Enzymes


This chapter explains how to use ETGEMs to identify bottleneck enzymes.



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

In the integrated enzyme-constrained model, for reactions where the enzyme cost is higher, it can be considered that they exert stronger control over the pathway, making them critical enzymes in the pathway. Therefore, with the determination of growth and MDF levels, the minimum enzyme cost required for each reaction can be calculated. A higher numerical value indicates a more critical role. Generally, in enzyme-constrained models, when enzyme quantity constraints come into play, the variability of enzyme levels in pathway reactions becomes small or even non-existent. In contrast, when enzyme quantity constraints are not met, there is more variability in the enzyme costs of reactions in the pathway. Thus, minimizing the enzyme cost is informative regardless of whether enzyme quantity constraints are met, whereas maximizing the enzyme cost may not be as meaningful.

Taking the production of arginine in Escherichia coli as an example, specific data sources can be found in Chapter 6.


```python
use_result = pd.read_csv("./use_result.csv",index_col=0)
norm_model.add_demand_reaction("ser__L_c")
norm_model.obj_name = "DM_ser__L_c"
B_value = 3.034507565524991
obj_enz_constraint = 16.78392460294908
use_result
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
      <td>NaN</td>
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
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>RNTR2c2_num3</th>
      <td>0.648538</td>
      <td>1.0</td>
      <td>-9999.000000</td>
      <td>3.343511e-07</td>
      <td>;dgtp_c : 0.02000000000856292;flxr_c : 4.99999...</td>
      <td>2.0 flxr_c + gtp_c + 2.0 h_c &lt;=&gt; dgtp_c + 2.0 ...</td>
      <td>b2895 and b4238</td>
    </tr>
    <tr>
      <th>NDPK1_num2</th>
      <td>0.648538</td>
      <td>1.0</td>
      <td>8.138367</td>
      <td>3.501028e-06</td>
      <td>;adp_c : 0.0020000000008562916;atp_c : 0.02000...</td>
      <td>atp_c + gdp_c &lt;=&gt; adp_c + gtp_c</td>
      <td>b2518</td>
    </tr>
    <tr>
      <th>FORtppi_num2</th>
      <td>0.541656</td>
      <td>1.0</td>
      <td>-9999.000000</td>
      <td>8.558936e-08</td>
      <td>;for_c : 0.02000000000856292;for_p : 1.0</td>
      <td>for_c &lt;=&gt; for_p</td>
      <td>b2492</td>
    </tr>
    <tr>
      <th>FDH4pp_num1</th>
      <td>0.541656</td>
      <td>1.0</td>
      <td>117.671279</td>
      <td>4.767550e-07</td>
      <td>;co2_p : 1.0;for_p : 1.0;h_c : 1.0;h_p : 1.0;q...</td>
      <td>for_p + 2.0 h_c + q8_c &lt;=&gt; co2_p + h_p + q8h2_c</td>
      <td>b3894 and b3893 and b3892</td>
    </tr>
    <tr>
      <th>PDH</th>
      <td>0.066070</td>
      <td>1.0</td>
      <td>46.713522</td>
      <td>1.048153e-04</td>
      <td>;accoa_c : 0.011820926251385226;co2_c : 9.9999...</td>
      <td>coa_c + nad_c + pyr_c &lt;=&gt; accoa_c + co2_c + na...</td>
      <td>b0116 and b0115 and b0114</td>
    </tr>
  </tbody>
</table>
<p>84 rows × 7 columns</p>
</div>




```python
use_result_sort = use_result.sort_values(by='enz',ascending = False)
use_result_sort.head()
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
      <th>ATPS4rpp_num2</th>
      <td>29.519810</td>
      <td>1.0</td>
      <td>-9999.000000</td>
      <td>0.103569</td>
      <td>;adp_c : 0.0020000000008562916;atp_c : 0.02000...</td>
      <td>adp_c + 4.0 h_p + pi_c &lt;=&gt; atp_c + h2o_c + 3.0...</td>
      <td>b3732 and b3735 and b3731 and b3733 and b3734 ...</td>
    </tr>
    <tr>
      <th>PGCD</th>
      <td>13.868163</td>
      <td>1.0</td>
      <td>3.034508</td>
      <td>0.024225</td>
      <td>;3pg_c : 0.015054490135661808;3php_c : 4.99999...</td>
      <td>3pg_c + nad_c &lt;=&gt; 3php_c + h_c + nadh_c</td>
      <td>b2913</td>
    </tr>
    <tr>
      <th>THRA_num2</th>
      <td>2.915762</td>
      <td>1.0</td>
      <td>22.629070</td>
      <td>0.019905</td>
      <td>;acald_c : 3.2016953742420415e-06;gly_c : 0.02...</td>
      <td>thr__L_c &lt;=&gt; acald_c + gly_c</td>
      <td>b0870</td>
    </tr>
    <tr>
      <th>RPE</th>
      <td>4.527349</td>
      <td>1.0</td>
      <td>3.400000</td>
      <td>0.015440</td>
      <td>;ru5p__D_c : 0.02000000000856292;xu5p__D_c : 0...</td>
      <td>ru5p__D_c &lt;=&gt; xu5p__D_c</td>
      <td>b3386</td>
    </tr>
    <tr>
      <th>G6PDH2r</th>
      <td>14.527349</td>
      <td>1.0</td>
      <td>28.455846</td>
      <td>0.013447</td>
      <td>;6pgl_c : 4.999999992621078e-07;g6p_c : 0.0162...</td>
      <td>g6p_c + nadp_c &lt;=&gt; 6pgl_c + h_c + nadph_c</td>
      <td>b1852</td>
    </tr>
  </tbody>
</table>
</div>




```python
e_threshold = norm_model.E_total*0.01
enz_use_reaction_list = list(use_result_sort[use_result_sort['enz'] > e_threshold].index)
enz_use_reaction_list
```




    ['ATPS4rpp_num2',
     'PGCD',
     'THRA_num2',
     'RPE',
     'G6PDH2r',
     'HSDy_reverse_num1',
     'PGL',
     'PFL_num2',
     'THRS']



 Determination of key enzymes.


```python
max_min_E_list_fixed = max_Max_Min_E_for_list(norm_model,B_value,norm_model.obj_name,obj_enz_constraint,enz_use_reaction_list,4)
max_min_E_list_fixed = max_min_E_list_fixed.sort_values(by=['min_E'],ascending=False)
max_min_E_list_fixed
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
      <th>max_E</th>
      <th>min_E</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>ATPS4rpp_num2</th>
      <td>0.103756</td>
      <td>0.102866</td>
    </tr>
    <tr>
      <th>PGCD</th>
      <td>0.024252</td>
      <td>0.024204</td>
    </tr>
    <tr>
      <th>THRA_num2</th>
      <td>0.019993</td>
      <td>0.019705</td>
    </tr>
    <tr>
      <th>RPE</th>
      <td>0.015616</td>
      <td>0.015196</td>
    </tr>
    <tr>
      <th>G6PDH2r</th>
      <td>0.013488</td>
      <td>0.013381</td>
    </tr>
    <tr>
      <th>PGL</th>
      <td>0.010283</td>
      <td>0.010197</td>
    </tr>
    <tr>
      <th>THRS</th>
      <td>0.006764</td>
      <td>0.006697</td>
    </tr>
    <tr>
      <th>PFL_num2</th>
      <td>0.009365</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>HSDy_reverse_num1</th>
      <td>0.010312</td>
      <td>0.000000</td>
    </tr>
  </tbody>
</table>
</div>




```python
Limiting_enzyme = max_min_E_list_fixed[(max_min_E_list_fixed['max_E'] - max_min_E_list_fixed['min_E']) <= 0.001]
for index,row in Limiting_enzyme.iterrows():
    Limiting_enzyme.loc[index,"gpr"] = norm_model.get_reaction_by_id(index).gene_reaction_rule
Limiting_enzyme
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
      <th>max_E</th>
      <th>min_E</th>
      <th>gpr</th>
    </tr>
    <tr>
      <th>Unnamed: 0</th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>ATPS4rpp_num2</th>
      <td>0.103756</td>
      <td>0.102866</td>
      <td>b3732 and b3735 and b3731 and b3733 and b3734 ...</td>
    </tr>
    <tr>
      <th>PGCD</th>
      <td>0.024252</td>
      <td>0.024204</td>
      <td>b2913</td>
    </tr>
    <tr>
      <th>THRA_num2</th>
      <td>0.019993</td>
      <td>0.019705</td>
      <td>b0870</td>
    </tr>
    <tr>
      <th>RPE</th>
      <td>0.015616</td>
      <td>0.015196</td>
      <td>b3386</td>
    </tr>
    <tr>
      <th>G6PDH2r</th>
      <td>0.013488</td>
      <td>0.013381</td>
      <td>b1852</td>
    </tr>
    <tr>
      <th>PGL</th>
      <td>0.010283</td>
      <td>0.010197</td>
      <td>b0767</td>
    </tr>
    <tr>
      <th>THRS</th>
      <td>0.006764</td>
      <td>0.006697</td>
      <td>b0004</td>
    </tr>
  </tbody>
</table>
</div>


