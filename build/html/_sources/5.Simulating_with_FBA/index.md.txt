# 5. Simulating with FBA

This chapter primarily introduces how several of the most common flux analysis algorithms for metabolic network models operate on ETGEMs.

## 5.1. Flux balance analysis (FBA)

First, read the model in JSON format.


```python
from autoETGEMs.io.json import *
model = read_json_model("./autoETGEMs/model/iML1515_etcmodel_20231108.json")
```

Check if the model's growth medium includes the substrate.


```python
model.medium
```




    {'EX_glc__D_e_reverse': 10}



Solutions were obtained for four models separately (GEM: basic metabolic network model, ECM: with enzyme constraints, TCM: with thermodynamic constraints, ETM: with both enzyme and thermodynamic constraints). The default objective function for the models is biomass synthesis reaction.

When calculating models with added thermodynamic constraints, it is necessary to provide a B_value to constrain the minimum thermodynamic driving force of the model. In the example, we set the B_value default to 0.


```python
from autoETGEMs.analysis.flux_balance_analysis import *
opt_gem = calculate_GEM(model,"maximize")
opt_ecm = calculate_ECM(model,"maximize")
opt_tcm = calculate_TCM(model,"maximize",0.0)
opt_etm = calculate_ETM(model,"maximize",0.0)
```

To view the flux value of the objective function, use opt.obj().


```python
print("gem_obj_flux: " + str(opt_gem.obj()))
print("ecm_obj_flux: " + str(opt_ecm.obj()))
print("tcm_obj_flux: " + str(opt_tcm.obj()))
print("etm_obj_flux: " + str(opt_etm.obj()))
```

    gem_obj_flux: 0.876996163704355
    ecm_obj_flux: 0.6772977358819029
    tcm_obj_flux: 0.872879761713119
    etm_obj_flux: 0.6608464145195769
    

If you want to view the flux value of a specific reaction, the example code would be as follows:


```python
opt_gem.reaction["PDH"].value
```




    9.689659737887645



If you want to check the concentration value (in ln form) of a specific metabolite in a model with added thermodynamic constraints, the code is as follows:


```python
opt_tcm.metabolite["atp_c"].value
```




    -12.206072647005957



If you want to optimize the maximum production rate of a simulated product, you need to first use the add_demand_reaction() function to add a demand reaction for the target product. 


```python
model.add_demand_reaction("arg__L_c")
```

Replace model.obj_name with the added demand reaction, and the reaction name is 'DM_' + the ID of the target product.


```python
model.obj_name = "DM_arg__L_c"
```

Solve using the four models.


```python
opt_gem = calculate_GEM(model,"maximize")
opt_ecm = calculate_ECM(model,"maximize")
opt_tcm = calculate_TCM(model,"maximize",0.0)
opt_etm = calculate_ETM(model,"maximize",0.0)
```


```python
print("gem_arg_yield: " + str(opt_gem.obj()))
print("ecm_arg_yield: " + str(opt_ecm.obj()))
print("tcm_arg_yield: " + str(opt_tcm.obj()))
print("etm_arg_yield: " + str(opt_etm.obj()))
```

    gem_arg_yield: 8.340480000007807
    ecm_arg_yield: 7.607841176597471
    tcm_arg_yield: 7.822809917363128
    etm_arg_yield: 5.438335308720897
    

## 5.2. parsimonious enzyme usage FBA (pFBA)

ETGEMs also support pFBA. Here, we take the calculation of the basic metabolic network model as an example.

Perform pFBA calculation using Min_Flux_Sum_Calculation(). In the parameters, besides specifying the model file, you also need to provide the model type, the target flux value for FBA, and the name of the target reaction.


```python
[min_v, ConcreteModel] = Min_Flux_Sum_Calculation(model,"GEM",8.340480000007807,"DM_arg__L_c",0.0)
opt_gem_pfba = Model_Solve(ConcreteModel,model.solver)
print("minimum flux sum : ", min_v)
print("arg field by pFBA : ", opt_gem_pfba.reaction["DM_arg__L_c"].value)
```

    minimum flux sum :  867.8473600009307
    arg field by pFBA :  8.340480000007807
    

## 5.3. Flux variability analysis

ETGEMs supports Flux Variability Analysis (FVA), using the basic metabolic network model to simulate growth as an example.

Change model.obj_name to the biomass synthesis reaction.


```python
model.obj_name = "BIOMASS_Ec_iML1515_core_75p37M"
```


```python
opt_gem = calculate_GEM(model,"maximize")
```


```python
opt_gem.obj()
```




    0.8769961637043371



To perform analysis using the flux_variability_analysis() function, parameters such as the model object, model type, target reaction flux value, target reaction name, etc., need to be provided. Additionally, fraction_of_optimum is a value less than 1. In the example, setting it to 0.95 signifies that the flux variability analysis is conducted for all reactions when the target reaction flux value is above 95% of the optimum


```python
from autoETGEMs.analysis.flux_variability_analysis import *
# 通量可变性分析
fva_gem = flux_variability_analysis(model, model_type="GEM", obj_value=opt_gem.obj(), 
                                    fraction_of_optimum=0.95, B_value=0.0, biomass_id="BIOMASS_Ec_iML1515_core_75p37M",max_workers=4)
```


```python
fva_gem
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
      <th>minimize</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>DHORD5</th>
      <td>1.065904</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>GLYCTO2</th>
      <td>15.542667</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>GTHOr</th>
      <td>1000.000000</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>G3PD5_num1</th>
      <td>11.657000</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>TRPS2</th>
      <td>7.982256</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>4ABZGLUtex_reverse_num3</th>
      <td>1000.000000</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>AI2tex_reverse_num4</th>
      <td>1000.000000</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>METGLCURtex_reverse_num2</th>
      <td>1000.000000</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>AI2tex_reverse_num2</th>
      <td>1000.000000</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>AI2tex_reverse_num3</th>
      <td>1000.000000</td>
      <td>0.0</td>
    </tr>
  </tbody>
</table>
<p>5918 rows × 2 columns</p>
</div>


