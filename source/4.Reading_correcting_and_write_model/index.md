# 4. Reading, Correcting, and Writing Models.

## 4.1. Reading models 

You can use the read_json_model() function to read ETGEMs in JSON format.


```python
from autoETGEMs.io.json import *
model = read_json_model("./autoETGEMs/model/iML1515_etcmodel_test.json")
```

## 4.2. Correcting models

When incorporating the thermodynamic parameters obtained in Chapter 3 into the model, thermodynamic over-constraints may sometimes occur. Specifically, thermodynamics-constraint metabolic network model simulates a growth solution of zero. Therefore, it is necessary to investigate and rectify thermodynamic over-constraints in the model.


```python
from autoETGEMs.analysis.flux_balance_analysis import *
opt_tcm = calculate_TCM(model,"maximize",0.0)
print("TCM optimize : ", opt_tcm.obj())
```

    TCM optimize :  0.0
    

The find_bottleneck_reaction() function can be used to check for thermodynamic over-constrained reactions in the model and adjust their standard Gibbs free energy.


```python
from autoETGEMs.tcmodel.construct import *
fbr_result = find_bottleneck_reaction(model = model, obj_target = "maximize", B_value = 0.0, max_workers = 4)
fbr_result
```




    {'S7PI': 1.8138139876344725e-11,
     'F6PA_num1': 1.6156061328368017e-11,
     'PMPK': 1.923621594122025e-11,
     'TALA_reverse_num1': 6.120614163115733e-12,
     'AIRC3_reverse': 0.8769934969068253,
     'DTMPK_reverse': 4.590460622336699e-12,
     'PGMT_reverse_num3': 5.100511802596291e-12}



In the given key-value pairs, the keys are the names of reactions with unreasonable standard Gibbs free energy, and the values are the objective reaction flux values of the model after removing the reaction from thermodynamic constraints. Since most values are close to zero, except for the value of AIRC3_reverse, which is more reasonable, we modified the standard Gibbs free energy of this reaction.


```python
model.get_reaction_by_id('AIRC3_reverse').set_g0(0)
opt_tcm = calculate_TCM(model,"maximize",0.0)
print("TCM optimize : ", opt_tcm.obj())
```

    TCM optimize :  0.8769934969068253
    

After correcting the standard Gibbs free energy for over-constrained reactions, the thermodynamic constrained model can be solved normally.

## Writing models

Save the ETGEMs object in JSON format at the specified path.


```python
write_json_model(model,"./autoETGEMs/model/iML1515_etcmodel_test.json")
```
