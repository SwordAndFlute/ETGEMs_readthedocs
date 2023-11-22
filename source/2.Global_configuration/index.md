# 2. Global Configuration

When converting the GEM structure to the ETGEMs format, it is necessary to assign initial values to various parameters in the model. This chapter primarily provides the initial values of parameters in ETGEMs before the automated acquisition of parameters.


The function trans_model2standard_json_etgem() can be used to convert a standard GEM into the standard ETGEMs format and save it at the specified path.


```python
import autoETGEMs
from autoETGEMs.io.convert import *
from autoETGEMs.io.json import *
trans_model2standard_json_etgem("./autoETGEMs/model/iML1515_github.json","./autoETGEMs/model/iML1515_standard.json")
```


The read_json_model() function is used to read a JSON-formatted ETGEMs model.


```python
model = read_json_model("./autoETGEMs/model/iML1515_standard.json")
```

## 2.1. Reaction bounds

The default value for the reaction flux upper limit in GEM is set to 1000, consistent with Cobra settings. However, in contrast, in ETGEMs where reversible reactions are not present, the default value for the reaction flux lower limit is set to 0.


```python
reaction = model.get_reaction_by_id("PDH")
print("default upper bound : ",reaction.upper_bound)
print("default lower bound : ",reaction.lower_bound)
```

    default upper bound :  1000.0
    default lower bound :  0.0
    

## 2.2. Metabolite bounds

In ETGEMs, the concentration range for most metabolites is set to 0.5 μM to 20 mM. The concentration limits are stored in the corresponding metabolites in ln form. 

Additionally, constraints are applied to the concentration ratios of five pairs of metabolites: ATP:ADP, ADP:AMP, NAD+:NADH, NADPH:NADP+, and HCO₃⁻: CO₂. These concentration ratios are respectively set to 10:1, 1:1, 10:1, 10:1, and 2:1.


```python
metabolite = model.get_metabolite_by_id("atp_c")
```


```python
print("default concentration upper bound : ",metabolite.concentration_ub)
print("default concentration lower bound : ",metabolite.concentration_lb)
```

    default concentration upper bound :  -3.912023005
    default concentration lower bound :  -14.50865774
    

## 2.3. Enzyme-related parameters

Due to the calculation of enzyme usage in terms of enzyme molecular weight (MW) divided by turnover number (kcat), at the initial state, to approach an enzyme usage close to 0, we set the default value of MW to 1 and the default value of kcat to 99999.


```python
print("default kcat : ", reaction.kcat)
print("default MW : ", reaction.mw)
```

    default kcat :  99999.0
    default MW :  1.0
    

In addition to the enzyme-related parameters for reactions, the calculation of enzyme constraints also requires defining the total enzyme weight. The total enzyme weight parameter can be queried using model.E_total with a default value of 1.


```python
print("default total enzyme weight : ", model.E_total)
```

    default total enzyme weight :  1.0
    

## 2.4. Thermodynamic-related parameters


At the initial state, by default, all reactions are thermodynamically feasible, and the default standard Gibbs free energy of reactions is set to -99999.


```python
print("default standard Gibbs free energy: ", reaction.g0)
```

    default standard Gibbs free energy:  -99999.0
    

In ETGEMs, when computing thermodynamic constraints, the calculation of the K_value involves the standard Gibbs free energy of reactions. The default value for K_value is set to 0.


```python
print("default K_value: ", model.K_value)
```

    default K_value:  0.0
    


In addition to the parameters required for the aforementioned calculations, in the process of automating the retrieval of standard Gibbs free energy of reactions using the equilibrator-api, specific conditions need to be provided. These parameters are stored in model.environment, and their default values are as follows.


```python
print("default ph: ", model.environment.p_h)
print("default pmg: ", model.environment.p_mg)
print("default temperature: ", model.environment.temperature)
print("default ionic strength: ", model.environment.ionic_strength)
```

    default ph:  7.5
    default pmg:  3.0
    default temperature:  310.15K
    default ionic strength:  0.1M
    

## 2.5. Medium

By default, in ETGEMs, the medium is empty, and attempting to directly compute fluxes will result in an error. Therefore, it is necessary to use model.add_medium() to add the components of the medium beforehand.


```python
model.medium
```




    {}


