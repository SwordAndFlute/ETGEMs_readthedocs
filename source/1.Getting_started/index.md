# 1. Getting Started

## 1.1. Loading a model and inspecting it
ETGEMs comes with pre-built Enzyme and Thermodynic constrained Genome-scale mEtabolic network Models (ETGEMs) for *Escherichia coli*, *Bacillus subtilis*, and *Corynebacterium glutamicum*. To load the ETGEMs of *Escherichia coli*, type


```python
import autoETGEMs
from autoETGEMs.io.json import *
model = read_json_model("../autoETGEMs/model/iML1515_etcmodel.json")
```

Reactions, metabolites, and genes are individually established as objects and stored in the form of dictionaries in model.reactions, model.metabolites, and model.genes. The keys of the dictionaries represent the names, and the values correspond to the respective objects.


```python
print(len(model.reactions))
print(len(model.metabolites))
print(len(model.genes))
```

    5918
    1877
    1516
    


```python
model.metabolites["coa_c"]
```




    <autoETGEMs.core.metabolite.Metabolite at 0x22276010fd0>




```python
model.reactions["PDH"]
```




    <autoETGEMs.core.reaction.Reaction at 0x222757c99b0>



In addition, the medium is also separately stored as a dictionary in the model.


```python
model.medium
```




    {'EX_glc__D_e_reverse': 10}



The construction and computation of ETGEMs require some additional parameters, which are stored in model.environment.


```python
model.environment
```




    <autoETGEMs.core.environment.Environment at 0x222760a9208>



## 1.2. Reactions

You can use get_reaction_by_id() to find the corresponding reaction by its name.


```python
reaction = model.get_reaction_by_id("PDH")
print(reaction)
```

    PDH: coa_c + nad_c + pyr_c <=> accoa_c + co2_c + nadh_c
    

You can directly view the information contained in a reaction by accessing its attributes or properties. For example, you might inspect properties such as reaction.id, reaction.name, reaction.metabolites, etc..


```python
print(reaction.id)
print(reaction.metabolites)
```

    PDH
    {'accoa_c': 1.0, 'co2_c': 1.0, 'coa_c': -1.0, 'nad_c': -1.0, 'nadh_c': 1.0, 'pyr_c': -1.0}
    

The reactions include upper and lower bounds for flux. Unlike cobrapy, where all reactions in the model are bidirectional, in this model, all reactions are split into irreversibility reactions, hence the upper and lower bounds for reactions are both greater than 0.


```python
print("lower_bound : ", reaction.lower_bound)
print("upper_bound : ", reaction.upper_bound)
```

    lower_bound :  0.0
    upper_bound :  1000.0
    

If the lower bound of a reaction is set to be greater than the upper bound, an error will be raised.


```python
reaction.set_lower_bound(1100)
```


    ---------------------------------------------------------------------------

    ValueError                                Traceback (most recent call last)

    <ipython-input-34-81bd5304f44f> in <module>
    ----> 1 reaction.set_lower_bound(1100)
    

    d:\autoETGEMs\etgems-master\autoETGEMs_test\autoETGEMs\core\reaction.py in set_lower_bound(self, value)
        246         """
        247         # Validate bounds before setting them.
    --> 248         self.check_bounds(value, self.upper_bound)
        249         self.lower_bound = value
        250 
    

    d:\autoETGEMs\etgems-master\autoETGEMs_test\autoETGEMs\core\reaction.py in check_bounds(self, lb, ub)
        266         if lb > ub:
        267             raise ValueError(
    --> 268                 f"The lower bound must be less than or equal to the upper bound "
        269                 f"({lb} <= {ub})."
        270             )    
    

    ValueError: The lower bound must be less than or equal to the upper bound (1100 <= 1000.0).


For ETGEMs, reactions also store enzyme turnover number (kcat), enzyme molecular weight (MW), and the standard Gibbs free energy of the reaction to meet the computational requirements for enzyme constraints and thermodynamic constraints.


```python
print("reaction kcat : ", reaction.kcat)
print("reaction MW : ", reaction.mw)
print("reaction g0 : ", reaction.g0)
```

    reaction kcat :  37.9
    reaction MW :  216452.0
    reaction g0 :  -34.4
    

## 1.3. Metabolites

Similar to reactions, you can find a specific metabolite by its name. 


```python
atp = model.get_metabolite_by_id("atp_c")
atp
```




    <autoETGEMs.core.metabolite.Metabolite at 0x22275fa2d30>



You can easily view information such as the metabolite's name, charge, chemical formula, and more.


```python
print(atp.name)
print(atp.charge)
print(atp.formula)
print(atp.compartment)
```

    ATP C10H12N5O13P3
    -4
    C10H12N5O13P3
    c
    

To meet the computational requirements of MDF, we have also added upper and lower bounds for metabolite concentrations, similar to setting upper and lower bounds for reaction fluxes. If the lower bound of a metabolite concentration is set higher than the upper bound, the code will automatically raise an error.


```python
print("concentration upperbound : ",atp.concentration_ub)
print("concentration lowerbound : ",atp.concentration_lb)
atp.set_concentration_lb(-3)
```

    concentration upperbound :  -3.912023005
    concentration lowerbound :  -14.50865774
    


    ---------------------------------------------------------------------------

    ValueError                                Traceback (most recent call last)

    <ipython-input-38-9efbbdfb0d1a> in <module>
          1 print("concentration upperbound : ",atp.concentration_ub)
          2 print("concentration lowerbound : ",atp.concentration_lb)
    ----> 3 atp.set_concentration_lb(-3)
    

    d:\autoETGEMs\etgems-master\autoETGEMs_test\autoETGEMs\core\metabolite.py in set_concentration_lb(self, concentration_lb)
         59             The lower bound of the concentration.
         60         """
    ---> 61         self.check_bounds(concentration_lb, self.concentration_ub)
         62         self.concentration_lb = concentration_lb
         63 
    

    d:\autoETGEMs\etgems-master\autoETGEMs_test\autoETGEMs\core\metabolite.py in check_bounds(self, lb, ub)
         79         if lb > ub:
         80             raise ValueError(
    ---> 81                 f"The lower bound must be less than or equal to the upper bound "
         82                 f"({lb} <= {ub})."
         83             )
    

    ValueError: The lower bound must be less than or equal to the upper bound (-3 <= -3.912023005).


## 1.4. Genes

Similarly to reactions and metabolites, genes can also be queried by name. 
Additionally, genes encoding enzymes catalyzing a reaction are saved in the GPR (Gene-Protein-Reaction) relationship of the reaction.


```python
print(reaction.gene_reaction_rule)
b0115 = model.get_gene_by_id("b0115")
```

    b0116 and b0115 and b0114
    

## 1.5. Medium

In the medium, the substrates used for model growth are saved in dictionary form. The keys of the dictionary are the names of the reactions for substrate uptake, and the values are the upper limits for those reactions. You can view information about the model's substrates through model.medium.


```python
model.medium
```




    {'EX_glc__D_e_reverse': 10}



You can add and remove components from the model's medium using add_medium() and remove_medium() functions, respectively.


```python
model.remove_medium("EX_glc__D_e_reverse")
model.medium
```




    {}




```python
model.add_medium("EX_glc__D_e_reverse",10)
model.medium
```




    {'EX_glc__D_e_reverse': 10}



## 1.6. Environment

Parameters used for obtaining standard Gibbs free energy, as well as parameters for enzyme constraints and thermodynamic constraints calculations, are stored in the environment of the model. By changing the parameter values, you can simulate bacterial growth under different conditions.


```python
print("ionic strength : ",model.environment.get_ionic_strength())
model.environment.set_ionic_strength("0.2M")
print("new ionic strength : ",model.environment.get_ionic_strength())
```

    ionic strength :  0.1M
    new ionic strength :  0.2M
    

Among them, K_value and E_total can be quickly queried.


```python
print("K_value : ", model.K_value)
print("E_total : ", model.E_total)
```

    K_value :  10353.64521565178
    E_total :  0.22680000000000003
    
