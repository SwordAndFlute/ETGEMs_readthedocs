
# 3. Build a Model


This chapter primarily encompasses the entire process of generating ETGEMs from a metabolic network model, including the conversion of model structure, automated retrieval of enzyme data, and automated retrieval of thermodynamic data.

## 3.1. Model structure conversion


```python
import datetime 
import autoETGEMs
from autoETGEMs.io.convert import *
from autoETGEMs.io.json import *
from autoETGEMs.ecmpy.AutoPACMEN_function import *
from autoETGEMs.ecmpy.ECMpy_function import *
from autoETGEMs.tcmodel.construct import *
```

First, it is necessary to convert the metabolic network model into the standard structure of ETGEMs, with default parameter values as described in Chapter 2.


```python
trans_model2standard_json_etgem("./autoETGEMs/model/iML1515_github.json","./autoETGEMs/model/iML1515_standard.json")
model = read_json_model("./autoETGEMs/model/iML1515_standard.json")
```

## 3.2. Automated Retrieval of Enzyme Data


Enzyme data retrieval primarily relies on the ECMpy toolkit, and here, only a brief overview of the process is provided. For the complete workflow, please refer to the ECMpy documentation.

Define the file paths for input and output; here, we use the Escherichia coli metabolic network model iML1515 as an example.


```python
# input files
autopacmen_folder = "./analysis/get_kcat_mw_by_AutoPACMEN/"
kcat_gap_fill= 'mean'
reaction_gap_fill='mean'
sbml_path = "./data/iML1515_new.xml"
organism = "Escherichia coli"
project_name = "iML1515_%s"%kcat_gap_fill
create_file(autopacmen_folder)
protein_kcat_database_path = "none"
bigg_metabolites_file = "./data/bigg_models_metabolites.txt"
brenda_textfile_path = "./data/brenda_2023_1.txt"
uniprot_data_file='./data/uniprot_data_accession_key.json'

#output files
brenda_json_path = "%skcat_database_brenda.json"%autopacmen_folder
brenda_json_path2 = "%ssa_database_brenda.json"%autopacmen_folder
sabio_rk_json_path = "%skcat_database_sabio_rk.json"%autopacmen_folder
bigg_id_name_mapping_path = "%sbigg_id_name_mapping.json"%autopacmen_folder
brenda_output_json_path = "%skcat_database_brenda_for_model.json"%autopacmen_folder
combined_output_path = "%skcat_database_combined.json"%autopacmen_folder
sub_description_path = '%sget_gene_subunitDescription.csv'%autopacmen_folder
gene_subnum_path = "%sgene_subnum.csv"%autopacmen_folder
reaction_mw_path = "%sreaction_mw.json"%autopacmen_folder
reaction_kcat_mw_path = '%sreaction_kcat_MW.csv'%autopacmen_folder
```

    Path exists
    

Based on the ECMpy toolkit, use AutoPacmen for the automated retrieval of enzyme data. The first step involves parsing metabolite BIGG-formatted IDs.


```python
print("Starting to deal BIGG metabolites text file...")
parse_bigg_metabolites_file(bigg_metabolites_file, autopacmen_folder)
print("BIGG metabolites text file done!")
print()
```

    Starting to deal BIGG metabolites text file...
    BIGG metabolites text file done!
    
    

### 3.2.1 Enzyme turnover number

ECMpy supports retrieving the required kcat values for the model from BRENDA.


```python
parse_brenda_textfile(brenda_textfile_path, autopacmen_folder, brenda_json_path, brenda_json_path2) 
parse_brenda_json_for_model(sbml_path, brenda_json_path, brenda_output_json_path)
```

ECMpy also supports retrieving the required kcat values for the model from SABIO-RK.


```python
print("Starting EC numbers kcat search in SABIO-RK...")
parse_sabio_rk_for_model_with_sbml(sbml_path, sabio_rk_json_path, bigg_id_name_mapping_path)
print("SABIO-RK done!")
```

    Starting EC numbers kcat search in SABIO-RK...
    SABIO-RK done!
    
    0:39:02.525710
    

### 3.2.2 Enzyme molecule weight

As the model's gene-protein-reaction (gpr) relationships do not include the subunit number of enzymes, it is necessary to obtain the subunit number of enzymes first when calculating the enzyme molecular weight.


```python
starttime=datetime.datetime.now()
print("Starting to fetch subunit number of each enzyme")
model=cobra.io.read_sbml_model(sbml_path)
get_gene_subunitDescription(sub_description_path,model)
subbnumdf = get_subunit_number(sub_description_path,gene_subnum_path)
print("Calculation done!")

endtime=datetime.datetime.now()
print(endtime-starttime)
```

    Starting to fetch subunit number of each enzyme
    Start downloading from UniProt...
    Calculation done!
    
    1:26:25.348115
    

Calculate the enzyme's molecular weight based on the obtained subunit number and the molecular weight of each subunit.


```python
starttime=datetime.datetime.now()

print("Starting UniProt ID<->Protein mass search using UniProt...")
get_protein_mass_mapping_from_local(sbml_path, autopacmen_folder, project_name, uniprot_data_file)
get_reaction_mw(sbml_path,autopacmen_folder, project_name, reaction_mw_path, gene_subnum_path)
print("Protein ID<->Mass mapping done!")

endtime=datetime.datetime.now()
print(endtime-starttime)
```

    Starting UniProt ID<->Protein mass search using UniProt...
    Protein ID<->Mass mapping done!
    
    0:00:20.203922
    


In addition to this, the total enzyme amount needs to be consulted from the literature and saved in model.E_total.

## 3.3. Automated Retrieval of Thermodynamic Data

By using the equilibrator-api, you can quickly retrieve the standard Gibbs free energy of reactions in the model and save it in ETGEMs.


```python
from autoETGEMs.tcmodel.construct import *
tcmodel = tcmodel_construction(model)
```

After the retrieval of the standard Gibbs free energy of reactions is completed, ETGEMs will automatically calculate the K value and save it in model.K_value.
