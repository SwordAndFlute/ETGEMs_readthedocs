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
    [0;38;2;66;227;35m100.00%[0;38;2;250;205;229m|█████████████████████████████|[0;38;2;85;85;85m 0:00:00|1:26:17 [0;38;2;146;52;247m ETC: 10-09 11:26:53[0m[K
    Success downloading! :-)
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


## 3.3. Machine learning for obtaining enzyme kinetic parameters.

### 3.3.1 DLKcat

ETGEMs supports one-step acquisition of reaction's kcat in the model using the DLKcat process in ECMpy. Input a model in cobrapy format and generate a CSV file containing kcat and molecular weight information for all reactions.


```python
DLKcat_folder = "./analysis/get_kcat_mw_by_DLKcat/"
get_reaction_kcatmw_onestop_by_DLKcat(DLKcat_folder,sbml_path)
```

The obtained results for kcat and mw are as follows. You can directly modify the kcat for the corresponding reactions in the model by assigning specific values, such as reaction.kcat = specific_value.


```python
pd.read_csv("./analysis/get_kcat_mw_by_DLKcat/reaction_kcat_MW.csv",index_col=0)
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
      <th>data_type</th>
      <th>kcat</th>
      <th>MW</th>
      <th>kcat_MW</th>
    </tr>
    <tr>
      <th>reactions</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>MDH_reverse</th>
      <td>DLkcat</td>
      <td>2343.8073</td>
      <td>34907.6530</td>
      <td>60428.771020</td>
    </tr>
    <tr>
      <th>FUM</th>
      <td>DLkcat</td>
      <td>1854.7952</td>
      <td>50243.3157</td>
      <td>33224.632108</td>
    </tr>
    <tr>
      <th>UREAAH</th>
      <td>DLkcat</td>
      <td>1843.6494</td>
      <td>61450.8998</td>
      <td>73530.485598</td>
    </tr>
    <tr>
      <th>PRFGS</th>
      <td>DLkcat</td>
      <td>972.9622</td>
      <td>8743.9328</td>
      <td>30962.982836</td>
    </tr>
    <tr>
      <th>AKGDH</th>
      <td>DLkcat</td>
      <td>765.9681</td>
      <td>50651.0874</td>
      <td>1210.495279</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>LEUt2rpp_reverse_num2</th>
      <td>DLkcat</td>
      <td>0.0090</td>
      <td>11479.8881</td>
      <td>2.822327</td>
    </tr>
    <tr>
      <th>LEUt2rpp_num2</th>
      <td>DLkcat</td>
      <td>0.0090</td>
      <td>11479.8881</td>
      <td>2.822327</td>
    </tr>
    <tr>
      <th>UPPRT_num2</th>
      <td>DLkcat</td>
      <td>0.0053</td>
      <td>20951.4081</td>
      <td>0.910679</td>
    </tr>
    <tr>
      <th>ALCD19_reverse_num3</th>
      <td>DLkcat</td>
      <td>0.0023</td>
      <td>12076.6632</td>
      <td>0.685620</td>
    </tr>
    <tr>
      <th>PPNDH</th>
      <td>DLkcat</td>
      <td>0.0017</td>
      <td>33739.4344</td>
      <td>0.181390</td>
    </tr>
  </tbody>
</table>
<p>738 rows × 4 columns</p>
</div>



### 3.3.2 TurNuP

The relevant files are saved in zhaojx/autoETGEMs/TurNup. Some portions of the code are from https://github.com/AlexanderKroll/kcat_prediction


```python
input_file = pd.read_csv("./analysis/get_kcat_mw_by_DLKcat/comdf.csv",index_col=0)
input_file = input_file.groupby(input_file.index).first()
unique_reaction_list = input_file.index.unique().tolist()
result_df = pd.DataFrame(index=unique_reaction_list)
```

Retrieve the metabolite InChI through the CKB for a reaction.


```python
import sqlite3

conn = sqlite3.connect('./../../etgems/autoETGEMs_test/data/compounds.sqlite') 

cursor = conn.cursor()
for index,row in result_df.iterrows():
    reaction = model.get_reaction_by_id(index)
    reactants = reaction.reactants
    reactants_inchi = ""
    for metabolite in reactants:
        accession_to_query = metabolite[:-2]
        cursor.execute(f"SELECT compound_id FROM compound_identifiers WHERE accession = '{accession_to_query}'")
        result = cursor.fetchone()
        if result:
            compound_id = result[0]

            cursor.execute(f"SELECT inchi FROM compounds WHERE id = '{compound_id}'")
            inchi_result = cursor.fetchone()

            if inchi_result is not None:
                inchi_result = inchi_result[0]
                if inchi_result:
                    reactants_inchi = reactants_inchi + inchi_result + ";"
    if reactants_inchi:
        reactants_inchi = reactants_inchi[:-1]
    products = reaction.products
    products_inchi = ""
    for metabolite in products:
        accession_to_query = metabolite[:-2]
        cursor.execute(f"SELECT compound_id FROM compound_identifiers WHERE accession = '{accession_to_query}'")
        result = cursor.fetchone()
        if result:
            compound_id = result[0]

            cursor.execute(f"SELECT inchi FROM compounds WHERE id = '{compound_id}'")
            inchi_result = cursor.fetchone()

            if inchi_result:
                inchi_result = inchi_result[0]
                if inchi_result:
                    products_inchi = products_inchi + inchi_result + ";"
    if products_inchi:
        products_inchi = products_inchi[:-1]
    result_df.loc[index,"substrates"] = reactants_inchi
    result_df.loc[index,"products"] = products_inchi
    result_df.loc[index,"enzyme"] = input_file.loc[index,"prosequence"]
conn.close()
```


```python
result_df
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
      <th>substrates</th>
      <th>products</th>
      <th>enzyme</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>2MAHMP</th>
      <td>InChI=1S/C6H11N3O7P2/c1-4-8-2-5(6(7)9-4)3-15-1...</td>
      <td>InChI=1S/C6H10N3O4P/c1-4-8-2-5(6(7)9-4)3-13-14...</td>
      <td>MKTWKTWGVVGASGLLIILSWLSSSSPMLADAFMIAAAIVAGWPIA...</td>
    </tr>
    <tr>
      <th>2METS</th>
      <td>InChI=1S/H2O/h1H2;InChI=1S/C4H4O5/c5-2(4(8)9)1...</td>
      <td>InChI=1S/C21H36N7O16P3S/c1-21(2,16(31)19(32)24...</td>
      <td>MSSATTTDVRKGLYGVIADYTAVSKVMPETNSLTYRGYAVEDLVEN...</td>
    </tr>
    <tr>
      <th>2METS_reverse</th>
      <td>InChI=1S/C21H36N7O16P3S/c1-21(2,16(31)19(32)24...</td>
      <td>InChI=1S/H2O/h1H2;InChI=1S/C4H4O5/c5-2(4(8)9)1...</td>
      <td>MSSATTTDVRKGLYGVIADYTAVSKVMPETNSLTYRGYAVEDLVEN...</td>
    </tr>
    <tr>
      <th>3NTD7pp</th>
      <td>InChI=1S/C10H14N5O7P/c11-8-5-9(13-2-12-8)15(3-...</td>
      <td>InChI=1S/C10H13N5O4/c11-8-5-9(13-2-12-8)15(3-1...</td>
      <td>MKRLSRAALAVVATTAVSFSALAVPAFADEASNVELNILGVTDFHG...</td>
    </tr>
    <tr>
      <th>3NTD9pp</th>
      <td>InChI=1S/C10H14N5O8P/c11-10-13-7-4(8(18)14-10)...</td>
      <td>InChI=1S/C10H13N5O5/c11-10-13-7-4(8(19)14-10)1...</td>
      <td>MKRLSRAALAVVATTAVSFSALAVPAFADEASNVELNILGVTDFHG...</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>XTSNH_num1</th>
      <td>InChI=1S/H2O/h1H2;InChI=1S/C10H12N4O6/c15-1-3-...</td>
      <td>InChI=1S/C5H10O5/c6-1-2-3(7)4(8)5(9)10-2/h2-9H...</td>
      <td>MAQRKLASVIGAALAASAVLVGLMTPATAQSSGSSSTDITRALTSS...</td>
    </tr>
    <tr>
      <th>XTSNH_num2</th>
      <td>InChI=1S/H2O/h1H2;InChI=1S/C10H12N4O6/c15-1-3-...</td>
      <td>InChI=1S/C5H10O5/c6-1-2-3(7)4(8)5(9)10-2/h2-9H...</td>
      <td>MPRSTIKRVVAVLAASTALSPFLVSMPTAAAQENIRWEECPPQVDI...</td>
    </tr>
    <tr>
      <th>XTSNH_num3</th>
      <td>InChI=1S/H2O/h1H2;InChI=1S/C10H12N4O6/c15-1-3-...</td>
      <td>InChI=1S/C5H10O5/c6-1-2-3(7)4(8)5(9)10-2/h2-9H...</td>
      <td>MRIALLQISTNSDKMDNFALLRDAAEKAAEQGARVLVFPEATSQSF...</td>
    </tr>
    <tr>
      <th>XYLK</th>
      <td>InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(...</td>
      <td>InChI=1S/C10H15N5O10P2/c11-8-5-9(13-2-12-8)15(...</td>
      <td>MKFWIPQRKRIFKMALVLGIDSSTQSCKALLVDAATGQVIDEGRAS...</td>
    </tr>
    <tr>
      <th>XYLUt2pp</th>
      <td>InChI=1S/p+1;InChI=1S/C5H10O5/c6-1-3(8)5(10)4(...</td>
      <td>InChI=1S/p+1;InChI=1S/C5H10O5/c6-1-3(8)5(10)4(...</td>
      <td>MKFWIPQRKRIFKMALVLGIDSSTQSCKALLVDAATGQVIDEGRAS...</td>
    </tr>
  </tbody>
</table>
<p>986 rows × 3 columns</p>
</div>




```python
from kcat_prediction import *
for index,row in result_df.iloc.iterrows():
    print(index)
    print(substrates)
    substrates = [result_df.loc[index,"substrates"]]
    products = [result_df.loc[index,"products"]]
    enzymes = [result_df.loc[index,"enzyme"]]
    print(len(substrates))
    if substrates != [np.nan] and products != [np.nan] and enzymes != [np.nan]:
        result = kcat_predicton(substrates = substrates,
                products = products,
                enzymes = enzymes)
        if "kcat [s^(-1)]" in result.columns:
            result_df.loc[index,"kcat"] = result["kcat [s^(-1)]"].tolist()[0]
result_df.to_csv("turnup_result.csv")
```


```python
pd.read_csv("turnup_result.csv",index_col=0)
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
      <th>substrates</th>
      <th>products</th>
      <th>enzyme</th>
      <th>kcat</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>2MAHMP</th>
      <td>InChI=1S/C6H11N3O7P2/c1-4-8-2-5(6(7)9-4)3-15-1...</td>
      <td>InChI=1S/C6H10N3O4P/c1-4-8-2-5(6(7)9-4)3-13-14...</td>
      <td>MKTWKTWGVVGASGLLIILSWLSSSSPMLADAFMIAAAIVAGWPIA...</td>
      <td>11.735885</td>
    </tr>
    <tr>
      <th>2METS</th>
      <td>InChI=1S/H2O/h1H2;InChI=1S/C4H4O5/c5-2(4(8)9)1...</td>
      <td>InChI=1S/C21H36N7O16P3S/c1-21(2,16(31)19(32)24...</td>
      <td>MSSATTTDVRKGLYGVIADYTAVSKVMPETNSLTYRGYAVEDLVEN...</td>
      <td>21.775869</td>
    </tr>
    <tr>
      <th>2METS_reverse</th>
      <td>InChI=1S/C21H36N7O16P3S/c1-21(2,16(31)19(32)24...</td>
      <td>InChI=1S/H2O/h1H2;InChI=1S/C4H4O5/c5-2(4(8)9)1...</td>
      <td>MSSATTTDVRKGLYGVIADYTAVSKVMPETNSLTYRGYAVEDLVEN...</td>
      <td>24.529543</td>
    </tr>
    <tr>
      <th>3NTD7pp</th>
      <td>InChI=1S/C10H14N5O7P/c11-8-5-9(13-2-12-8)15(3-...</td>
      <td>InChI=1S/C10H13N5O4/c11-8-5-9(13-2-12-8)15(3-1...</td>
      <td>MKRLSRAALAVVATTAVSFSALAVPAFADEASNVELNILGVTDFHG...</td>
      <td>47.222221</td>
    </tr>
    <tr>
      <th>3NTD9pp</th>
      <td>InChI=1S/C10H14N5O8P/c11-10-13-7-4(8(18)14-10)...</td>
      <td>InChI=1S/C10H13N5O5/c11-10-13-7-4(8(19)14-10)1...</td>
      <td>MKRLSRAALAVVATTAVSFSALAVPAFADEASNVELNILGVTDFHG...</td>
      <td>36.123520</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>XTSNH_num1</th>
      <td>InChI=1S/H2O/h1H2;InChI=1S/C10H12N4O6/c15-1-3-...</td>
      <td>InChI=1S/C5H10O5/c6-1-2-3(7)4(8)5(9)10-2/h2-9H...</td>
      <td>MAQRKLASVIGAALAASAVLVGLMTPATAQSSGSSSTDITRALTSS...</td>
      <td>176.085449</td>
    </tr>
    <tr>
      <th>XTSNH_num2</th>
      <td>InChI=1S/H2O/h1H2;InChI=1S/C10H12N4O6/c15-1-3-...</td>
      <td>InChI=1S/C5H10O5/c6-1-2-3(7)4(8)5(9)10-2/h2-9H...</td>
      <td>MPRSTIKRVVAVLAASTALSPFLVSMPTAAAQENIRWEECPPQVDI...</td>
      <td>15.612193</td>
    </tr>
    <tr>
      <th>XTSNH_num3</th>
      <td>InChI=1S/H2O/h1H2;InChI=1S/C10H12N4O6/c15-1-3-...</td>
      <td>InChI=1S/C5H10O5/c6-1-2-3(7)4(8)5(9)10-2/h2-9H...</td>
      <td>MRIALLQISTNSDKMDNFALLRDAAEKAAEQGARVLVFPEATSQSF...</td>
      <td>39.495121</td>
    </tr>
    <tr>
      <th>XYLK</th>
      <td>InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(...</td>
      <td>InChI=1S/C10H15N5O10P2/c11-8-5-9(13-2-12-8)15(...</td>
      <td>MKFWIPQRKRIFKMALVLGIDSSTQSCKALLVDAATGQVIDEGRAS...</td>
      <td>103.108665</td>
    </tr>
    <tr>
      <th>XYLUt2pp</th>
      <td>InChI=1S/p+1;InChI=1S/C5H10O5/c6-1-3(8)5(10)4(...</td>
      <td>InChI=1S/p+1;InChI=1S/C5H10O5/c6-1-3(8)5(10)4(...</td>
      <td>MKFWIPQRKRIFKMALVLGIDSSTQSCKALLVDAATGQVIDEGRAS...</td>
      <td>116.249046</td>
    </tr>
  </tbody>
</table>
<p>986 rows × 4 columns</p>
</div>



### 3.3.3 UniKP

The relevant files are saved in zhaojx/autoETGEMs/UniKP. Some portions of the code are from https://github.com/HanselYu/UniKP.


```python
def smiles_to_vec(Smiles):
    pad_index = 0
    unk_index = 1
    eos_index = 2
    sos_index = 3
    mask_index = 4
    vocab = WordVocab.load_vocab('./vocab.pkl')
    def get_inputs(sm):
        seq_len = 220
        sm = sm.split()
        if len(sm)>218:
            print('SMILES is too long ({:d})'.format(len(sm)))
            sm = sm[:109]+sm[-109:]
        ids = [vocab.stoi.get(token, unk_index) for token in sm]
        ids = [sos_index] + ids + [eos_index]
        seg = [1]*len(ids)
        padding = [pad_index]*(seq_len - len(ids))
        ids.extend(padding), seg.extend(padding)
        return ids, seg
    def get_array(smiles):
        x_id, x_seg = [], []
        for sm in smiles:
            a,b = get_inputs(sm)
            x_id.append(a)
            x_seg.append(b)
        return torch.tensor(x_id), torch.tensor(x_seg)
    trfm = TrfmSeq2seq(len(vocab), 256, len(vocab), 4)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    trfm.load_state_dict(torch.load('./trfm_12_23000.pkl', map_location=device))
    # trfm.load_state_dict(torch.load('./trfm_12_23000.pkl'))
    trfm.eval()
    x_split = [split(sm) for sm in Smiles]
    xid, xseg = get_array(x_split)
    X = trfm.encode(torch.t(xid))
    return X


def Seq_to_vec(Sequence):
    for i in range(len(Sequence)):
        if len(Sequence[i]) > 1000:
            Sequence[i] = Sequence[i][:500] + Sequence[i][-500:]
    sequences_Example = []
    for i in range(len(Sequence)):
        zj = ''
        for j in range(len(Sequence[i]) - 1):
            zj += Sequence[i][j] + ' '
        zj += Sequence[i][-1]
        sequences_Example.append(zj)
    tokenizer = T5Tokenizer.from_pretrained("./Prot-T5-XL-UniRef50/prot_t5_xl_uniref50/", do_lower_case=False)
    model = T5EncoderModel.from_pretrained("./Prot-T5-XL-UniRef50/prot_t5_xl_uniref50/")
    gc.collect()
    print(torch.cuda.is_available())
    # 'cuda:0' if torch.cuda.is_available() else
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    model = model.to(device)
    model = model.eval()
    features = []
    for i in range(len(sequences_Example)):
        print('For sequence ', str(i+1))
        sequences_Example_i = sequences_Example[i]
        sequences_Example_i = [re.sub(r"[UZOB]", "X", sequences_Example_i)]
        ids = tokenizer.batch_encode_plus(sequences_Example_i, add_special_tokens=True, padding=True)
        input_ids = torch.tensor(ids['input_ids']).to(device)
        attention_mask = torch.tensor(ids['attention_mask']).to(device)
        with torch.no_grad():
            embedding = model(input_ids=input_ids, attention_mask=attention_mask)
        embedding = embedding.last_hidden_state.cpu().numpy()
        for seq_num in range(len(embedding)):
            seq_len = (attention_mask[seq_num] == 1).sum()
            seq_emd = embedding[seq_num][:seq_len - 1]
            features.append(seq_emd)
    features_normalize = np.zeros([len(features), len(features[0][0])], dtype=float)
    for i in range(len(features)):
        for k in range(len(features[0][0])):
            for j in range(len(features[i])):
                features_normalize[i][k] += features[i][j][k]
            features_normalize[i][k] /= len(features[i])
    return features_normalize
```


```python
input_file = pd.read_csv("./analysis/get_kcat_mw_by_DLKcat/DLinput.tsv",sep='\t',index_col=0)
```


```python
sequences = input_file['prosequence']
Smiles = input_file['similes']
seq_vec = Seq_to_vec(sequences)
smiles_vec = smiles_to_vec(Smiles)
fused_vector = np.concatenate((smiles_vec, seq_vec), axis=1)

# For kcat
with open('./UniKP for kcat.pkl', "rb") as f:
    model = pickle.load(f)
# For Km
# with open('UniKP/UniKP for Km.pkl', "rb") as f:
#     model = pickle.load(f)
# For kcat/Km
# with open('UniKP/UniKP for kcat_Km.pkl', "rb") as f:
#     model = pickle.load(f)

Pre_label = model.predict(fused_vector)
Pre_label_pow = [math.pow(10, Pre_label[i]) for i in range(len(Pre_label))]
print(len(Pre_label))
res = pd.DataFrame({'sequences': sequences, 'Smiles': Smiles, 'Pre_label': Pre_label})
res.to_csv('Kinetic_parameters_predicted_label.csv')
```


```python
for index,row in res.iterrows():
    res.loc[index,"Kcat value (1/s)"] = 10**row["Pre_label"]
res.to_csv("./Unikp_result.csv")
```


```python
from ECMpy_function import *
comdf = pd.read_csv("./analysis/get_kcat_mw_by_DLKcat/comdf.csv")
DL_reaction_kcat_mw = DL_kcat_mw_calculation(res, comdf)
DL_reaction_kcat_mw.to_csv("unikp_kcat_mw.csv")
```


```python
pd.read_csv("unikp_kcat_mw.csv",index_col=0)
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
      <th>reactions</th>
      <th>data_type</th>
      <th>kcat</th>
      <th>MW</th>
      <th>kcat_MW</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>3OAR60_reverse</td>
      <td>DLkcat</td>
      <td>21.104117</td>
      <td>32579.5582</td>
      <td>1527.343066</td>
    </tr>
    <tr>
      <th>1</th>
      <td>GLYCK2</td>
      <td>DLkcat</td>
      <td>21.104117</td>
      <td>48197.8658</td>
      <td>1576.310903</td>
    </tr>
    <tr>
      <th>2</th>
      <td>ME1_num1</td>
      <td>DLkcat</td>
      <td>20.866252</td>
      <td>37729.4565</td>
      <td>1990.977701</td>
    </tr>
    <tr>
      <th>3</th>
      <td>ME2</td>
      <td>DLkcat</td>
      <td>20.866252</td>
      <td>40932.3755</td>
      <td>1835.185611</td>
    </tr>
    <tr>
      <th>4</th>
      <td>MDH</td>
      <td>DLkcat</td>
      <td>20.866252</td>
      <td>34907.6530</td>
      <td>537.980214</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>981</th>
      <td>LEUt2rpp_reverse_num3</td>
      <td>DLkcat</td>
      <td>NaN</td>
      <td>44831.8483</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>982</th>
      <td>URAt2rpp_reverse_num2</td>
      <td>DLkcat</td>
      <td>NaN</td>
      <td>66934.4117</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>983</th>
      <td>PHETA1_reverse_num2</td>
      <td>DLkcat</td>
      <td>NaN</td>
      <td>39917.7769</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>984</th>
      <td>SHK3Dr_reverse_num2</td>
      <td>DLkcat</td>
      <td>NaN</td>
      <td>29304.7230</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>985</th>
      <td>TYRTA_reverse_num2</td>
      <td>DLkcat</td>
      <td>NaN</td>
      <td>39917.7769</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
<p>986 rows × 5 columns</p>
</div>




In addition to this, the total enzyme amount needs to be consulted from the literature and saved in model.E_total.

## 3.4. Automated Retrieval of Thermodynamic Data

By using the equilibrator-api, you can quickly retrieve the standard Gibbs free energy of reactions in the model and save it in ETGEMs.


```python
from autoETGEMs.tcmodel.construct import *
# 自动化获取热力学数据
tcmodel = tcmodel_construction(model)
```

After the retrieval of the standard Gibbs free energy of reactions is completed, ETGEMs will automatically calculate the K value and save it in model.K_value.
