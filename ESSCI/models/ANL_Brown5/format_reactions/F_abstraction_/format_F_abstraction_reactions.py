import cantera as ct
import os 
import re
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.data.kinetics.family import KineticsFamily
from rmgpy.data.rmg import RMGDatabase
from rmgpy import settings
from rmgpy.data.base import Database
from rmgpy.chemkin import load_species_dictionary

#make the species dictionary

unneeded_species = ['CHF(T)',
'CF2(T)',
'CH3CF(T)',
'CH2FCH(T)',
'CH2FCF(T)',
'CHF2CH(S)',
'CHF2CF(T)',
'CF3CH(S)',
'CF3CF(T)',
'CF3CO']


species_dictionary = {}


file = './hcof_anl0_with_c2_20231126.yaml'

with open(file, 'r') as f: 
    lines = f.readlines()
#species are in lines 12-103
for line in lines[11:103]:
    split_species = line.split(', #')
    [label, smiles] = split_species
    
    #handle the label first
    label = label.replace(' ', '')
        
    #now the smiles
    smiles = smiles.replace(' ', '').replace('\n', '').replace(',','')
    match=re.search('([a-z]+)', smiles)
    if match:
        also_remove = match.group(1)
        smiles = smiles.replace(also_remove, '')
        
    #now save to species dictionary
    if label in unneeded_species: 
        continue
    else:
        species_dictionary[label] = smiles

#add the last couple of stragglers

species_dictionary['H']= "[H]"
species_dictionary['H2']= "[H][H]"
species_dictionary['CH']= "[CH]"
species_dictionary['CH2']= "[CH2]"
species_dictionary['CH3']= "[CH3]"
species_dictionary['CH4']= "C"
species_dictionary['O']= "[O]"
species_dictionary['OH']= "[OH]"
species_dictionary['H2O']= "O"
species_dictionary['O2']= "[O][O]"
species_dictionary['HO2']= "[O]O"
species_dictionary['H2O2']= "OO"
species_dictionary['CO']= "[C-]#[O+]"
species_dictionary['CO2']= "O=C=O"
species_dictionary['HCO']= "[CH]=O"
species_dictionary['CH2O'] = "C=O"


########### now let's make the reaction.py reactions 

gas = ct.Solution('/work/westgroup/nora/Code/projects/PFAS/ESSCI/models/ANL_Brown5/format_reactions/F_abstraction/hcof_anl0_with_c2_20231126.yaml')

#load in the database
database = RMGDatabase()
database.load(
            path = settings['database.directory'],
            thermo_libraries = ['Klippenstein_Glarborg2016', 'BurkeH2O2', 'thermo_DFT_CCSDTF12_BAC', 
                               'DFT_QCI_thermo',
                           'primaryThermoLibrary', 'primaryNS', 'NitrogenCurran', 'NOx2018', 'FFCM1(-)',
'SulfurLibrary', 'SulfurGlarborgH2S','SABIC_aromatics'],
            transport_libraries = [],
            reaction_libraries = [],
            seed_mechanisms = [],#['BurkeH2O2inN2','ERC-FoundationFuelv0.9'],
            kinetics_families = 'all',
            kinetics_depositories = ['training'],
            #frequenciesLibraries = self.statmechLibraries,
            depository = False, # Don't bother loading the depository information, as we don't use it
        )

#existing F abstraction training dictionary 
e_path = '/home/khalil.nor/Code/RMG-database/input/kinetics/families/F_Abstraction/training/dictionary.txt'
existing_dict = load_species_dictionary(e_path)


#function definitions

def get_labeled_adj_lists(rxn): 
    """
    given a reaction, use RMG's .get_labeled_reactants_and_products to get the actual adjacency lists with all the *#s 

    """
    
    labeled_adjacencies = {}

    reactants = list(rxn.reactants.keys())
    products = list(rxn.products.keys())

    #make all Molecule objects
    reactant_Molecule = [Molecule(smiles=species_dictionary[x]) for x in reactants]
    product_Molecule = [Molecule(smiles=species_dictionary[x]) for x in products]

    reactant_Species = [Species(molecule=[Molecule(smiles=species_dictionary[x])], label=x) for x in reactants]
    product_Species = [Species(molecule=[Molecule(smiles=species_dictionary[x])], label=x) for x in products]


    generated_labels = database.kinetics.families[rmg_family].get_labeled_reactants_and_products(reactant_Molecule,product_Molecule, relabel_atoms=True)

    #reactants
    for reactant_s in reactant_Species: 
        for reactant in generated_labels[0]: #generated_labels[0] is list of reactants from .get_labeled_reactants_and_products
            if reactant.is_isomorphic(reactant_s.molecule[0]):  #let's match which generated label goes with which reactant
                reactant_label = reactant_s.label
                if reactant_label not in labeled_adjacencies.keys():
                    labeled_adjacencies[reactant_label]=[reactant.to_adjacency_list(), reactant_s]

    #products 
    for product_s in product_Species: 
        for product in generated_labels[1]: 
            if product.is_isomorphic(product_s.molecule[0]):
                product_label = product_s.label
                if product_label not in labeled_adjacencies.keys():
                    labeled_adjacencies[product_label]=[product.to_adjacency_list(), product_s]
    
    return labeled_adjacencies #labeled_adjacencies is a dictionary with key = species label, value = [generated adj list, Species object]

def get_matching_species_label(labeled_adjacencies):   

    switch_labels = {}
    for key, value in labeled_adjacencies.items():
        [adj_list, spec] = value
        starred = re.findall('(\*[0-9] [A-Z])', adj_list)
        for label, species in existing_dict.items(): #let's see if theres an existing species name that precisely matches this adj list
            if species.is_isomorphic(spec):
                flag = 0 
                for starred_str in starred:
                    test_str = species.to_adjacency_list()
                    if starred_str in test_str: #that there is an existing species in dictionary.txt that has a "*# Atom_symbol"
                        pass
                    else: 
                        flag+=1
                        
                    #make sure theres no extra stars in the existing species adj list that we're looking at 
                    already_present = re.findall('(\*[0-9] [A-Z])', test_str)
                    if len(already_present)>len(starred):
                        #this means theres an extra
                        flag+=1
                if flag == 0:
                    switch_labels[key]=label #if everything has worked out...
                    # print(species.to_adjacency_list())
                    # print('matched to')
                    # print(adj_list)

    return switch_labels




rmg_family = 'F_Abstraction'
H_abstraction_lines = []
continued_count = 90 #90 was last index in H_abstraction reactions.py
for index, rxn in enumerate(gas.reactions()[41:81]):
    continued_count+=1
 
    #get rate parameters
    A = "{:e}".format(rxn.rate.pre_exponential_factor*(100*100*100)/1000) #convert from m3 / kmol s to cm3 / mol s
    n = rxn.rate.temperature_exponent
    Ea = rxn.rate.activation_energy/(1000*1000) #convert from J/kmol to kJ/mol
    
    #get the degeneracy using RMG
    reactants = list(rxn.reactants.keys())
    products = list(rxn.products.keys())

    reactant_Molecule = [Molecule(smiles=species_dictionary[x]) for x in reactants]
    product_Molecule = [Molecule(smiles=species_dictionary[x]) for x in products]

    try: 
        template_reaction = database.kinetics.families[rmg_family].generate_reactions(reactants=reactant_Molecule, products=product_Molecule)[0] #just choose the first
        degeneracy = database.kinetics.families[rmg_family].calculate_degeneracy(template_reaction)
    except IndexError as e:
        print(index, rxn)
        print(e)
        continue
    
    #now match to species already existing in reactions.py
    labeled_adjacencies = get_labeled_adj_lists(rxn)
    switch_labels = get_matching_species_label(labeled_adjacencies)

    
    #and you should change the labels to existing labels in the reactions.py file if the species already exists there 
    final_reactants = []
    for reactant in reactants: 
        try:
            final_label = switch_labels[reactant]
        except KeyError as e: 
            print(index, rxn)
            print(e)
            final_label = reactant
        final_reactants.append(final_label)

    equation_string = ''
    for reactant in final_reactants: 
        equation_string+=f'{reactant} + '
    equation_string=equation_string[:-3] + ' <=> ' #remove the last bit when the loop finishes and add a '<=>'


    final_products = []
    for product in products: 
        try:
            final_label = switch_labels[product]
        except KeyError as e: 
            final_label = product
        final_products.append(final_label)

    #now add to the equation string
    for product in final_products: 
        equation_string+=f'{product} + '
    equation_string=equation_string[:-3] #remove the last bit when the loop finishes and add a '<=>'

    f_string = f'''
entry(
    index = {continued_count},
    label = "{equation_string}",
    degeneracy = {degeneracy},
    kinetics = Arrhenius(A=({A},'cm^3/(mol*s)'), n={n}, Ea=({Ea},'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K')),
    rank = 5,
    shortDesc = """ANL0 method from Brown""",
    longDesc = 
"""
Electronic structures done using ANL0 compound method. 
Torsional scans with  M06-2X/cc-pVTZ. 
""",
)
    '''
    print(f_string)