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

#load in the database
database = RMGDatabase()
database.load(
            path = settings['database.directory'],
            thermo_libraries = [],
            transport_libraries = [],
            reaction_libraries = [],
            seed_mechanisms = [],#['BurkeH2O2inN2','ERC-FoundationFuelv0.9'],
            kinetics_families = 'all',
            kinetics_depositories = ['training'],
            #frequenciesLibraries = self.statmechLibraries,
            depository = False, # Don't bother loading the depository information, as we don't use it
        )

reaction_smiles = [(['[H]', 'O=CF'], ['F', '[CH]=O']),
 (['[H]', 'FC=C(F)F'], ['F', 'F[C]=CF']),
 (['[H]', 'FC=C(F)F'], ['F', '[CH]=C(F)F']),
 (['F[C]=CF', 'FC=C(F)F'], ['[CH]=C(F)F', 'FC=C(F)F']),
 (['O=CF', '[CH]=C(F)F'], ['[CH]=O', 'FC=C(F)F']),
 (['[H]', 'FC(F)F'], ['F', 'F[CH]F']),
 (['F[CH]F', 'O=CF'], ['FC(F)F', '[CH]=O']),
 (['F[CH]F', 'FC=C(F)F'], ['FC(F)F', '[CH]=C(F)F']),
 (['FC(F)F', 'F[C]=CF'], ['F[CH]F', 'FC=C(F)F']),
 (['F[CH]F', 'C#CF'], ['FC(F)F', '[C]#C'])]




def get_kinetics(reactant_smiles, product_smiles): 
    
    [react_1, react_2] = reactant_smiles
    [prod_1, prod_2] = product_smiles
    
    reactant_Molecule = [Molecule(smiles=react_1), Molecule(smiles=react_2)]
    product_Molecule = [Molecule(smiles=prod_1), Molecule(smiles=prod_2),] 


    template_reaction = database.kinetics.families['F_Abstraction'].generate_reactions(reactants=reactant_Molecule, products=product_Molecule) 


    rxn = template_reaction[0]

    family = database.kinetics.families['F_Abstraction']

    list_of_kinetics = family.get_kinetics(rxn, rxn.template)

    Arr = list_of_kinetics[0][0]

    rxn.kinetics = Arr

    rate_coefficient = rxn.get_rate_coefficient(1400,101325)

    return rate_coefficient 


for tuple_ in reaction_smiles: 
    (reactant_smiles, product_smiles)= tuple_
    new_rate_coeff = get_kinetics(reactant_smiles, product_smiles)
    print(new_rate_coeff)