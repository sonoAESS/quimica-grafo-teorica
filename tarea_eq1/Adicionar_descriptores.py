import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, MolSurf, GraphDescriptors
from rdkit.Chem import rdMolDescriptors
import warnings
warnings.filterwarnings('ignore')

def calcular_descriptores_1d_2d(smiles):
    """
HERRAMIENTAS USADAS:
    IDE: Spyder (Python 3.88)
    BIBLIOTECAS: rdkit, pandas, numpy, warnings
    

******************************************************************************
    DESCRIPTORES 1D (Físico-Químicos y de Composición)
******************************************************************************
Composición Atómica:

PesoMolecular - Peso molecular total
NumAtomos - Número total de átomos
NumEnlaces - Número total de enlaces
NumAtomosCarbono - Número de átomos de carbono
NumAtomosOxigeno - Número de átomos de oxígeno
NumAtomosNitrogeno - Número de átomos de nitrógeno
NumAtomosAzufre - Número de átomos de azufre
NumAtomosHidrogeno - Número de átomos de hidrógeno
NumAtomosHalogenos - Número de átomos de halógenos (F, Cl, Br, I)
NumHeteroatomos - Número total de heteroátomos (no C ni H)
NumAtomosPesados - Número de átomos no hidrógeno



Regla de los 5 de Lipinski:
    
HDonadores - Número de donadores de enlaces H
HAceptores - Número de aceptores de enlaces H
LogP - Coeficiente de partición octanol-agua
TPSA - Área de superficie polar topológica
NumRotables - Número de enlaces rotables



Fracciones Atómicas:

FraccionC - Fracción de átomos de carbono
FraccionO - Fracción de átomos de oxígeno
FraccionN - Fracción de átomos de nitrógeno
FraccionS - Fracción de átomos de azufre
FraccionH - Fracción de átomos de hidrógeno
FraccionHeteroatomos - Fracción de heteroátomos







******************************************************************************
DESCRIPTORES 2D (Topológicos y Estructurales)
******************************************************************************

Estructura y Anillos:
    
NumAnillosAromaticos - Número de anillos aromáticos
NumAnillos - Número total de anillos
FraccionSP3 - Fracción de carbonos con hibridación SP3
Forma y Tamaño Molecular:
RadioGiro - Radio de giro molecular
AreaSuperficie - Área de superficie molecular (Labute ASA)



Propiedades Físicas:

IndiceRefractivo - Índice de refractividad molar

******************************************************************************
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    descriptores = {}
    
    # DESCRIPTORES 1D (Básicos)
    descriptores['PesoMolecular'] = Descriptors.MolWt(mol)
    descriptores['NumAtomos'] = mol.GetNumAtoms()
    descriptores['NumEnlaces'] = mol.GetNumBonds()
    descriptores['NumAtomosCarbono'] = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    descriptores['NumAtomosOxigeno'] = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8])
    descriptores['NumAtomosNitrogeno'] = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7])
    descriptores['NumAtomosAzufre'] = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16])
    descriptores['NumAtomosHidrogeno'] = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1])
    descriptores['NumHeteroatomos'] = Lipinski.NumHeteroatoms(mol)
    
    # DESCRIPTORES DE LIPINSKI
    descriptores['HDonadores'] = Lipinski.NumHDonors(mol)
    descriptores['HAceptores'] = Lipinski.NumHAcceptors(mol)
    descriptores['LogP'] = Descriptors.MolLogP(mol)
    descriptores['TPSA'] = Descriptors.TPSA(mol)
    descriptores['NumRotables'] = Lipinski.NumRotatableBonds(mol)
    
    # DESCRIPTORES 2D (Topológicos)
    descriptores['NumAnillosAromaticos'] = Lipinski.NumAromaticRings(mol)
    
    # Fracción SP3 - versión compatible
    try:
        carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
        if carbon_atoms:
            sp3_carbons = len([atom for atom in carbon_atoms if atom.GetHybridization() == Chem.HybridizationType.SP3])
            descriptores['FraccionSP3'] = sp3_carbons / len(carbon_atoms)
        else:
            descriptores['FraccionSP3'] = 0
    except:
        descriptores['FraccionSP3'] = 0
    
    # Número de anillos
    try:
        ssr = Chem.GetSymmSSSR(mol)
        descriptores['NumAnillos'] = len(ssr)
    except:
        descriptores['NumAnillos'] = 0
    
    # DESCRIPTORES DE FORMA Y TAMAÑO
    try:
        descriptores['RadioGiro'] = Descriptors.RadiusOfGyration(mol)
    except:
        descriptores['RadioGiro'] = 0
    
    try:
        descriptores['AreaSuperficie'] = MolSurf.LabuteASA(mol)
    except:
        descriptores['AreaSuperficie'] = 0
    
    # DESCRIPTORES ADICIONALES
    descriptores['NumAtomosHalogenos'] = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() in [9, 17, 35, 53]])
    
    try:
        descriptores['IndiceRefractivo'] = Descriptors.MolMR(mol)
    except:
        descriptores['IndiceRefractivo'] = 0
    
    # DESCRIPTORES DE HIDRÓGENO
    try:
        descriptores['NumAtomosPesados'] = Descriptors.HeavyAtomCount(mol)
    except:
        descriptores['NumAtomosPesados'] = descriptores['NumAtomos'] - descriptores['NumAtomosHidrogeno']
    
    # DESCRIPTORES DE FRACCIONES ATÓMICAS
    total_atoms = mol.GetNumAtoms()
    if total_atoms > 0:
        descriptores['FraccionC'] = descriptores['NumAtomosCarbono'] / total_atoms
        descriptores['FraccionO'] = descriptores['NumAtomosOxigeno'] / total_atoms
        descriptores['FraccionN'] = descriptores['NumAtomosNitrogeno'] / total_atoms
        descriptores['FraccionS'] = descriptores['NumAtomosAzufre'] / total_atoms
        descriptores['FraccionH'] = descriptores['NumAtomosHidrogeno'] / total_atoms
        descriptores['FraccionHeteroatomos'] = descriptores['NumHeteroatomos'] / total_atoms
    else:
        descriptores['FraccionC'] = descriptores['FraccionO'] = descriptores['FraccionN'] = 0
        descriptores['FraccionS'] = descriptores['FraccionH'] = descriptores['FraccionHeteroatomos'] = 0
    
    return descriptores

# Procesar el dataset
df = pd.read_csv('dataset_qsar_tiofeno.csv')

print(f"Procesando {len(df)} moléculas...")

# Calcular descriptores para cada molécula
descriptores_lista = []
smiles_validos = []
ids_validos = []
ic50_means = []
ic50_stds = []

for idx, row in df.iterrows():
    smiles = row['SMILES']
    descriptores = calcular_descriptores_1d_2d(smiles)
    
    if descriptores is not None:
        descriptores_lista.append(descriptores)
        smiles_validos.append(smiles)
        ids_validos.append(row['id'])
        ic50_means.append(row['IC50_mean'])
        ic50_stds.append(row['IC50_std'])
    else:
        print(f"SMILES inválido: {smiles}")

# Crear DataFrame con los descriptores
if descriptores_lista:
    df_descriptores = pd.DataFrame(descriptores_lista)
    df_descriptores.insert(0, 'SMILES', smiles_validos)
    df_descriptores.insert(0, 'id', ids_validos)
    df_descriptores['IC50_mean'] = ic50_means
    df_descriptores['IC50_std'] = ic50_stds
    
    # Guardar resultados
    archivo_salida = 'dataset_con_descriptores.csv'
    df_descriptores.to_csv(archivo_salida, index=False)
    print(f"✅ Archivo guardado: {archivo_salida}")
    print(f"✅ Moléculas procesadas: {len(descriptores_lista)}")
    print(f"✅ Descriptores calculados: {len(df_descriptores.columns) - 4}")  # Excluyendo id, SMILES, IC50_mean, IC50_std
else:
    print("❌ No se pudieron calcular descriptores")