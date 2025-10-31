# app_aspirina.py
import warnings
warnings.filterwarnings('ignore')

import logging
logging.getLogger('streamlit').setLevel(logging.ERROR)

import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Configuración de la página DEBE SER LA PRIMERA instrucción de Streamlit
st.set_page_config(
    page_title="Análisis Computacional de la Aspirina",
    page_icon="💊",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Ahora importamos RDKit después de la configuración de Streamlit
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem, MACCSkeys, DataStructs, Lipinski, QED
    from rdkit.Chem import Draw, ChemicalFeatures
    from rdkit.Chem.Draw import rdMolDraw2D
    from rdkit import RDConfig
    RDKIT_AVAILABLE = True
except ImportError as e:
    st.error(f"Error importando RDKit: {e}")
    RDKIT_AVAILABLE = False

# Configuración de matplotlib para evitar advertencias
plt.rcParams['figure.figsize'] = (10, 6)
plt.rcParams['font.size'] = 12
sns.set_style("whitegrid")

# Estilo CSS personalizado
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        color: #1f77b4;
        text-align: center;
        margin-bottom: 2rem;
    }
    .section-header {
        font-size: 1.8rem;
        color: #2e86ab;
        border-bottom: 2px solid #2e86ab;
        padding-bottom: 0.5rem;
        margin-top: 2rem;
    }
    .code-block {
        background-color: #f5f5f5;
        border-left: 4px solid #2e86ab;
        padding: 1rem;
        margin: 1rem 0;
        border-radius: 0.25rem;
    }
    .descriptor-card {
        background-color: #f8f9fa;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #2e86ab;
        margin: 0.5rem 0;
    }
    .metric-card {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        padding: 1rem;
        border-radius: 0.5rem;
        text-align: center;
    }
</style>
""", unsafe_allow_html=True)

def inicializar_molecula():
    """Inicializa y retorna la molécula de aspirina"""
    try:
        aspirin_smiles = 'CC(=O)OC1=CC=CC=C1C(=O)O'
        mol = Chem.MolFromSmiles(aspirin_smiles)
        if mol is None:
            st.error("Error: No se pudo generar la molécula desde el SMILES proporcionado")
            return None
        return mol
    except Exception as e:
        st.error(f"Error inicializando la molécula: {e}")
        return None

def calcular_descriptores_basicos(mol):
    """Calcula descriptores básicos de constitución"""
    try:
        descriptores = {
            'Masa Molecular (Da)': Descriptors.MolWt(mol),
            'Átomos Totales': mol.GetNumAtoms(),
            'Átomos Pesados': Descriptors.HeavyAtomCount(mol),
            'Enlaces Totales': mol.GetNumBonds(),
            'Anillos Totales': Descriptors.RingCount(mol),
            'Anillos Aromáticos': Descriptors.NumAromaticRings(mol)
        }
        return descriptores
    except Exception as e:
        st.error(f"Error calculando descriptores básicos: {e}")
        return {}

def calcular_composicion_atomica(mol):
    """Calcula la composición atómica de la molécula"""
    try:
        atom_counts = {}
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
        return atom_counts
    except Exception as e:
        st.error(f"Error calculando composición atómica: {e}")
        return {}

def calcular_descriptores_topologicos(mol):
    """Calcula descriptores topológicos"""
    try:
        descriptores = {
            'Chi0n (Orden 0)': Descriptors.Chi0n(mol),
            'Chi1n (Orden 1)': Descriptors.Chi1n(mol),
            'Chi2n (Orden 2)': Descriptors.Chi2n(mol),
            'Kappa1': Descriptors.Kappa1(mol),
            'Kappa2': Descriptors.Kappa2(mol),
            'Kappa3': Descriptors.Kappa3(mol),
            'BertzCT (Complejidad)': Descriptors.BertzCT(mol)
        }
        return descriptores
    except Exception as e:
        st.error(f"Error calculando descriptores topológicos: {e}")
        return {}

def calcular_propiedades_farmacologicas(mol):
    """Calcula propiedades farmacológicas y verifica Regla de Lipinski"""
    try:
        descriptores = {
            'LogP (Partición octanol/agua)': Descriptors.MolLogP(mol),
            'TPSA (Å²)': Descriptors.TPSA(mol),
            'Donadores de H': Lipinski.NumHDonors(mol),
            'Aceptores de H': Lipinski.NumHAcceptors(mol),
            'Enlaces Rotables': Descriptors.NumRotatableBonds(mol),
            'QED (Drug-likeness)': QED.qed(mol)
        }
        
        # Verificación de Regla de Lipinski
        lipinski_limits = {
            'Masa Molecular': 500,
            'LogP': 5,
            'Donadores H': 5,
            'Aceptores H': 10
        }
        
        compliance_data = []
        for prop, limit in lipinski_limits.items():
            if prop == 'Masa Molecular':
                value = Descriptors.MolWt(mol)
            elif prop == 'LogP':
                value = Descriptors.MolLogP(mol)
            elif prop == 'Donadores H':
                value = Lipinski.NumHDonors(mol)
            elif prop == 'Aceptores H':
                value = Lipinski.NumHAcceptors(mol)
            
            compliant = value <= limit
            compliance_data.append({
                'Propiedad': prop,
                'Valor': f"{value:.2f}",
                'Límite': limit,
                'Cumple': '✅' if compliant else '❌'
            })
        
        return descriptores, compliance_data
    except Exception as e:
        st.error(f"Error calculando propiedades farmacológicas: {e}")
        return {}, []

def generar_fingerprints(mol):
    """Genera diferentes tipos de fingerprints moleculares"""
    try:
        from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
        
        # Morgan Fingerprints
        gen = GetMorganGenerator(radius=2, fpSize=1024)
        fp_morgan = gen.GetFingerprint(mol)
        arr_morgan = np.zeros((1024,), dtype=int)
        DataStructs.ConvertToNumpyArray(fp_morgan, arr_morgan)
        
        # MACCS Keys
        fp_maccs = MACCSkeys.GenMACCSKeys(mol)
        arr_maccs = np.zeros((167,), dtype=int)
        DataStructs.ConvertToNumpyArray(fp_maccs, arr_maccs)
        
        # Topological Fingerprints
        fp_topo = Chem.RDKFingerprint(mol, maxPath=7, fpSize=2048)
        arr_topo = np.zeros((2048,), dtype=int)
        DataStructs.ConvertToNumpyArray(fp_topo, arr_topo)
        
        fingerprint_stats = {
            'Morgan': (np.sum(arr_morgan), len(arr_morgan)),
            'MACCS': (np.sum(arr_maccs), len(arr_maccs)),
            'Topological': (np.sum(arr_topo), len(arr_topo))
        }
        
        return fingerprint_stats
    except Exception as e:
        st.error(f"Error generando fingerprints: {e}")
        return {}

def analizar_farmacoforo(mol):
    """Realiza análisis farmacóforo de la molécula"""
    try:
        fdefPath = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
        factory = ChemicalFeatures.BuildFeatureFactory(fdefPath)
        feats = factory.GetFeaturesForMol(mol)
        
        farmacoforo_data = {}
        for f in feats:
            feat_type = f.GetType()
            farmacoforo_data[feat_type] = farmacoforo_data.get(feat_type, 0) + 1
        
        return farmacoforo_data
    except Exception as e:
        st.error(f"Error en análisis farmacóforo: {e}")
        return {}

def main():
    """Función principal de la aplicación"""
    
    # Verificar disponibilidad de RDKit
    if not RDKIT_AVAILABLE:
        st.error("""
        ❌ RDKit no está disponible. Por favor, instálelo usando:
        ```bash
        conda install -c conda-forge rdkit
        ```
        """)
        return

    # Título principal
    st.markdown('<h1 class="main-header">💊 Análisis Computacional de Descriptores Moleculares: Aspirina</h1>', 
                unsafe_allow_html=True)
    
    st.markdown('<h3 class="main-header">Antonio Elias Sánchez Soto</h1>', 
                unsafe_allow_html=True)

    # Introducción
    st.markdown("""
    Esta aplicación presenta un análisis computacional completo de los descriptores moleculares 
    del ácido acetilsalicílico (aspirina) utilizando herramientas de quimioinformática.
    """)

    # Sidebar con información
    with st.sidebar:
        st.header("🔬 Información del Proyecto")
        st.markdown("""
        **Objetivo Académico:**
        - Demostrar el cálculo de descriptores moleculares
        - Analizar propiedades farmacológicas
        - Visualizar características estructurales
        
        **Molécula:** Ácido Acetilsalicílico
        **SMILES:** `CC(=O)OC1=CC=CC=C1C(=O)O`
        **Fórmula:** C₉H₈O₄
        """)
        
        st.header("📚 Librerías Utilizadas")
        st.markdown("""
        - **RDKit**: Quimioinformática y ML
        - **Pandas/NumPy**: Análisis de datos
        - **Matplotlib/Seaborn**: Visualización
        - **Streamlit**: Interfaz web
        """)

    # Inicializar molécula
    mol = inicializar_molecula()
    if mol is None:
        return

    # =============================================================================
    # SECCIÓN 1: CONFIGURACIÓN INICIAL Y ESTRUCTURA MOLECULAR
    # =============================================================================

    st.markdown('<h2 class="section-header">1. Configuración y Estructura Molecular</h2>', 
                unsafe_allow_html=True)

    col1, col2 = st.columns([1, 1])

    with col1:
        st.subheader("🔧 Configuración Inicial")
        
        st.markdown("""
        **Librerías principales utilizadas:**
        - `rdkit.Chem`: Manipulación de estructuras moleculares
        - `rdkit.Chem.Descriptors`: Cálculo de descriptores
        - `rdkit.Chem.AllChem`: Fingerprints y conformaciones
        """)
        
        with st.expander("📝 Código: Importación de librerías"):
            st.code("""
import streamlit as st
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import AllChem, MACCSkeys, DataStructs
from rdkit.Chem import Lipinski, QED, Draw
            """, language='python')

    with col2:
        st.subheader("🧬 Estructura Molecular")
        
        try:
            img = Draw.MolToImage(mol, size=(300, 300))
            st.image(img, caption='Estructura 2D de la Aspirina')
            
            formula = rdMolDescriptors.CalcMolFormula(mol)
            st.markdown(f"**Fórmula Molecular:** {formula}")
            st.markdown(f"**SMILES:** `CC(=O)OC1=CC=CC=C1C(=O)O`")
        except Exception as e:
            st.error(f"Error mostrando estructura molecular: {e}")

    # =============================================================================
    # SECCIÓN 2: DESCRIPTORES BÁSICOS Y DE CONSTITUCIÓN
    # =============================================================================

    st.markdown('<h2 class="section-header">2. Descriptores Básicos y de Constitución</h2>', 
                unsafe_allow_html=True)

    st.markdown("""
    Los descriptores de constitución describen la composición atómica y conectividad básica 
    de la molécula sin considerar la disposición espacial.
    """)

    with st.expander("📝 Código: Cálculo de descriptores básicos"):
        st.code("""
# SMILES de aspirina
aspirin_smiles = 'CC(=O)OC1=CC=CC=C1C(=O)O'
mol = Chem.MolFromSmiles(aspirin_smiles)

# Cálculo de descriptores básicos
mw = Descriptors.MolWt(mol)                    # Masa molecular
n_atoms = mol.GetNumAtoms()                    # Número total de átomos
n_bonds = mol.GetNumBonds()                    # Número total de enlaces
heavy_atoms = Descriptors.HeavyAtomCount(mol)  # Átomos pesados (no H)
rings = Descriptors.RingCount(mol)             # Número de anillos
formula = rdMolDescriptors.CalcMolFormula(mol) # Fórmula química
        """, language='python')

    # Calcular y mostrar descriptores básicos
    descriptores_basicos = calcular_descriptores_basicos(mol)
    atom_counts = calcular_composicion_atomica(mol)
    
    if descriptores_basicos and atom_counts:
        col1, col2 = st.columns([1, 1])
        
        with col1:
            st.subheader("📊 Descriptores Numéricos")
            df_basicos = pd.DataFrame(list(descriptores_basicos.items()), 
                                    columns=['Descriptor', 'Valor'])
            st.dataframe(df_basicos, use_container_width=True)
        
        with col2:
            st.subheader("🧪 Composición Atómica")
            df_composicion = pd.DataFrame(list(atom_counts.items()), 
                                        columns=['Elemento', 'Cantidad'])
            st.dataframe(df_composicion, use_container_width=True)
            
            # Gráfico de composición
            fig, ax = plt.subplots(figsize=(8, 6))
            colors = plt.cm.Set3(np.linspace(0, 1, len(atom_counts)))
            ax.pie(atom_counts.values(), labels=atom_counts.keys(), autopct='%1.1f%%',
                   colors=colors)
            ax.set_title('Composición Atómica', fontweight='bold')
            st.pyplot(fig)

    # =============================================================================
    # SECCIÓN 3: DESCRIPTORES TOPOLÓGICOS
    # =============================================================================

    st.markdown('<h2 class="section-header">3. Descriptores Topológicos</h2>', 
                unsafe_allow_html=True)

    st.markdown("""
    Los índices topológicos describen la conectividad molecular basándose en la teoría de grafos. 
    Estos descriptores son invariantes a la conformación molecular y capturan información 
    sobre la complejidad estructural.
    """)

    with st.expander("📝 Código: Cálculo de descriptores topológicos"):
        st.code("""
# Índices de conectividad de Kier-Hall
chi0 = Descriptors.Chi0n(mol)  # Índice de orden 0
chi1 = Descriptors.Chi1n(mol)  # Índice de orden 1  
chi2 = Descriptors.Chi2n(mol)  # Índice de orden 2

# Índices de forma kappa
kappa1 = Descriptors.Kappa1(mol)  # Índice de forma 1
kappa2 = Descriptors.Kappa2(mol)  # Índice de forma 2
kappa3 = Descriptors.Kappa3(mol)  # Índice de forma 3

# Complejidad molecular
bertz_ct = Descriptors.BertzCT(mol)  # Índice de complejidad de Bertz
        """, language='python')

    descriptores_topo = calcular_descriptores_topologicos(mol)
    
    if descriptores_topo:
        col1, col2 = st.columns([1, 1])
        
        with col1:
            st.subheader("📐 Descriptores Topológicos")
            df_topo = pd.DataFrame(list(descriptores_topo.items()), 
                                columns=['Descriptor', 'Valor'])
            st.dataframe(df_topo, use_container_width=True)
        
        with col2:
            st.subheader("📈 Visualización de Índices")
            
            # Gráfico de índices topológicos
            fig, ax = plt.subplots(figsize=(10, 6))
            
            # Seleccionar solo algunos para el gráfico
            topo_plot = {k: v for k, v in list(descriptores_topo.items())[:4]}
            
            bars = ax.bar(topo_plot.keys(), topo_plot.values(), 
                        color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'])
            ax.set_ylabel('Valor del Índice')
            ax.set_title('Principales Índices Topológicos', fontweight='bold')
            plt.xticks(rotation=45, ha='right')
            
            # Añadir valores en las barras
            for bar in bars:
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                    f'{height:.2f}', ha='center', va='bottom')
            
            st.pyplot(fig)

    # =============================================================================
    # SECCIÓN 4: PROPIEDADES FARMACOLÓGICAS
    # =============================================================================

    st.markdown('<h2 class="section-header">4. Propiedades Farmacológicas y Regla de Lipinski</h2>', 
                unsafe_allow_html=True)

    st.markdown("""
    La **Regla de Lipinski** (Rule of Five) es un conjunto de criterios para evaluar 
    la "drug-likeness" de un compuesto. Un compuesto tiene alta probabilidad de ser 
    un fármaco oralmente activo si cumple con estos criterios:
    """)

    st.markdown("""
    - **Masa Molecular** ≤ 500 Da
    - **LogP** ≤ 5
    - **Donadores de H** ≤ 5  
    - **Aceptores de H** ≤ 10
    - **Enlaces rotables** ≤ 10 (variante)
    """)

    with st.expander("📝 Código: Cálculo de propiedades farmacológicas"):
        st.code("""
# Propiedades farmacológicas clave
logp = Descriptors.MolLogP(mol)                    # Coeficiente de partición octanol/agua
tpsa = Descriptors.TPSA(mol)                       # Superficie polar accesible
h_donors = Lipinski.NumHDonors(mol)                # Donadores de hidrógeno
h_acceptors = Lipinski.NumHAcceptors(mol)          # Aceptores de hidrógeno
rotatable_bonds = Descriptors.NumRotatableBonds(mol) # Enlaces rotables

# Métrica QED (Drug-likeness)
qed_score = QED.qed(mol)                           # Puntuación QED

# Verificación de la Regla de Lipinski
lipinski_violations = 0
if Descriptors.MolWt(mol) > 500: lipinski_violations += 1
if Descriptors.MolLogP(mol) > 5: lipinski_violations += 1  
if Lipinski.NumHDonors(mol) > 5: lipinski_violations += 1
if Lipinski.NumHAcceptors(mol) > 10: lipinski_violations += 1
        """, language='python')

    descriptores_farma, compliance_data = calcular_propiedades_farmacologicas(mol)
    
    if descriptores_farma and compliance_data:
        col1, col2 = st.columns([1, 1])
        
        with col1:
            st.subheader("💊 Propiedades Farmacológicas")
            df_farma = pd.DataFrame(list(descriptores_farma.items()), 
                                columns=['Descriptor', 'Valor'])
            st.dataframe(df_farma, use_container_width=True)
        
        with col2:
            st.subheader("📋 Cumplimiento Regla de Lipinski")
            df_compliance = pd.DataFrame(compliance_data)
            st.dataframe(df_compliance, use_container_width=True)
            
            # Calcular violaciones
            violations = sum(1 for item in compliance_data if item['Cumple'] == '❌')
            st.metric("Violaciones de la Regla de Lipinski", violations)
            
            if violations == 0:
                st.success("✅ La aspirina cumple perfectamente con los criterios de drug-likeness")
            elif violations <= 1:
                st.info("ℹ️ La aspirina presenta una violación menor de la regla de Lipinski")
            else:
                st.warning("⚠️ La aspirina presenta múltiples violaciones de la regla de Lipinski")

    # =============================================================================
    # SECCIÓN 5: FINGERPRINTS MOLECULARES
    # =============================================================================

    st.markdown('<h2 class="section-header">5. Fingerprints Moleculares</h2>', 
                unsafe_allow_html=True)

    st.markdown("""
    Los fingerprints moleculares son representaciones vectoriales de características 
    estructurales que permiten comparaciones cuantitativas entre moléculas. 
    Se utilizan en búsqueda de similitud, QSAR y machine learning.
    """)

    with st.expander("📝 Código: Generación de fingerprints"):
        st.code("""
# Morgan Fingerprints (circular fingerprints)
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
gen = GetMorganGenerator(radius=2, fpSize=1024)
fp_morgan = gen.GetFingerprint(mol)
arr_morgan = np.zeros((1024,), dtype=int)
DataStructs.ConvertToNumpyArray(fp_morgan, arr_morgan)

# MACCS Keys (166 bits estructurales)
fp_maccs = MACCSkeys.GenMACCSKeys(mol)
arr_maccs = np.zeros((167,), dtype=int)
DataStructs.ConvertToNumpyArray(fp_maccs, arr_maccs)

# Topological Fingerprints
fp_topo = Chem.RDKFingerprint(mol, maxPath=7, fpSize=2048)
arr_topo = np.zeros((2048,), dtype=int)
DataStructs.ConvertToNumpyArray(fp_topo, arr_topo)
        """, language='python')

    fingerprint_stats = generar_fingerprints(mol)
    
    if fingerprint_stats:
        col1, col2 = st.columns([1, 1])
        
        with col1:
            st.subheader("🔢 Estadísticas de Fingerprints")
            
            fp_data = []
            for name, (active, total) in fingerprint_stats.items():
                density = active / total
                fp_data.append({
                    'Tipo': name,
                    'Bits Activos': active,
                    'Total Bits': total,
                    'Densidad (%)': f"{density*100:.2f}%"
                })
            
            df_fp = pd.DataFrame(fp_data)
            st.dataframe(df_fp, use_container_width=True)
        
        with col2:
            st.subheader("📊 Densidad de Bits")
            
            fig, ax = plt.subplots(figsize=(10, 6))
            
            names = list(fingerprint_stats.keys())
            densities = [stats[0]/stats[1] for stats in fingerprint_stats.values()]
            
            bars = ax.bar(names, densities, color=['#1f77b4', '#ff7f0e', '#2ca02c'])
            ax.set_ylabel('Densidad de Bits Activos')
            ax.set_title('Comparación de Densidad entre Fingerprints', fontweight='bold')
            ax.set_ylim(0, 0.5)
            
            for bar, density in zip(bars, densities):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                    f'{density:.3f}', ha='center', va='bottom')
            
            st.pyplot(fig)

    # =============================================================================
    # SECCIÓN 6: CARACTERÍSTICAS FARMACÓFORAS
    # =============================================================================

    st.markdown('<h2 class="section-header">6. Características Farmacóforas</h2>', 
                unsafe_allow_html=True)

    st.markdown("""
    El análisis farmacóforo identifica grupos funcionales clave que pueden participar 
    en interacciones moleculares con receptores biológicos. Estas características 
    definen el "farmacóforo" - el conjunto de features estructurales responsables 
    de la actividad biológica.
    """)

    with st.expander("📝 Código: Análisis farmacóforo"):
        st.code("""
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig
import os

# Cargar definición de características
fdefPath = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
factory = ChemicalFeatures.BuildFeatureFactory(fdefPath)

# Extraer características farmacóforas
feats = factory.GetFeaturesForMol(mol)

# Contar características por tipo
farmacoforo_data = {}
for f in feats:
    feat_type = f.GetType()
    farmacoforo_data[feat_type] = farmacoforo_data.get(feat_type, 0) + 1
        """, language='python')

    farmacoforo_data = analizar_farmacoforo(mol)
    
    if farmacoforo_data:
        col1, col2 = st.columns([1, 1])
        
        with col1:
            st.subheader("🎯 Características Identificadas")
            
            df_farmacoforo = pd.DataFrame(list(farmacoforo_data.items()), 
                                        columns=['Característica', 'Cantidad'])
            st.dataframe(df_farmacoforo, use_container_width=True)
        
        with col2:
            st.subheader("📈 Distribución de Features")
            
            fig, ax = plt.subplots(figsize=(10, 6))
            
            features = list(farmacoforo_data.keys())
            counts = list(farmacoforo_data.values())
            
            bars = ax.bar(features, counts, color=sns.color_palette("Set2", len(features)))
            ax.set_ylabel('Número de Características')
            ax.set_title('Distribución de Características Farmacóforas', fontweight='bold')
            plt.xticks(rotation=45, ha='right')
            
            for bar, count in zip(bars, counts):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height + 0.05,
                    str(count), ha='center', va='bottom')
            
            st.pyplot(fig)
    else:
        st.info("No se detectaron características farmacóforas específicas")

    # =============================================================================
    # SECCIÓN 7: CONCLUSIONES
    # =============================================================================

    st.markdown('<h2 class="section-header">7. Conclusiones</h2>', 
                unsafe_allow_html=True)
    
    st.markdown("""
    **Hallazgos Principales:**
    
    1. **Estructura y Composición**: La aspirina presenta una estructura aromática 
       con grupos funcionales carboxilo y éster que definen sus propiedades.
    
    2. **Drug-likeness**: Cumple favorablemente con la Regla de Lipinski, 
       explicando su buena biodisponibilidad oral.
    
    3. **Propiedades Físico-Químicas**: LogP moderado y TPSA adecuado sugieren 
       buen balance hidrofilia/lipofilia.
    
    4. **Fingerprints Moleculares**: Los diferentes tipos de fingerprints capturan 
       aspectos complementarios de la estructura molecular.
    
    **Aplicaciones en Investigación:**
    - Diseño de análogos estructurales
    - Estudios QSAR y modelado predictivo
    - Búsqueda de similitud molecular en bases de datos
    - Optimización de propiedades ADMET
    """)

    # =============================================================================
    # Pie de página
    # =============================================================================

    st.markdown("---")
    st.markdown("""
    **🔬 Herramientas Computacionales en Química Farmacéutica**  
    *Análisis realizado con RDKit, Streamlit y Python para fines académicos*
    """)

if __name__ == "__main__":
    main()