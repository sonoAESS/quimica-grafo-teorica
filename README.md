# 🧬 Proyectos de Química Grafo-Teórica

Este repositorio recopila proyectos dedicados al estudio computacional de moléculas mediante técnicas de teoría de grafos, descriptores químicos y aprendizaje automático. Su estructura modular facilita la investigación, el aprendizaje y el análisis reproducible en bioinformática y química computacional.

---

## 📂 Estructura de Carpetas

```
quimica-grafo-teorica/
├── scripts/
│ └── [Ficheros de conversión, extracción y procesamiento automatizado]
├── notebooks/
│ └── [Jupyter notebooks con actividades investigativas y ML]
├── results/
│ └── [Imágenes, gráficos y archivos generados en los estudios]
├── README.md
```

---

## ✨ Instalación de requisitos y preparación del entorno

### 1. Clonación del repositorio

```
git clone https://tu-repo-url/quimica-grafo-teorica.git
cd quimica-grafo-teorica
```

### 2. Creación y activación de un entorno virtual

> Recomendado para aislar las dependencias y facilitar reproducibilidad

```
python -m venv env
source env/bin/activate # En Linux/Mac
env\Scripts\activate # En Windows
```

### 3. Instalación de dependencias principales

#### Requisitos:

- **Python ≥ 3.8** 🐍
- **RDKit** 🔬 - Quimioinformática y manipulación molecular
- **matplotlib** 📊 - Visualización de datos
- **seaborn** 🌊 - Visualización estadística
- **pandas** 🐼 - Manipulación y análisis de datos estructurados
- **numpy** 🔢 - Cálculo científico

#### Instalación rápida (con `pip`):

```
pip install rdkit matplotlib seaborn pandas numpy
```

> RDKit puede requerir instalación desde fuentes oficiales si tu sistema operativo es poco común  
> Para instrucciones específicas:  
> [RDKit Installation Guide](https://www.rdkit.org/docs/Install.html)

---

## 🏗️ Descripción de carpetas

- **`scripts/`**  
  Scripts funcionales para conversión de formatos químicos (SMILES, MOL, SDF), extracción de descriptores, generación de fingerprints, etc.
- **`notebooks/`**  
  Actividades prácticas computacionales, análisis exploratorios y estudios de machine learning en química molecular.
- **`results/`**  
  Imágenes moleculares, gráficos estadísticos y archivos resultado de los análisis desarrollados.

---

## 📚 Ejemplo de uso

Para ejecutar un notebook, activa el entorno y abre Jupyter Lab o Notebook:

```
source env/bin/activate
jupyter lab
```

---

## ⚗️ Recomendaciones académicas

- Mantén los scripts bien documentados y modularizados para facilitar la extensión.
- Utiliza los notebooks para registrar prácticas de laboratorio computacional y experimentos reproducibles.
- Guarda siempre los resultados relevantes en la carpeta `results/` para trazabilidad.

---

## 💬 Contacto y contribuciones

¿Quieres colaborar, reportar un error o sugerir nuevas actividades?  
Abre un *Issue* o envía un *Pull Request* 👉

---

¡Que la química computacional y la grafo-teoría te acompañen! 🧪✨
