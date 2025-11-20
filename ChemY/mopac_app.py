import streamlit as st
import tempfile
import subprocess
import os
import re

# --- Funciones auxiliares para MOPAC ---

def run_mopac(mop_file_path, mopac_executable="MOPAC6J.exe"):
    """
    Ejecuta MOPAC6J.exe sobre un archivo .mop.
    """
    try:
        # Añadimos un mensaje de depuración
        st.write(f"🔍 Intentando ejecutar: {mopac_executable} {mop_file_path}")
        # Verificamos si el archivo .mop existe
        if not os.path.exists(mop_file_path):
            st.error(f"❌ El archivo .mop temporal no existe: {mop_file_path}")
            return -1, "", "Archivo .mop no encontrado"
        # Verificamos si el ejecutable existe (opcional, pero útil)
        # import shutil
        # if not shutil.which(mopac_executable):
        #     st.error(f"❌ El ejecutable {mopac_executable} no está en el PATH.")
        #     return -1, "", "Executable not found in PATH"

        result = subprocess.run([mopac_executable, mop_file_path], capture_output=True, text=True, timeout=300) # Añadido timeout
        # Añadimos otro mensaje
        st.write(f"✅ MOPAC ejecutado. Código de retorno: {result.returncode}")
        return result.returncode, result.stdout, result.stderr
    except FileNotFoundError:
        st.error(f"❌ No se encontró el ejecutable {mopac_executable}. Asegúrate de que esté instalado y en el PATH o en la ruta correcta.")
        return -1, "", "Executable not found"
    except subprocess.TimeoutExpired:
        st.error("❌ La ejecución de MOPAC tardó demasiado y se interrumpió.")
        return -1, "", "Timeout"
    except Exception as e:
        st.error(f"❌ Ocurrió un error inesperado al ejecutar MOPAC: {e}")
        return -1, "", str(e)

def parse_mopac_output(out_content):
    """
    Extrae información específica del archivo .out de MOPAC.
    """
    results = {
        "charges": None,
        "homo_lumo": None,
        "error": None
    }

    try:
        # --- Cargas Netas ---
        # Buscar la sección de cargas atómicas
        # El patrón puede ser más robusto, pero este es un inicio
        # Busca la línea que dice "NET ATOMIC CHARGES" y captura las siguientes líneas de átomos
        # Patrón para la sección de cargas
        charge_pattern = r"NET ATOMIC CHARGES AND DIPOLE CONTRIBUTIONS\s*\n\s*ATOM NO\.\s+TYPE\s+CHARGE.*?\n((?:\s*\d+\s+[A-Z]\s+[-+]?\d+\.\d+\s+\d+\.\d+\n?)+)"
        charge_match = re.search(charge_pattern, out_content)
        if charge_match:
            results["charges"] = charge_match.group(1).strip()
            # st.write("DEBUG: Se encontraron cargas.") # Mensaje de depuración
        else:
            # st.warning("DEBUG: No se encontraron cargas netas con el patrón estándar.") # Mensaje de depuración
            results["error"] = "No se encontraron cargas netas."

        # --- Energías HOMO/LUMO ---
        # Buscar el número de orbitales llenos
        filled_levels_match = re.search(r"NO\. OF FILLED LEVELS = (\d+)", out_content)
        if filled_levels_match:
            num_filled = int(filled_levels_match.group(1))
            homo_orbital = num_filled
            lumo_orbital = num_filled + 1

            # Buscar la tabla de EIGENVECTORS y extraer energías
            # El patrón busca bloques de 6 columnas de orbitales (ROOT NO. X X X X X X)
            eigen_pattern = r"ROOT NO\.\s+(.+?)\n((?:\s*[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\n?)+)"
            all_eigenblocks = re.findall(eigen_pattern, out_content)

            if all_eigenblocks:
                # Recorrer bloques para encontrar HOMO y LUMO
                found_homo = None
                found_lumo = None
                current_orbital = 1
                for block_header, block_data in all_eigenblocks:
                    orbital_numbers = list(map(int, block_header.split()))
                    num_orbitals_in_block = len(orbital_numbers)
                    if current_orbital <= homo_orbital < current_orbital + num_orbitals_in_block:
                        # HOMO está en este bloque
                        energies = list(map(float, block_data.split()))
                        idx_in_block = homo_orbital - current_orbital
                        found_homo = energies[idx_in_block]
                    if current_orbital <= lumo_orbital < current_orbital + num_orbitals_in_block:
                        # LUMO está en este bloque
                        energies = list(map(float, block_data.split()))
                        idx_in_block = lumo_orbital - current_orbital
                        found_lumo = energies[idx_in_block]
                    # Avanzar al siguiente bloque
                    current_orbital += num_orbitals_in_block
                    if found_homo is not None and found_lumo is not None:
                        break # Ya encontramos ambos

                if found_homo is not None and found_lumo is not None:
                    results["homo_lumo"] = {
                        "HOMO": f"{found_homo:.5f}",
                        "LUMO": f"{found_lumo:.5f}",
                        "Orbital HOMO": homo_orbital,
                        "Orbital LUMO": lumo_orbital
                    }
                elif found_homo is not None:
                    results["homo_lumo"] = {
                        "HOMO": f"{found_homo:.5f}",
                        "LUMO": "N/A (No encontrado)",
                        "Orbital HOMO": homo_orbital,
                        "Orbital LUMO": lumo_orbital
                    }
                elif found_lumo is not None:
                    results["homo_lumo"] = {
                        "HOMO": "N/A (No encontrado)",
                        "LUMO": f"{found_lumo:.5f}",
                        "Orbital HOMO": homo_orbital,
                        "Orbital LUMO": lumo_orbital
                    }
                else:
                    results["homo_lumo"] = {"HOMO": "N/A", "LUMO": "N/A"}
                    if not results["error"]:
                        results["error"] = f"No se encontraron energías para los orbitales HOMO ({homo_orbital}) o LUMO ({lumo_orbital})."
            else:
                results["homo_lumo"] = {"HOMO": "N/A", "LUMO": "N/A"}
                if not results["error"]:
                    results["error"] = "No se encontró la tabla de EIGENVECTORS."
        else:
            results["homo_lumo"] = {"HOMO": "N/A", "LUMO": "N/A"}
            if not results["error"]:
                results["error"] = "No se pudo determinar el número de orbitales llenos (NO. OF FILLED LEVELS)."


    except Exception as e:
        results["error"] = f"Error al parsear el archivo .out: {e}"

    return results

# --- Interfaz de la página de MOPAC ---
def show_mopac_page():
    st.title("🧪 Procesamiento con MOPAC6J")
    st.markdown("**Creado por:** Antonio Elias Sánchez Soto")

    with st.expander("📖 Guía de Usuario para MOPAC", expanded=False):
        st.markdown("""
        ### 📥 Cómo usar esta herramienta

        1. **Cargar archivo .zmt**:
           - Sube un archivo con extensión `.zmt` que contenga una geometría molecular optimizada (por ejemplo, desde HyperChem).
           - Este archivo se usará como base para generar la entrada de MOPAC.

        2. **Seleccionar método semiempírico**:
           - Elige entre `PM3`, `AM1` o `RM1`.

        3. **Añadir keywords opcionales**:
           - Puedes agregar otras keywords para el cálculo (por ejemplo, `PRECISE`, `VECTORS`, `CHARGE=0`, `UHF`, etc.).

        4. **Generar archivo .mop**:
           - La aplicación creará el archivo de entrada `.mop` con las keywords seleccionadas.

        5. **Ejecutar MOPAC**:
           - Haz clic en "Ejecutar MOPAC". La aplicación llamará al ejecutable `MOPAC6J.exe`.
           - Se generará un archivo de salida `.out`.

        6. **Ver resultados**:
           - Tras la ejecución, la aplicación leerá el archivo `.out` y mostrará:
             - Cargas netas sobre átomos.
             - Energías del HOMO y LUMO.
        """)

    # Subir archivo .zmt
    zmt_file = st.file_uploader("Sube tu archivo .zmt", type=["zmt"])
    if zmt_file is not None:
        st.success("✅ Archivo .zmt cargado.")

        # Seleccionar método
        method = st.selectbox("Selecciona el método semiempírico", ["PM3", "AM1", "RM1"])

        # Keywords adicionales
        keywords_input = st.text_input("Keywords adicionales (opcional)", value="PRECISE VECTORS")

        # Keywords finales
        keywords = f"{method} {keywords_input}".strip()

        # Botón para generar archivo .mop
        if st.button("⚙️ Generar archivo .mop"):
            # Leer contenido del archivo .zmt
            zmt_content = zmt_file.getvalue().decode("utf-8")
            # Generar archivo .mop
            mop_content = f"{keywords}\n\n" + zmt_content

            st.download_button(
                label="📥 Descargar archivo .mop generado",
                data=mop_content.encode("utf-8"),
                file_name="entrada.mop",
                mime="chemical/x-mopac-input"
            )

            # Guardar archivo .mop temporalmente para ejecución
            with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".mop") as f:
                f.write(mop_content)
                temp_mop_path = f.name
            st.write(f"📁 Archivo .mop temporal guardado en: {temp_mop_path}") # Mensaje de depuración

            # Botón para ejecutar MOPAC
            if st.button("▶️ Ejecutar MOPAC6J"):
                st.write("⏳ Ejecutando MOPAC...") # Mensaje de estado
                # Añadimos un mensaje *antes* de llamar a la función
                st.write("🔍 Llamando a la función run_mopac...")
                return_code, stdout, stderr = run_mopac(temp_mop_path)

                # Este mensaje se mostrará siempre que se llame a run_mopac
                st.write("🔍 Función run_mopac finalizada.")

                if return_code == 0:
                    st.success("✅ MOPAC ejecutado correctamente.")
                    # st.write("STDOUT de MOPAC (opcional):") # Mostrar salida estándar
                    # st.text(stdout)
                    if stderr: # Mostrar errores estándar si existen
                        st.warning("STDERR de MOPAC (opcional):")
                        st.text(stderr)

                    # Leer archivo .out
                    out_file_path = temp_mop_path.replace(".mop", ".out")
                    st.write(f"🔍 Buscando archivo .out en: {out_file_path}") # Mensaje de depuración
                    if os.path.exists(out_file_path):
                        st.write("✅ Archivo .out encontrado.")
                        with open(out_file_path, "r") as f:
                            out_content = f.read()

                        # Parsear resultados
                        parsed_results = parse_mopac_output(out_content)

                        if parsed_results["error"]:
                            st.warning(f"⚠️ {parsed_results['error']}")

                        if parsed_results["charges"]:
                            st.subheader("Cargas Netas sobre Átomos")
                            st.text(parsed_results["charges"])

                        if parsed_results["homo_lumo"]:
                            st.subheader("Energías Orbitales (HOMO / LUMO)")
                            homo_val = parsed_results["homo_lumo"].get("HOMO", "N/A")
                            lumo_val = parsed_results["homo_lumo"].get("LUMO", "N/A")
                            orb_homo = parsed_results["homo_lumo"].get("Orbital HOMO", "N/A")
                            orb_lumo = parsed_results["homo_lumo"].get("Orbital LUMO", "N/A")
                            st.write(f"**Orbital HOMO:** {orb_homo} (Energía: {homo_val} eV)")
                            st.write(f"**Orbital LUMO:** {orb_lumo} (Energía: {lumo_val} eV)")


                        # Mostrar el archivo .out completo (opcional, útil para debug)
                        # st.subheader("Contenido completo del archivo .out")
                        # st.text(out_content)

                    else:
                        st.error("❌ No se encontró el archivo .out generado por MOPAC.")
                        st.write(f"Verifique si MOPAC creó el archivo en: {out_file_path}")
                        # Opcional: listar archivos en el directorio temporal
                        # import pathlib
                        # st.write("Archivos en el directorio actual:")
                        # st.write(os.listdir(pathlib.Path(temp_mop_path).parent))

                else:
                    st.error(f"❌ MOPAC falló con código: {return_code}")
                    if stderr:
                        st.text("stderr de MOPAC:")
                        st.text(stderr)
                    if stdout:
                        st.text("stdout de MOPAC:")
                        st.text(stdout)

            # Opcional: Borrar archivos temporales después de la ejecución
            # try:
            #     os.remove(temp_mop_path)
            #     out_file_path = temp_mop_path.replace(".mop", ".out")
            #     if os.path.exists(out_file_path):
            #         os.remove(out_file_path)
            # except:
            #     pass # Ignorar errores al borrar

    else:
        st.warning("⚠️ Sube un archivo .zmt para comenzar.")
