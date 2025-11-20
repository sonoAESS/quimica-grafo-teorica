import os
import subprocess
import shutil
from pathlib import Path

class MOPACRunner:
    """
    Clase para automatizar la ejecución de cálculos con MOPAC v6
    """
    def __init__(self, mopac_executable_path="mopac6j.exe"):
        """
        Inicializa el runner de MOPAC
        :param mopac_executable_path: Ruta al ejecutable de MOPAC (mopac6j.exe)
        """
        self.mopac_executable = mopac_executable_path
        self._validate_executable()

    def _validate_executable(self):
        """
        Valida que el ejecutable de MOPAC esté disponible
        """
        # Verificar si el ejecutable existe en la carpeta actual
        if not Path(self.mopac_executable).exists():
            raise FileNotFoundError(f"No se encontró el ejecutable {self.mopac_executable} en la carpeta actual")

    def create_mop_file(self, zmt_file_path, method="PM3", keywords="PRECISE VECTORS", comment=""):
        """
        Crea un archivo .mop a partir de un archivo .zmt
        :param zmt_file_path: Ruta al archivo .zmt
        :param method: Método semiempírico (PM3, AM1, RM1, etc.)
        :param keywords: Palabras clave adicionales
        :param comment: Comentario para el archivo
        :return: Ruta al archivo .mop creado
        """
        zmt_path = Path(zmt_file_path)
        
        if not zmt_path.exists():
            raise FileNotFoundError(f"Archivo .zmt no encontrado: {zmt_file_path}")

        # Leer el contenido del archivo .zmt
        with open(zmt_path, 'r') as f:
            zmt_content = f.read()

        # Crear el contenido del archivo .mop
        mop_content = f"{method} {keywords}\n"
        mop_content += f"{comment}\n"
        mop_content += "\n"  # Línea en blanco
        mop_content += zmt_content

        # Guardar archivo .mop
        mop_file_path = zmt_path.with_suffix('.mop')
        with open(mop_file_path, 'w') as f:
            f.write(mop_content)

        print(f"Archivo .mop creado: {mop_file_path}")
        return mop_file_path

    def run_mopac(self, mop_file_path):
        """
        Ejecuta MOPAC con el archivo .mop especificado
        :param mop_file_path: Ruta al archivo .mop
        :return: Ruta al archivo .out generado
        """
        mop_path = Path(mop_file_path)
        
        if not mop_path.exists():
            raise FileNotFoundError(f"Archivo .mop no encontrado: {mop_file_path}")

        # Verificar que el ejecutable existe antes de intentar ejecutar
        if not Path(self.mopac_executable).exists():
            raise FileNotFoundError(f"El ejecutable de MOPAC no se encuentra: {self.mopac_executable}")

        # El archivo de salida de MOPAC tendrá el mismo nombre base que el .mop
        expected_out_path = mop_path.with_suffix('.out')
        
        # Eliminar archivo de salida anterior si existe
        if expected_out_path.exists():
            expected_out_path.unlink()
            print(f"Archivo .out anterior eliminado: {expected_out_path}")
        
        print(f"Ejecutando MOPAC con: {mop_path}")
        print(f"Esperando archivo de salida: {expected_out_path}")

        # Ejecutar MOPAC en el directorio actual
        try:
            # Cambiar al directorio del archivo .mop
            original_dir = os.getcwd()
            os.chdir(mop_path.parent)
            
            # Ejecutar MOPAC
            cmd = [self.mopac_executable, str(mop_path.name)]
            print(f"Comando ejecutado: {' '.join(cmd)}")
            
            # Ejecutar MOPAC y capturar salida
            result = subprocess.run(cmd, 
                                  capture_output=True, 
                                  text=True, 
                                  timeout=300,  # Timeout de 5 minutos
                                  check=False)  # No lanzar excepción si hay error
            
            # Volver al directorio original
            os.chdir(original_dir)
            
            # Mostrar la salida de MOPAC
            print("Salida de MOPAC (stdout):")
            if result.stdout:
                print(result.stdout)
            else:
                print("(stdout vacío)")
            
            print("Errores de MOPAC (stderr):")
            if result.stderr:
                print(result.stderr)
            else:
                print("(stderr vacío)")
            
            print(f"Código de retorno de MOPAC: {result.returncode}")
            
            # Verificar si se generó el archivo de salida
            if expected_out_path.exists():
                # Verificar si el archivo está vacío
                if expected_out_path.stat().st_size == 0:
                    print(f"ADVERTENCIA: El archivo de salida {expected_out_path} se generó pero está vacío.")
                    print("Esto puede indicar un problema con el archivo de entrada .mop o con la ejecución de MOPAC.")
                    return expected_out_path
                else:
                    print(f"MOPAC ejecutado exitosamente. Archivo generado: {expected_out_path}")
                    print(f"Tamaño del archivo: {expected_out_path.stat().st_size} bytes")
                    return expected_out_path
            else:
                print(f"ERROR: No se generó el archivo de salida esperado: {expected_out_path}")
                print(f"Contenido del directorio {mop_path.parent}:")
                for file in mop_path.parent.iterdir():
                    print(f"  - {file.name}")
                return None
                
        except subprocess.TimeoutExpired:
            print(f"ERROR: La ejecución de MOPAC excedió el tiempo límite de 5 minutos.")
            os.chdir(original_dir)
            return None
        except Exception as e:
            print(f"ERROR: Error ejecutando MOPAC: {e}")
            import traceback
            traceback.print_exc()
            os.chdir(original_dir)
            return None

    def extract_key_data(self, out_file_path):
        """
        Extrae datos clave del archivo .out basado en el formato real de MOPAC v6
        :param out_file_path: Ruta al archivo .out
        :return: Diccionario con datos extraídos
        """
        if out_file_path is None:
            print("Advertencia: No hay archivo de salida para extraer datos.")
            return {
                'homo_energy': None,
                'homo_level': None,
                'lumo_energy': None,
                'lumo_level': None,
                'net_charges': [],
                'atomic_coordinates': [],
                'all_orbital_energies': [],
                'ionization_potential': None,
                'heat_of_formation': None,
                'total_energy': None,
                'electronic_energy': None,
                'core_core_repulsion': None,
                'molecular_weight': None,
                'scf_calculations': None
            }
        
        out_path = Path(out_file_path)
        
        if not out_path.exists():
            print(f"Advertencia: El archivo .out no se generó: {out_file_path}")
            return {
                'homo_energy': None,
                'homo_level': None,
                'lumo_energy': None,
                'lumo_level': None,
                'net_charges': [],
                'atomic_coordinates': [],
                'all_orbital_energies': [],
                'ionization_potential': None,
                'heat_of_formation': None,
                'total_energy': None,
                'electronic_energy': None,
                'core_core_repulsion': None,
                'molecular_weight': None,
                'scf_calculations': None
            }

        # Verificar si el archivo está vacío
        if out_path.stat().st_size == 0:
            print(f"Advertencia: El archivo .out está vacío: {out_file_path}")
            print("No se pueden extraer datos de un archivo vacío.")
            return {
                'homo_energy': None,
                'homo_level': None,
                'lumo_energy': None,
                'lumo_level': None,
                'net_charges': [],
                'atomic_coordinates': [],
                'all_orbital_energies': [],
                'ionization_potential': None,
                'heat_of_formation': None,
                'total_energy': None,
                'electronic_energy': None,
                'core_core_repulsion': None,
                'molecular_weight': None,
                'scf_calculations': None
            }

        data = {
            'homo_energy': None,
            'homo_level': None,
            'lumo_energy': None,
            'lumo_level': None,
            'net_charges': [],
            'atomic_coordinates': [],
            'all_orbital_energies': [],
            'ionization_potential': None,
            'heat_of_formation': None,
            'total_energy': None,
            'electronic_energy': None,
            'core_core_repulsion': None,
            'molecular_weight': None,
            'scf_calculations': None
        }

        with open(out_path, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()

        # Extraer valores generales
        for line in lines:
            if 'IONIZATION POTENTIAL' in line and '=' in line:
                try:
                    value = float(line.split('=')[1].strip().split()[0])
                    data['ionization_potential'] = value
                except (ValueError, IndexError):
                    continue
            elif 'NO. OF FILLED LEVELS' in line and '=' in line:
                try:
                    filled_levels = int(line.split('=')[1].strip())
                    print(f"Número de niveles llenos detectado: {filled_levels}")
                except (ValueError, IndexError):
                    filled_levels = 12  # Valor por defecto
                    continue
            elif 'FINAL HEAT OF FORMATION' in line and '=' in line:
                try:
                    value = float(line.split('=')[1].strip().split()[0])
                    data['heat_of_formation'] = value
                except (ValueError, IndexError):
                    continue
            elif 'TOTAL ENERGY' in line and '=' in line:
                try:
                    value = float(line.split('=')[1].strip().split()[0])
                    data['total_energy'] = value
                except (ValueError, IndexError):
                    continue
            elif 'ELECTRONIC ENERGY' in line and '=' in line:
                try:
                    value = float(line.split('=')[1].strip().split()[0])
                    data['electronic_energy'] = value
                except (ValueError, IndexError):
                    continue
            elif 'CORE-CORE REPULSION' in line and '=' in line:
                try:
                    value = float(line.split('=')[1].strip().split()[0])
                    data['core_core_repulsion'] = value
                except (ValueError, IndexError):
                    continue
            elif 'MOLECULAR WEIGHT' in line and '=' in line:
                try:
                    value = float(line.split('=')[1].strip())
                    data['molecular_weight'] = value
                except (ValueError, IndexError):
                    continue
            elif 'SCF CALCULATIONS' in line and '=' in line:
                try:
                    value = int(line.split('=')[1].strip())
                    data['scf_calculations'] = value
                except (ValueError, IndexError):
                    continue

        # Buscar energías orbitales en la sección de EIGENVECTORS
        orbital_energies = []
        current_line_idx = 0
        
        while current_line_idx < len(lines):
            line = lines[current_line_idx]
            if 'EIGENVECTORS' in line:
                # Buscar bloques de ROOT NO. y las energías en la línea siguiente
                j = current_line_idx + 1
                while j < len(lines):
                    if 'ROOT NO.' in lines[j]:
                        # Extraer energías de la línea siguiente
                        energies_line = lines[j + 1].strip()
                        if energies_line:
                            # Dividir por espacios y convertir a float
                            parts = energies_line.split()
                            for part in parts:
                                try:
                                    energy = float(part)
                                    orbital_energies.append(energy)
                                except ValueError:
                                    continue
                    elif 'S  C' in lines[j]:  # Comienza la parte de coeficientes
                        break
                    j += 1
                break  # Solo procesar la primera sección de EIGENVECTORS
            current_line_idx += 1

        # Asignar HOMO y LUMO basado en el número de niveles llenos (por defecto 12 según el archivo de ejemplo)
        filled_levels = 12
        if len(orbital_energies) >= filled_levels:
            # En el archivo de ejemplo, los orbitales están numerados del 1 al 24
            # El orbital 12 es el HOMO (índice 11 en la lista 0-indexed)
            # El orbital 13 es el LUMO (índice 12 en la lista 0-indexed)
            if len(orbital_energies) >= 12:  # Aseguramos tener al menos 12 orbitales
                data['homo_energy'] = orbital_energies[11]  # Orbital 12
                data['homo_level'] = 12
                if len(orbital_energies) > 12:
                    data['lumo_energy'] = orbital_energies[12]  # Orbital 13
                    data['lumo_level'] = 13
            data['all_orbital_energies'] = orbital_energies

        # Buscar cargas netas (en el formato específico de MOPAC)
        for i, line in enumerate(lines):
            if 'NET ATOMIC CHARGES AND DIPOLE CONTRIBUTIONS' in line:
                # Las cargas están en las líneas siguientes
                j = i + 1
                while j < len(lines):
                    charge_line = lines[j].strip()
                    if charge_line and charge_line.startswith(' '):
                        # Formato: número átomo, tipo, carga, electrones
                        parts = charge_line.split()
                        if len(parts) >= 4:
                            try:
                                atom_num = int(parts[0])
                                atom_type = parts[1]
                                charge = float(parts[2])
                                data['net_charges'].append((atom_num, atom_type, charge))
                            except (ValueError, IndexError):
                                pass
                    elif 'DIPOLE' in lines[j] or j > i + 20:  # Finaliza cuando encuentra DIPOLE o después de 20 líneas
                        break
                    j += 1
                break  # Solo tomar la primera aparición

        # Buscar coordenadas atómicas (la segunda aparición después de la optimización)
        found_coords_section = 0
        for i, line in enumerate(lines):
            if 'CARTESIAN COORDINATES' in line:
                found_coords_section += 1
                if found_coords_section == 2:  # Tomar la segunda aparición (después de la optimización)
                    j = i + 1
                    while j < len(lines):
                        coord_line = lines[j].strip()
                        if coord_line and coord_line.startswith(' '):
                            parts = coord_line.split()
                            if len(parts) >= 5 and parts[0].isdigit():
                                try:
                                    atom_num = int(parts[0])
                                    atom_type = parts[1]
                                    x = float(parts[2])
                                    y = float(parts[3])
                                    z = float(parts[4])
                                    data['atomic_coordinates'].append((atom_num, atom_type, x, y, z))
                                except (ValueError, IndexError):
                                    pass
                        elif 'TOTAL' in lines[j] or j > i + 50:  # Finaliza cuando encuentra TOTAL o después de 50 líneas
                            break
                        j += 1
                    break  # Solo tomar la segunda aparición de coordenadas

        return data

    def process_molecule(self, zmt_file_path, method="PM3", keywords="PRECISE VECTORS", comment=""):
        """
        Procesa completamente una molécula: crear .mop, ejecutar MOPAC, extraer datos
        :param zmt_file_path: Ruta al archivo .zmt
        :param method: Método semiempírico
        :param keywords: Palabras clave
        :param comment: Comentario
        :return: Diccionario con resultados
        """
        print(f"Procesando molécula: {zmt_file_path}")
        
        # Crear archivo .mop
        mop_file = self.create_mop_file(zmt_file_path, method, keywords, comment)
        
        # Ejecutar MOPAC
        out_file = self.run_mopac(mop_file)
        
        # Extraer datos
        results = self.extract_key_data(out_file)
        
        results['input_file'] = str(zmt_file_path)
        results['output_file'] = str(out_file) if out_file else None
        
        print(f"Procesamiento completado para: {zmt_file_path}")
        return results

def main():
    """
    Ejemplo de uso del script - ahora pide el archivo .zmt
    """
    print("=== AUTOMATIZADOR DE CÁLCULOS MOPAC ===")
    
    # Verificar que el ejecutable existe antes de continuar
    if not Path("mopac6j.exe").exists():
        print("ERROR: No se encontró mopac6j.exe en la carpeta actual.")
        print("Asegúrate de que el archivo mopac6j.exe esté en la misma carpeta que este script.")
        return
    
    # Pedir al usuario el archivo .zmt a procesar
    zmt_file = input("Introduce la ruta al archivo .zmt: ").strip()
    
    # Verificar que el archivo existe
    if not Path(zmt_file).exists():
        print(f"ERROR: El archivo {zmt_file} no existe.")
        return
    
    # Pedir método semiempírico
    print("\nMétodos disponibles: PM3, AM1, RM1")
    method = input("Introduce el método semiempírico (PM3 por defecto): ").strip()
    if not method:
        method = "PM3"
    
    # Pedir palabras clave adicionales
    keywords = input("Introduce palabras clave adicionales (PRECISE VECTORS por defecto): ").strip()
    if not keywords:
        keywords = "PRECISE VECTORS"
    
    # Pedir comentario
    comment = input("Introduce un comentario para el archivo .mop: ").strip()
    if not comment:
        comment = f"Cálculo {method} para {Path(zmt_file).stem}"
    
    # Inicializar el runner con el ejecutable correcto
    runner = MOPACRunner(mopac_executable_path="mopac6j.exe")
    
    try:
        # Procesar la molécula
        results = runner.process_molecule(zmt_file, method, keywords, comment)
        
        # Mostrar resultados
        print("\n=== RESULTADOS ===")
        print(f"Molécula: {results['input_file']}")
        print(f"Archivo de salida: {results['output_file']}")
        
        # Si el archivo de salida existe y no está vacío, mostrar los resultados
        if results['output_file'] and Path(results['output_file']).stat().st_size > 0:
            # Datos generales
            if results['ionization_potential'] is not None:
                print(f"Potencial de ionización: {results['ionization_potential']:.5f} eV")
            if results['heat_of_formation'] is not None:
                print(f"Heat of formation: {results['heat_of_formation']:.5f} kcal/mol")
            if results['total_energy'] is not None:
                print(f"Total energy: {results['total_energy']:.5f} eV")
            if results['electronic_energy'] is not None:
                print(f"Electronic energy: {results['electronic_energy']:.5f} eV")
            if results['core_core_repulsion'] is not None:
                print(f"Core-core repulsion: {results['core_core_repulsion']:.5f} eV")
            if results['molecular_weight'] is not None:
                print(f"Molecular weight: {results['molecular_weight']:.3f}")
            if results['scf_calculations'] is not None:
                print(f"SCF calculations: {results['scf_calculations']}")
            
            # Energías HOMO/LUMO
            if results['homo_energy'] is not None:
                print(f"Energía HOMO (nivel {results['homo_level']}): {results['homo_energy']:.5f} eV")
            if results['lumo_energy'] is not None:
                print(f"Energía LUMO (nivel {results['lumo_level']}): {results['lumo_energy']:.5f} eV")
            
            if results['homo_energy'] is not None and results['lumo_energy'] is not None:
                gap = results['lumo_energy'] - results['homo_energy']
                print(f"Gap HOMO-LUMO: {gap:.5f} eV")
            
            if results['all_orbital_energies']:
                print(f"Número total de energías orbitales encontradas: {len(results['all_orbital_energies'])}")
                if results['all_orbital_energies']:
                    print(f"Primera energía: {results['all_orbital_energies'][0]:.5f} eV")
                    print(f"Última energía: {results['all_orbital_energies'][-1]:.5f} eV")
            
            print(f"Número de cargas atómicas: {len(results['net_charges'])}")
            
            # Mostrar las cargas atómicas si existen
            if results['net_charges']:
                print("\nCargas atómicas:")
                for atom_num, atom_type, charge in results['net_charges']:
                    print(f"  Átomo {atom_num} ({atom_type}): {charge:.4f}")
            
            print(f"Número de coordenadas atómicas: {len(results['atomic_coordinates'])}")
            
            # Mostrar las coordenadas atómicas si existen
            if results['atomic_coordinates']:
                print("\nCoordenadas atómicas (después de optimización):")
                for atom_num, atom_type, x, y, z in results['atomic_coordinates']:
                    print(f"  Átomo {atom_num} ({atom_type}): X={x:.4f}, Y={y:.4f}, Z={z:.4f}")
        else:
            print("\nNo se pudieron extraer datos porque el archivo de salida está vacío o no se generó.")
            print("Esto puede deberse a:")
            print("1. Problemas con el archivo .mop de entrada (formato incorrecto)")
            print("2. Problemas con la instalación o ejecución de MOPAC")
            print("3. Restricciones del modelo (máximo 50 átomos pesados, 50 hidrógenos)")
            print("4. Problemas de permisos o rutas de acceso")
    
    except Exception as e:
        print(f"Error procesando el archivo: {e}")
        import traceback
        print("Detalles del error:")
        traceback.print_exc()

if __name__ == "__main__":
    main()