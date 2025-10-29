import tkinter as tk
from tkinter import filedialog, messagebox
from rdkit import Chem
from rdkit.Chem import AllChem

class MolecularConverterApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Molecular Data Converter")

        self.create_widgets()

    def create_widgets(self):
        self.label = tk.Label(self.root, text="Select a file to convert:")
        self.label.pack(pady=10)

        self.select_button = tk.Button(self.root, text="Select File", command=self.select_file)
        self.select_button.pack(pady=5)

        self.convert_button = tk.Button(self.root, text="Convert", command=self.convert_file, state=tk.DISABLED)
        self.convert_button.pack(pady=5)

        self.status_label = tk.Label(self.root, text="")
        self.status_label.pack(pady=10)

    def select_file(self):
        self.file_path = filedialog.askopenfilename(filetypes=[("All Files", "*.*")])
        if self.file_path:
            self.status_label.config(text=f"Selected file: {self.file_path}")
            self.convert_button.config(state=tk.NORMAL)
        else:
            self.status_label.config(text="No file selected")
            self.convert_button.config(state=tk.DISABLED)

    def convert_file(self):
        try:
            with open(self.file_path, 'r') as file:
                data = file.read()

            if self.file_path.endswith('.smi'):
                mol = Chem.MolFromSmiles(data)
                if mol:
                    AllChem.Compute2DCoords(mol)
                    Chem.MolToMolFile(mol, 'output.mol2')
                    self.status_label.config(text="Conversion successful: output.mol2")
                else:
                    self.status_label.config(text="Invalid SMILES string")
            elif self.file_path.endswith('.mol2'):
                mol = Chem.MolFromMol2File(self.file_path)
                if mol:
                    smiles = Chem.MolToSmiles(mol)
                    with open('output.smi', 'w') as out_file:
                        out_file.write(smiles)
                    self.status_label.config(text="Conversion successful: output.smi")
                else:
                    self.status_label.config(text="Invalid MOL2 file")
            else:
                self.status_label.config(text="Unsupported file format")
        except Exception as e:
            messagebox.showerror("Error", str(e))

if __name__ == "__main__":
    root = tk.Tk()
    app = MolecularConverterApp(root)
    root.mainloop()
