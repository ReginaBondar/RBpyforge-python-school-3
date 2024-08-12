from fastapi import FastAPI, status, HTTPException
from rdkit import Chem


app = FastAPI()


molecules_db = [
    {"molecule_id": 1, "molecule": "CCO"},
    {"molecule_id": 2, "molecule": "Cc1ccccc1"},
    {"molecule_id": 3, "molecule": "CC(=O)O"},
    {"molecule_id": 4, "molecule": "CC(=O)Oc1ccccc1C(=O)O"}
]
appended_molecules = []
molecules_values = [molecules_db[i]["molecule"] for i in range(len(molecules_db))]


@app.get("/molecules", tags=["MOLECULES"], summary="Retrieve all molecules")
def retrieve_molecule():
    return molecules_db


@app.get("/molecules/{molecule_id}", tags=["MOLECULES"], summary="Retrieve molecule by id")
def retrieve_molecule_by_id(molecule_id: int):
    for index, molecule in enumerate(molecules_db):
        if molecule["molecule_id"] == molecule_id:
            return molecule
    raise HTTPException(status_code=404, detail="Molecule is not found")


@app.put("/molecules/{molecule_id}", tags=["MOLECULES"], summary="Update molecule by id", status_code=201,
         response_description="Molecule is updated")
def update_molecule(molecule_id: int, update_molecule: dict):
    '''
       Update a molecule with all information:

       - **molecule_id**: each molecule must have a identifier
       - **molecule**: each molecule must have a name
    '''
    for index, molecule in enumerate(molecules_db):
        if molecule["molecule_id"] == molecule_id:
            molecules_db[index] = update_molecule
            return molecules_db
    raise HTTPException(status_code=404, detail="Molecule is not found")


@app.delete("/molecules/{molecule_id}", tags=["MOLECULES"], summary="Delete molecule by id",
            response_description="Molecule is deleted")
def delete_molecule_by_id(molecule_id: int):
    for index, molecule in enumerate(molecules_db):
        if molecule["molecule_id"] == molecule_id:
            deleted_molecule = molecules_db.pop(index)
            return molecule
    raise HTTPException(status_code=404, detail="Molecule is not found")


@app.post("/molecules/{molecule_id}", tags=["MOLECULES"], summary="Create new molecule", status_code=201,
          response_description="Molecule is created")
def add_molecule(molecule_id: int, molecule: str):
    '''
    Create a molecule with all information:

    - **molecule_id**: each molecule must have a identifier
    - **molecule**: each molecule must have a name
    '''
    global appended_molecules
    if {"molecule_id": molecule_id, "molecule": molecule} in molecules_db:
        raise HTTPException(status_code=404, detail="The molecule is already in the molecule database")
    else:
        appended_molecules.append({"molecule_id": molecule_id, "molecule": molecule})
        molecules_db.append({"molecule_id": molecule_id, "molecule": molecule})
        return appended_molecules


@app.get("/molecules-search", tags=["SEARCH"], summary="Substructure_search")
def substructure_search(molecule_id: int, molecule: str):
    '''
        Enter the data of the molecule that was added:

        - **molecule_id**: each molecule must have a identifier
        - **molecule**: each molecule must have a name
    '''
    flag = False
    if len(appended_molecules) != 0:
        if {"molecule_id": molecule_id, "molecule": molecule} in appended_molecules:
            new_molecule = Chem.MolFromSmiles(molecule)
            flag = True
    if flag:
        molecules_container = [Chem.MolFromSmiles(smiles) for smiles in molecules_values]
        substructure_rezult = [molecules_db[index] for index, molec in enumerate(molecules_container)
                               if molec.HasSubstructMatch(new_molecule)]
        return substructure_rezult
    else:
        raise HTTPException(status_code=404, detail="Molecule is not added")
