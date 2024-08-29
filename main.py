from fastapi import FastAPI, HTTPException
from rdkit import Chem
from iterator import IteratorMolecules
import logging

logging.basicConfig(level=logging.DEBUG, filename="py_log.log", filemode="w",
                    format="%(asctime)s %(levelname)s %(message)s")
logging.debug("A DEBUG Message")
logging.info("An INFO")
logging.warning("A WARNING")
logging.critical("A message of CRITICAL severity")

app = FastAPI()

molecules_db = [
    {"molecule_id": 1, "molecule": "CCO"},
    {"molecule_id": 2, "molecule": "Cc1ccccc1"},
    {"molecule_id": 3, "molecule": "CC(=O)O"},
    {"molecule_id": 4, "molecule": "CC(=O)Oc1ccccc1C(=O)O"}
]
num_all_molecules = len(molecules_db)
molecules_values = [molecules_db[i]["molecule"]
                    for i in range(len(molecules_db))]


@app.get("/molecules", tags=["MOLECULES"], summary="Retrieve all molecules")
def retrieve_molecule():
    logging.info(f"Molecules database successfully displayed.\n"
                 f"Response:{molecules_db}")
    return molecules_db


@app.get("/molecules/{molecule_id}", tags=["MOLECULES"],
         summary="Retrieve molecule by id")
def retrieve_molecule_by_id(molecule_id: int):
    for index, molecule in enumerate(molecules_db):
        if molecule["molecule_id"] == molecule_id:
            logging.info(f"Molecule with id={molecule_id} "
                         f"successfully released.\n"
                         f"Response:{molecule}")
            return molecule
    logging.info(f"Impossible to display molecule.\n"
                 f"Molecule with id={molecule_id} is not found",
                 exc_info=HTTPException(status_code=404))
    raise HTTPException(status_code=404, detail="Molecule is not found")


@app.put("/molecules/{molecule_id}", tags=["MOLECULES"],
         summary="Update molecule by id", status_code=201,
         response_description="Molecule is updated")
def update_molecule(molecule_id: int, updated_molecule: dict):
    """
       Update a molecule with all information:

       - **molecule_id**: each molecule must have a identifier
       - **molecule**: each molecule must have a name
    """
    for index, molecule in enumerate(molecules_db):
        if molecule["molecule_id"] == molecule_id:
            molecules_db[index] = updated_molecule
            logging.info(f"Molecule with id={molecule_id} "
                         f"successfully updated.\n"
                         f"Response:{molecules_db}")
            return molecules_db
    logging.info(f"Impossible to update molecule.\n"
                 f"Molecule with id={molecule_id} is not found",
                 exc_info=HTTPException(status_code=404))
    raise HTTPException(status_code=404, detail="Molecule is not found")


@app.delete("/molecules/{molecule_id}", tags=["MOLECULES"],
            summary="Delete molecule by id",
            response_description="Molecule is deleted")
def delete_molecule_by_id(molecule_id: int):
    for index, molecule in enumerate(molecules_db):
        if molecule["molecule_id"] == molecule_id:
            molecules_db.pop(index)
            logging.info(f"Molecule with id={molecule_id} "
                         f"successfully deleted.\n"
                         f"Response:{molecule}")
            return molecule
    logging.info(f"Impossible to delete molecule.\n"
                 f"Molecule with id={molecule_id} is not found",
                 exc_info=HTTPException(status_code=404))
    raise HTTPException(status_code=404, detail="Molecule is not found")


@app.post("/molecules/{molecule_id}", tags=["MOLECULES"],
          summary="Create new molecule", status_code=201,
          response_description="Molecule is created")
def add_molecule(molecule_id: int, molecule: str):
    """
    Create a molecule with all information:

    - **molecule_id**: each molecule must have an identifier
    - **molecule**: each molecule must have a name
    """
    molecule_to_add = {"molecule_id": molecule_id, "molecule": molecule}
    if molecule_to_add in molecules_db:
        logging.info(f"Impossible to add molecule.\n"
                     f"Molecule with id={molecule_id} "
                     f"is already in the molecule database",
                     exc_info=HTTPException(status_code=404))
        raise HTTPException(status_code=404,
                            detail="The molecule is already "
                                   "in the molecule database")
    else:
        molecules_db.append(molecule_to_add)
        logging.info(f"Molecule successfully added.\n"
                     f"Response:{molecule_to_add}")
        return molecule_to_add


@app.get("/molecules-search", tags=["SEARCH"], summary="Substructure_search")
def substructure_search(molecule_id: int, molecule: str,
                        number_of_molecules: int = 1):
    """
        Enter the data of the molecule that was added:

        - **molecule_id**: each molecule must have a identifier
        - **molecule**: each molecule must have a name
        - **number_of_molecules**: number of molecules in answer
    """
    flag = False
    if len(molecules_db) >= num_all_molecules:
        if {"molecule_id": molecule_id, "molecule": molecule} in molecules_db:
            new_molecule = Chem.MolFromSmiles(molecule)
            flag = True
    if flag:
        molecules_container = [Chem.MolFromSmiles(smiles)
                               for smiles in molecules_values]
        smiles_db = IteratorMolecules(molecules_container, number_of_molecules)
        substructure_result = [molecules_db[index] for index, molec
                               in enumerate(smiles_db)
                               if molec.HasSubstructMatch(new_molecule)]

        logging.info(f"The search for the  substructure of the molecules"
                     f" is successfully completed.\n"
                     f"Molecule smiles:{molecule}.\n"
                     f"Response:{substructure_result}")
        return substructure_result
    else:
        logging.info(f"Impossible to search for "
                     f"the  substructure of the molecules.\n"
                     f"Molecule with id={molecule_id} "
                     f"and smiles: {molecule} is not added",
                     exc_info=HTTPException(status_code=404))
        raise HTTPException(status_code=404, detail="Molecule is not added")
