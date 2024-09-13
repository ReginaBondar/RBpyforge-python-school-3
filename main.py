from fastapi import FastAPI, HTTPException
from rdkit import Chem
from iterator import SearchMolecules
import logging
import redis
import json

logging.basicConfig(level=logging.DEBUG, filename="py_log.log", filemode="w",
                    format="%(asctime)s %(levelname)s %(message)s")
logging.debug("A DEBUG Message")
logging.info("An INFO")
logging.warning("A WARNING")
logging.critical("A message of CRITICAL severity")

app = FastAPI()
redis_client = redis.Redis(host='redis', port=6379, db=0)

molecules_db = [
    {"molecule_id": 1, "molecule": "CCO"},
    {"molecule_id": 2, "molecule": "Cc1ccccc1"},
    {"molecule_id": 3, "molecule": "CC(=O)O"},
    {"molecule_id": 4, "molecule": "CC(=O)Oc1ccccc1C(=O)O"}
]


def get_cached_result(key: str):
    result = redis_client.get(key)
    if result:
        return json.loads(result)
    return None


def set_cache(key: str, value: dict, expiration: int = 60):
    redis_client.setex(key, expiration, json.dumps(value))


def set_cache_search(key: str, value: list, expiration: int = 30):
    redis_client.setex(key, expiration, json.dumps(value))


@app.get("/molecules", tags=["MOLECULES"], summary="Retrieve all molecules")
def retrieve_molecule():
    logging.info(f"Molecules database successfully displayed.\n"
                 f"Response:{molecules_db}")
    return molecules_db


@app.get("/molecules/{molecule_id}", tags=["MOLECULES"],
         summary="Retrieve molecule by id")
def retrieve_molecule_by_id(molecule_id: int):
    cache_key = f"search:{molecule_id}"
    cached_result = get_cached_result(cache_key)

    if cached_result:
        return f'source: cache, data: {cached_result}'

    for index, molecule in enumerate(molecules_db):
        if molecule["molecule_id"] == molecule_id:
            logging.info(f"Molecule with id={molecule_id} "
                         f"successfully released.\n"
                         f"Response:{molecule}")
            set_cache(cache_key, molecule)
            return f'source: database, data: {molecule}'

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
                        limit: int = 1):
    """
        Enter the data of the molecule that was added:

        - **molecule_id**: each molecule must have a identifier
        - **molecule**: each molecule must have a name
        - **limit**: number of molecules in answer
    """
    cache_key = f"search:{molecule_id}"
    cached_result = get_cached_result(cache_key)

    if cached_result:
        return f'source: cache, data: {cached_result}'

    if {"molecule_id": molecule_id, "molecule": molecule} in molecules_db:
        molecule_add = Chem.MolFromSmiles(molecule)
        molecules_values = [molecules_db[i]["molecule"]
                            for i in range(len(molecules_db))]
        molecules_container = [Chem.MolFromSmiles(smiles)
                               for smiles in molecules_values]
        index_db = SearchMolecules(molecules_container, molecule_add, limit)
        substructure_result = [molecules_db[i]
                               for i in index_db if i is not None]

        logging.info(f"The search for the  substructure of the molecules"
                     f" is successfully completed.\n"
                     f"Molecule smiles:{molecule}.\n"
                     f"Response:{substructure_result}")
        set_cache_search(cache_key, substructure_result)
        return f'source: database, data: {substructure_result}'
    else:
        logging.info(f"Impossible to search for "
                     f"the  substructure of the molecules.\n"
                     f"Molecule with id={molecule_id} "
                     f"and smiles: {molecule} is not added",
                     exc_info=HTTPException(status_code=404))
        raise HTTPException(status_code=404, detail="Molecule is not added")
