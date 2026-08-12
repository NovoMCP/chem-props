"""
ChemProps Service - STATELESS VERSION
Chemical property calculation service without database dependencies

This service is completely stateless:
- Receives molecules via REST API
- Calculates chemical properties using RDKit
- Returns results to caller
- Does NOT store anything in database
"""

import sys
import os
from dotenv import load_dotenv
from fastapi import FastAPI, HTTPException, Depends, Header
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field
from typing import List, Dict, Any, Optional
from datetime import datetime
from contextlib import asynccontextmanager
import uvicorn
import logging
import io
import base64

# RDKit imports
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, Lipinski, QED, Draw, AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import rdMolDescriptors, Crippen
from qed_v2 import qed_v2
from rdkit.Contrib.SA_Score import sascorer

load_dotenv()

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    stream=sys.stderr,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Environment variables
API_KEY = os.getenv("API_KEY", None)  # Optional API key for authentication
SERVICE_PORT = int(os.getenv("PORT", "8003"))

# Application lifecycle handler
@asynccontextmanager
async def lifespan(app: FastAPI):
    logger.info(f"Starting ChemProps Service (STATELESS mode)...")
    logger.info(f"ChemProps ready for requests on port {SERVICE_PORT}")
    yield
    logger.info("ChemProps Service shutting down...")

# Initialize FastAPI
app = FastAPI(
    title="ChemProps Service (Stateless)",
    description="Stateless chemical property calculation service",
    version="2.0.0",
    lifespan=lifespan
)

# CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Models
class Molecule(BaseModel):
    id: str
    smiles: str
    properties: Optional[Dict[str, Any]] = None

class MoleculesList(BaseModel):
    molecules: List[Molecule]

class PropertyResponse(BaseModel):
    id: str
    smiles: str
    properties: Dict[str, Any]
    calculated_at: str

# Authentication
async def verify_api_key(x_api_key: Optional[str] = Header(None)):
    """Verify API key if one is configured."""
    if API_KEY is None:
        return True  # No API key configured, allow all
    
    if x_api_key is None:
        raise HTTPException(
            status_code=401,
            detail="API key required",
            headers={"WWW-Authenticate": "ApiKey"}
        )
    
    if x_api_key != API_KEY:
        raise HTTPException(status_code=401, detail="Invalid API key")
    
    return True

# Core calculation functions
def calculate_properties(smiles: str) -> Dict[str, Any]:
    """Calculate chemical properties for a molecule."""
    try:
        # Validate SMILES first - log errors clearly
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            error_msg = f"Invalid SMILES string: {smiles}"
            logger.warning(error_msg)
            return {"error": error_msg, "error_type": "validation"}

        # Add hydrogens for more accurate calculations
        mol_with_h = Chem.AddHs(mol)
        
        properties = {}
        
        # Basic properties
        properties["molecular_weight"] = round(Descriptors.MolWt(mol), 2)
        properties["exact_mass"] = round(Descriptors.ExactMolWt(mol), 4)
        properties["logp"] = round(Descriptors.MolLogP(mol), 2)
        properties["tpsa"] = round(Descriptors.TPSA(mol), 2)
        properties["qed"] = round(qed_v2(mol), 3)
        
        # Lipinski's Rule of Five
        properties["h_bond_donors"] = Lipinski.NumHDonors(mol)
        properties["h_bond_acceptors"] = Lipinski.NumHAcceptors(mol)
        properties["rotatable_bonds"] = Descriptors.NumRotatableBonds(mol)
        
        # Lipinski violations
        lipinski_violations = 0
        if properties["molecular_weight"] > 500:
            lipinski_violations += 1
        if properties["logp"] > 5:
            lipinski_violations += 1
        if properties["h_bond_donors"] > 5:
            lipinski_violations += 1
        if properties["h_bond_acceptors"] > 10:
            lipinski_violations += 1
        properties["lipinski_violations"] = lipinski_violations
        
        # Veber rules
        veber_violations = 0
        if properties["rotatable_bonds"] > 10:
            veber_violations += 1
        if properties["tpsa"] > 140:
            veber_violations += 1
        properties["veber_violations"] = veber_violations
        
        # Additional descriptors
        properties["num_atoms"] = mol.GetNumAtoms()
        properties["num_heavy_atoms"] = Lipinski.HeavyAtomCount(mol)
        properties["num_rings"] = Descriptors.RingCount(mol)
        properties["num_aromatic_rings"] = Descriptors.NumAromaticRings(mol)
        properties["fraction_sp3"] = round(Descriptors.FractionCSP3(mol), 3)
        properties["bertz_ct"] = round(Descriptors.BertzCT(mol), 2)  # Complexity
        
        # Calculate logD (pH 7.4) - approximation using logP
        # logD = logP for neutral molecules, adjusted for ionizable groups
        logD = properties["logp"]
        # Simple adjustment for ionizable groups
        num_basic_nitrogens = sum(1 for atom in mol.GetAtoms() 
                                  if atom.GetAtomicNum() == 7 and atom.GetTotalDegree() < 4)
        num_acidic_oxygens = sum(1 for atom in mol.GetAtoms() 
                                 if atom.GetAtomicNum() == 8 and atom.GetTotalDegree() == 1)
        
        # Rough approximation: basic groups increase logD at pH 7.4, acids decrease it
        if num_basic_nitrogens > 0:
            logD -= 0.5 * num_basic_nitrogens  # Basic groups are protonated at pH 7.4
        if num_acidic_oxygens > 0:
            logD -= 0.7 * num_acidic_oxygens  # Acidic groups are deprotonated at pH 7.4
        properties["logD"] = round(logD, 2)
        
        # Calculate molecular volume
        try:
            # Use LabuteASA (Accessible Surface Area) as a proxy for volume
            asa = Descriptors.LabuteASA(mol)
            # Convert ASA to approximate volume (empirical relationship)
            properties["volume"] = round(asa * 0.65, 2)
        except:
            # Fallback: estimate from molecular weight
            properties["volume"] = round(properties["molecular_weight"] * 0.7, 2)
        
        # Calculate molar refractivity
        properties["refractivity"] = round(Crippen.MolMR(mol), 2)
        
        # Drug-likeness scores
        properties["drug_likeness"] = round(1.0 - (lipinski_violations / 4.0), 2)
        properties["veber_score"] = round(1.0 - (veber_violations / 2.0), 2)
        
        # Synthetic accessibility
        try:
            properties["synthetic_accessibility"] = round(
                sascorer.calculateScore(mol), 2
            )
        except Exception as e:
            logger.warning(f"SA score calculation failed: {e}")
            properties["synthetic_accessibility"] = None
        
        # Solubility estimate
        solubility_score = 0.5
        if properties["logp"] < 3:
            solubility_score += 0.25
        if properties["tpsa"] > 75:
            solubility_score += 0.25
        properties["solubility_estimate"] = min(1.0, solubility_score)
        
        # Formal charge
        formal_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
        properties["formal_charge"] = formal_charge
        
        # Scaffolds
        try:
            scaffold = MurckoScaffold.GetScaffoldForMol(mol)
            scaffold_smiles = Chem.MolToSmiles(scaffold) if scaffold else ""
            properties["murcko_scaffold"] = scaffold_smiles
        except:
            properties["murcko_scaffold"] = ""
        
        return properties
        
    except Exception as e:
        error_msg = f"Error calculating properties for {smiles}: {str(e)}"
        logger.error(error_msg)
        return {"error": str(e), "error_type": "calculation"}

def generate_molecule_image(smiles: str, size=(300, 300)) -> str:
    """Generate a base64-encoded PNG image of the molecule."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        # Generate 2D coordinates
        AllChem.Compute2DCoords(mol)
        
        # Create image
        img = Draw.MolToImage(mol, size=size)
        
        # Convert to base64
        buffered = io.BytesIO()
        img.save(buffered, format="PNG")
        img_base64 = base64.b64encode(buffered.getvalue()).decode()
        
        return img_base64
        
    except Exception as e:
        logger.error(f"Error generating image: {str(e)}")
        return None

# API Endpoints

@app.get("/health")
async def health():
    """Health check endpoint."""
    return {
        "status": "healthy",
        "service": "chem-props",
        "version": "2.0.0",
        "mode": "stateless"
    }

@app.get("/chem-props/health")
async def service_health():
    """Service-specific health check."""
    return await health()

@app.post("/chem-props/calculate", response_model=List[PropertyResponse])
async def calculate_batch(
    request: MoleculesList,
    authenticated: bool = Depends(verify_api_key)
):
    """
    Calculate properties for multiple molecules - STATELESS operation.

    This endpoint:
    1. Receives molecules
    2. Calculates chemical properties
    3. Returns results immediately
    4. Does NOT store anything in database
    """
    results = []

    for molecule in request.molecules:
        properties = calculate_properties(molecule.smiles)

        # Log errors for monitoring/debugging
        if "error" in properties:
            if properties.get("error_type") == "validation":
                logger.warning(f"Molecule {molecule.id} failed validation: {properties['error']}")
            else:
                logger.error(f"Molecule {molecule.id} calculation error: {properties['error']}")

        results.append(PropertyResponse(
            id=molecule.id,
            smiles=molecule.smiles,
            properties=properties,
            calculated_at=datetime.utcnow().isoformat()
        ))

    return results

@app.post("/chem-props/calculate_single")
async def calculate_single(
    smiles: str,
    authenticated: bool = Depends(verify_api_key)
):
    """Calculate properties for a single SMILES string."""
    properties = calculate_properties(smiles)

    # Log errors for monitoring/debugging
    if "error" in properties:
        if properties.get("error_type") == "validation":
            logger.warning(f"SMILES '{smiles}' failed validation: {properties['error']}")
        else:
            logger.error(f"SMILES '{smiles}' calculation error: {properties['error']}")

    return {
        "smiles": smiles,
        "properties": properties,
        "calculated_at": datetime.utcnow().isoformat()
    }

@app.get("/chem-props/image/{smiles}")
async def get_molecule_image(
    smiles: str,
    size: int = 300,
    authenticated: bool = Depends(verify_api_key)
):
    """Generate molecular structure image."""
    img_base64 = generate_molecule_image(smiles, size=(size, size))
    
    if img_base64 is None:
        raise HTTPException(status_code=400, detail="Invalid SMILES or image generation failed")
    
    return {
        "smiles": smiles,
        "image": f"data:image/png;base64,{img_base64}",
        "size": size
    }

@app.get("/chem-props/druglikeness/{smiles}")
async def get_druglikeness(
    smiles: str,
    authenticated: bool = Depends(verify_api_key)
):
    """Calculate drug-likeness metrics for a molecule."""
    properties = calculate_properties(smiles)

    if "error" in properties:
        error_detail = f"{properties['error']} (type: {properties.get('error_type', 'unknown')})"
        logger.warning(f"Drug-likeness request failed for '{smiles}': {error_detail}")
        raise HTTPException(status_code=400, detail=error_detail)

    return {
        "smiles": smiles,
        "lipinski_violations": properties.get("lipinski_violations", 0),
        "veber_violations": properties.get("veber_violations", 0),
        "drug_likeness_score": properties.get("drug_likeness", 0),
        "veber_score": properties.get("veber_score", 0),
        "qed": properties.get("qed", 0),
        "synthetic_accessibility": properties.get("synthetic_accessibility", 5.0),
        "details": {
            "molecular_weight": properties.get("molecular_weight"),
            "logp": properties.get("logp"),
            "h_bond_donors": properties.get("h_bond_donors"),
            "h_bond_acceptors": properties.get("h_bond_acceptors"),
            "tpsa": properties.get("tpsa"),
            "rotatable_bonds": properties.get("rotatable_bonds")
        }
    }

@app.post("/chem-props/tanimoto")
async def pairwise_tanimoto(
    smiles_a: str,
    smiles_b: List[str],
    authenticated: bool = Depends(verify_api_key)
):
    """Compute pairwise Tanimoto similarity (Morgan radius-2, 2048-bit) between a seed and candidates."""
    seed_mol = Chem.MolFromSmiles(smiles_a)
    if seed_mol is None:
        raise HTTPException(status_code=400, detail=f"Invalid seed SMILES: {smiles_a}")

    seed_fp = AllChem.GetMorganFingerprintAsBitVect(seed_mol, 2, nBits=2048)

    results = []
    for smi in smiles_b:
        mol_b = Chem.MolFromSmiles(smi)
        if mol_b is None:
            results.append({"smiles": smi, "tanimoto": None, "patent_risk": "error", "note": "Invalid SMILES"})
            continue

        fp_b = AllChem.GetMorganFingerprintAsBitVect(mol_b, 2, nBits=2048)
        tc = round(DataStructs.TanimotoSimilarity(seed_fp, fp_b), 4)

        if tc > 0.7:
            risk, note = "high", "Tc > 0.7 — too similar to seed, likely in same patent family"
        elif tc > 0.4:
            risk, note = "low", "Tc 0.4-0.7 — true scaffold hop, likely patentable"
        else:
            risk, note = "novel", "Tc < 0.4 — highly novel, verify pharmacophore retention"

        results.append({"smiles": smi, "tanimoto": tc, "patent_risk": risk, "note": note})

    return {"seed_smiles": smiles_a, "comparisons": results, "count": len(results)}


@app.get("/chem-props/status")
async def service_status():
    """Get detailed service status."""
    return {
        "service": "chem-props",
        "version": "2.0.0",
        "mode": "stateless",
        "capabilities": {
            "property_calculation": True,
            "image_generation": True,
            "druglikeness_scoring": True,
            "scaffold_analysis": True,
            "lipinski_evaluation": True,
            "veber_evaluation": True
        },
        "supported_properties": [
            "molecular_weight", "logp", "tpsa", "qed",
            "h_bond_donors", "h_bond_acceptors",
            "rotatable_bonds", "num_rings",
            "drug_likeness", "synthetic_accessibility"
        ]
    }

# Main
if __name__ == "__main__":
    try:
        print(f"=== Starting ChemProps Service (STATELESS) ===")
        print(f"Mode: STATELESS - No database connections")
        print(f"Port: {SERVICE_PORT}")
        uvicorn.run(app, host="0.0.0.0", port=SERVICE_PORT, log_level="info")
    except Exception as e:
        print(f"FATAL ERROR: Failed to start service: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)