# ChemProps Service (Stateless)

## Overview
The ChemProps Service is a **stateless microservice** that calculates chemical properties for molecular structures using RDKit. This service receives molecules via REST API, computes their properties, and returns results without any database dependencies.

## Architecture
```
Frontend/Client → Quanta MCP → ChemProps Service → Results back to Quanta MCP → DB Manager → Database
```

**Key Design Principles:**
- **Stateless**: No database connections or storage
- **Single Responsibility**: Only calculates molecular properties
- **Scalable**: Can run multiple instances without coordination
- **Fast**: Sub-second response times for property calculations

## Core Functionality
- **Property Calculation**: 20+ molecular descriptors (MW, LogP, TPSA, QED, etc.)
- **Drug-likeness Scoring**: Lipinski's Rule of Five and Veber rules
- **Scaffold Analysis**: Murcko scaffold extraction
- **Molecular Imaging**: Generate 2D structure images
- **Synthetic Accessibility**: Estimate synthesis difficulty

## Technology Stack
- **Language**: Python 3.9
- **Framework**: FastAPI
- **Chemistry Library**: RDKit 2023.9.5
- **Deployment**: AWS ECS Fargate
- **Container**: Docker

## API Endpoints

### Calculate Properties (Batch)
```
POST /chem-props/calculate
Headers:
  Content-Type: application/json
  X-API-Key: your-api-key-here  # Required

Body:
{
  "molecules": [
    {
      "id": "mol-123",
      "smiles": "CC(=O)Oc1ccccc1C(=O)O"
    }
  ]
}

Response:
[
  {
    "id": "mol-123",
    "smiles": "CC(=O)Oc1ccccc1C(=O)O",
    "properties": {
      "molecular_weight": 180.16,
      "logp": 1.31,
      "tpsa": 63.6,
      "qed": 0.55,
      "h_bond_donors": 1,
      "h_bond_acceptors": 3,
      "lipinski_violations": 0,
      "drug_likeness": 1.0,
      "synthetic_accessibility": 5.0,
      // ... 15+ more properties
    },
    "calculated_at": "2024-01-15T10:30:00Z"
  }
]
```

### Calculate Single Molecule
```
POST /chem-props/calculate_single?smiles=CCO
Headers:
  X-API-Key: your-api-key-here
```

### Get Drug-likeness Score
```
GET /chem-props/druglikeness/{smiles}
Headers:
  X-API-Key: your-api-key-here

Response:
{
  "smiles": "...",
  "lipinski_violations": 0,
  "veber_violations": 0,
  "drug_likeness_score": 1.0,
  "qed": 0.822,
  "details": {...}
}
```

### Generate Molecular Image
```
GET /chem-props/image/{smiles}?size=300
Headers:
  X-API-Key: your-api-key-here

Response:
{
  "smiles": "...",
  "image": "data:image/png;base64,..."
}
```

### Health Check
```
GET /health
GET /chem-props/health
# No authentication required
```

### Service Status
```
GET /chem-props/status
Headers:
  X-API-Key: your-api-key-here
```

## Properties Calculated

| Property | Description |
|----------|-------------|
| molecular_weight | Molecular weight in g/mol |
| exact_mass | Exact molecular mass |
| logp | Octanol-water partition coefficient |
| tpsa | Topological polar surface area |
| qed | Quantitative estimate of drug-likeness |
| h_bond_donors | Number of hydrogen bond donors |
| h_bond_acceptors | Number of hydrogen bond acceptors |
| rotatable_bonds | Number of rotatable bonds |
| lipinski_violations | Count of Lipinski's Rule of Five violations |
| veber_violations | Count of Veber rule violations |
| num_rings | Total number of rings |
| num_aromatic_rings | Number of aromatic rings |
| fraction_sp3 | Fraction of sp3 carbons |
| bertz_ct | Molecular complexity |
| drug_likeness | Drug-likeness score (0-1) |
| synthetic_accessibility | Synthetic accessibility score (1-10, lower is easier) |
| murcko_scaffold | SMILES of Murcko scaffold |
| solubility_estimate | Estimated aqueous solubility |

## Environment Variables
```bash
PORT=8005
API_KEY=your-secure-api-key  # Required for authentication
```

## Security

### API Authentication
The service requires API key authentication for all endpoints except health checks:
- Header: `X-API-Key: your-api-key`
- Returns 401 Unauthorized without valid key

### Example with Authentication
```bash
curl -X POST http://your-alb-url/chem-props/calculate \
  -H "Content-Type: application/json" \
  -H "X-API-Key: your-api-key" \
  -d '{
    "molecules": [
      {"id": "aspirin", "smiles": "CC(=O)Oc1ccccc1C(=O)O"}
    ]
  }'
```

## Deployment

The service is a single stateless container. Run the published image, or build it yourself.

```bash
# Pull and run the published image
docker run -p 8003:8003 ghcr.io/novomcp/chem-props:latest

# Or build from source
docker build -t chem-props .
docker run -p 8003:8003 chem-props

# Test
curl http://localhost:8003/health
```

Set `API_KEY` to require an `X-API-Key` header; leave it unset to run open (useful for local development).

Point the NovoMCP engine at this service by setting `CHEM_PROPS_URL` to its URL.

## Performance
- **Response Time**: < 100ms per molecule for property calculation
- **Batch Size**: Up to 100 molecules per request
- **Concurrent Requests**: Handles multiple requests in parallel
- **Memory Usage**: ~500MB per container

## Error Handling
- Invalid SMILES returns error message with 400 status
- Missing API key returns 401 Unauthorized
- Server errors return 500 with error details

## Monitoring
- **Health Endpoint**: `/health` for load-balancer health checks
- **Logs**: structured logs to stdout (collect with your platform's log driver)

## Integration Example (Python)

```python
import requests

API_KEY = "your-api-key"
URL = "http://chem-props-service.com/chem-props/calculate"

molecules = [
    {"id": "mol1", "smiles": "CC(=O)Oc1ccccc1C(=O)O"},
    {"id": "mol2", "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"}
]

response = requests.post(
    URL,
    headers={
        "Content-Type": "application/json",
        "X-API-Key": API_KEY
    },
    json={"molecules": molecules}
)

results = response.json()
for mol in results:
    print(f"{mol['id']}: MW={mol['properties']['molecular_weight']}")
```

## Troubleshooting

### Common Issues
1. **401 Unauthorized**: Check API key in X-API-Key header
2. **Invalid SMILES**: Verify SMILES string syntax
3. **Timeout**: Reduce batch size for complex molecules
4. **500 Error**: Check CloudWatch logs for details

## Support
For issues or questions, create an issue in the GitHub repository.

## License

Code is licensed under the Apache License 2.0 (see `LICENSE`).

The QED v2 desirability coefficients (`qed_v2_coefficients.json`) are refit on ChEMBL 35 approved oral drugs and are licensed CC-BY-SA-3.0, with attribution to ChEMBL (EMBL-EBI). See `NOTICE`.