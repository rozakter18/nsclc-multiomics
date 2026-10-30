PY := python

.PHONY: all inventory download verify preprocess integrate model knowledge validate clean-inventory

all:
	$(PY) scripts/01_sample_inventory.py
	$(PY) scripts/02_download.py
	$(PY) scripts/03_preprocess.py
	$(PY) scripts/04_integrate_mofa.py
	$(PY) scripts/05_stage_model.py
	$(PY) scripts/06_knowledge_mining.py
	$(PY) scripts/07_validate_geo.py --gse GSE37745

inventory: results/inventory/inventory_report.md

results/inventory/inventory_report.md: scripts/01_sample_inventory.py
	$(PY) scripts/01_sample_inventory.py

download:
	$(PY) scripts/02_download.py

verify:
	$(PY) scripts/02_download.py --verify

clean-inventory:
	$(PY) -c "from pathlib import Path; import shutil; path = Path('results/inventory'); shutil.rmtree(path) if path.exists() else None"

preprocess:
	$(PY) scripts/03_preprocess.py

integrate:
	$(PY) scripts/04_integrate_mofa.py

model:
	$(PY) scripts/05_stage_model.py

knowledge:
	$(PY) scripts/06_knowledge_mining.py

validate:
	$(PY) scripts/07_validate_geo.py --gse GSE37745