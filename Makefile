PY      ?= .venv/bin/python
GOLDENS ?= ../lulcc-goldens/goldens

.PHONY: help env test check ensemble compare maps audit audit-power clean

help:
	@echo "make env     create .venv with the pinned packages"
	@echo "make test    unit tests (seconds)"
	@echo "make check   compare with the deterministic golden, stage by stage (~1 min)"
	@echo "             GOLDENS=path/to/lulcc-goldens/goldens"
	@echo "make ensemble  our 30 seeds (~20 min)"
	@echo "make compare   apply acceptance.toml to both ensembles, draw the figure"
	@echo "make maps      map figure of the article (observed and simulated 1999)"
	@echo "make audit     checks that the comparison is what it claims (~3 min)"
	@echo "make audit-power  audit + ensemble negative control (~15 min more)"

env:
	python3 -m venv .venv
	.venv/bin/pip install -q -r requirements.txt

test:
	$(PY) -m pytest -q tests

check:
	$(PY) scripts/check_deterministic.py $(GOLDENS) --out results/deterministic.json

ensemble:
	$(PY) scripts/run_ensemble.py $(GOLDENS)/data --seeds 30 --out results/ensemble_python.csv

compare:
	$(PY) scripts/compare_ensemble.py $(GOLDENS)/ensemble results/ensemble_python.csv --outdir results

maps:
	$(PY) scripts/make_maps.py $(GOLDENS) --out results/fig_maps.png

audit:
	$(PY) scripts/audit.py $(GOLDENS) --out results/audit.json

audit-power:
	$(PY) scripts/audit.py $(GOLDENS) --power --out results/audit.json

clean:
	rm -rf .pytest_cache
