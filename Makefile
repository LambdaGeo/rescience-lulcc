PY      ?= .venv/bin/python
GOLDENS ?= ../lulcc-goldens/goldens

.PHONY: help env test check clean

help:
	@echo "make env     create .venv with the pinned packages"
	@echo "make test    unit tests (seconds)"
	@echo "make check   compare with the deterministic golden, stage by stage (~1 min)"
	@echo "             GOLDENS=path/to/lulcc-goldens/goldens"

env:
	python3 -m venv .venv
	.venv/bin/pip install -q -r requirements.txt

test:
	$(PY) -m pytest -q tests

check:
	$(PY) scripts/check_deterministic.py $(GOLDENS) --out results/deterministic.json

clean:
	rm -rf .pytest_cache
