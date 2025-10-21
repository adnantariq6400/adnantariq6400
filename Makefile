PY=python
PIP=pip
OUT?=outputs/demo

.PHONY: install run-demo

install:
	$(PIP) install -r requirements.txt

run-demo:
	$(PY) -m biomarker_pipeline.cli demo --output-dir $(OUT)
