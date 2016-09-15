# Instrution to use scripts

assay:
	python assay.py --download
	python assay.py --build
	python assay.py --pivot

subset: assay
	python assay.py --subset
