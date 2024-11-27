# Enrichproqt 

Enrichproqt is a Python package for pathway enrichment analysis.

**Designed by** `qinti`

## Installation 

You can install it using pip:

```shell
pip install Enrichproqt
```

## Example Usage

```python
from Enrichproqt import PathwayEnrichment

protein_list = ["P12345", "Q67890"]
enrichment_results = PathwayEnrichment.analyze("ALL", protein_list)
print(enrichment_results)
```

