import pickle
import pandas as pd
import scipy.stats as stats
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import pkg_resources

class PathwayEnrichment:
    def __init__(self, pathway_type="ALL"):
        self.valid_types = {"MF", "CC", "BP", "ALL"}
        if pathway_type not in self.valid_types:
            raise ValueError(f"Invalid pathway type. Choose from {self.valid_types}.")
        self.pathway_type = pathway_type
        self.pathway_dict = self._load_pathway_dict()

    def _load_pathway_dict(self):
        pathway_dict = {}
        try:
            if self.pathway_type == "ALL":
                for pt in ["MF", "CC", "BP"]:
                    file_path = pkg_resources.resource_filename(
                        __name__, f"data/{pt}_pathway_dict.pkl")
                    with open(file_path, "rb") as f:
                        pathway_dict.update(pickle.load(f))
            else:
                file_path = pkg_resources.resource_filename(
                    __name__, f"data/{self.pathway_type}_pathway_dict.pkl")
                with open(file_path, "rb") as f:
                    pathway_dict = pickle.load(f)
        except FileNotFoundError as e:
            print(f"file not found: {e}")
            raise
        except Exception as e:
            print(f"loading error: {e}")
            raise
        return pathway_dict

    def _calculate_enrichment(self, pathway, proteins, protein_set, total_proteins, N, n):
        """
        Helper function to calculate enrichment for a single pathway.
        """
        K = len(proteins)
        k = len(protein_set & proteins)

        if k > 0:
            p_value = 1 - stats.hypergeom.sf(k - 1, N, K, n)
            return {
                "Pathway": pathway,
                "P_value": p_value,
                "Subset_Proteins": k,
                "Pathway_Proteins": K
            }
        return None

    def enrich(self, protein_list, num_threads=4):
        """
        Calculate enrichment results based on protein list for the selected pathway type.
        This function uses multi-threading to parallelize the enrichment calculation.
        """
        total_proteins = set(protein for proteins in self.pathway_dict.values() for protein in proteins)
        N = len(total_proteins) 
        n = len(protein_list)  
        protein_set = frozenset(protein_list)

        results = []
        with ThreadPoolExecutor(max_workers=num_threads) as executor:

            future_to_pathway = {
                executor.submit(self._calculate_enrichment, pathway, proteins, protein_set, total_proteins, N, n):
                pathway for pathway, proteins in self.pathway_dict.items()
            }

            for future in as_completed(future_to_pathway):
                result = future.result()
                if result:
                    results.append(result)

        return pd.DataFrame(results).sort_values(by="P_value").reset_index(drop=True)

    @classmethod
    def analyze(cls, pathway_type, protein_list, num_threads=4):
        """
        Class method to analyze pathway enrichment for a given pathway type using multi-threading.
        :param pathway_type: One of 'MF', 'CC', 'BP', or 'ALL'.
        :param protein_list: List of proteins for analysis.
        :param num_threads: Number of threads to use for parallel computation.
        :return: DataFrame with enrichment results.
        """
        instance = cls(pathway_type)
        return instance.enrich(protein_list, num_threads)
