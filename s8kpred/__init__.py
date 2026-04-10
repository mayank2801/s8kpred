"""
S8kPred: Protein Secondary Structure Prediction Tool
=====================================================
Predicts 3-state (H/E/C) and 8-state (H/G/I/E/B/T/S/L) secondary structure
from protein sequences using XGBoost models trained on PSSM and tripeptide
propensity features.

Usage
-----
    from s8kpred import predict
    results = predict("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL")
"""

__version__ = "0.1.0"
__author__ = "Mayank"
__license__ = "MIT"

from s8kpred.api import predict, predict_file, predict_fasta_string

__all__ = ["predict", "predict_file", "predict_fasta_string", "__version__"]
