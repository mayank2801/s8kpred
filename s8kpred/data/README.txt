# S8kPred data directory
#
# This directory must contain the following files before running predictions:
# data.zip contains 
#   TriPeptidePropensityThreeStateSecStructure2AND.csv
#   TriPeptidePropensityEightStateSecStructure.csv
#   TripeptideBinaryTable_60.csv
#   model_3state.json     (3-state XGBoost model)
#   model_8state.ubj      (8-state XGBoost model)
#
# Download them with:
#   python scripts/download_models.py
#
# Or from the GitHub Releases page:
#   https://github.com/mayank2801/s8kpred/releases
