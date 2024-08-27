# PolyUniverse
Generation of a Large-scale Polymer Library Using Rule-Based Polymerization Reactions for Polymer Informatics.

All generation results and small compounds datasets can be found at https://zenodo.org/records/12585902

The generation codes for vitrimer, epoxy,polybenzimidazole (PBI), and copolymer poly(imide-imine) (PI-PIM), along with a small-scale example dataset (due to GitHub's file size limitations), can be found in the Generation folder. For generating large-scale hypothetical datasets, please use the extensive small molecule dataset provided at https://zenodo.org/records/12585902. For the generation of Polyimide, Polyolefin, Polyester, Polyamide, and Polyurethane, please refer to: https://github.com/PEJpOhno/SMiPoly

The pre-trained models needed for prediction can be found at: https://zenodo.org/records/12587825 or https://huggingface.co/ytl0410/PolyUniverse/tree/main. To predict polymer properties, place the dataframe containing the polymer's SMILES information into a .csv file, then use a command such as:

> python Thermal_Property.py PBI.csv
