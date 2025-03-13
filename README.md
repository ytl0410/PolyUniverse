# PolyUniverse
Generation of a Large-scale Polymer Library Using Rule-Based Polymerization Reactions for Polymer Informatics.

All generation results and small compounds datasets can be found at https://zenodo.org/records/12585902

The generation codes for vitrimer, epoxy,polybenzimidazole (PBI), and copolymer poly(imide-imine) (PI-PIM), along with a small-scale example dataset (due to GitHub's file size limitations), can be found in the Generation folder. For generating large-scale hypothetical datasets, please use the extensive small molecule dataset provided at https://zenodo.org/records/12585902. For the generation of Polyimide, Polyolefin, Polyester, Polyamide, and Polyurethane, please refer to: https://github.com/PEJpOhno/SMiPoly

Before running, please set up your development environment: ensure that all necessary libraries and dependencies are installed. You can then run:
```
conda create -n py39 python=3.9
conda install -c conda-forge cudatoolkit=11.2 cudnn=8.1.0
python -m pip install "tensorflow<2.11"
pip install numpy==1.26.4
pip install matplotlib==3.8.0
pip install pandas==1.5.3
pip install rdkit==2023.9.1
pip install scikit-learn==1.3.0
pip install keras-tuner==1.4.5
pip install seaborn==0.13.0
```

The pre-trained models needed for prediction can be found at: https://zenodo.org/records/12587825 or https://huggingface.co/ytl0410/PolyUniverse/tree/main. To predict polymer properties, place the dataframe containing the polymer's SMILES information into a .csv file, then use a command such as:

```
python Thermal_Property.py PBI.csv
```
