# Beyond Affinity: A Benchmark of 1D, 2D, and 3D Methods Reveals Critical Trade-offs in Structure-Based Drug Design

This repository hosts an open-source benchmark for Structure-based Drug Design, to facilitate the transparent and reproducible evaluation of algorithmic advances in molecular optimization. This repository supports 15 Structure-based Drug Design algorithms on 4 tasks categories.


## 15 Methods


|                  `Model`                                                                                        | `Dimension` | `Tested molcules (Avg) `   | `requires_gpu` |
|-----------------------------------------------------------------------------------------------------------|------------|-----------------------------------------------|---------|
| **[3DSBDD](https://arxiv.org/abs/2203.10446)**                              | 3D      |  771          |    yes     |
| **[Pocket2mol](https://pubs.acs.org/doi/10.1021/acs.jcim.8b00839)**                                       | 3D     | 928          |    yes     |
| **[PocketFlow](https://arxiv.org/abs/2205.07249)**               | 3D    | 1000          |    yes    |
| **[RenGen](https://www.nature.com/articles/s42256-023-00712-7)**                                             | 3D    |  631          |    yes     |
| **[DST](https://arxiv.org/abs/2109.10469)**                      | 2D       |  1001          |    no     |
| **[Graph GA](https://pubs.rsc.org/en/content/articlelanding/2019/sc/c8sc05372c)**                                 | 2D     | 643          |    no    |
| **[MIMOSA](https://arxiv.org/abs/2010.02318)**                                                   | 2D     |  1001          |    yes     |
| **[MolDQN](https://arxiv.org/abs/1810.08678)**                                | 2D    |  501          |    yes    |
| **[Pasithea](https://arxiv.org/abs/2012.09712)**                                                   | 1D     |  914          |    yes     |
| **[REINVENT](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-017-0235-x)**                                                  | 1D    |  1000          |    yes     |
| **[SELFIES-VAE-BO](https://arxiv.org/abs/1610.02415)**                                                   | 1D   | 200          |    yes     |
| **[SMILES-GA](https://arxiv.org/abs/1804.02134)**                                                          | 1D     |  584          |    no    |
| **[SMILES-LSTM-HC](https://arxiv.org/abs/1811.09621)**                                             | 1D     | 501          |    no    |
| **[SMILES-VAE-BO](https://arxiv.org/abs/1610.02415)**                                                 | 1D    | 200          |    yes     | 
| **[TargetDiff](https://arxiv.org/abs/2303.03543)**                                                    | 3D     | 75          |    yes     |



## Receptor information

All the receptors used in our benchmark can be found in receptor/ folder

## Sampling and evaluating

For 3DSBDD and Pocket2mol, we use this command to generate:

```bash
python sample_for_pdb.py --pdb_path [your pdb] --center=[centers] --bbox_size [box size] --outdir [your outdir]
```

Also need to change the num_samples in the sample_for_pdb.yml

For PocketFlow, we use this command to generate:

```bash
python main_generate.py -pkt [your pdb] --ckpt ckpt/ZINC-pretrained-255000.pt -n 1000 -d cuda:0 --root_path [your outdir] --name [pdb name] -at 1.0 -bt 1.0 --max_atom_num 35 -ft 0.5 -cm True --with_print True
```

For ResGen, we first convert our pdb file to sdf file and use this command to generate:

```bash
python gen.py --pdb_file [your pdb] --sdf_file [correspond sdf] --outdir [your outdir]
```

For the rest of models that are under **[PMO](https://github.com/wenhao-gao/mol_opt)**, we use the following command to generate

```bash
python -u run.py --method [method to test] --task production --n_runs 3 --wandb offline --max_oracle_calls 5 --custom_docking [your pdb] --custom_docking_config [your pdb config]
```

After getting molecules, we could start evaluate. First we calculate the docking score using vina_test_model_pmo.py or vina_test_model_pmo.py. Then we could run Heuristic Oracles and Molecule Generation Oracles using heuristic_merics.py and generation_metrics.py. At last we could run pose test with posebuster_test.py and clash_test.py

## Environment

To run vina, please use vina_environment

To run Heuristic Oracles and Molecule Generation Oracles, please use pyscreener_environment

To run pose test, please use posebusters_environment