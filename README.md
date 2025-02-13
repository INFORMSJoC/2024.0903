[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# An Exact Solution Approach for Hierarchical Clustering

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data
that were used in the research reported on in the paper 
[An Exact Solution Approach for Hierarchical Clustering](https://pubsonline.informs.org/doi/10.1287/ijoc.2024.0903) by Rick Willemsen, Carlo Cavicchia, Wilco van den Heuvel, and Michel van de Velden.

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2024.0903

https://doi.org/10.1287/ijoc.2024.0903.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{Willemsen2024,
  author =        {Willemsen, Rick and Cavicchia, Carlo and Van den Heuvel, Wilco and Van de Velden, Michel},
  publisher =     {INFORMS Journal on Computing},
  title =         {An Exact Solution Approach for Hierarchical Clustering},
  year =          {2024},
  doi =           {10.1287/ijoc.2024.0903.cd},
  url =           {https://github.com/INFORMSJoC/2024.0903},
  note =          {Available for download at https://github.com/INFORMSJoC/2024.0903},
}  
```

## Description

This repository contains datasets for hierarchical clustering problems and several solution methods.

## Software

The `CPLEX` library needs to be installed from [IBM CPLEX](https://www.ibm.com/docs/nl/icos/22.1.2?topic=cplex-installing).

## Datasets

The `data` folder contains the datasets.

The first line of a file may contain two input arguments: `n` `p`. This means that the file contains `n` observations in `p`-dimensional space. Each row corresponds to a coordinate.

The first line of a file may contain three input arguments: `n` `p` `k`. This means that the file contains `n` observations in `p`-dimensional space, which are partitioned into `k` clusters. Each row corresponds to a coordinate and a cluster assignment.

## Replicating

The `src` folder contains all the source code. To replicate the results of the paper, please ensure to import the correct data and run the appropriate `BatchRun` script (e.g. `BatchRunPDCGM.java`). The algorithm to be run can be changed within the script, the parameters can be adjusted within the `cl.data.GlobalParam.java` file.