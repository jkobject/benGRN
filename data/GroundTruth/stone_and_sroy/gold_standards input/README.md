# Experimental "gold standard" networks

This file describes the experimentally derived "gold standard" networks we used
as our ground truth to measure the accuracy of the inferred networks.

## Yeast

| Network | Assay | Description |
| --- | --- | --- |
| hu | Perturbation | Interactions from Hu |
| mac2 | ChIP | Interactions from MacIsaac |
| yeastract_c3 | ChIP,Perturbation | Interactions from the YEASTRACT database that appear in at least three studies |
| yeastract_t2 | ChIP,Perturbation | Interactions from the YEASTRACT database that appear in at least two types of study |

## Mouse embryonic stem cell (mESC)

| Network | Assay | Description |
| --- | --- | --- |
| mESC_chip | ChIP | Interactions from ENCODE database |
| mESC_escapechip | ChIP | Interactions from ESCAPE database |
| mESC_chip_union | ChIP | Union of the Encode and Escape networks |
| mESC_logof | Perturbation | Interactions from LOGOF |
| mESC_nishiyama | Perturbation | Interactions from Nishiyama |
| mESC_KDUnion | Perturbation | Union of the LOGOF and Nishiyama networks |
| mESC_chip_KDUnion_intersect | ChIP,Perturbation | Intersection of the ChIP and knockdown union networks |
| mESC_chip_KDUnion_c2 | ChIP,Perturbation | Edges from at least two of the above networks (2x ChIP, 2x perturb, or 1 of each) |
| mESC_ATAC_sridharan | ATAC | Network generated from Rupa's scATAC-seq experiment |
| mESC_lit | Multiple | Network of interactions reported in the literature |

## Mouse dendritic cell (mDC)

| Network | Assay | Description |
| --- | --- | --- |
| mDC_chip | ChIP | Interactions from Garber |
| mDC_KDUnion | Perturbation | Union of interactions from Parnas, Dixit, and Chevrier |
| mDC_chip_KDUnion_intersect | ChIP,Perturbation | Intersection of the chip and KDUnion networks |

## Human embryonic stem cell (hESC)

| Network | Assay | Description |
| --- | --- | --- |
| hESC_chip | ChIP | Interactions from ENCODE database |
| hESC_escapechip | ChIP | Interactions from ESCAPE database (Mouse orthologs) |
| hESC_chip_union | ChIP | Union of the Encode and Escape networks |
| hESC_logof | Perturbation | Interactions from LOGOF (Mouse orthologs) |
| hESC_nishiyama | Perturbation | Interactions from Nishiyama (Mouse orthologs) |
| hESC_KDUnion | Perturbation | Union of the LOGOF and Nishiyama networks |
| hESC_chip_KDUnion_intersect | ChIP,Perturbation | Intersection of the ChIP and knockdown union networks |
| hESC_chip_KDUnion_c2 | ChIP,Perturbation | Edges from at least two of the above networks (2x ChIP, 2x perturb, or 1 of each) |
| hESC_lit | Multiple | Network of interactions reported in the literature (Mouse orthologs?) |
