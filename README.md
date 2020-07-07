# CaribbeanExtinctions-WTWTW ("Where the wild things were")

### README
#### Files for joining the VertLife mammal tree w/ Caribbean taxa.
#### code: Nathan S. Upham (author line: Samuel T. Turvey, Clare Duncan, Nathan S. Upham, Xavier Harrison, Liliana M. Dávalos)

### FILES

_pruningCode_MamPhy-to-CaribbeanTaxa.R_
- Full code for the pruning of the Mammalia trees, addition of unsampled Caribbean taxa, and resolution of polytomies, with visualization steps for checking the work along the way.

_plottedTree_CaribbeanMam_77taxa_MamPhy_added42taxa_all100_resolvedUltra_FINAL.pdf_
- PDF plot of 100 trees of 'CaribbeanMam_77taxa_MamPhy_added42taxa_random100trees_resolvedUltra_FINAL.nex' for visualizing how the polytomy resolutions influenced the tree topology.

_CaribbeanMam_77taxa_MamPhy_added42taxa_random100trees_resolvedUltra_FINAL.nex_
- Mammal trees (100) with 77 total taxa, now completely bifurcating with the addition of short branch lengths (0.0001 million years) for those nodes that were previously unresolved of length zero (i.e., polytomies). 

_plottedTree_CaribbeanMam_77taxa_MamPhy_added42taxa_all100_wPolytomies.pdf_
- PDF plot of 100 trees of 'CaribbeanMam_77taxa_MamPhy_added42taxa_random100trees_wPolytomies.nex' for visualizing how the polytomies look relative to the addition of 42 taxa.

_CaribbeanMam_77taxa_MamPhy_added42taxa_random100trees_wPolytomies.nex_
- Mammal trees (100) with 77 total taxa, adding 42 and pruning out 3 that were placeholders. Note that the polytomies are unresolved at this point.

_plottedTree_CaribbeanMam_35taxa-plus3_MamPhy_1of100.pdf_
- PDF plot of 1 tree of the 100 in 'CaribbeanMam_35taxa-plus3_MamPhy_random100trees.nex' to see what the starting tree looks like, with renamed taxa from the VertLife MamPhy.

_CaribbeanMam_35taxa-plus3_MamPhy_random100trees.nex_
- Mammal trees (100) pruned to 38 taxa for use in adding unsampled Caribbean taxa in their known clades, probabilistically across the sample of trees.

_taxonList_WTWTW_dataset_July2020.xlsx_
- Annotated listing of how taxon-matching decisions were made.  See the 'mamPhy-taxonList' tab and 'Reason_ifManual' field for justifications.

_nameToAdd_mrca1_mrca2_MamPhy.txt_
- Text file; columns from 'mamPhy-taxonList' tab of 'taxonList_WTWTW_dataset_July2020.xlsx' file. For use in adding the unsampled Caribbean taxa within phylogenetic constraints along the sample of 100 trees (i.e., between mrca1 and mrca2).

_MamPhy_CaribbeanTaxaToPruneOut.txt_
- Text file; Two-column list of correspondences among tips in the mammal tree ("MamPhy_tip") and taxa to add or change their name labels ("nameToAdd"). See the 'taxonList_WTWTW_dataset_July2020.xlsx' file for more details.

_MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_v2_sample100_nexus.trees_
- Raw sample of 100 trees downloaded from http://vertlife.org/phylosubsets/


