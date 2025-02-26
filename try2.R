library(infercnv)

#1
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="../oligodendroglioma_expression_downsampled.counts.matrix",
                                    annotations_file="../oligodendroglioma_annotations_downsampled.txt",
                                    delim="\t",
                                    gene_order_file="../gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt",
                                    ref_group_names=c("Microglia/Macrophage","Oligodendrocytes (non-malignant)")
                                    )

#2
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1, 
                             out_dir="try2",
                             cluster_by_groups=F, 
                             analysis_mode="subclusters",
                             denoise=TRUE,
                             HMM=TRUE,
                             num_threads=1)