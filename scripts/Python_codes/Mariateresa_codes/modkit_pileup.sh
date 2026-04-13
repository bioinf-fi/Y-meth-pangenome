for file in Work/Data/BAM/ont_nucflag_alignments_Pille_hallast/results_ceph/*_chrY.bam; do   tag=$(basename "${file}" _chrY.bam);    fasta_ref="Work/Data/FASTA/pangenome_sequences_Jax/ceph/${tag}_chrY.fa";    modkit pileup "${file}" "${file}.bedgraph" --filter-threshold C:0.8 --force-allow-implicit --cpg --ref "${fasta_ref}" --combine-strands --bedgraph; done

