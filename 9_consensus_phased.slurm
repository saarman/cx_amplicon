phased_dir="../cx_amplicon_phased"

for gene in cqm1 ace2 COi; do
  echo "Processing phased samples for gene: $gene"
  region_file="../cx_amplicon_scripts/ref_list_${gene}.txt"
  region=$(cat "$region_file")
  mkdir -p "${out_dir}/${gene}"
  mkdir -p "${vcf_dir}/${gene}"

  for vcf in ${phased_dir}/*.phased.vcf.gz; do
    sample=$(basename "$vcf" .phased.vcf.gz)

    echo "Generating phased consensus for $sample in $gene..."
    bcftools view -r "$region" "$vcf" -Ov -o "${vcf_dir}/${gene}/${sample}_${gene}_phased.vcf"
    bgzip -f "${vcf_dir}/${gene}/${sample}_${gene}_phased.vcf"
    tabix -f -p vcf "${vcf_dir}/${gene}/${sample}_${gene}_phased.vcf.gz"
    rm -f "${vcf_dir}/${gene}/${sample}_${gene}_phased.vcf"

    name=$(bcftools query -l "${vcf_dir}/${gene}/phased/${sample}_${gene}.vcf.gz")

    # Check for phased GT
    phased=$(bcftools query -f '[%GT]\n' "${vcf_dir}/${gene}/${sample}_${gene}_phased.vcf.gz" | grep "|" | wc -l)

    if [ "$phased" -gt 0 ]; then
      echo "Phased genotypes found for $sample. Creating hapA and hapB."

      samtools faidx "$ref" "$region" | \
      bcftools consensus --missing N --absent "-" --sample "$name" --haplotype 1 "${vcf_dir}/${gene}/${sample}_${gene}_phased.vcf.gz" | \
      awk -v sample="$sample" '/^>/{print ">"sample"|hapA"} !/^>/' > "${out_dir}/${gene}/${sample}_${gene}_phased_consensus_hapA.fa"

      samtools faidx "$ref" "$region" | \
      bcftools consensus --missing N --absent "-" --sample "$name" --haplotype 2 "${vcf_dir}/${gene}/${sample}_${gene}_phased.vcf.gz" | \
      awk -v sample="$sample" '/^>/{print ">"sample"|hapB"} !/^>/' > "${out_dir}/${gene}/phased/${sample}_${gene}_phased_consensus_hapB.fa"

    else
      echo "No phased genotypes for $sample. Creating single consensus."
      samtools faidx "$ref" "$region" | \
      bcftools consensus --missing N --absent "-" --iupac-codes --sample "$name" "${vcf_dir}/${gene}/${sample}_${gene}_phased.vcf.gz" | \
      awk -v sample="$sample" '/^>/{print ">"sample"} !/^>/' > "${out_dir}/${gene}/${sample}_${gene}_phased_consensus.fa"
    fi
  done

  # Combine phased outputs
  cat ${out_dir}/${gene}/*_phased_consensus*.fa > ${out_dir}/all_${gene}_phased_consensus.fa
done
