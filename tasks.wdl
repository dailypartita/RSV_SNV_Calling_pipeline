version 1.0

task QualityControl {
    input {
        File r1_raw
        File r2_raw
        Int cpu = 4
        Int mem = 20
        Int disk = 20
    }
    command {
        fastp -i ${r1_raw} -I ${r2_raw} -o R1.fq -O R2.fq --thread ~{cpu} --detect_adapter_for_pe --qualified_quality_phred 30 --length_required 36
    }
    runtime {
        cpu: "~{cpu}"
        memory: "~{mem} GB"
        disk: "~{disk} GB"
        continueOnReturnCode: [0, 2]
        docker: "registry-vpc.miracle.ac.cn/biocontainers/fastp:0.23.2--h79da9fb_0"
    }
    output {
        File qc_r1 = "R1.fq"
        File qc_r2 = "R2.fq"
    }
}

task HostRemove {
    input {
        File r1
        File r2
        Int cpu = 4
        Int mem = 20
        Int disk = 20
    }
    command {
        bwa mem -t ${cpu} /opt/hg38/hg38.fasta ${r1} ${r2} > host.sam
        samtools view -@ ${cpu} -b -h -f 4 host.sam > host.bam
        samtools fastq -@ ${cpu} host.bam -1 R1.fq -2 R2.fq -0 /dev/null -s /dev/null
    }
    runtime {
        cpu: "~{cpu}"
        memory: "~{mem} GB"
        disk: "~{disk} GB"
        docker: "registry-vpc.miracle.ac.cn/gznl/rmhost:1.0.0"
    }
    output {
        File rh_r1 = "R1.fq"
        File rh_r2 = "R2.fq"
    }
}

task RefGb2Fa{
    input {
        File ref_gb
        Int cpu = 4
        Int mem = 20
        Int disk = 20
    }
    command {
        any2fasta -u ${ref_gb} > ref.fa
        samtools faidx ref.fa
        picard CreateSequenceDictionary -R ref.fa -O ref.dict
        grep -o '^>.*' ref.fa | sed 's/^>//' > ref.fa.id
    }
    runtime{
        cpu: "~{cpu}"
        memory: "~{mem} GB"
        disk: "~{disk} GB"
        docker: "registry-vpc.miracle.ac.cn/gznl/gb2fa:1.0.0"
    }
    output {
        File ref_fa = "ref.fa"
        File ref_fai = "ref.fa.fai"
        File ref_dic = "ref.dict"
        String fa_id = read_string("ref.fa.id")
    }
}

task BwaMap {
    input {
        File ref_fa
        File r1
        File r2
        Int cpu = 4
        Int mem = 20
        Int disk = 20
    }
    command {
        bwa index ${ref_fa}
        bwa mem -t ~{cpu} ${ref_fa} ${r1} ${r2} > bwamap.sam
    }
    runtime{
        cpu: "~{cpu}"
        memory: "~{mem} GB"
        disk: "~{disk} GB"
        docker: "registry-vpc.miracle.ac.cn/biocontainers/bwa:0.7.8--hed695b0_5"
    }
    output { File sam = "bwamap.sam" }
}

task SortBam {
    input {
        File sam
        Int cpu = 4
        Int mem = 20
        Int disk = 20
    }
    command {
        samtools view -@ ~{cpu} -bS ${sam} > raw.bam
        samtools sort -@ ~{cpu} -o sort.bam raw.bam
        samtools index -@ ~{cpu} sort.bam
    }
    runtime {
        cpu: "~{cpu}"
        memory: "~{mem} GB"
        disk: "~{disk} GB"
        docker: "registry-vpc.miracle.ac.cn/biocontainers/samtools:1.15.1--h1170115_0"
    }
    output {
        File bam = "sort.bam"
        File bai = "sort.bam.bai"
    }
}

task MarkDup {
    input {
        File bam
        Int cpu = 4
        Int mem = 20
        Int disk = 20
    }
    command {
        picard MarkDuplicates -I ${bam} -O dedup.bam -M marked_dup_metrics.txt --REMOVE_DUPLICATES true
    }
    runtime{
        cpu: "~{cpu}"
        memory: "~{mem} GB"
        disk: "~{disk} GB"
        docker: "registry-vpc.miracle.ac.cn/biocontainers/picard:2.26.9--hdfd78af_0"
    }
    output { File dedup_bam = "dedup.bam" }
}

task AddHead {
    input {
        String sample_id
        File bam
        Int cpu = 4
        Int mem = 20
        Int disk = 20
    }
    command {
        samtools addreplacerg -r '@RG\tID:${sample_id}\tSM:${sample_id}' ${bam} -o addreplacerg.bam
        samtools index -@ "~{cpu}" addreplacerg.bam
    }
    runtime {
        cpu: "~{cpu}"
        memory: "~{mem} GB"
        disk: "~{disk} GB"
        docker: "registry-vpc.miracle.ac.cn/biocontainers/samtools:1.15.1--h1170115_0"
    }
    output {
        File headbam = "addreplacerg.bam"
        File headbai = "addreplacerg.bam.bai"
    }
}

task StatDepth {
    input {
        String sample_id
        File bam
        File bai
        Int cpu = 4
        Int mem = 20
        Int disk = 20
    }
    command {
        mv ${bam} map.bam
        mv ${bai} map.bam.bai
        samtools depth -a map.bam > depth_${sample_id}.txt
        python /opt/Draw_SequencingDepth.py depth_${sample_id}.txt depth_${sample_id}
    }
    runtime {
        cpu: "~{cpu}"
        memory: "~{mem} GB"
        disk: "~{disk} GB"
        docker: "registry-vpc.miracle.ac.cn/gznl/statdepth:1.0.0"
    }
    output {
        File depth = "depth_~{sample_id}.txt"
        File depth_png = "depth_~{sample_id}.jpg"
    }
}

task SNV_bcftools {
    input {
        String sample_id
        File ref_fa
        File ref_fai
        File bam
        File bai
        Int cpu = 4
        Int mem = 20
        Int disk = 20
    }
    command <<<
        mv ~{bam} map.bam
        mv ~{bai} map.bam.bai
        mv ~{ref_fa} ref.fa
        mv ~{ref_fai} ref.fa.fai
        bcftools mpileup --threads ~{cpu} --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR -f ref.fa map.bam | bcftools call -mv -o bcftools_raw.vcf
        bcftools +fill-tags bcftools_raw.vcf -Ov -o bcftools_ft.vcf -- -t INFO/AF
        
        python <<CODE
        import vcf
        vcf_in = 'bcftools_ft.vcf'
        vcf_out = 'bcftools.vcf'
        vcf_in_data = vcf.Reader(open(vcf_in, 'r'))
        with open(vcf_out, 'w') as f:
            writer = vcf.Writer(f, vcf_in_data)
            for rec in vcf_in_data:
                ad = rec.genotype("~{sample_id}").data.AD
                rec.INFO['AF'] = round(ad[1] / sum(ad), 4)
                writer.write_record(rec)
        CODE

        bgzip -c bcftools.vcf > bcftools.vcf.gz
        bcftools index -t bcftools.vcf.gz
    >>>
    runtime {
        cpu: "~{cpu}"
        memory: "~{mem} GB"
        disk: "~{disk} GB"
        docker: "registry-vpc.miracle.ac.cn/gznl/bcftools:1.0.1"
    }
    output {
        File vcf_gz = "bcftools.vcf.gz"
        File vcf_index = "bcftools.vcf.gz.tbi"
    }
}

task SNV_freebayes {
    input {
        String sample_id
        File ref_fa
        File ref_fai
        File bam
        File bai
        Int cpu = 4
        Int mem = 20
        Int disk = 20
    }
    command {
        mv ~{bam} map.bam
        mv ~{bai} map.bam.bai
        mv ~{ref_fa} ref.fa
        mv ~{ref_fai} ref.fa.fai
        freebayes -f  ref.fa map.bam > freebayes_raw.vcf
        bcftools +fill-tags freebayes_raw.vcf -Ov -o freebayes_ft.vcf -- -t INFO/AF
        
        python <<CODE
        import vcf
        vcf_in = 'freebayes_ft.vcf'
        vcf_out = 'freebayes.vcf'
        vcf_in_data = vcf.Reader(open(vcf_in, 'r'))
        with open(vcf_out, 'w') as f:
            writer = vcf.Writer(f, vcf_in_data)
            for rec in vcf_in_data:
                ad = rec.genotype("~{sample_id}").data.AD
                rec.INFO['AF'] = round(ad[1] / sum(ad), 4)
                writer.write_record(rec)
        CODE

        bgzip -c freebayes.vcf > freebayes.vcf.gz
        bcftools index -t freebayes.vcf.gz
    }
    runtime {
        cpu: "~{cpu}"
        memory: "~{mem} GB"
        disk: "~{disk} GB"
        docker: "registry-vpc.miracle.ac.cn/gznl/freebayes136pyvcf:1.0.0"
    }
    output {
        File vcf_gz = "freebayes.vcf.gz"
        File vcf_index = "freebayes.vcf.gz.tbi"
    }
}

task SNV_mutect2 {
    input {
        File ref_fa
        File ref_fai
        File ref_dic
        File bam
        File bai
        Int cpu = 4
        Int mem = 20
        Int disk = 20
    }
    command {
        mv ${ref_fa} ref.fa
        mv ${ref_fai} ref.fa.fai
        mv ${ref_dic} ref.dict
        mv ${bam} map.bam
        mv ${bai} map.bam.bai
        gatk Mutect2 -R ref.fa -I map.bam -O mutect2.vcf.gz
    }
    runtime {
        cpu: "~{cpu}"
        memory: "~{mem} GB"
        disk: "~{disk} GB"
        continueOnReturnCode: [0, 2]
        docker: "registry-vpc.miracle.ac.cn/biocontainers/gatk4:4.2.5.0--hdfd78af_0"
    }
    output {
        File vcf_gz = "mutect2.vcf.gz"
        File vcf_index = "mutect2.vcf.gz.tbi"
    }
}

task SNV_varscan {
    input {
        File ref_fa
        File ref_fai
        File bam
        File bai
        Int cpu = 4
        Int mem = 20
        Int disk = 20
    }
    command {
        samtools mpileup -f ~{ref_fa} ~{bam} | varscan pileup2snp --p-value 0.01 > varscan.tsv
    }
    runtime {
        cpu: "~{cpu}"
        memory: "~{mem} GB"
        disk: "~{disk} GB"
        docker: "registry-vpc.miracle.ac.cn/gznl/varscan:2.4.6"
    }
    output {
        File vcf_tsv = "varscan.tsv"
    }
}

task SNV_deepvariant {
    input {
        String sample_id
        File ref_fa
        File ref_fai
        File bam
        File bai
        Int cpu = 4
        Int mem = 20
        Int disk = 20
    }
    command {
        /opt/deepvariant/bin/run_deepvariant \
            --model_type=WGS \
            --ref=~{ref_fa} \
            --reads=~{bam} \
            --output_vcf=deepvariant_raw.vcf
        python <<CODE
        import vcf
        vcf_in = 'deepvariant_raw.vcf'
        vcf_out = 'deepvariant.vcf'
        vcf_in_data = vcf.Reader(open(vcf_in, 'r'))
        with open(vcf_out, 'w') as f:
            writer = vcf.Writer(f, vcf_in_data)
            for rec in vcf_in_data:
                af = float(rec.genotype("~{sample_id}").data.VAF)
                rec.INFO['AF'] = round(af, 4)
                writer.write_record(rec)
        CODE
        bgzip -c deepvariant.vcf > deepvariant.vcf.gz
        bcftools index -t deepvariant.vcf.gz
    }
    runtime {
        cpu: "~{cpu}"
        memory: "~{mem} GB"
        disk: "~{disk} GB"
        docker: "registry-vpc.miracle.ac.cn/gznl/deepvariant:1.6.1"
    }
    output {
        File vcf_gz = "deepvariant.vcf.gz"
        File vcf_index = "deepvariant.vcf.gz.tbi"
    }
}

task CombineVcfs {
    input {
        String sample_id
        File ref_fa
        Array[File] vcf_files
        Array[File] vcf_index
        Int depth_cutoff
        Float freq_cutoff
        Int cpu = 4
        Int mem = 20
        Int disk = 20
    }
    command <<<
        picard MergeVcfs -I ~{sep = " -I " vcf_files} -O combine_raw.vcf
        sed -n 's/^>//p' ~{ref_fa} | sed 's/$/ ./' > chrome_name.txt
        bcftools annotate --rename-chrs chrome_name.txt combine_raw.vcf > combine_rename.vcf
        
        python <<CODE
        import vcf
        sample_id = "~{sample_id}"
        vcf_in = 'combine_rename.vcf'
        vcf_out = 'combine.vcf'
        depth_cutoff = int("~{depth_cutoff}")
        freq_cutoff = float("~{freq_cutoff}")
        vcf_in_data = vcf.Reader(open(vcf_in, 'r'))
        with open(vcf_out, 'w') as f:
            writer = vcf.Writer(f, vcf_in_data)
            for rec in vcf_in_data:
                if 'DP' in rec.INFO:
                    if rec.INFO['DP'] < depth_cutoff:
                        continue
                else:
                    if int(rec.genotype(sample_id).data.DP) < depth_cutoff:
                        continue
                if 'AF' in rec.INFO:
                    if type(rec.INFO['AF']) == type(0.1):
                        if rec.INFO['AF'] > freq_cutoff:
                            writer.write_record(rec)
                    elif type(rec.INFO['AF']) == type(1):
                        continue
                    elif type(rec.INFO['AF']) == type([]):
                        print(rec.INFO['AF'])
                        continue
                    if rec.INFO['AF'] < freq_cutoff or rec.INFO['AF'] == 1:
                        continue
                else:
                    ad = rec.genotype(sample_id).data.AD
                    if type(ad) == type(1):
                        writer.write_record(rec)
                        continue
                    elif type(ad) == type([]):
                        rec.INFO['AF'] = round(ad[1] / sum(ad), 4)
                    if rec.INFO['AF'] < freq_cutoff:
                        continue
                writer.write_record(rec)
        CODE
    >>>
    runtime {
        cpu: "~{cpu}"
        memory: "~{mem} GB"
        disk: "~{disk} GB"
        continueOnReturnCode: [0, 2]
        docker: "registry-vpc.miracle.ac.cn/gznl/combinevcf:1.2"
    }
    output {
        File combine_vcf = "combine.vcf"
        File chrome_name = "chrome_name.txt"
    }
}

task Annotation {
    input {
        String sample_id
        File vcf
        File ref_gb
        File chrome_name
        Int cpu = 4
        Int mem = 20
        Int disk = 20
    }
    command <<<
        # build annotation database
        mkdir -p data/new
        mv ~{ref_gb} data/new/genes.gbk
        echo "new.genome : new" >> snpEff.config
        snpEff build -v new
        # annotate vcf
        mv ~{vcf} raw.vcf
        snpEff -v new raw.vcf > annotation_~{sample_id}_raw.vcf
        awk '{print $2, $1}' ~{chrome_name} > chrome_name_re.txt
        bcftools annotate --rename-chrs chrome_name_re.txt annotation_~{sample_id}_raw.vcf > annotation_~{sample_id}.vcf
    >>>
    runtime {
        cpu: "~{cpu}"
        memory: "~{mem} GB"
        disk: "~{disk} GB"
        docker: "registry-vpc.miracle.ac.cn/gznl/snpeffandbcftools:5.2"
    }
    output {
        File ann_vcf = "annotation_~{sample_id}.vcf"
    }
}