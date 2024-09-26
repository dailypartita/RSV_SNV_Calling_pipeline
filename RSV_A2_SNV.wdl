version 1.0

import "tasks.wdl" as tasks
import "report.wdl" as report

workflow ViralSNV {
    input {
        String sample_id
        String snp_only = 'no'
        String remove_host = 'no'
        File ref_gb = "s3://bioos-wconlua5eig466kt0jlt0/reference_gb_or_gbk/RSV_A2_NCBI.gb"
        File r1_raw
        File r2_raw
        Int depth_cutoff = 100
        Float freq_cutoff = 0.01
        Int cpu = 32
        Int mem = 128
        Int disk = 256
    }
    call tasks.QualityControl as qc {
        input: r1_raw = r1_raw, r2_raw = r2_raw, cpu = cpu, mem = mem, disk = disk
    }
    call tasks.RefGb2Fa as gb2fa {
        input: ref_gb = ref_gb, cpu = cpu, mem = mem, disk = disk
    }
    if (remove_host == 'yes') {
        call tasks.HostRemove as rmhost {
            input: r1 = qc.qc_r1, r2 = qc.qc_r2, cpu = cpu, mem = mem, disk = disk
        }
    }
    call tasks.BwaMap as map {
        input:
            ref_fa = gb2fa.ref_fa,
            r1 = select_first([rmhost.rh_r1, qc.qc_r1]),
            r2 = select_first([rmhost.rh_r2, qc.qc_r2]),
            cpu = cpu, mem = mem, disk = disk
    }
    call tasks.SortBam as sort {
        input: sam = map.sam, cpu = cpu, mem = mem, disk = disk
    }
    call tasks.MarkDup as dedup {
        input: bam = sort.bam, cpu = cpu, mem = mem, disk = disk
    }
    call tasks.AddHead as addhead {
        input: bam = dedup.dedup_bam, sample_id = sample_id, cpu = cpu, mem = mem, disk = disk
    }
    call tasks.StatDepth as getdepth {
        input: sample_id = sample_id, bam = addhead.headbam, bai = addhead.headbai, cpu = cpu, mem = mem, disk = disk
    }
    call tasks.SNV_bcftools as snv_bcftools {
        input: sample_id = sample_id, ref_fa = gb2fa.ref_fa, ref_fai = gb2fa.ref_fai, bam = addhead.headbam, bai = addhead.headbai, cpu = cpu, mem = mem, disk = disk
    }
    call tasks.SNV_freebayes as snv_freebayes {
        input: sample_id = sample_id, ref_fa = gb2fa.ref_fa, ref_fai = gb2fa.ref_fai, bam = addhead.headbam, bai = addhead.headbai, cpu = cpu, mem = mem, disk = disk
    }
    call tasks.SNV_mutect2 as snv_mutect2 {
        input: ref_fa = gb2fa.ref_fa, ref_fai = gb2fa.ref_fai, ref_dic = gb2fa.ref_dic, bam = addhead.headbam, bai = addhead.headbai, cpu = cpu, mem = mem, disk = disk
    }
    call tasks.SNV_deepvariant as snv_deepvariant {
        input: sample_id = sample_id, ref_fa = gb2fa.ref_fa, ref_fai = gb2fa.ref_fai, bam = addhead.headbam, bai = addhead.headbai, cpu = cpu, mem = mem, disk = disk
    }
    # call tasks.SNV_varscan as snv_varscan {
    #     input: ref_fa = gb2fa.ref_fa, ref_fai = gb2fa.ref_fai, bam = addhead.headbam, bai = addhead.headbai, cpu = cpu, mem = mem, disk = disk
    # }
    call tasks.CombineVcfs as combine {
        input:
            sample_id = sample_id,
            ref_fa = gb2fa.ref_fa,
            vcf_files = [snv_bcftools.vcf_gz, snv_freebayes.vcf_gz, snv_mutect2.vcf_gz, snv_deepvariant.vcf_gz],
            vcf_index = [snv_bcftools.vcf_index, snv_freebayes.vcf_index, snv_mutect2.vcf_index, snv_deepvariant.vcf_index],
            depth_cutoff = depth_cutoff,
            freq_cutoff = freq_cutoff,
            cpu = cpu, mem = mem, disk = disk
    }
    call tasks.Annotation as ann {
        input: sample_id = sample_id, vcf = combine.combine_vcf, ref_gb = ref_gb, chrome_name = combine.chrome_name, cpu = cpu, mem = mem, disk = disk
    }
    call report.MakeXlsx as mkxlsx {
        input:
            sample_id = sample_id,
            snp_only = snp_only,
            ref_gb = ref_gb,
            depth_txt = getdepth.depth,
            ann_vcf = ann.ann_vcf,
            vcfgz_list = [snv_bcftools.vcf_gz, snv_freebayes.vcf_gz, snv_mutect2.vcf_gz, snv_deepvariant.vcf_gz],
            # vcf_tsv = snv_varscan.vcf_tsv,
            cpu = cpu, mem = mem, disk = disk
    }
    call report.MakePlot as plot {
        input:
            sample_id = sample_id,
            fa_id = gb2fa.fa_id,
            bam = addhead.headbam,
            bai = addhead.headbai,
            vcf = ann.ann_vcf,
            ref_gb = ref_gb,
            cpu = cpu, mem = mem, disk = disk
    }
    output {
        File depth_txt = getdepth.depth
        File depth_png = getdepth.depth_png
        File ann_vcf = ann.ann_vcf
        File xlsx = mkxlsx.xlsx
        File html = plot.html
        File pdf = plot.pdf
    }
}
