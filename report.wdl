version 1.0

task MakeXlsx {
    input {
        String sample_id
        String snp_only
        File ref_gb
        File depth_txt
        File ann_vcf
        Array[File] vcfgz_list
        File? vcf_tsv
        Int cpu = 4
        Int mem = 20
        Int disk = 20
    }
    command <<<
        mv ~{sep=" " vcfgz_list} ./
        
        python <<CODE
        import vcf
        import numpy as np
        import pandas as pd
        from Bio import SeqIO
        from pathlib import Path

        sample_id = "~{sample_id}"
        ref_gb = "~{ref_gb}"
        depth_file = "~{depth_txt}"
        ann_vcf = "~{ann_vcf}"
        vcf_tsv = "~{vcf_tsv}" if "~{vcf_tsv}" != "" else None
        snp_only = "~{snp_only}"

        vcfgz_list = [str(i) for i in list(Path(".").glob("*.vcf.gz"))]
        depth_data = pd.read_csv(depth_file, sep="\t", header=None)
        vcf_data = vcf.Reader(open(ann_vcf, "r"))

        def get_depth(pos) -> int:
            return int(depth_data[depth_data[1] == pos][2])

        def get_gene(pos) -> str:
            for rec in SeqIO.parse(open(ref_gb, "r"), "genbank"):
                for feature in rec.features:
                    if feature.type == "CDS" and feature.location.start <= pos <= feature.location.end:
                        try:
                            return str(feature.qualifiers["label"][0])
                        except:
                            return str(feature.qualifiers["gene"][0])
            return None

        def get_AF_tsv(pos) -> float:
            vcf_tsv_data = pd.read_csv(vcf_tsv, sep="\t")
            for i in range(len(vcf_tsv_data)):
                if int(vcf_tsv_data['Position'][i]) == int(pos):
                    return round(float(vcf_tsv_data['VarFreq'][i][:-1]) / 100, 4)
            return np.nan

        def get_AF(pos, vcf_gz) -> float:
            print("reading: ", vcf_gz)
            vcf_reader = vcf.Reader(open(vcf_gz, "rb"))
            for rec in vcf_reader:
                if rec.POS == pos:
                    if "AF" in rec.genotype(sample_id).data._fields:
                        af = rec.genotype(sample_id).data.AF
                    if "VAF" in rec.genotype(sample_id).data._fields:
                        af = rec.genotype(sample_id).data.VAF
                    else:
                        if "AF" in rec.INFO.keys():
                            af = rec.INFO['AF']
                    if type(af) == list:
                        return np.average(list(map(float, af))) # 如何处理indel的VAF？暂时取均值
                    else:
                        return float(af)
            return np.nan

        res_list = []
        pos_dic = {}
        for rec in vcf_data:
            if pd.isna(get_gene(rec.POS)):
                continue
            if rec.POS in pos_dic:
                continue
            if snp_only == 'yes':
                if rec.var_type != "snp":
                    continue
            pos_dic[rec.POS] = 1
            tmp_dic = {
                '基因': get_gene(rec.POS),
                '位点': rec.POS,
                '核酸变异': rec.REF + '->' + str(rec.ALT[0]),
                '核酸变异类型': rec.var_type,
                '氨基酸变异': rec.INFO['ANN'][0].split('|')[10].replace('p.', ''),
                '氨基酸变异类型': rec.INFO['ANN'][0].split('|')[1],
                '测序深度': get_depth(rec.POS),
                '平均突变频率': np.nan}
            af_list = []
            for vcf_gz in vcfgz_list:
                sw = Path(vcf_gz).stem.split('.')[0]
                af = get_AF(rec.POS, vcf_gz)
                if not np.isnan(af):
                    af_list.append(af)
                    tmp_dic['突变频率_' + sw] = af
            if vcf_tsv:
                af_vs2 = get_AF_tsv(rec.POS)
                if not np.isnan(af_vs2):
                    af_list.append(af_vs2)
                    tmp_dic['突变频率_varscan2'] = af_vs2
            tmp_dic['平均突变频率'] = np.average(af_list)
            print(tmp_dic, af_list)
            res_list.append(tmp_dic)
        res_df = pd.DataFrame(res_list)
        res_df.to_excel(f"report_~{sample_id}.xlsx", index=False)
        CODE
    >>>
    runtime {
        cpu: "~{cpu}"
        memory: "~{mem} GB"
        disk: "~{disk} GB"
        continueOnReturnCode: [0, 2]
        docker: "registry-vpc.miracle.ac.cn/gznl/viralsnvpyreport:1.0"
    }
    output {
        File xlsx = "report_~{sample_id}.xlsx"
    }
}

task MakePlot {
    input {
        String sample_id
        String fa_id
        File bam
        File bai
        File vcf
        File ref_gb
        Int cpu = 4
        Int mem = 20
        Int disk = 20
    }
    command <<<
        cp ~{bam} map.bam
        cp ~{bai} map.bam.bai
        bcftools reheader -h /opt/head.vcf ~{vcf} > reheader.vcf
        mv ~{ref_gb} ref.gb
        bamdash -b map.bam -r ~{fa_id} -t reheader.vcf ref.gb -e pdf
        mv ~{fa_id}_plot.html report_~{sample_id}.html
        mv ~{fa_id}_plot.pdf report_~{sample_id}.pdf
    >>>
    runtime {
        cpu: "~{cpu}"
        memory: "~{mem} GB"
        disk: "~{disk} GB"
        continueOnReturnCode: [0, 1, 2]
        docker: "registry-vpc.miracle.ac.cn/gznl/bamdash:1.0.2"
    }
    output {
        File html = "report_~{sample_id}.html"
        File pdf = "report_~{sample_id}.pdf"
    }
}
