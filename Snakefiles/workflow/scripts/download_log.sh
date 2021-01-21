#################################################################################
## Get list of files from ENCODE Portal for ABC predictions
# review of fastq and bam 

# download all experiments, only file that possesses 
curl -L -u '3WWN6YZV:qynzeuuoauil6utx' -o all_experiments_20201224.tsv "https://www.encodeproject.org/report.tsv?type=Experiment"
grep -E "H3K27ac|DNase-seq|ATAC-seq" all_experiments_20201224.tsv > rel_experiments_20201224.tsv
grep "hg19" rel_experiments_20201224.tsv > hg19_rel_experiments_20201224.tsv
grep "GRCh28" rel_experiments_20201224.tsv > hg38_rel_experiments_20201224.tsv

## Download metadata files from ENCODE portal
wget --quiet -O bam_hg19_DHS.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=hg19&assay_title=DNase-seq&audit.ERROR.category%21=extremely+low+spot+score&aud%5Cit.ERROR.category%21=extremely+low+read+depth&audit.NOT_COMPLIANT.category%21=insufficient+read+depth%5C&files.status=released&files.assembly=hg19" 
wget --quiet -O bam_hg19_ATAC.tsv "https://www.encodeproject.org/metadata/?type=Experiment&files.file_type=bam&assembly=hg19&assay_title=ATAC-seq&files.assembly=hg19&&audit.ERROR.category%21=extremely+low+spot+score&aud%5Cit.ERROR.category%21=extremely+low+read+depth&audit.NOT_COMPLIANT.category%21=insufficient+read+depth%5C&files.status=released" 
wget --quiet -O bam_hg19_H3K27ac.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=hg19&assay_title=Histone+ChIP-seq&target.label=H3K27ac&audit.ERROR.category%21=extremely+low+spot+score&aud%5Cit.ERROR.category%21=extremely+low+read+depth&audit.NOT_COMPLIANT.category%21=insufficient+read+depth%5C&files.status=released&files.assembly=hg19"

# downloading fastq files 
wget --quiet -O fastq_hg19_DHS.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&assembly=hg19&assay_title=DNase-seq&files.file_type=fastq"
wget --quiet -O fastq_hg19_ATAC.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&assembly=hg19&assay_title=ATAC-seq&files.file_type=fastq"
wget --quiet -O fastq_hg19_H3K27ac.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&assembly=hg19&assay_title=Histone+ChIP-seq&target.label=H3K27ac&files.file_type=fastq"

wget --quiet -O bam_GRCh38_DHS.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=GRCh38&assay_title=DNase-seq&audit.ERROR.category%21=extremely+low+spot+score&aud%5Cit.ERROR.category%21=extremely+low+read+depth&audit.NOT_COMPLIANT.category%21=insufficient+read+depth%5C&files.status=released" 
wget --quiet -O bam_GRCh38_ATAC.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=GRCh38&assay_title=ATAC-seq&audit.ERROR.category%21=extremely+low+spot+score&aud%5Cit.ERROR.category%21=extremely+low+read+depth&audit.NOT_COMPLIANT.category%21=insufficient+read+depth%5C&files.status=released" 
wget --quiet -O bam_GRCh38_H3K27ac.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&files.file_type=bam&assembly=GRCh38&assay_title=Histone+ChIP-seq&target.label=H3K27ac&audit.ERROR.category%21=extremely+low+spot+score&aud%5Cit.ERROR.category%21=extremely+low+read+depth&audit.NOT_COMPLIANT.category%21=insufficient+read+depth%5C&files.status=released"

# downloading fastq files
wget --quiet -O fastq_GRCh38_DHS.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&assembly=GRCh38&assay_title=DNase-seq&files.file_type=fastq"
wget --quiet -O fastq_GRCh38_ATAC.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&assembly=GRCh38&assay_title=ATAC-seq&files.file_type=fastq"
wget --quiet -O fastq_GRCh38_H3K27ac.tsv "https://www.encodeproject.org/metadata/?type=Experiment&status=released&assembly=GRCh38&assay_title=Histone+ChIP-seq&target.label=H3K27ac&files.file_type=fastq"
