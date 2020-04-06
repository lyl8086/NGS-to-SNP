#!/bin/bash

BWA="bwa-mem2"
SAMTOOLS="samtools"
BCFTOOLS="bcftools"
BEDTOOLS="bedtools"
FASTP="fastp"
MULTIQC="multiqc"
QUALIMAP="/opt/bio/qualimap_v2.2.1/qualimap"
PICARD="/opt/bio/picard/picard.jar"
in=$1
out=$2
popmap=$3
T=$4
ref=$5
stage=$6
typ=$7
[ $# -lt 5 ] && echo $0 [inpath] [outpath] [popmap] [threads] [reference] && exit 1
###### initial ######
[ -z "$typ" ] && typ='old'
let t=$T/2
base_dir=$out
mkdir -p $out/log $out/bam $out/fq $out/multi_qc $out/SNP $out/bam_qc
mkdir -p $out/log/fastp $out/log/bamqc $out/log/bwa $out/log/samtools $out/log/bcftools 
sed -i -E "s/\s+/\t/g" $popmap
sed -i -E "/^$/d" $popmap
tot=`awk '{c++}END{print c}' $popmap`
######         ######

function fastp_process {
    # fastp processing...
    local j
    echo "[+] Processing raw reads..."
    in=$in
    out=$base_dir/fq
    :>$base_dir/log/fastp/fastp.logs
    for i in `cat $popmap|awk '{print $1}'`; do 
        R1=${i}_1.fq.gz
        R2=${i}_2.fq.gz
        let j=j+1
        # check
        [ ! -f "$in/$R1" ] || [ ! -f "$in/$R2" ] && echo -e "\033[33m[!] no fastq \"$R1 and $R2\" in \"$in\"!\033[0m" && exit 1
        
        echo -n "    Running fastp on $i, $j of $tot..."
        
        $FASTP -i $in/$R1 -I $in/$R2 -o $out/$R1 -O $out/$R2 \
        -l 38 -w 6 -h $base_dir/log/fastp/$i.fastp.html -j $base_dir/log/fastp/$i.fastp.json 2>>$base_dir/log/fastp/fastp.logs
        echo 'done.'
    done
    echo
}

function multi_qc {
    # fastq qc
    echo "[+] Generating multiqc reports..."
    in="$base_dir/log"
    out="$base_dir/multi_qc"
    $MULTIQC -f -x '*.log' -x '*.bed' -o $out $in
    cp -f $out/multiqc_report.html $base_dir/NGSpipeline.qc.html
    echo
}

function bwa_mem {
    # bwa mapping...
    local j
    [ ! -f $ref ] && echo "[!] no reference file in $ref!" && exit 1
    :>$base_dir/log/bwa/bwa.logs 
    :>$base_dir/log/samtools/samtools.logs
    :>$base_dir/log/samtools/picard.logs
    [ ! -f "$ref.bwt.2bit.64" ] && \
    echo -n "[+] Generating bwa index..." && \
    $BWA index $ref >&$base_dir/log/bwa/bwa.logs&& \
    echo "done."
    echo "[+] Running bwa mapping on clean reads..."
    in="$base_dir/fq"
    out="$base_dir/bam"

    for i in `cat $popmap|awk '{print $1}'`; do 
        R1=${i}_1.fq.gz
        R2=${i}_2.fq.gz
        
        # check
        [ ! -f "$in/$R1" ] || [ ! -f "$in/$R2" ] && echo -e "\033[33m[!] no fastq \"$R1 and $R2\" in \"$in\"!\033[0m" && exit 1
        let j=j+1
        # bwa mem and sort by name.
        echo -n "    Processing $i, $j of $tot...Mapping...";
        
        $BWA mem -R "@RG\tID:$i\tPL:ILLUMINA\tSM:$i" -t $T -M $ref $in/$R1 $in/$R2 \
        2>>$base_dir/log/bwa/bwa.logs | \
        $SAMTOOLS sort -n -m10G -@4 -o $out/$i.sort.bam 2>>$base_dir/log/samtools/samtools.logs;
        
        markdup_samtools
        
        echo "Generating stats...(done in sub pid)"
        # bed
        $BEDTOOLS bamtobed -i $out/$i.bam | \
        $BEDTOOLS merge -i - > $base_dir/log/samtools/$i.bed && \
        # flag
        $SAMTOOLS flagstat $out/$i.bam >$base_dir/log/samtools/$i && \
        # stats
        $SAMTOOLS stats $out/$i.bam >$base_dir/log/samtools/$i.bam.stat && \
        # idxstats
        $SAMTOOLS idxstats $out/$i.bam >$base_dir/log/samtools/$i.bam.idxstat &
    done
    
    for i in `cat $popmap|awk '{print $1}'`;
    do
        log="$i.bam.idxstat"
        while [ ! -e "$base_dir/log/samtools/$log" ];
        do
            printf "    Waiting the Job done [$i]...\r"
        done
        wait
    done
    echo
}

function markdup_samtools {

    echo -n "Samtools markdup..."
    # fixmate.
    $SAMTOOLS fixmate -m -@$t $out/$i.sort.bam - \
    2>>$base_dir/log/samtools/samtools.logs | \
    # sort by coordinates, best compress.
    $SAMTOOLS sort -m 10G -@$t \
    2>>$base_dir/log/samtools/samtools.logs | \
    
    # mark duplications, create index.
    $SAMTOOLS markdup -@$t - $out/$i.bam \
    --write-index \
    2>>$base_dir/log/samtools/samtools.logs \
    -f $base_dir/log/samtools/$i.markdup.stat && \
    rm $out/$i.sort.bam
}

function markdup_picard {
            
    # markdup using picard tools
    echo -n "Picard markdup..."
    $SAMTOOLS sort -m 10G -@$T \
    -o $out/$i.sort.bam $out/$i.sort.bam 2>> $base_dir/log/samtools/samtools.logs
    #
    java -Xmx60G \
    -Djava.io.tmpdir=temp/ \
    -XX:ParallelGCThreads=$T \
    -jar $PICARD MarkDuplicates \
    TMP_DIR=temp/ \
    I=$out/$i.sort.bam \
    O=$out/$i.bam \
    ASSUME_SORT_ORDER=coordinate \
    M=$base_dir/log/samtools/${i}_dup_metrics.txt \
    2>>$base_dir/log/samtools/picard.logs
    rm $out/$i.sort.bam
    $SAMTOOLS index -@$T $out/$i.bam   
}

function bam_qc {
    # bam qc
    echo -n "[+] Generating qc on bam files..."
    in="$base_dir/bam"
    out="$base_dir/bam_qc"
    ls $in/*.bam >bam.lst
    
    # check
    [ ! -s bam.lst ] && echo -e "\033[33m[!] no bam files in bam.lst!\033[0m" && exit 1
    
    awk -v p=$in '{print $1"\t"p"/"$1".bam"}' $popmap >bamfile
    out="$base_dir/bam_qc" # out dir for bamqc
    $QUALIMAP multi-bamqc -r -c -d bamfile -outdir $out \
    --java-mem-size=4G -outformat PDF:HTML >& $base_dir/log/bamqc/qualimap.logs
    mv bamfile bam.lst $out
    cp --remove-destination -r -l $base_dir/bam/*_stats $base_dir/log/bamqc/
    echo 'done.'
    echo
}

function call_snp {
    # bcftools call SNP...
    echo "[+] Calling SNP on bam files..."
    in="$base_dir/bam"
    out="$base_dir/SNP"
    awk -v p=$in '{print p"/"$1".bam"}' $popmap >bam.lst
    
    # check.
    [ ! -s bam.lst ] && echo -e "\033[33m[!] no bam files in bam.lst!\033[0m" && exit 1
    [ ! -f $ref ] && echo -e "\033[33m[!] no reference file in $ref!\033[0m" && exit 1
    :>$base_dir/log/bcftools/bcftools.logs
    # new bcf.
    bcf_v=`$BCFTOOLS --version |head -1|awk '{print $2}'`
    if [ `awk -v v=$bcf_v 'BEGIN {if (v<1.10) {print 0} else {print 1}}'` -ne 0 ];then
        flg="-G $popmap"
    else
        flg=''
    fi
    # call snp.
    if [ $typ == 'old' ];then
        # get fasta names.
        $SAMTOOLS faidx $ref
        cut -f1 $ref.fai | shuf > tmp.ref.lst
        # split into 8 chunks.
        split -d -n r/8 tmp.ref.lst ref_bed
        :>call_snp.cmd
        :>call_snp.cmd.completed
        # save cmd into one file.
        for f in ref_bed*; do
            let i=i+1
            if [ $i -eq 1 ];then
                echo "$BCFTOOLS mpileup -Ou -b bam.lst -d 1000 -f $ref -R $f -q 20 -a DP,AD,INFO/AD -p -P ILLUMINA --ff 0x400 --ff 0x100 2>>$base_dir/log/bcftools/bcftools.logs | $BCFTOOLS call -M -Oz -f GQ -v -m -o $i.snp.vcf.gz $flg 2>>$base_dir/log/bcftools/bcftools.logs" >> call_snp.cmd
            else
                echo "$BCFTOOLS mpileup -Ou -b bam.lst -d 1000 -f $ref -R $f -q 20 -a DP,AD,INFO/AD -p -P ILLUMINA --ff 0x400 --ff 0x100 2>>$base_dir/log/bcftools/bcftools.logs | $BCFTOOLS call -M -Ov -f GQ -v -m $flg 2>>$base_dir/log/bcftools/bcftools.logs | grep -v '^#' | bgzip -c >$i.snp.vcf.gz" >> call_snp.cmd
            fi
        done
        # ParaFly...
        ParaFly -c call_snp.cmd -CPU 8 -v
        # combine the bgzipped file into one.
        cat {1..8}.snp.vcf.gz > $base_dir/SNP/raw.snp.vcf.gz
        # rm tmp files.
        mv tmp.ref.lst call_snp.cmd call_snp.cmd.completed \
        bam.lst ref_bed* $base_dir/SNP/
        rm {1..8}.snp.vcf.gz
    else    
        # new bcftools
        [ -z "$flg" ] && echo "Please use bcftools >= 1.10" && exit 1
        $BCFTOOLS mpileup -Ou -b bam.lst -d 1000 -f $ref \
        --threads $t \
        -q 20 -a DP,AD,INFO/AD -p -P ILLUMINA \
        --ff 0x400 --ff 0x100 \
        2>>$base_dir/log/bcftools/bcftools.logs | \
        $BCFTOOLS call -M -Oz -f GQ -v -m -o $base_dir/SNP/raw.snp.vcf.gz \
        -G $popmap --threads $t \
        2>>$base_dir/log/bcftools/bcftools.logs
    fi
    
    # generating stats.
    echo -n "[+] Generating vcf stats..."
    tabix -p vcf $base_dir/SNP/raw.snp.vcf.gz
    for i in `cat $popmap|awk '{print $1}'`; do
        $BCFTOOLS stats --threads $T --af-bins 0.05,0.1,0.5,1 -s $i \
        -F $ref $base_dir/SNP/raw.snp.vcf.gz \
        >$base_dir/log/bcftools/$i.stat \
        2>>$base_dir/log/bcftools/bcftools.logs
        sed -i "s/raw.snp.vcf.gz/$i.vcf.gz/" $base_dir/log/bcftools/$i.stat
    done
    $BCFTOOLS stats --threads $T --af-bins 0.05,0.1,0.5,1 \
        -s - \
        -F $ref $base_dir/SNP/raw.snp.vcf.gz \
        >$base_dir/log/bcftools/all.stat \
        2>>$base_dir/log/bcftools/bcftools.logs
    plot-vcfstats -p $base_dir/SNP/stats/ \
    -s $base_dir/log/bcftools/all.stat \
    2>>$base_dir/log/bcftools/bcftools.logs
    echo "done."
    echo
}

case $stage in
    qc)
    fastp_process
    multi_qc
    ;;
    
    map)
    bwa_mem
    bam_qc
    multi_qc
    ;;
    
    snp)
    call_snp
    multi_qc
    ;;
    mapsnp)
    bwa_mem
    call_snp
    bam_qc
    multi_qc
    ;;
    
    *)
    fastp_process
    bwa_mem
    call_snp
    bam_qc
    multi_qc
    ;;
esac  
echo ":-) NGS pipeline done."
