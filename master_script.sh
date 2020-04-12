#!/bin/bash

# this script is used to write individual running scripts and steps.

YJI_SPACE=/run/media/qu/7de37f23-d1d5-4201-b381-28d1d2810570/yji/

cd $YJI_SPACE/bgi/01
ls -d */ | sed 's/\///' > all_folders_01.list

while read folder; do
    echo '#!/bin/bash' > 0_trimmomatic_${folder}.sh
    echo "cd $folder" >> 0_trimmomatic_${folder}.sh
    echo 'for file in `ls *.1.fq.gz`; do
        base=${file%%.1.fq.gz}
        gunzip -c $file  > ../raw_fq/${base}.1.cp.fq
        gunzip -c ${base}.2.fq.gz > ../raw_fq/${base}.2.cp.fq
        cd .. ' >> 0_trimmomatic_${folder}.sh
    echo '
        java -jar /run/media/qu/7de37f23-d1d5-4201-b381-28d1d2810570/share/Trimmomatic-0.39/Trimmomatic-0.39/trimmomatic-0.39.jar \
            PE -threads 4 -phred33 \
            raw_fq/${base}.1.cp.fq \
            raw_fq/${base}.2.cp.fq \
            trimmed_reads/${base}_p1.fq.gz trimmed_reads/${base}_u1.fq.gz \
            trimmed_reads/${base}_p2.fq.gz trimmed_reads/${base}_u2.fq.gz \
            LEADING:20 TRAILING:20 MINLEN:10' >>  0_trimmomatic_${folder}.sh

    echo '
        echo $base done' >> 0_trimmomatic_${folder}.sh
    echo "cd $folder" >> 0_trimmomatic_${folder}.sh
    echo "done" >> 0_trimmomatic_${folder}.sh
done < all_folders_01.list 

i=1
for job in 0_trimmomatic_*.sh; do
    if [ $i -gt 10 -a $i -lt 21 ]; then
        bash $job > ${job}.out 2>&1 &
        echo $job
        ((i=i+1))
    else
        ((i=i+1))
    fi
done

mkdir fa
cd trimmed_reads
for file in `ls *.fq.gz`; do
    base=${file%%.fq.gz}
    id=${base%%_*}
    species=`grep $id ../id-species.txt | cut -f2`
    gunzip $file
    paste - - - - < ${base}.fq | cut -f1,2 | sed  's/^@/>/' | tr "\t" "\n" | sed 's/ /-/; s/#/-/' >> ../fa/${species}.fa
    echo "$base done"
done

rm raw_fq/*.cp.fq

mkdir done_trimmed-reads
mv trimmed_reads/*.fq done_trimmed-reads
cd done_trimmed-reads
for fq in *.fq; do
    gzip $fq
    echo "$fq gzip trimmed fq done"
done
cd ..

mkdir trf
cd trf
for fasta in `ls ../fa/*.fa`; do
    species0=${fasta%%.fa}
    species=${species0##*/}
    trf_base=${species}.fa.2.5.7.80.10.20.7 
    echo '#!/bin/bash' > ${species}_trf.sh
    echo "trf $fasta 2 5 7 80 10 20 7 -f -h > /dev/null 2>&1" >> ${species}_trf.sh
    #to be run
    echo 'dat_nseq=`grep -c Sequence: ' "${species}.fa.2.5.7.80.10.20.7.dat" '`' >> ${species}_trf.sh
    echo 'nseq=`grep -c ">"' " $fasta " '`' >> ${species}_trf.sh 
    echo 'if [ $nseq -eq $dat_nseq ]; then' >> ${species}_trf.sh
    echo '    echo number of sequence match for '"$species" >> ${species}_trf.sh
    echo "    trf-sort.pl ${trf_base}.dat > ${trf_base}.sorted.dat" >> ${species}_trf.sh
    echo "    trf-redundancy_v2.pl -i ${trf_base}.sorted.dat -f" >> ${species}_trf.sh
    echo '    awk '"'"'$5>=5'"'"' \' >> ${species}_trf.sh
    echo "    ${trf_base}.sorted.mdat > ${trf_base}.sorted.als5.mdat" >> ${species}_trf.sh
    
    echo "    column-transfer.pl -i ${trf_base}.sorted.als5.mdat,15,17 \
                  -r $YJI_SPACE/bgi/msat1-7_L110_collapsed.txt,3,2 \
                  > ${trf_base}.sorted.als5.c17.mdat " >> ${species}_trf.sh
    echo "    column-transfer.pl -i ${trf_base}.sorted.als5.c17.mdat,17,18 \
                  -r $YJI_SPACE/bgi/msat1-7_L110_collapsed.txt,3,2 \
                  > ${trf_base}.sorted.als5.c18.mdat " >> ${species}_trf.sh

    echo "    rm ${trf_base}.sorted.dat" >> ${species}_trf.sh
    echo "    rm ${trf_base}.sorted.mdat" >> ${species}_trf.sh
    echo "    rm ${trf_base}.sorted.als5.mdat" >> ${species}_trf.sh
    echo "    rm ${trf_base}.sorted.als5.c17.mdat" >> ${species}_trf.sh
    echo "    nohup gzip ${trf_base}.dat &" >> ${species}_trf.sh
    echo "fi" >> ${species}_trf.sh
done

## 1. prepare 15-bp flanks
mkdir 1_prepare-flank-15bp
cd 1_prepare-flank-15bp
while read motif; do
    mkdir $motif
done < ../../all_motif.list

for file in `ls ../trf/*.sorted.als5.c18.mdat`; do
    file_name=${file##*/}
    species=${file_name%%.fa.2.5.7.80.10.20.7.sorted.als5.c18.mdat}
    echo '#!/bin/bash' > ${species}_1-prepare-flank-15bp.sh
    while read motif; do
        echo "cd $motif" >> ${species}_1-prepare-flank-15bp.sh
        echo "awk -v motif="'"'$motif'"' "'" '$18==motif' "'" "../../trf/${species}.fa.2.5.7.80.10.20.7.sorted.als5.c18.mdat \
         >> ${species}_${motif}-only_re-fmt.sorted.als5.c18.mdat" >> ${species}_1-prepare-flank-15bp.sh
        echo "cut -f1 ${species}_${motif}-only_re-fmt.sorted.als5.c18.mdat | sort | uniq > ${species}_${motif}-only.list" >> ${species}_1-prepare-flank-15bp.sh 
        echo "faSomeRecords ../../fa/${species}.fa ${species}_${motif}-only.list ${species}_${motif}-only.fasta" >> ${species}_1-prepare-flank-15bp.sh
        echo "bioawk -c fastx '" '{print $name, length($seq) }' "'" "< ${species}_${motif}-only.fasta > ${species}_${motif}-only_length.txt" >> ${species}_1-prepare-flank-15bp.sh
        echo "mdat_flank_filter.pl -i ${species}_${motif}-only_re-fmt.sorted.als5.c18.mdat \
                -l 15 -t ${species}_${motif}-only_length.txt \
                -o ${species}_${motif}-only " >> ${species}_1-prepare-flank-15bp.sh
        echo "mdat_flank2bed.pl -i ${species}_${motif}-only_w-both-flanks.mdat -l 15 \
                -t ${species}_${motif}-only_length.txt \
                > ${species}_${motif}-only_both-flanks.bed" >> ${species}_1-prepare-flank-15bp.sh
        echo "clip_bed.pl -i ${species}_${motif}-only_both-flanks.bed \
                 -s ${species}_${motif}-only.fasta \
                 > ${species}_${motif}-only_both-flanks.fasta" >> ${species}_1-prepare-flank-15bp.sh
        echo "cat_flanks.pl ${species}_${motif}-only_both-flanks.fasta \
                 > ${species}_${motif}-only_cat-flanks.fasta" >> ${species}_1-prepare-flank-15bp.sh
        echo "rm ${species}_${motif}-only.list " >> ${species}_1-prepare-flank-15bp.sh
        echo "rm ${species}_${motif}-only_both-flanks.bed " >> ${species}_1-prepare-flank-15bp.sh
        echo "rm ${species}_${motif}-only_both-flanks.fasta " >> ${species}_1-prepare-flank-15bp.sh
        echo "cd .." >> ${species}_1-prepare-flank-15bp.sh
        echo "echo $species $motif done" >>  ${species}_1-prepare-flank-15bp.sh
    done < ../../all_motif.list
done 

mkdir 2_self-blast
cd 2_self-blast
while read motif; do
    mkdir $motif
done < ../../all_motif.list
          
for file in `ls ../trf/*.sorted.als5.c18.mdat`; do
    file_name=${file##*/}
    species=${file_name%%.fa.2.5.7.80.10.20.7.sorted.als5.c18.mdat}
    echo '#!/bin/bash' > ${species}_2-self-blast.sh
    while read motif; do
        echo "cd $motif" >> ${species}_2-self-blast.sh
        echo "ln -s ../../1_prepare-flank-15bp/${motif}/${species}_${motif}-only_cat-flanks.fasta ." >> ${species}_2-self-blast.sh
        echo "filterN.pl -s ${species}_${motif}-only_cat-flanks.fasta -n 2" >> ${species}_2-self-blast.sh 
        echo "makeblastdb -in ${species}_${motif}-only_cat-flanks_filtered-N.fasta -dbtype nucl" >> ${species}_2-self-blast.sh
        echo "blastn -task blastn-short \
        -query ${species}_${motif}-only_cat-flanks_filtered-N.fasta \
        -db ${species}_${motif}-only_cat-flanks_filtered-N.fasta \
        -outfmt 6 \
        -qcov_hsp_perc 100 \
        -out ${species}_${motif}-only_cat-flanks_all-v-all.fmt6 " >> ${species}_2-self-blast.sh
        echo "silix -i 0.93 -r 0.99 -l 29 -m 0.99 \
                  ${species}_${motif}-only_cat-flanks_filtered-N.fasta \
                  ${species}_${motif}-only_cat-flanks_all-v-all.fmt6 \
                  > ${species}_${motif}-only_cat-flanks_silix.fnodes" >> ${species}_2-self-blast.sh
        echo "mkdir ${species}_${motif}-only_silix_clusters" >> ${species}_2-self-blast.sh
        echo "cd ${species}_${motif}-only_silix_clusters" >> ${species}_2-self-blast.sh
        echo "silix-split ../${species}_${motif}-only_cat-flanks_filtered-N.fasta \
                  ../${species}_${motif}-only_cat-flanks_silix.fnodes" >> ${species}_2-self-blast.sh
        echo "for i in *.fasta; do " >> ${species}_2-self-blast.sh
        echo '  num=`grep -c ">" $i`' >> ${species}_2-self-blast.sh
        echo '  if [ "$num" -le 100 ]; then ' >> ${species}_2-self-blast.sh
        echo '      head -n2 $i >> ' "../${species}_${motif}-only_cat-flanks_reduced-100.fasta" >> ${species}_2-self-blast.sh
        echo '  fi' >> ${species}_2-self-blast.sh
        echo 'done' >> ${species}_2-self-blast.sh
        echo "cd .." >> ${species}_2-self-blast.sh
        echo "tar -cf ${species}_${motif}-only_silix_clusters.tar ${species}_${motif}-only_silix_clusters" >> ${species}_2-self-blast.sh
        echo 'if [ $? -eq 0 ]; ' "then rm -r ${species}_${motif}-only_silix_clusters; fi" >> ${species}_2-self-blast.sh 
        echo "rm ${species}_${motif}-only_cat-flanks_filtered-N.fasta.*" >> ${species}_2-self-blast.sh
        echo "rm ${species}_${motif}-only_cat-flanks.fasta" >> ${species}_2-self-blast.sh
        echo "cd .." >> ${species}_2-self-blast.sh
        echo "echo $species $motif self blast done" >>  ${species}_2-self-blast.sh
    done < ../../all_motif.list
done


cd trf


