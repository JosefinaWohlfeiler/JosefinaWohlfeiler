###############    STAR ######################

mkdir star

mkdir star/starIndex

############################    generating genome indexes files


STAR --runThreadN 28 --runMode genomeGenerate \
--genomeFastaFiles \
../../../../../vcru_share_s/carrot/LNRQ01/data/130.DHv2.fna \
--sjdbGTFfile ../../51.geneswithphase.gtf \
--genomeDir starIndex \
--sjdbOverhang 150




#dio este error


Jun 03 19:33:41 ..... started STAR run
Jun 03 19:33:41 ... starting to generate Genome files
Jun 03 19:34:09 ... starting to sort Suffix Array. This may take a long time...
Jun 03 19:34:26 ... sorting Suffix Array chunks and saving them to disk...
Jun 03 19:36:53 ... loading chunks from disk, packing SA...
Jun 03 19:37:18 ... finished generating suffix array
Jun 03 19:37:18 ... generating Suffix Array index
Jun 03 19:39:12 ... completed Suffix Array index
Jun 03 19:39:12 ..... processing annotations GTF

Fatal INPUT FILE error, no exon lines in the GTF file: ../../51.geneswithphase.gtf
Solution: check the formatting of the GTF file, it must contain some lines with exon in the 3rd column.
          Make sure the GTF file is unzipped.
          If exons are marked with a different word, use --sjdbGTFfeatureExon .

Jun 03 19:39:12 ...... FATAL ERROR, exiting





QUIZAS HACR UN SED REEMPLAZANDO TODOS LOS CDS POR EXON EN EL GTF FILE Y USAR ESE

POR LAS DUDAS ELIMINAR TODA LA CARPETA CON SUS SUBCARPETAS ANTES DE VOLVER A CORRER STAR (NO LO HICE ASI QUE CREO QUE SOBREESCRIBIO)



#hice esto
sed -i 's/CDS/exon/g' "../../51.geneswithphase.gtf" | more

#y lo volvi a correr

Jun 05 11:20:40 ..... started STAR run
Jun 05 11:20:40 ... starting to generate Genome files
Jun 05 11:21:06 ... starting to sort Suffix Array. This may take a long time...
Jun 05 11:21:22 ... sorting Suffix Array chunks and saving them to disk...
Jun 05 11:23:58 ... loading chunks from disk, packing SA...
Jun 05 11:24:24 ... finished generating suffix array
Jun 05 11:24:24 ... generating Suffix Array index
Jun 05 11:26:25 ... completed Suffix Array index
Jun 05 11:26:25 ..... processing annotations GTF
Jun 05 11:26:28 ..... inserting junctions into the genome indices
Jun 05 11:28:41 ... writing Genome to disk ...
Jun 05 11:28:45 ... writing Suffix Array to disk ...
Jun 05 11:28:50 ... writing SAindex to disk
Jun 05 11:28:53 ..... finished successfully
 
 



####################        Running mapping jobs


mkdir star/starMapping



cd star/starMapping


STAR --outSAMstrandField intronMotif \
--outFilterIntronMotifs RemoveNoncanonical \
--runThreadN 28 \
--genomeDir ../starIndex \
--sjdbGTFfile ../../../51.geneswithphase.gtf \
--outSAMtype BAM Unsorted \
--readFilesCommand zcat ../../trimmed_fastq/B51-1_S42_L003_R1_001.trim.paired.fq.gz ../../trimmed_fastq/B51-1_S42_L003_R2_001.trim.paired.fq.gz

#me dio

Jun 06 13:39:53 ..... started STAR run
Jun 06 13:39:53 ..... loading genome
Jun 06 13:40:00 ..... processing annotations GTF
Jun 06 13:40:03 ..... inserting junctions into the genome indices
Jun 06 13:40:38 ..... started mapping
gzip: Read2.gz: No such file or directory
gzip: Read1.gz: No such file or directory
Jun 06 14:18:32 ..... finished successfully


###########   por las dudas voy a borrar todo lo que esta dentro de starMapping y lo voy a correr de nuevo

# me dio

Jun 06 19:55:12 ..... started STAR run
Jun 06 19:55:12 ..... loading genome
Jun 06 19:55:19 ..... processing annotations GTF
Jun 06 19:55:22 ..... inserting junctions into the genome indices
Jun 06 19:55:53 ..... started mapping
gzip: gzip: Read2.gz: No such file or directory
Read1.gz: No such file or directory
Jun 06 20:33:48 ..... finished successfully

##### al ver este output

more Log.final.out 

# me da un % altisimo de unmapped reads     % of reads unmapped: too short  86.16%

 Started job on |	Jun 06 19:55:12
                             Started mapping on |	Jun 06 19:55:53
                                    Finished on |	Jun 06 20:33:48
       Mapping speed, Million of reads per hour |	147.63

                          Number of input reads |	93296530
                      Average input read length |	281
                                    UNIQUE READS:
                   Uniquely mapped reads number |	18014
                        Uniquely mapped reads % |	0.02%
                          Average mapped length |	202.46
                       Number of splices: Total |	199
            Number of splices: Annotated (sjdb) |	56
                       Number of splices: GT/AG |	180
                       Number of splices: GC/AG |	19
                       Number of splices: AT/AC |	0
               Number of splices: Non-canonical |	0
                      Mismatch rate per base, % |	4.92%
                         Deletion rate per base |	0.01%
                        Deletion average length |	1.84
                        Insertion rate per base |	0.02%
                       Insertion average length |	1.08
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |	12828337
             % of reads mapped to multiple loci |	13.75%
        Number of reads mapped to too many loci |	465
             % of reads mapped to too many loci |	0.00%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |	0.00%
                 % of reads unmapped: too short |	86.16%
                     % of reads unmapped: other |	0.07%
                                  CHIMERIC READS:
                       Number of chimeric reads |	0
                            % of chimeric reads |	0.00%


# entonces busque el google y sugirieron agregar estos parametros


--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0

# borro todo y lo vuelvo a correr

STAR --outSAMstrandField intronMotif \
--outFilterIntronMotifs RemoveNoncanonical \
--runThreadN 28 \
--genomeDir ../starIndex \
--sjdbGTFfile ../../../51.geneswithphase.gtf \
--outSAMtype BAM Unsorted \
--readFilesCommand zcat ../../trimmed_fastq/B51-1_S42_L003_R1_001.trim.paired.fq.gz ../../trimmed_fastq/B51-1_S42_L003_R2_001.trim.paired.fq.gz \
--outFilterScoreMinOverLread 0 \
--outFilterMatchNminOverLread 0 \
--outFilterMatchNmin 0

# dio lo mismo y tambien me da 

gzip: Read1.gz: No such file or directory
gzip: Read2.gz: No such file or directory

--outFilterScoreMinOverLread: command not found

#esto creo q es porque me falto el enter   


# va de nuevo
#me dio
Jun 07 11:07:56 ..... started STAR run
Jun 07 11:07:56 ..... loading genome
Jun 07 11:08:04 ..... processing annotations GTF
Jun 07 11:08:07 ..... inserting junctions into the genome indices
Jun 07 11:08:43 ..... started mapping
gzip: gzip: Read1.gz: No such file or directory
Read2.gz: No such file or directory
Jun 07 11:56:05 ..... finished successfully



######### ahora solucionar este error
gzip: gzip: Read2.gz: No such file or directory
Read1.gz: No such file or directory
Jun 06 20:33:48 ..... finished successfully

#lo busque en google, voy a ver si me sale con lo que encontre

mkdir starMapping2
cd StarMapping2


R1=/dcdata2/josefina/rnaseq/trimmed_fastq/B51-1_S42_L003_R1_001.trim.paired.fq.gz
R2=/dcdata2/josefina/rnaseq/trimmed_fastq/B51-1_S42_L003_R2_001.trim.paired.fq.gz

STAR --outSAMstrandField intronMotif \
--outFilterIntronMotifs RemoveNoncanonical \
--runThreadN 28 \
--genomeDir ../starIndex \
--sjdbGTFfile ../../../51.geneswithphase.gtf \
--outSAMtype BAM Unsorted \
--readFilesCommand zcat $R1 $R2 \
--outFilterScoreMinOverLread 0 \
--outFilterMatchNminOverLread 0 \
--outFilterMatchNmin 0

# me dio lo mismo q el otro

Jun 07 11:07:56 ..... started STAR run
Jun 07 11:07:56 ..... loading genome
Jun 07 11:08:04 ..... processing annotations GTF
Jun 07 11:08:07 ..... inserting junctions into the genome indices
Jun 07 11:08:43 ..... started mapping
gzip: gzip: Read1.gz: No such file or directory
Read2.gz: No such file or directory
Jun 07 11:56:05 ..... finished successfully

# al parecer si reconoce los files igualmente, a pesar de el error

#hacer una prueba deszipeando la muestra para ver si no me da ese error











###########################################################################

##############                    Mapping multiple files in one run

cd star

mkdir mapped

cd mapped


for file in ../../trimmed_fastq/*R1*trim.paired.fq.gz
do STAR --outSAMstrandField intronMotif \
--outFilterIntronMotifs RemoveNoncanonical \
--runThreadN 28 \
--genomeDir ../starIndex \
--sjdbGTFfile ../../../51.geneswithphase.gtf \
--outSAMtype BAM Unsorted \
--readFilesCommand zcat $file ${file/R1/R2} \
--outFilterScoreMinOverLread 0 \
--outFilterMatchNminOverLread 0 \
--outFilterMatchNmin 0 \
--outFileNamePrefix ./${file%R1*}
done

#deje corriendo esto domingo 14.30, termino lunes



###################### LO QUE PASO CON ESTE LOOP FUE QUE ME FALTO PONER UN _ ANTES DE R1   for file in ../../trimmed_fastq/*_R1*trim.paired.fq.gz

#################################    ENTONCES TOMO A LAS CR1 COMO UN ANALISIS DENTRO DEL LOOP ASI QUE TIENE ERRORES





#########   El input de cuffdiff requires the BAM files to be sorted by genomic coordinates, let's do this using samtools

#lo hago con una sola muestra primero B51-1 que es la que esta en este directorio

cd star/starMapping-prueba


#recordar que esta muestra prueba es B51-1_S42_L003

samtools sort -O bam -o Aligned.out.sorted.bam -T tmp Aligned.out.bam

#me dio

[bam_sort_core] merging from 224 files and 1 in-memory blocks...



#Ahora hacer un sorted de B54 para tener dos muestras que comparar con cuffdiff

#primero copiar y pegar los files de B54 outputs de STAR en star/starMapping

cp B54_S12_L001*.out* ../star/starMapping-prueba

cp -r B54_S12_L001__STARgenome ../star/starMapping-prueba


#Ahora si hago el sorted

cd star/starMapping-prueba

samtools sort -O bam -o B54_S12_L001_Aligned.out.sorted.bam -T tmp B54_S12_L001_Aligned.out.bam



#Entonces ahora ya tengo 2 muestras para comparar su expresion con cuffdiff

mkdir star/cuffdiff_star_prueba

cd star/cuffdiff_star_prueba


cuffdiff -o star/cuffdiff_star_prueba \
-L B51,B54 \
-b ../../../../../vcru_share_s/carrot/LNRQ01/data/130.DHv2.fna \
-u --library-type fr-unstranded \
../../../51.geneswithphase.gtf \
../starMapping/Aligned.out.sorted.bam \
../starMapping/B54_S12_L001_Aligned.out.sorted.bam


#dio
>                     Estimated Mean: 447.27
>                  Estimated Std Dev: 133.44
> Map Properties:
>       Normalized Map Mass: 18892975.93
>       Raw Map Mass: 58339961.81
>       Number of Multi-Reads: 44902626 (with 105836675 total hits)
>       Fragment Length Distribution: Empirical (learned)
>                     Estimated Mean: 312.17
>                  Estimated Std Dev: 178.14
[09:56:48] Calculating preliminary abundance estimates
> Processed 32045 loci.                        [*************************] 100%
[11:35:54] Learning bias parameters.
[12:29:58] Testing for differential expression and regulation in locus.
> Processed 32045 loci.                        [*************************] 100%
Performed 644 isoform-level transcription difference tests
Performed 0 tss-level transcription difference tests
Performed 644 gene-level transcription difference tests
Performed 0 CDS-level transcription difference tests
Performed 0 splicing tests
Performed 0 promoter preference tests
Performing 0 relative CDS output tests
Writing isoform-level FPKM tracking
Writing TSS group-level FPKM tracking
Writing gene-level FPKM tracking
Writing CDS-level FPKM tracking
Writing isoform-level count tracking
Writing TSS group-level count tracking
Writing gene-level count tracking
Writing CDS-level count tracking
Writing isoform-level read group tracking
Writing TSS group-level read group tracking
Writing gene-level read group tracking
Writing CDS-level read group tracking
Writing read group info
Writing run info



#revisar el script anterior. Al parecer creo una carpeta dentro de esa con esos nombres, asi que no volver a poner los nombres.






#Despues para poder ver los genes diferencialmente expresados
#to see which are the most significantly differentially expressed genes
#use a sort command to sort the file and write the sorted file in a different one called gene_exp_sorted.diff


sort -t$'\t' -g -k 13 star/cuffdiff_star_prueba/gene_exp.diff > star/cuffdiff_star_prueba/gene_exp_sorted.diff

#aca me ordeno por p-value de menor a mayor asi que las primeras son las significativas


#entonces este fue el de prueba, comparacion de B51 con B54 (esta en /dcdata2/josefina/rnaseq/star/cuffdiff_star_prueba/star/cuffdiff_star_prueba)




#hacer loop de BAM files sorted by genomic coordinates para todas las muestras asi tengo los input de cuffdiff
#primero mover todos los output files de STAR a la carpeta mapped


mv *out* ../star/mapped
mv *STAR* ../star/mapped/

cd star/mapped

for file in *Aligned.out.bam
do samtools sort -O bam -o ${file%bam}sorted.bam -T tmp $file
done
#deje corriendo martes 15 hs

#miercoles 16:30 hs van 14 muestras (de 31)
#miercoles 20 hs van 16
#miercoles 22.30 hs van 18
#jueves 11.30 hs van 26
*jueves 16.45 hs van 29

#nota: quizas tendria que haber puesto un threads 28 y demoraba menos

#copio la salida
[bam_sort_core] merging from 224 files and 1 in-memory blocks...
[bam_sort_core] merging from 237 files and 1 in-memory blocks...
[bam_sort_core] merging from 485 files and 1 in-memory blocks...
[bam_sort_core] merging from 250 files and 1 in-memory blocks...
[bam_sort_core] merging from 284 files and 1 in-memory blocks...
[bam_sort_core] merging from 166 files and 1 in-memory blocks...
[bam_sort_core] merging from 218 files and 1 in-memory blocks...
[bam_sort_core] merging from 258 files and 1 in-memory blocks...
[bam_sort_core] merging from 117 files and 1 in-memory blocks...
[bam_sort_core] merging from 174 files and 1 in-memory blocks...
[bam_sort_core] merging from 137 files and 1 in-memory blocks...
[bam_sort_core] merging from 104 files and 1 in-memory blocks...
[bam_sort_core] merging from 178 files and 1 in-memory blocks...
[bam_sort_core] merging from 178 files and 1 in-memory blocks...
[bam_sort_core] merging from 153 files and 1 in-memory blocks...
[bam_sort_core] merging from 226 files and 1 in-memory blocks...
[bam_sort_core] merging from 149 files and 1 in-memory blocks...
[bam_sort_core] merging from 70 files and 1 in-memory blocks...
[bam_sort_core] merging from 66 files and 1 in-memory blocks...
[bam_sort_core] merging from 137 files and 1 in-memory blocks...
[bam_sort_core] merging from 122 files and 1 in-memory blocks...
[bam_sort_core] merging from 70 files and 1 in-memory blocks...
[bam_sort_core] merging from 335 files and 1 in-memory blocks...
[bam_sort_core] merging from 170 files and 1 in-memory blocks...
[bam_sort_core] merging from 172 files and 1 in-memory blocks...
[bam_sort_core] merging from 375 files and 1 in-memory blocks...
[bam_sort_core] merging from 386 files and 1 in-memory blocks...
[bam_sort_core] merging from 224 files and 1 in-memory blocks...
[bam_sort_core] merging from 473 files and 1 in-memory blocks...
[bam_sort_core] merging from 248 files and 1 in-memory blocks...
[bam_sort_core] merging from 380 files and 1 in-memory blocks...
[bam_sort_core] merging from 291 files and 1 in-memory blocks...





##############################   me di cuenta de que me falto un parametro asi que volver a correr el loop

#########  primero hacer una prueba

mkdir mapped-quant

cd mapped-quant

R1=/dcdata2/josefina/rnaseq/trimmed_fastq/B51-1_S42_L003_R1_001.trim.paired.fq.gz
R2=/dcdata2/josefina/rnaseq/trimmed_fastq/B51-1_S42_L003_R2_001.trim.paired.fq.gz

STAR --outSAMstrandField intronMotif \
--outFilterIntronMotifs RemoveNoncanonical \
--runThreadN 28 \
--genomeDir ../starIndex \
--sjdbGTFfile ../../../51.geneswithphase.gtf \
--outSAMtype BAM Unsorted \
--readFilesCommand zcat $R1 $R2 \
--outFilterScoreMinOverLread 0 \
--outFilterMatchNminOverLread 0 \
--outFilterMatchNmin 0 \
--quantMode TranscriptomeSAM GeneCounts

#deje corriendo 11 jul 18.15 hs


#en Reads.per.Gene me dio esto

N_unmapped	71704184	71704184	71704184
N_multimapping	81877861	81877861	81877861
N_noFeature	49944	57743	58042
N_ambiguous	0	0	0


#entonces por las dudas probar lo mismo pero con outSAMtype BAM SortedByCoordinate

mkdir mapped-quant2

cd mapped-quant2

R1=/dcdata2/josefina/rnaseq/trimmed_fastq/B51-1_S42_L003_R1_001.trim.paired.fq.gz
R2=/dcdata2/josefina/rnaseq/trimmed_fastq/B51-1_S42_L003_R2_001.trim.paired.fq.gz

STAR --outSAMstrandField intronMotif \
--outFilterIntronMotifs RemoveNoncanonical \
--runThreadN 28 \
--genomeDir ../starIndex \
--sjdbGTFfile ../../../51.geneswithphase.gtf \
--outSAMtype BAM SortedByCoordinate \
--readFilesCommand zcat $R1 $R2 \
--outFilterScoreMinOverLread 0 \
--outFilterMatchNminOverLread 0 \
--outFilterMatchNmin 0 \
--quantMode TranscriptomeSAM GeneCounts

#deje corriendo 11 jul 22.30 hs


mkdir mapped2

cd mapped2


for file in ../../trimmed_fastq/*R1*trim.paired.fq.gz
do STAR --outSAMstrandField intronMotif \
--outFilterIntronMotifs RemoveNoncanonical \
--runThreadN 28 \
--genomeDir ../starIndex \
--sjdbGTFfile ../../../51.geneswithphase.gtf \
--outSAMtype BAM Unsorted \
--readFilesCommand zcat $file ${file/R1/R2} \
--outFilterScoreMinOverLread 0 \
--outFilterMatchNminOverLread 0 \
--outFilterMatchNmin 0 \
--outFileNamePrefix ./${file%R1*}
done













