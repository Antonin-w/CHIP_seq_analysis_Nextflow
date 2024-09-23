# <a name="rmd-show-all-code"></a><a name="rmd-hide-all-code"></a>UE ADG - TP CHIP-seq - Année 2021-2022
#### **Antonin Weber**
#### **31 March, 2022**

# Echantillons et séquençage

## Les échantillons
Les échantillons biologiques sont des issus de trois ASC (ASC28, ASC42 et ASC63) sur lesquelles différentes conditions ont été testées :

- H3K27ac-CTR : Marque histone H3K27ac contrôle
- H3K27ac-TLT : Marque histone H3K27ac traité TLT
- H3K27me3-CTR : Marque histone H3K27me3 contrôle
- H3K27me3-TLT : Marque histone H3K27me3 traité TLT
- INPUT-CTR : Input contrôle
- INPUT-TLT : Input traité TLT

## Annotations des échantillons

```r
library(knitr)
annot <- read.table("~/Desktop/S2/GN3/Ressources/Annotation_samples_Chipseq.txt", header = TRUE, sep="\t")
kable(annot, caption = "Table 1 : Tableau d'annotation des échantillons")
```

<details>
  <summary>
    Voir le tableau
  </summary>
<table><tr><th colspan="7">Table 1 : Tableau d’annotation des échantillons</th></tr>
<tr><td><b>Samples</b></td><td><b>Samble.renamed</b></td><td><b>Conditions</b></td><td><b>Ctrl_or_TLT</b></td><td><b>Type</b></td><td><b>Marque.histone</b></td><td><b>Cond_name</b></td></tr>
<tr><td>A1140C01</td><td>C01</td><td>H3K27me3-CTR</td><td>CTR</td><td>ASC28</td><td>Broad</td><td>H3K27me3</td></tr>
<tr><td>A1141C20</td><td>C20</td><td>H3K27me3-TLT</td><td>TLT</td><td>ASC28</td><td>Broad</td><td>H3K27me3</td></tr>
<tr><td>A1140C05</td><td>C05</td><td>H3K4me1-CTR</td><td>CTR</td><td>ASC28</td><td>Broad</td><td>H3K4me1</td></tr>
<tr><td>A1140C08</td><td>C08</td><td>INPUT-TLT</td><td>TLT</td><td>ASC28</td><td>Broad</td><td>INPUT</td></tr>
<tr><td>A1141C19</td><td>C19</td><td>INPUT-CTR</td><td>CTR</td><td>ASC28</td><td>Broad</td><td>INPUT</td></tr>
<tr><td>A1141C24</td><td>C24</td><td>H3K4me1-TLT</td><td>TLT</td><td>ASC28</td><td>Broad</td><td>H3K4me1</td></tr>
<tr><td>A1140C09</td><td>C09</td><td>H3K27ac-CTR</td><td>CTR</td><td>ASC28</td><td>Narrow</td><td>H3K27ac</td></tr>
<tr><td>A1141C16</td><td>C16</td><td>H3K27ac-TLT</td><td>TLT</td><td>ASC28</td><td>Narrow</td><td>H3K27ac</td></tr>
<tr><td>A1140C05</td><td>C05</td><td>H3K4me1-CTR</td><td>CTR</td><td>ASC28</td><td>Narrow</td><td>H3K4me1</td></tr>
<tr><td>A1140C08</td><td>C08</td><td>INPUT-TLT</td><td>TLT</td><td>ASC28</td><td>Narrow</td><td>INPUT</td></tr>
<tr><td>A1141C19</td><td>C19</td><td>INPUT-CTR</td><td>CTR</td><td>ASC28</td><td>Narrow</td><td>INPUT</td></tr>
<tr><td>A1141C24</td><td>C24</td><td>H3K4me1-TLT</td><td>TLT</td><td>ASC28</td><td>Narrow</td><td>H3K4me1</td></tr>
<tr><td>A1140C07</td><td>C07</td><td>H3K27me3-CTR</td><td>CTR</td><td>ASC42</td><td>Broad</td><td>H3K27me3</td></tr>
<tr><td>A1141C14</td><td>C14</td><td>H3K27me3-TLT</td><td>TLT</td><td>ASC42</td><td>Broad</td><td>H3K27me3</td></tr>
<tr><td>A1141C15</td><td>C15</td><td>INPUT-CTR</td><td>CTR</td><td>ASC42</td><td>Broad</td><td>INPUT</td></tr>
<tr><td>A1141C17</td><td>C17</td><td>INPUT-TLT</td><td>TLT</td><td>ASC42</td><td>Broad</td><td>INPUT</td></tr>
<tr><td>A1141C18</td><td>C18</td><td>H3K4me1-TLT</td><td>TLT</td><td>ASC42</td><td>Broad</td><td>H3K4me1</td></tr>
<tr><td>A1141C23</td><td>C23</td><td>H3K4me1-CTR</td><td>CTR</td><td>ASC42</td><td>Broad</td><td>H3K4me1</td></tr>
<tr><td>A1140C03</td><td>C03</td><td>H3K27ac-CTR</td><td>CTR</td><td>ASC42</td><td>Narrow</td><td>H3K27ac</td></tr>
<tr><td>A1140C10</td><td>C10</td><td>H3K27ac-TLT</td><td>TLT</td><td>ASC42</td><td>Narrow</td><td>H3K27ac</td></tr>
<tr><td>A1141C15</td><td>C15</td><td>INPUT-CTR</td><td>CTR</td><td>ASC42</td><td>Narrow</td><td>INPUT</td></tr>
<tr><td>A1141C17</td><td>C17</td><td>INPUT-TLT</td><td>TLT</td><td>ASC42</td><td>Narrow</td><td>INPUT</td></tr>
<tr><td>A1141C18</td><td>C18</td><td>H3K4me1-TLT</td><td>TLT</td><td>ASC42</td><td>Narrow</td><td>H3K4me1</td></tr>
<tr><td>A1141C23</td><td>C23</td><td>H3K4me1-CTR</td><td>CTR</td><td>ASC42</td><td>Narrow</td><td>H3K4me1</td></tr>
<tr><td>A1140C02</td><td>C02</td><td>H3K27me3-TLT</td><td>TLT</td><td>ASC63</td><td>Broad</td><td>H3K27me3</td></tr>
<tr><td>A1141C13</td><td>C13</td><td>H3K27me3-CTR</td><td>CTR</td><td>ASC63</td><td>Broad</td><td>H3K27me3</td></tr>
<tr><td>A1140C06</td><td>C06</td><td>INPUT-TLT</td><td>TLT</td><td>ASC63</td><td>Broad</td><td>INPUT</td></tr>
<tr><td>A1140C11</td><td>C11</td><td>H3K4me1-CTR</td><td>CTR</td><td>ASC63</td><td>Broad</td><td>H3K4me1</td></tr>
<tr><td>A1140C12</td><td>C12</td><td>H3K4me1-TLT</td><td>TLT</td><td>ASC63</td><td>Broad</td><td>H3K4me1</td></tr>
<tr><td>A1141C22</td><td>C22</td><td>INPUT-CTR</td><td>CTR</td><td>ASC63</td><td>Broad</td><td>INPUT</td></tr>
<tr><td>A1140C04</td><td>C04</td><td>H3K27ac-TLT</td><td>TLT</td><td>ASC63</td><td>Narrow</td><td>H3K27ac</td></tr>
<tr><td>A1141C21</td><td>C21</td><td>H3K27ac-CTR</td><td>CTR</td><td>ASC63</td><td>Narrow</td><td>H3K27ac</td></tr>
<tr><td>A1140C06</td><td>C06</td><td>INPUT-TLT</td><td>TLT</td><td>ASC63</td><td>Narrow</td><td>INPUT</td></tr>
<tr><td>A1140C11</td><td>C11</td><td>H3K4me1-CTR</td><td>CTR</td><td>ASC63</td><td>Narrow</td><td>H3K4me1</td></tr>
<tr><td>A1140C12</td><td>C12</td><td>H3K4me1-TLT</td><td>TLT</td><td>ASC63</td><td>Narrow</td><td>H3K4me1</td></tr>
<tr><td>A1141C22</td><td>C22</td><td>INPUT-CTR</td><td>CTR</td><td>ASC63</td><td>Narrow</td><td>INPUT</td></tr>
</table>
</details>


<br>

# Expérience CHiP-seq
## Séquençage
Les échantillons ont été séquencés sur la plateforme de séquençage de l’Institut Curie à l’aide d’un Illumina HiSeq -Rapid run en paired-end avec une longueur de reads de 100 pb.

# Pipeline
*pipeline\_netflow\_diagram*
pipeline_workflow.png

Cette analyse a été faite grâce à `Nextflow`. Le principal avantage de cet outil est sa reproductbilité, extrêmement importante en Science, ainsi que sa paralélistion permettant une optimisation de la rapidité des analyses.
Toutes les versions de logiciels utilisées pour chacune des étapes de cette analyse sont indiqués dans le fichier `nextflow.config` suivant, cela permet aux autres chercheurs de pouvoir réitérer l'analyse avec exactement les mêmes configurations. Conda à été utililé pour l'installation automatique des logiciels dans l'environnement de travail. 

```bash
/*
 * -------------------------------------------------
 *  Nextflow base config file
 * -------------------------------------------------
 */

// CHIP-seq Pipeline
// Automatic pipeline | CHIPseq Analysis | Illumina HiSeq- Rapid run | paired-end | reads : 100pb

fastqc = "0.11.9"
bowtie2 = "2.4.5"
samtools = "1.9"
picard = "2.18.14"
macs2 = "2.2.7.1"
bedtools = "2.30.0"

process {
  withName:fastqc{
    conda = "fastqc=${fastqc}"
  }
  withName:alignement_human {
  conda = "bowtie2=${bowtie2} samtools=${samtools}"
  cpus = 12
  memory = 80.GB
  }
  withName:remove_pcr_duplicate {
    conda = "picard=${picard} samtools=${samtools}"
    memory = 32.GB
  }
  withName:samtools_view_mito {
    conda = "samtools=${samtools}"
    memory = 32.GB
  }
  withName:macs2 {
    conda = "macs2=${macs2}"
  }
  withName:macs2_filt_bed {
    conda = "bedtools=${bedtools}"
  }
  withName:macs2_region {
    conda = "bedtools=${bedtools} samtools=${samtools}"
  }
}

params {
  max_memory = 250.GB
  max_cpus = 16
  max_time = 240.h
}
```

# Analyses bioinformatiques
## Paramètres utilisés lors de l’analyse

L'une des principales force de `Nextflow`, et plus généralement, d'une pipeline d'analyse, est qu'il est très facile de pouvoir relancer une analyse complète très rapidement et facilement en changeant simplement les options de départ. Par exemple, ré-utiliser cette pipeline sur d'autres jeux de données, ou encore sur les mêmes données de départ, mais en spécifiant certains arguments induisant un traitement différent de l'information, ou encore en mettant à jour les versions de logiciel.

Des arguments peuvent être rentrés lors de l'éxécution du script pour directement changer certains paramètres, notamment les fichiers d'entrés, de références etc...
Les paramètres par défaut sont spécifiés au début de la pipeline :

```bash
#!/usr/bin/env nextflow

/*
* SET UP CONFIGURATION VARIABLES
*/
params.executor = 'slurm'
params.outdir = 'results'

params.reads = "1_Brut/*R{1,2}.fastq.gz"
params.index_human = "Annotations/Bowtie2_indexes_Hg19/bt2/hg19"
params.index_droso = "Annotations/Bowtie2_indexes_dm6/dm6"
params.reads_filtered = "Bed_*_WO_Norm/*_filt.bed"
params.fileBlacklistInit = "Annotations/hg38.blacklist.bed"
params.chrNameLength = "Annotations/chrNameLength_chr.txt"
```

L'**executor** que nous utiliserons sera toujours le même pour tous les **process**, le Worfload Manager de `Genouest` : `Slurm`. 
Les fichiers sortant seront enregistrés dans un dossier `results`.
Les paramètres suivants indiquent le chemin menant aux fichiers nécéssaires à l'analyse. Si un utilisateur souhaite changer un de ces arguments, il suffit de les indiquer comme dans l'exemple ci-dessous:

```bash
nextflow run ChIP_pipeline.nf -reads "/path/*R{1,2}.fastq.gz" -executor 'local'
```

De plus, `Nextflow` dispose de l'option **-resume** qui permet de reprendre l'analyse ou elle s'est arrêter la dernière fois. Pratique si une erreur ou un problème est apparu après un certain nombre d'étape, pour éviter de tout recommencer depuis le début.

## Contrôle qualité sur les données brutes

Avant de réaliser l’alignement des fichiers bruts obtenus en sortie de séquençage, il est nécessaire d’effectuer un contrôle qualité de ces fichiers. En effet, les premières paires de bases d’un read sont séquencées avec beaucoup de fiabilité, mais plus nous avancons dans la séquence, moins c’est précis. Un programme comme FastQC produira divers graphiques permettant de savoir s’il faut nettoyer les reads, et de quelle longueur le “trimming” doit alors se faire. Cette étape est trés importante puisque des reads de mauvaise qualité s’aligneront difficilement sur le génome, et rendront toute l’analyse inutile.


```bash
Channel
    .fromFilePairs( params.reads )
    .set { raw_reads_fastq }
```

Nextflow utilise principalement un système de **Channel** afin de gérer les fichiers d’entrées et de sortie. Cela va permettre d’indiquer directement d’un **process** à l’autre quels fichiers sortant doivent servir de fichier d’entrées lors de la prochaine étape. Le chemin menant aux fichiers à analyser est présent dans le paramètre **param.reads**. Comme enoncé précedemment, nos données sont paired-end, nous allons donc utiliser une fonction que possède Nextflow afin de grouper par paires les fichiers de départ. Ils sont ensuite placés dans **raw\_reads\_fastq**.

```bash
process fastqc {
    executor params.executor
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'symlink'

    input:
    set val(name), file(reads) from raw_reads_fastqc

    output:
    file "*.zip"
    file "*.html"

    script:
    """
    fastqc -q $reads
    """
}
```

Le **process fastqc** prendra en input tous les fichiers présents dans le channel `raw_reads_fastq`. Les fichiers d'outputs seront des *.zip* et des *.html* permettant d'évaluer la qualité de nos échantillons. Tous les fichiers générés par `Nextflow` sont présent dans un répertoire créer spécialement : **work**. Cependant, si nous souhaitons consutler ces fichiers, il est plutôt compliqué de s'y retrouver dedans. Nous allons donc ajouter un lien symbolique pour tous les fichiers que nous souhaitons avoir pour plus tard dans un dossier grâce à : *publishDir "${params.outdir}/fastqc", mode: 'symlink'*.

Les fichiers de sortie du *fastqc* sont les suivants :  

* *.fastqc.html  
* *.fastqc.zip

Afin de visualiser l’ensemble des résultats du contrôle qualité des fastq sur un seul graphique, il faut appliquer le script suivant, qui va générer un tableau récapitulatif contenant le nombre de reads séquencés, le pourcentage de reads dupliqués, le pourcentage de GC et la longueur des reads pour les séquences R1 et R2 de chaque échantillon.

```r
library(fastqcr)

setwd("2_fastQC")
qc.dir <- getwd()

list.files(qc.dir)
qc_agr<-qc_aggregate(qc.dir , progressbar = F)
e<-summary(qc_agr)
d<-qc_stats(qc_agr)

write.table(d,"Result_fastQC.txt", sep="\t", row.names =T)
```

```r
fastqc <- read.table("Result_fastQC.txt", sep="\t", header = T)
library(tidyverse)

my_plot <- ggplot(fastqc, aes (x = sample, group=1, fill=tot.seq)) +
  geom_bar(aes(y = tot.seq/1000000), fill = "light blue", stat = "identity") +
  geom_line(aes(y = pct.gc, colour = "GC%")) +  
  geom_line(aes(y = pct.dup, colour = "Duplicats%")) +
  geom_line(aes(y = 101-35, colour = "Longueur des reads"))

my_plot <- my_plot +
  scale_color_manual(values = c('GC%' = 'red', 'Duplicats%' = 'green', 'Longueur des reads' = 'blue')) +
  labs(title = "Contrôle qualité des fastq",
       y = "Pourcentage de GC et des duplicats",
       x = "Echantillons",
       color = "Légende") +
  theme_classic() +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        axis.text.x = element_text(angle = 60, hjust = c(1,2), size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13, face = "bold"))

my_plot
```

![](Aspose.Words.4dd4bbe9-9463-44b2-8e50-5fef52002936.002.png)

Tout d'abord, nous pouvons observer que pour tous les échantillons, il y a le même nombre de reads entre les échantillons `R1` et `R2`. Cela montre à première vue qu'il n'y a pas eu d'erreurs flagrante de séquencage. De plus, la longueur des reads est similaire sur l'ensemble des samples. Le pourcentage de GC ne dérive pas vers des valeurs extremes, mais est bien contenu entre 40 et 60 pourcents, assez cohérent pour des échantillons provenant d'Humains. 
Pour terminer, le pourcentage de duplicats est compris entre 10 et 20%, cela est plutôt standard pour ce type de séquençage, il ne faudra cependant pas oublié de supprimer les duplicats de PCR après l'étape d'alignement. 

## Alignements sur génome de référence avec Bowtie 2 et nettoyage des fichiers bam

### Alignement via Bowtie2

Nous allons maintenant déterminer l’emplacement de nos reads sur le génome humain, pour cela nous utiliserons Bowtie2 sur le génome de référence GRCh38 release 90.

Justification des options utilisées pour `Bowtie2`:  

* `sensitive` : optimisation de l’alignement A DETAILLER 
* `p` : nombre maximum de mismatches autorisés: **20**  
* `1` : fichier R1 (fastq)
* `2` : fichier R2 (fastq)  
* `S` : fichier de sortie (SAM)

Justification des options utilisées pour `samtools`: 

* `m` : mémoire maximale requise exprimée en Giga 
* `bhS` : transforme le fichier sam en fichier bam

```bash
Channel.fromFilePairs( params.reads ).set { read_pairs_ch }

process alignement_human {
    executor params.executor

    input:
    set sampleId, file(reads) from read_pairs_ch

    output:
    set sampleId, file("${sampleId}_sorted.bam") into sorted_bam_human

    script:
    """
    bowtie2 --sensitive -p 20 -x ${params.index_human} -1 ${reads[0]} -2 ${reads[1]} -S H${sampleId}_bowtie2.sam
    samtools view -m 12g -bhS -o C${sampleId}_samview.bam H${sampleId}_bowtie2.sam
    rm H${sampleId}_bowtie2.sam
    samtools sort -m 12g -o ${sampleId}_sorted.bam C${sampleId}_samview.bam
    rm C${sampleId}_samview.bam
    """
}
```

## Suppression des duplicats de PCR

Le processus de `Polymerase Chain Reaction PCR` induit inévitablement des duplicats qui doivent être supprimés afin de ne pas fausser nos analyses. En effet, si ces derniers sont laissés, alors les reads en question pourraient être

```bash
process remove_pcr_duplicate {
    executor params.executor

    input:
    set sampleId, file(reads) from sorted_bam_human

    output:
    set sampleId, file("${sampleId}_sorted_rmdup.bam") into bam_sorted_rmdup
    file "${sampleId}_Bowtie2_dups.txt"

    script:
    """
    picard -Xmx12g MarkDuplicates \\
        INPUT=${reads} \\
        OUTPUT=${sampleId}_sorted_rmdup.bam \\
        METRICS_FILE=${sampleId}_Bowtie2_dups.txt \\
        REMOVE_DUPLICATES=true \\
        CREATE_INDEX=true
    rm ${reads}
    """
}
```

## Suppression des chromosomes mitochondriaux

```bash
process samtools_view_mito {
    executor params.executor
    publishDir "$params.outdir/3_Align_Bowtie2/Human", mode: 'symlink'

    input:
    set sampleId, file(reads) from bam_sorted_rmdup

    output:
    file("${sampleId}_wo_chrM") into bam_wo_chrM

    script:
    """
    samtools index ${reads}
    samtools idxstats ${reads} | cut -f 1 | grep -v M | xargs samtools view -b ${reads} > ${sampleId}_wo_chrM
    rm ${reads}
    """
}
```

## Peak calling

Il est maintenant nécessaire de détecter les peaks grâce au calcul du chevauchement maximum des tags décalés selon un score minimum significatif défini. Pour cela, nous savons les types de marques pour chaque histone :

- H3K27me3 : **Broad**
- H3K27ac : **Narrow**
- H3K4me1 : **Broad et Narrow**

*Table d’annotation contenant le nom des échantillons, la marque d’histone, et l’échantillon contrôle*

```r
annot_Peak <- read.csv("Annot_Peak_Calling.csv")
kable(annot_Peak)
```

|**Name**|**Sample**|**Histone**|**Control**|
| :- | :- | :- | :- |
|A1140C01|C01|Broad|C19|
|A1140C02|C02|Broad|C06|
|A1140C05|C05|Broad|C19|
|A1140C06|C06|Broad|C06|
|A1140C07|C07|Broad|C15|
|A1140C08|C08|Broad|C08|
|A1140C11|C11|Broad|C22|
|A1140C12|C12|Broad|C06|
|A1141C13|C13|Broad|C22|
|A1141C14|C14|Broad|C17|
|A1141C15|C15|Broad|C15|
|A1141C17|C17|Broad|C17|
|A1141C18|C18|Broad|C17|
|A1141C19|C19|Broad|C19|
|A1141C20|C20|Broad|C08|
|A1141C22|C22|Broad|C22|
|A1141C23|C23|Broad|C15|
|A1141C24|C24|Broad|C08|
|A1140C03|C03|Narrow|C15|
|A1140C04|C04|Narrow|C06|
|A1140C05|C05|Narrow|C19|
|A1140C06|C06|Narrow|C06|
|A1140C08|C08|Narrow|C08|
|A1140C09|C09|Narrow|C19|
|A1140C10|C10|Narrow|C17|
|A1140C11|C11|Narrow|C22|
|A1140C12|C12|Narrow|C06|
|A1141C15|C15|Narrow|C15|
|A1141C16|C16|Narrow|C08|
|A1141C17|C17|Narrow|C17|
|A1141C18|C18|Narrow|C17|
|A1141C19|C19|Narrow|C19|
|A1141C21|C21|Narrow|C22|
|A1141C22|C22|Narrow|C22|
|A1141C23|C23|Narrow|C15|
|A1141C24|C24|Narrow|C08|

Afin d’automatiser l’étape de Peak Calling, un script bash à été incorporé au process macs2, et permet à partir d’un tableau contenant 4 colonnes : Le nom de l’échantillon, la marque d’histone, le chemin menant à l’échantillon ainsi que le chemin menant à l’échantillon contrôle de réaliser le Peak Calling.

Paramètres `macs2` utilisé : 

* `t` : fichier d'entré 
* `n` : préfixe pour les fichiers de sortie
* `c` : fichier contrôle 
* `f BAMPE` : format du fichier d'entré : BAM Paired-End
* `g` : taille du génome alignable 
* `q` : p-value limite accepté 
* `--nomodel` : Construction d'un modèle de décalage
* `--outdir` : répertoire ou les fichiers de sortie seront enregistrés 
* `--broad` : broad peak calling


```bash
process macs2 {
    executor params.executor
    publishDir "$params.outdir/4_Peak_Calling", mode: 'symlink'

    input:
    file "*" from bam_wo_chrM.toList()

    output:
    file "*.broadPeak" into broad_peak 
    file "*.narrowPeak" into narrow_peak 
    file "*peaks.xls"
    file "*peaks.gappedPeak"
    
    script:
    """
    {
    read line ; IFS=,; while read -r Name Marque Sample_File Control_File; 
    do
    if [[ \$Marque = "Narrow" ]]
    then
        if [[ \$Control_File == *"None"* ]]
        then
            macs2 callpeak -t \$Sample_File -n Narrow\$Name -f BAMPE -g 2.7e9 -q 0.05 --nomodel --outdir 4_Peak_Calling/Narrow
        else
            macs2 callpeak -t \$Sample_File -c \$Control_File -n Narrow\$Name -f BAMPE -g 2.7e9 -q 0.05 --nomodel --outdir /4_Peak_Calling/Narrow
        fi
    else 
        if [[ \$Control_File == *"None"* ]]
        then
            macs2 callpeak -t \$Sample_File -n Broad\$Name -f BAMPE -g 2.7e9 --broad -q 0.05 --nomodel --outdir 4_Peak_Calling/Broad
        else
            macs2 callpeak -t \$Sample_File -c \$Control_File -n Broad\$Name -f BAMPE -g 2.7e9 --broad -q 0.05 --nomodel --outdir 4_Peak_Calling/Broad
        fi
    fi
    done
    } <Annotation_Peak.csv
    """
}
```

Cette étape va produire les fichiers de sortie suivants :

- `{Broad|Narrow}{NAME}\_peaks.broadPeak` : liste de l’ensemble des peaks indentifiés
- `{Broad|Narrow}{NAME}\_peaks.xls` : liste des peaks identifiés pour chaque sample
- `{Broad|Narrow}{NAME}\_peaks.gappedPeak` : représentation de pics étroits sous forme de blocs sur un pic large. Pour tromper les outils de visualisation, ils utilisent le même format que les modèles de gènes, mais utilisent les coordonnées de pic étroites comme coordonnées d’exons et les coordonnées de pic larges comme coordonnées de région codante.

## Génération des fichiers BED

Une fois le peak calling réalisé, il faut générer les fichiers .bed pour chaque échantillons mais aussi un fichier bed qui regroupe l’ensemble des peaks identifiés pour tous les échantillons afin de créer l’univers de référence pour la suite de l’analyse.

En entré, les fichiers filtrés avec une q-value de 10^3. Vérification de l’ordre des chromosomes entre blacklist et chrSize du génome de référence.

```bash
params.reads_filtered = "Bed_*_WO_Norm/*_filt.bed"
params.fileBlacklistInit = "hg38.blacklist.bed"
params.chrNameLength = "chrNameLength_bam.txt"

Channel. fromPath(params.reads_filtered). set { bed_filtered }

process macs2_filt_bed {
    executor params.executor
    publishDir "$params.outdir/6_filted_bed", mode: 'symlink'

    input:
    file(bed_filtered)

    output:
    file "*ChIP_Region.bed"
    file "*ChIP_Region_sorted.bed"
    file "*ChIP.bed"
    file "*ChIP_sorted.bed" into sorted_bed

    script:
    """
    cat ${bed_filtered} | sort -k1V -k2,2n | bedtools intersect -a - -b ${params.fileBlacklistInit} -v | bedtools merge -i - | awk '{OFS="\t"; print(\$0, "Region_"NR)}' >${bed_filtered.baseName}ChIP_Region.bed
    bedtools sort -g ${params.chrNameLength} -i ${bed_filtered.baseName}ChIP_Region.bed >${bed_filtered.baseName}ChIP_Region_sorted.bed
    
    cat ${bed_filtered} | sort -k1V -k2,2n | bedtools intersect -a - -b ${params.fileBlacklistInit} -v | bedtools merge -i - -c 4 -o distinct >${bed_filtered.baseName}ChIP.bed
    bedtools sort -g ${params.chrNameLength} -i ${bed_filtered.baseName}ChIP.bed >${bed_filtered.baseName}ChIP_sorted.bed
    """
}
```

Les fichiers d’ouputs sont les suivants :

- \***filtChIP.bed** : fichier bed correspond au merge de tous les fichiers xx\_filt.bed. Dans ce fichier les peaks sont nommés de la façon suivante : sample\_peak\_1
- \***filtChIP\_sorted.bed** : fichier ChIP\_Bowtie2.bed trié sur les chromosomes afin que l’ordre corresponde à l’ordre des fichiers .bam
- \***filtChIP\_Region.bed** : fichier ChIP\_Bowtie2.bed pour lequel les peaks sont nommés de Region\_1 à Region\_N
- \***filtChIP\_Region\_sorted.bed** : fichier ChIP\_Region\_Bowtie2.bed trié sur les chromosomes afin que l’ordre corresponde à l’ordre des fichiers .bam


## Contrôle qualité des fichiers bed

En utilisant une fenêtre de +/- 8000 pb autour du TSS (site de démarrage de la transcription), nous pouvons plotter les fréquences de counts pour voir où les reads se fixent par rapport au TSS. Pour ce projet, les échantillons sont plottés 3 par 3 en fonction des conditions (H3K27me3-CTR, H3K27me3-TLT, H3K27ac-CTR, H3K27ac-TLT, H3K4me1-CTR et H3K4me1-TLT).

```r
#! /usr/bin/env Rscript

require(org.Hs.eg.db)
require(ChIPseeker)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(clusterProfiler)
require(rtracklayer)

setwd("6_filted_bed/")
dir <- getwd()
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

file_H3K27me3_CTR <- list(H01="01_filtChIP.bed", H07="07_filtChIP.bed", H13="13_filtChIP.bed")
file_H3K27me3_TLT <- list(H02="02_filtChIP.bed", H14="14_filtChIP.bed", H20="20_filtChIP.bed")

file_H3K27ac_CTR <- list(H03="03_filtChIP.bed", H09="09_filtChIP.bed", H21="21_filtChIP.bed")
file_H3K27ac_TLT <- list(H04="04_filtChIP.bed", H10="10_filtChIP.bed", H16="16_filtChIP.bed")

file_H3K4me1_CTR <- list(H05="05_filtChIP.bed", H11="11_filtChIP.bed", H23="23_filtChIP.bed")
file_H3K4me1_TLT <- list(H12="12_filtChIP.bed", H18="18_filtChIP.bed", H24="24_filtChIP.bed")


samplefiles_H3K27me3_CTR <- as.list(file_H3K27me3_CTR)
names(samplefiles_H3K27me3_CTR) <- c("H01", "H07", "H13")
samplefiles_H3K27me3_TLT <- as.list(file_H3K27me3_TLT)
names(samplefiles_H3K27me3_TLT) <- c("H02", "H14", "H20")

samplefiles_H3K27ac_CTR <- as.list(file_H3K27ac_CTR)
names(samplefiles_H3K27ac_CTR) <- c("H03", "H09", "H21")
samplefiles_H3K27ac_TLT <- as.list(file_H3K27ac_TLT)
names(samplefiles_H3K27ac_TLT) <- c("H04", "H10", "H16")

samplefiles_H3K4me1_CTR <- as.list(file_H3K4me1_CTR)
names(samplefiles_H3K4me1_CTR) <- c("H05", "H11","H23")
samplefiles_H3K4me1_TLT <- as.list(file_H3K4me1_TLT)
names(samplefiles_H3K4me1_TLT) <- c("H12", "H18","H24")

for (ech in names(samplefiles_H3K27me3_CTR)) {
	name <- ech
	assign(name, readPeakFile(file_H3K27me3_CTR[[ech]]))
}
for (ech in names(samplefiles_H3K27me3_TLT)) {
	name <- ech
	assign(name, readPeakFile(file_H3K27me3_TLT[[ech]]))
}

for (ech in names(samplefiles_H3K27ac_CTR)) {
	name <- ech
	assign(name, readPeakFile(file_H3K27ac_CTR[[ech]]))
}
for (ech in names(samplefiles_H3K27ac_TLT)) {
	name <- ech
	assign(name, readPeakFile(file_H3K27ac_TLT[[ech]]))
}

for (ech in names(samplefiles_H3K4me1_CTR)) {
	name <- ech
	assign(name, readPeakFile(file_H3K4me1_CTR[[ech]]))
}
for (ech in names(samplefiles_H3K4me1_TLT)) {
	name <- ech
	assign(name, readPeakFile(file_H3K4me1_TLT[[ech]]))
}


# Prépartion des régions des promoteurs
promoter <- getPromoters(TxDb=txdb, upstream=8000, downstream=8000)

# Calcul de la matrice de tags
tagMatrixList_H3K27me3_CTR <- lapply(as.list(samplefiles_H3K27me3_CTR), getTagMatrix, windows=promoter)
pdf(file = "H3K27me3_CTR.pdf")
tagHeatmap(tagMatrixList_H3K27me3_CTR, xlim=c(-8000, 8000), color=NULL)
dev.off()

tagMatrixList_H3K27me3_TLT <- lapply(as.list(samplefiles_H3K27me3_TLT), getTagMatrix, windows=promoter)
pdf(file = "H3K27me3_TLT.pdf")
tagHeatmap(tagMatrixList_H3K27me3_TLT, xlim=c(-8000, 8000), color=NULL)
dev.off()

tagMatrixList_H3K27ac_CTR <- lapply(as.list(samplefiles_H3K27ac_CTR), getTagMatrix, windows=promoter)
pdf(file = "H3K27ac_CTR.pdf")
tagHeatmap(tagMatrixList_H3K27ac_CTR, xlim=c(-8000, 8000), color=NULL)
dev.off()

tagMatrixList_H3K27ac_TLT <- lapply(as.list(samplefiles_H3K27ac_TLT), getTagMatrix, windows=promoter)
pdf(file = "H3K27ac_TLT.pdf")
tagHeatmap(tagMatrixList_H3K27ac_TLT, xlim=c(-8000, 8000), color=NULL)
dev.off()

tagMatrixList_H3K4me1_CTR <- lapply(as.list(samplefiles_H3K4me1_CTR), getTagMatrix, windows=promoter)
pdf(file = "H3K4me1_CTR.pdf")
tagHeatmap(tagMatrixList_H3K4me1_CTR, xlim=c(-8000, 8000), color=NULL)
dev.off()

tagMatrixList_H3K4me1_TLT <- lapply(as.list(samplefiles_H3K4me1_TLT), getTagMatrix, windows=promoter)
pdf(file = "H3K4me1_TLT.pdf")
tagHeatmap(tagMatrixList_H3K4me1_TLT, xlim=c(-8000, 8000), color=NULL)
dev.off()
```

## Annotation des pics

L’annotation des pics est effectuée à l’aide la fonction annotatePeak. Il est possible de définir la région TSS (site de démarrage de la transcription). Par défaut, TSS est défini entre -8 Ko et + 8 Ko. le fichier de sortie de cette fonction peut être exporté au format txt.

L’objet TxDb contient les caractéristiques liées à la transcription d’un génome, il s’agit ici du génome de référence humain Hg38.
Toutes les informations sur les pics contenues dans peakfile seront conservées dans la sortie de annotatePeak. Les informations sur la position et le brin des gènes les plus proches sont rapportées. La distance du pic au TSS de son gène le plus proche est également rapportée. La région génomique du pic est indiquée dans la colonne d’annotation. Certaines annotations pouvant se chevaucher, ChIPseeker a adopté la priorité suivante en annotation génomique :

* Promoter
* 5’ UTR
* 3’ UTR
* Exon
* Intron
* Downstream
* Intergenic

**Les tableaux ci-dessous reprennent les fréquences de chaque critère d’annotation pour chaque condition**

Feature_H3K3.jpg
Feature_H3K4.jpg

**Distribution des features d’annotation autour du TSS**

merge_from_ofoct-2.jpg

## MACS2 Region Count

Le script MACS2_Region.sh permet de compter le nombre de reads à partir des fichiers bam nettoyés de leurs chromosomes mitochondriaux et de leurs duplicats de PCR et ce par région identifiées à l’étape précédente.

```bash

Channel.fromPath("Bowtie2_Indexes/Hg38_R90/chrNameLength.txt"). set { fileAnnot }
process macs2_region {
    executor params.executor
    publishDir "$params.outdir/7_region", mode: 'symlink'

    input:
    file(region_sorted_bed)
    file(bam_wo_chrM)
    file(fileAnnot)

    output:
    file "*"

    script:
    """
    # Comptage des reads associes a chaque region
    bedtools intersect -a ${region_sorted_bed} -b ${bam_wo_chrM} -c -sorted -g ${fileAnnot} > ${region_sorted_bed}_count.txt
    
    # calcul du count normalise = count brut / taille de la region / nb reads total du .bam
    awk -v '{OFS="/t"; print($0, ($5/($3-$2))/($(samtools view -c ${region_sorted_bed}_count.txt)/1e6))}'${region_sorted_bed}_count.txt >${region_sorted_bed}_count2.txt

    awk '{print ($4, $5)}' ${region_sorted_bed}_count2.txt > ${region_sorted_bed}_DESeq2.txt
    """
}
```

# Conclusion

Pour conclure, ce projet à permis une initiation à l''outil `Nextflow`, son utilisation est largement en train de se démocratiser pour les raisons qui ont pu être évoquée précedemment : Sa reproductibilité, Nextflow permet de renseigner de manière précise l'environnement de travail et la version de tous les logiciels. La paralélisation est aussi un point important, en effet, plutôt que d'exécuter un programme sur un échantillon à la fois. Nextflow va lancer un processus sur tous les échantillons en même temps.

Par exemple, si une analyse tel que FastQC prend 10 minutes par sample, avec 24 samples:

* `Script bash classique avec boucle` : 24 samples * 10 minutes = 240 minutes = 4 heures au total
* `Nextflow` : Tous les échantillons font l'étape en même temps = 10 minutes au total 

Cela nécéssite par contre un serveur assez puissant pour faire tourner en même temps le tout, en l'occurence, `Genouest` à la capacité de le faire.`sbatch` gère également automatiquement la paralélisation, donc cela reste très simple à implémenter. 

De plus, plutôt que de faire étapes par étapes, Nextflow va automatiquement créer des "liens" entre les process pour savoir l'ordre d'éxécution. Si par exemple il voit qu'il ne faut pas attendre d'obtenir les fichiers de sortie **fastqc** pour faire l'alignement, alors il exécutera **fastqc** et **alignement** simultanément, augmentant de nouveau la rapidité de l'analyse. 

Pour finir, la syntaxe de `Nextflow` est extremement clair et lisible pour déterminer aisément les différentes étapes d'une analyse. Une piste d'amélioration aurait été d'utiliser la syntaxe DSL-2. Cela permet d'avoir une liste de fichier contenant chacun un **process**, puis un fichier principale contenant le workflow, comme l'exemple ci-bas :


```bash
workflow {
  reads = Channel.fromFilePairs(params.reads) // Enregistrement des fichiers fastqc dans le Channel "reads"
  fastqc(reads) // Exécution du process fastqc sur ces reads
  trimming(reads) // Exécution du process trimming sur ces reads 
  ...
}
```

Afin d'améliorer cette pipeline, nous pourrions aussi utiliser l'outil `docker` (https://www.docker.com), permettant la création de `containers` retranscrivant l'environnement d'exécution de la pipeline de la manière la plus précise possible. (https://datascientest.com/docker-guide-complet)
Je n'ai malheureusement pas réussis à automatiser l'éxécution des `Scripts R`, mais cela est normalement possible au sein des process, sûrement plus facile à l'aide de `docker`.

N'ayant que peu d'information sur l'aspect biologique de l'experience, j'ai fait le choix de ne pas pousser les interpretations jusqu'a l'extrapolation, et plutôt profiter de ce projet pour décourvir et présenter une manière plus reproductible et rapide de réaliser une analyse ChIP-seq.

