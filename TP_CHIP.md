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

![pipeline_workflow.png](https://github.com/Antonin-w/CHIP_seq_analysis_Nextflow/blob/fd2b3daa3627f3a21c1c455fe1bee3a2acadc2c9/Images/pipeline_workflow.png?raw=true)


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

![qc_quality](https://github.com/Antonin-w/CHIP_seq_analysis_Nextflow/blob/fd2b3daa3627f3a21c1c455fe1bee3a2acadc2c9/Images/qc_fastq?raw=true)

Tout d'abord, nous pouvons observer que pour tous les échantillons, il y a le même nombre de reads entre les échantillons `R1` et `R2`. Cela montre à première vue qu'il n'y a pas eu d'erreurs flagrante de séquencage. De plus, la longueur des reads est similaire sur l'ensemble des samples. Le pourcentage de GC ne dérive pas vers des valeurs extremes, mais est bien contenu entre 40 et 60 pourcents, assez cohérent pour des échantillons provenant d'Humains. 
Pour terminer, le pourcentage de duplicats est compris entre 10 et 20%, cela est plutôt standard pour ce type de séquençage, il ne faudra cependant pas oublié de supprimer les duplicats de PCR après l'étape d'alignement. 

## Alignements sur génome de référence avec Bowtie 2 et nettoyage des fichiers bam

### Alignement via Bowtie2

Nous allons maintenant déterminer l’emplacement de nos reads sur le génome humain, pour cela nous utiliserons Bowtie2 sur le génome de référence GRCh38 release 90.

Justification des options utilisées pour `Bowtie2`:  

* `sensitive` : optimisation de l’alignement pré-reglage de Bowtie2 équivalent aux options suivantes : -D 15 -R 2 -N 0 -L 22 -i S,1,1.15
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

<details>
  <summary>
    Voir le tableau d'annotation
  </summary>

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
</details>

<br>

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

**H3K27me3 CTR Broad**
Des peaks sont facilement observables, pas de grande différence entre données normalisé ou non

![H3K27me3 CTR Broad](https://github.com/Antonin-w/CHIP_seq_analysis_Nextflow/blob/fd2b3daa3627f3a21c1c455fe1bee3a2acadc2c9/Images/HeatMap_H3K27me3_CTR.png?raw=true)

**H3K27me3 TLT Broad**
Même observation que pour le CTR, pas de différence entre les deux conditions

![H3K27me3 TLT Broad](https://github.com/Antonin-w/CHIP_seq_analysis_Nextflow/blob/fd2b3daa3627f3a21c1c455fe1bee3a2acadc2c9/Images/HeatMap_H3K27_TLT.png?raw=true)

**H3K27ac CTR Broad**
Peaks facilement observables

![H3K27ac CTR Broad](https://github.com/Antonin-w/CHIP_seq_analysis_Nextflow/blob/fd2b3daa3627f3a21c1c455fe1bee3a2acadc2c9/Images/HeatMapH3K27ac_CTR.png?raw=true)

**H3K27ac TLT Broad**
Pas de Peaks observables, différence avec la condition contrôle
![H3K27ac TLT Broad](https://github.com/Antonin-w/CHIP_seq_analysis_Nextflow/blob/fd2b3daa3627f3a21c1c455fe1bee3a2acadc2c9/Images/HeatMap_H3K27ac_TLT.png?raw=true)

**H3K4me1 CTR Broad**
Pas de Peaks observables
Pas de Peaks observables, différence avec la condition contrôle
![H3K4me1 CTR Broad](https://github.com/Antonin-w/CHIP_seq_analysis_Nextflow/blob/fd2b3daa3627f3a21c1c455fe1bee3a2acadc2c9/Images/HeatMap_H3K4me1CTR_B.png?raw=true)

**H3K4me1 TLT Broad**
Pas de Peaks observables non plus, comme dans la condition contrôle
![H3K4me1 TLT Broad](https://github.com/Antonin-w/CHIP_seq_analysis_Nextflow/blob/fd2b3daa3627f3a21c1c455fe1bee3a2acadc2c9/Images/HeatMap_H3K4me1TLT_B.png?raw=true)

**H3K4me1 TLT Narrow**
Peaks similaire à la condition contrôle de la marque d'histone Narrow
![H3K4me1 TLT Narrow](https://github.com/Antonin-w/CHIP_seq_analysis_Nextflow/blob/fd2b3daa3627f3a21c1c455fe1bee3a2acadc2c9/Images/HeatMap_H3K4me1_Narrow_TLT.png?raw=true)


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


<details>
  <summary>
    Voir le tableau H3K27me3
  </summary>


|**Sample**|**A1140C01**|**A1140C07**|**A1141C13**|**A1140C01\_Norm**|**A1140C07\_Norm**|**A1141C13\_Norm**|
| :-: | :-: | :-: | :-: | :-: | :-: | :-: |
|Nb Annoted Peaks|4629\.0000000|2825\.0000000|698\.0000000|3534\.0000000|1848\.0000000|698\.0000000|
|Promoter (<=1kb)|25\.2754375|25\.0265487|22\.0630372|25\.2122241|25\.9740260|22\.0630372|
|Promoter (1-2kb)|4\.0829553|4\.3539823|5\.3008596|3\.9898132|4\.0043290|5\.3008596|
|Promoter (2-3kb)|3\.4996760|3\.8230088|4\.4412607|3\.5370685|3\.7878788|4\.4412607|
|Promoter (3-4kb)|3\.1324260|4\.0353982|3\.5816619|2\.8862479|3\.8419913|3\.5816619|
|Promoter (4-5kb)|2\.8299849|2\.5840708|3\.0085960|2\.5466893|2\.6515152|3\.0085960|
|Promoter (5-6kb)|2\.8299849|2\.7256637|2\.5787966|3\.1126203|3\.0844156|2\.5787966|
|Promoter (6-7kb)|2\.2899114|2\.6902655|2\.4355301|2\.2354273|2\.8138528|2\.4355301|
|Promoter (7-8kb)|2\.1170879|1\.6991150|2\.1489971|2\.4335031|1\.6233766|2\.1489971|
|Promoter (8-9kb)|1\.9442644|1\.6991150|2\.0057307|1\.7260894|1\.4069264|2\.0057307|
|Promoter (9-10kb)|1\.7498380|1\.5575221|2\.2922636|1\.6411998|1\.2987013|2\.2922636|
|X5UTR|0\.5184705|0\.2477876|0\.1432665|0\.5093379|0\.2705628|0\.1432665|
|X3UTR|2\.1386909|2\.1238938|1\.7191977|2\.0373514|1\.9480519|1\.7191977|
|X1st Exon|3\.7157053|3\.1150442|2\.8653295|3\.0277306|2\.6515152|2\.8653295|
|Other Exon|4\.4718082|4\.1415929|2\.8653295|4\.0181098|3\.6255411|2\.8653295|
|X1st Intron|5\.7031756|6\.1238938|7\.0200573|5\.8856819|6\.8181818|7\.0200573|
|Other Intron|11\.5575718|11\.3274336|11\.3180516|11\.6015846|11\.5800866|11\.3180516|
|Downnstream (<=300)|0\.1728235|0\.3893805|NA|0\.2829655|0\.3246753|NA|
|Distal Intergenic|21\.9701879|22\.3362832|24\.2120344|23\.3163554|22\.2943723|24\.2120344|

|||
| :- | :- |

|**Sample**|**A1141C20**|**A114C14**|**A1140C02**|**A1141C20\_Norm**|**A1141C14\_Norm**|**A1140C02\_Norm**|
| :-: | :-: | :-: | :-: | :-: | :-: | :-: |
|Nb Annoted Peaks|2976\.0000000|3277\.0000000|1937\.000000|2885\.0000000|3317\.0000000|1747\.0000000|
|Promoter (<=1kb)|27\.4865591|26\.8538297|26\.277749|26\.9670711|24\.0277359|23\.9267315|
|Promoter (1-2kb)|3\.5954301|4\.0280745|4\.543108|3\.5355286|4\.0397950|4\.6365197|
|Promoter (2-3kb)|3\.5954301|3\.0820873|3\.975219|3\.8821490|4\.0699427|5\.2089296|
|Promoter (3-4kb)|2\.6881720|3\.0515716|3\.407331|2\.5649913|3\.2861019|2\.8048082|
|Promoter (4-5kb)|2\.9233871|3\.1126030|3\.149200|2\.6689775|3\.1655110|3\.2054951|
|Promoter (5-6kb)|2\.6545699|2\.4717730|2\.994321|2\.8769497|2\.4721134|2\.8048082|
|Promoter (6-7kb)|2\.4193548|2\.2886787|2\.736190|2\.3223570|2\.6529997|3\.0910132|
|Promoter (7-8kb)|1\.8481183|1\.8309429|2\.219928|1\.9410745|2\.2912270|2\.3468804|
|Promoter (8-9kb)|1\.6129032|1\.7088801|1\.961797|1\.8024263|1\.6279771|2\.3468804|
|Promoter (9-10kb)|1\.8817204|1\.7088801|1\.548787|1\.6984402|2\.1404884|2\.0034345|
|X5UTR|0\.4032258|0\.6103143|0\.464636|0\.4159445|0\.4823636|0\.3434459|
|X3UTR|2\.1169355|2\.1971315|1\.600413|2\.1837088|2\.1103407|1\.4310246|
|X1st Exon|3\.0577957|3\.6313702|3\.355705|3\.1889081|2\.9846247|3\.2627361|
|Other Exon|4\.5698925|4\.7299359|3\.252452|4\.7140381|4\.2206813|2\.8620492|
|X1st Intron|6\.0819892|6\.0726274|5\.472380|6\.1351820|6\.6023515|6\.0675444|
|Other Intron|11\.3239247|10\.9551419|11\.615901|11\.4731369|11\.6068737|11\.6199199|
|Downnstream (<=300)|0\.2016129|0\.2136100|0\.464636|0\.1733102|0\.2713295|0\.4006869|
|Distal Intergenic|21\.5389785|21\.4525481|20\.960248|21\.4558059|21\.9475430|21\.6370922|

</details>

<details>
  <summary>
    Voir le tableau H3K27ac
  </summary>

|**Sample**|**A1140C09**|**A1140C03**|**A1141C21**|**A1140C09\_Norm**|**A1140C03\_Norm**|**A1141C21\_Norm**|
| :-: | :-: | :-: | :-: | :-: | :-: | :-: |
|Nb Annoted Peaks|5\.630500e+04|5\.165600e+04|4\.176500e+04|4\.670000e+04|4\.707900e+04|3\.273000e+04|
|Promoter (<=1kb)|1\.116775e+01|1\.093774e+01|1\.170119e+01|1\.119589e+01|1\.102843e+01|1\.232223e+01|
|Promoter (1-2kb)|4\.733150e+00|4\.553198e+00|4\.855740e+00|4\.963566e+00|4\.737100e+00|5\.058568e+00|
|Promoter (2-3kb)|3\.939259e+00|4\.015022e+00|4\.158985e+00|4\.194171e+00|4\.052474e+00|4\.171636e+00|
|X5UTR|4\.209218e-01|3\.929844e-01|4\.142224e-01|3\.857694e-01|3\.465652e-01|4\.067651e-01|
|X3UTR|2\.788385e+00|2\.702493e+00|2\.789417e+00|2\.756108e+00|2\.706610e+00|2\.807597e+00|
|X1st Exon|1\.479442e+00|1\.403516e+00|1\.441398e+00|1\.369481e+00|1\.320350e+00|1\.397682e+00|
|Other Exon|4\.631915e+00|4\.454468e+00|4\.597151e+00|4\.312045e+00|4\.367146e+00|4\.621219e+00|
|X1st Intron|1\.320309e+01|1\.340793e+01|1\.322639e+01|1\.327261e+01|1\.326090e+01|1\.311741e+01|
|Other Intron|2\.782701e+01|2\.840328e+01|2\.777924e+01|2\.818903e+01|2\.857568e+01|2\.749488e+01|
|Downnstream (<=300)|6\.997602e-01|7\.220846e-01|7\.494313e-01|7\.329619e-01|7\.505369e-01|8\.288222e-01|
|Distal Intergenic|2\.910932e+01|2\.900728e+01|2\.828684e+01|2\.862838e+01|2\.885421e+01|2\.777319e+01|

|||
| :- | :- |

|**Sample**|**A1141C16**|**A1140C10**|**A1140C04**|**A1141C16\_Norm**|**A1140C10\_Norm**|**A1140C04\_Norm**|
| :-: | :-: | :-: | :-: | :-: | :-: | :-: |
|Nb Annoted Peaks|2\.884400e+04|5\.207700e+04|3\.544400e+04|2\.862300e+04|4\.440100e+04|2\.887800e+04|
|Promoter (<=1kb)|1\.172168e+01|1\.086276e+01|1\.224749e+01|1\.177707e+01|1\.105498e+01|1\.236434e+01|
|Promoter (1-2kb)|5\.259326e+00|4\.508708e+00|4\.627017e+00|5\.243984e+00|4\.717875e+00|5\.044901e+00|
|Promoter (2-3kb)|4\.191513e+00|3\.957601e+00|4\.059926e+00|4\.217606e+00|4\.048081e+00|4\.188482e+00|
|X5UTR|2\.842879e-01|3\.994086e-01|4\.373096e-01|2\.802396e-01|3\.788733e-01|4\.334108e-01|
|X3UTR|2\.607128e+00|2\.803541e+00|2\.920099e+00|2\.581707e+00|2\.773894e+00|2\.940259e+00|
|X1st Exon|1\.144085e+00|1\.349924e+00|1\.500959e+00|1\.117455e+00|1\.298994e+00|1\.376513e+00|
|Other Exon|3\.896824e+00|4\.366611e+00|4\.745514e+00|3\.891828e+00|4\.300663e+00|4\.628827e+00|
|X1st Intron|1\.317432e+01|1\.332450e+01|1\.302618e+01|1\.323081e+01|1\.318840e+01|1\.283936e+01|
|Other Intron|2\.863334e+01|2\.845210e+01|2\.726837e+01|2\.851788e+01|2\.859591e+01|2\.719046e+01|
|Downnstream (<=300)|8\.389960e-01|7\.258483e-01|7\.250875e-01|8\.617368e-01|7\.712778e-01|8\.182795e-01|
|Distal Intergenic|2\.824851e+01|2\.924900e+01|2\.844205e+01|2\.827968e+01|2\.887105e+01|2\.817517e+01|

</details>



<details>
  <summary>
    Voir le tableau H3K4me1 Broad
  </summary>

|**Sample**|**A1140C05**|**A1141C23**|**A1140C11**|**A1140C05\_Norm**|**A1141C23\_Norm**|**A1140C11\_Norm**|
| :-: | :-: | :-: | :-: | :-: | :-: | :-: |
|Nb Annoted Peaks|3\.964500e+04|4\.246300e+04|3\.725200e+04|2\.500200e+04|3\.545900e+04|3\.029400e+04|
|Promoter (<=1kb)|1\.495523e+01|1\.426654e+01|1\.493611e+01|1\.739461e+01|1\.509913e+01|1\.587113e+01|
|Promoter (1-2kb)|4\.431832e+00|4\.441514e+00|4\.378288e+00|4\.591633e+00|4\.495333e+00|4\.495940e+00|
|Promoter (2-3kb)|3\.942490e+00|3\.857476e+00|3\.785032e+00|4\.023678e+00|3\.897459e+00|3\.954579e+00|
|Promoter (3-4kb)|3\.513684e+00|3\.614912e+00|3\.548803e+00|3\.687705e+00|3\.629544e+00|3\.475936e+00|
|Promoter (4-5kb)|3\.238744e+00|3\.238113e+00|3\.175668e+00|3\.391729e+00|3\.288305e+00|3\.390110e+00|
|Promoter (5-6kb)|2\.966326e+00|2\.922544e+00|2\.821325e+00|2\.931766e+00|2\.941425e+00|2\.832244e+00|
|Promoter (6-7kb)|2\.590491e+00|2\.691755e+00|2\.660260e+00|2\.603792e+00|2\.619927e+00|2\.720011e+00|
|Promoter (7-8kb)|2\.302939e+00|2\.286697e+00|2\.313970e+00|2\.495800e+00|2\.318170e+00|2\.310689e+00|
|Promoter (8-9kb)|2\.262580e+00|2\.204272e+00|2\.179749e+00|2\.247820e+00|2\.208184e+00|2\.205057e+00|
|Promoter (9-10kb)|2\.108715e+00|2\.190142e+00|2\.206593e+00|2\.131829e+00|2\.199724e+00|2\.158843e+00|
|X5UTR|2\.572834e-01|2\.378541e-01|2\.523354e-01|2\.559795e-01|2\.256127e-01|2\.772826e-01|
|X3UTR|1\.722790e+00|1\.601394e+00|1\.629443e+00|1\.859851e+00|1\.655433e+00|1\.700007e+00|
|X1st Exon|1\.526044e+00|1\.441255e+00|1\.559648e+00|1\.627870e+00|1\.452382e+00|1\.577870e+00|
|Other Exon|3\.400177e+00|3\.329958e+00|3\.350156e+00|3\.663707e+00|3\.353168e+00|3\.456130e+00|
|X1st Intron|8\.326397e+00|8\.510939e+00|8\.504241e+00|7\.887369e+00|8\.401252e+00|8\.232653e+00|
|Other Intron|1\.818136e+01|1\.870805e+01|1\.819231e+01|1\.637069e+01|1\.822668e+01|1\.749851e+01|
|Downnstream (<=300)|3\.304326e-01|3\.249888e-01|3\.167615e-01|3\.799696e-01|3\.299585e-01|3\.102925e-01|
|Distal Intergenic|2\.394249e+01|2\.413160e+01|2\.418931e+01|2\.245420e+01|2\.365831e+01|2\.353271e+01|

|||
| :- | :- |

|**Sample**|**A1141C24**|**A1141C18**|**A1140C12**|**A1141C24\_Norm**|**A1141C18\_Norm**|**A1140C12\_Norm**|
| :-: | :-: | :-: | :-: | :-: | :-: | :-: |
|Nb Annoted Peaks|3\.498500e+04|4\.257800e+04|3\.961800e+04|3\.450800e+04|3\.569200e+04|3\.261000e+04|
|Promoter (<=1kb)|1\.500643e+01|1\.423270e+01|1\.438235e+01|1\.508056e+01|1\.475120e+01|1\.539405e+01|
|Promoter (1-2kb)|4\.664856e+00|4\.368453e+00|4\.482811e+00|4\.659789e+00|4\.538832e+00|4\.504753e+00|
|Promoter (2-3kb)|3\.944548e+00|3\.879938e+00|3\.861881e+00|3\.880260e+00|3\.902835e+00|3\.836247e+00|
|Promoter (3-4kb)|3\.581535e+00|3\.485368e+00|3\.548892e+00|3\.628144e+00|3\.653480e+00|3\.593990e+00|
|Promoter (4-5kb)|3\.169930e+00|3\.210578e+00|3\.268716e+00|3\.210850e+00|3\.191191e+00|3\.287335e+00|
|Promoter (5-6kb)|2\.978419e+00|2\.907605e+00|2\.910293e+00|2\.950041e+00|2\.958646e+00|2\.980681e+00|
|Promoter (6-7kb)|2\.663999e+00|2\.623420e+00|2\.604877e+00|2\.625478e+00|2\.656057e+00|2\.649494e+00|
|Promoter (7-8kb)|2\.352437e+00|2\.285218e+00|2\.319653e+00|2\.382056e+00|2\.325451e+00|2\.324440e+00|
|Promoter (8-9kb)|2\.163785e+00|2\.247640e+00|2\.201020e+00|2\.158920e+00|2\.252606e+00|2\.171113e+00|
|Promoter (9-10kb)|2\.152351e+00|2\.184227e+00|2\.226261e+00|2\.158920e+00|2\.202174e+00|2\.125115e+00|
|X5UTR|2\.343862e-01|2\.489549e-01|2\.372659e-01|2\.405239e-01|2\.381486e-01|2\.361239e-01|
|X3UTR|1\.657853e+00|1\.646390e+00|1\.552325e+00|1\.663382e+00|1\.644626e+00|1\.603803e+00|
|X1st Exon|1\.486351e+00|1\.444408e+00|1\.456409e+00|1\.477918e+00|1\.400874e+00|1\.521006e+00|
|Other Exon|3\.355724e+00|3\.386725e+00|3\.364632e+00|3\.349948e+00|3\.423736e+00|3\.351733e+00|
|X1st Intron|8\.437902e+00|8\.574851e+00|8\.559241e+00|8\.400951e+00|8\.444469e+00|8\.313401e+00|
|Other Intron|1\.806774e+01|1\.859881e+01|1\.849159e+01|1\.810015e+01|1\.821977e+01|1\.788408e+01|
|Downnstream (<=300)|3\.401458e-01|3\.170652e-01|3\.407542e-01|3\.274603e-01|3\.390115e-01|3\.710518e-01|
|Distal Intergenic|2\.374160e+01|2\.435765e+01|2\.419102e+01|2\.370465e+01|2\.385689e+01|2\.385158e+01|

</details>


<details>
  <summary>
    Voir le tableau H3K4me1 Narrow
  </summary>

|**Sample**|**A1140C05**|**A1141C23**|**A1140C11**|**A1140C05\_Norm**|**A1141C23\_Norm**|**A1140C11\_Norm**|
| :-: | :-: | :-: | :-: | :-: | :-: | :-: |
|Nb Annoted Peaks|3\.964500e+04|4\.246300e+04|3\.725200e+04|2\.500200e+04|3\.545900e+04|3\.029400e+04|
|Promoter (<=1kb)|1\.495523e+01|1\.426654e+01|1\.493611e+01|1\.739461e+01|1\.509913e+01|1\.587113e+01|
|Promoter (1-2kb)|4\.431832e+00|4\.441514e+00|4\.378288e+00|4\.591633e+00|4\.495333e+00|4\.495940e+00|
|Promoter (2-3kb)|3\.942490e+00|3\.857476e+00|3\.785032e+00|4\.023678e+00|3\.897459e+00|3\.954579e+00|
|Promoter (3-4kb)|3\.513684e+00|3\.614912e+00|3\.548803e+00|3\.687705e+00|3\.629544e+00|3\.475936e+00|
|Promoter (4-5kb)|3\.238744e+00|3\.238113e+00|3\.175668e+00|3\.391729e+00|3\.288305e+00|3\.390110e+00|
|Promoter (5-6kb)|2\.966326e+00|2\.922544e+00|2\.821325e+00|2\.931766e+00|2\.941425e+00|2\.832244e+00|
|Promoter (6-7kb)|2\.590491e+00|2\.691755e+00|2\.660260e+00|2\.603792e+00|2\.619927e+00|2\.720011e+00|
|Promoter (7-8kb)|2\.302939e+00|2\.286697e+00|2\.313970e+00|2\.495800e+00|2\.318170e+00|2\.310689e+00|
|Promoter (8-9kb)|2\.262580e+00|2\.204272e+00|2\.179749e+00|2\.247820e+00|2\.208184e+00|2\.205057e+00|
|Promoter (9-10kb)|2\.108715e+00|2\.190142e+00|2\.206593e+00|2\.131829e+00|2\.199724e+00|2\.158843e+00|
|X5UTR|2\.572834e-01|2\.378541e-01|2\.523354e-01|2\.559795e-01|2\.256127e-01|2\.772826e-01|
|X3UTR|1\.722790e+00|1\.601394e+00|1\.629443e+00|1\.859851e+00|1\.655433e+00|1\.700007e+00|
|X1st Exon|1\.526044e+00|1\.441255e+00|1\.559648e+00|1\.627870e+00|1\.452382e+00|1\.577870e+00|
|Other Exon|3\.400177e+00|3\.329958e+00|3\.350156e+00|3\.663707e+00|3\.353168e+00|3\.456130e+00|
|X1st Intron|8\.326397e+00|8\.510939e+00|8\.504241e+00|7\.887369e+00|8\.401252e+00|8\.232653e+00|
|Other Intron|1\.818136e+01|1\.870805e+01|1\.819231e+01|1\.637069e+01|1\.822668e+01|1\.749851e+01|
|Downnstream (<=300)|3\.304326e-01|3\.249888e-01|3\.167615e-01|3\.799696e-01|3\.299585e-01|3\.102925e-01|
|Distal Intergenic|2\.394249e+01|2\.413160e+01|2\.418931e+01|2\.245420e+01|2\.365831e+01|2\.353271e+01|

|||
| :- | :- |

|**Sample**|**A1141C24**|**A1141C18**|**A1140C12**|**A1141C24\_Norm**|**A1141C18\_Norm**|**A1140C12\_Norm**|
| :-: | :-: | :-: | :-: | :-: | :-: | :-: |
|Nb Annoted Peaks|3\.498500e+04|4\.257800e+04|3\.961800e+04|3\.450800e+04|3\.569200e+04|3\.261000e+04|
|Promoter (<=1kb)|1\.500643e+01|1\.423270e+01|1\.438235e+01|1\.508056e+01|1\.475120e+01|1\.539405e+01|
|Promoter (1-2kb)|4\.664856e+00|4\.368453e+00|4\.482811e+00|4\.659789e+00|4\.538832e+00|4\.504753e+00|
|Promoter (2-3kb)|3\.944548e+00|3\.879938e+00|3\.861881e+00|3\.880260e+00|3\.902835e+00|3\.836247e+00|
|Promoter (3-4kb)|3\.581535e+00|3\.485368e+00|3\.548892e+00|3\.628144e+00|3\.653480e+00|3\.593990e+00|
|Promoter (4-5kb)|3\.169930e+00|3\.210578e+00|3\.268716e+00|3\.210850e+00|3\.191191e+00|3\.287335e+00|
|Promoter (5-6kb)|2\.978419e+00|2\.907605e+00|2\.910293e+00|2\.950041e+00|2\.958646e+00|2\.980681e+00|
|Promoter (6-7kb)|2\.663999e+00|2\.623420e+00|2\.604877e+00|2\.625478e+00|2\.656057e+00|2\.649494e+00|
|Promoter (7-8kb)|2\.352437e+00|2\.285218e+00|2\.319653e+00|2\.382056e+00|2\.325451e+00|2\.324440e+00|
|Promoter (8-9kb)|2\.163785e+00|2\.247640e+00|2\.201020e+00|2\.158920e+00|2\.252606e+00|2\.171113e+00|
|Promoter (9-10kb)|2\.152351e+00|2\.184227e+00|2\.226261e+00|2\.158920e+00|2\.202174e+00|2\.125115e+00|
|X5UTR|2\.343862e-01|2\.489549e-01|2\.372659e-01|2\.405239e-01|2\.381486e-01|2\.361239e-01|
|X3UTR|1\.657853e+00|1\.646390e+00|1\.552325e+00|1\.663382e+00|1\.644626e+00|1\.603803e+00|
|X1st Exon|1\.486351e+00|1\.444408e+00|1\.456409e+00|1\.477918e+00|1\.400874e+00|1\.521006e+00|
|Other Exon|3\.355724e+00|3\.386725e+00|3\.364632e+00|3\.349948e+00|3\.423736e+00|3\.351733e+00|
|X1st Intron|8\.437902e+00|8\.574851e+00|8\.559241e+00|8\.400951e+00|8\.444469e+00|8\.313401e+00|
|Other Intron|1\.806774e+01|1\.859881e+01|1\.849159e+01|1\.810015e+01|1\.821977e+01|1\.788408e+01|
|Downnstream (<=300)|3\.401458e-01|3\.170652e-01|3\.407542e-01|3\.274603e-01|3\.390115e-01|3\.710518e-01|
|Distal Intergenic|2\.374160e+01|2\.435765e+01|2\.419102e+01|2\.370465e+01|2\.385689e+01|2\.385158e+01|

</details>



**Représentation graphique en barplot de l’annotation pour chaque population**

Feature_H3K3

![Feature_H3K3](https://github.com/Antonin-w/CHIP_seq_analysis_Nextflow/blob/fd2b3daa3627f3a21c1c455fe1bee3a2acadc2c9/Images/Feature_H3K3.jpg?raw=true)


Feature_H3K4
![Feature_H3K4](https://github.com/Antonin-w/CHIP_seq_analysis_Nextflow/blob/fd2b3daa3627f3a21c1c455fe1bee3a2acadc2c9/Images/feature_H3K4.jpg?raw=true)

**Distribution des features d’annotation autour du TSS**

![merge_from_ofoct-2.jpg](https://github.com/Antonin-w/CHIP_seq_analysis_Nextflow/blob/fd2b3daa3627f3a21c1c455fe1bee3a2acadc2c9/Images/merge_from_ofoct-2.jpg?raw=true)


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

