#!/bin/r env

library(argparse)
library(stringr)

parser = ArgumentParser()

parser$add_argument("-wd", "--wd", 
                    help="Path to output files")

parser$add_argument("-i", "--in", 
                    help="Path to input files")

parser$add_argument("-pattern1", "--RDataPattern", 
                    help="Pattern to recognize input files")

parser$add_argument("-fd", "--fd", 
                    help="Path to function files and reference files")

parser$add_argument("-gzvcf", "--gzvcf", 
                    help="name of the gzvcf file")

parser$add_argument("-o", "--out", 
                    help="Prefix of output files")

parser$add_argument("-rkeep", "--keep_regexp", 
                    help="Regular expression to identify samples of interest")

parser$add_argument("-ebed", "--exclude_bed", 
                    help="name of .bed file containing no core genomic regions to exclude from the VCF file")

parser$add_argument("-gff", "--ref_gff", 
                    help="name of .gff file containing coordinates of genomic regions")

parser$add_argument("-n", "--n", 
                    help="Number of iterations to split the vcf file")

parser$add_argument("-thres", "--thres", 
                    help="Minimun read depth to call an allele")

args = parser$parse_args()

wd = args$wd
setwd(wd)

fd = args$fd

# First filter arguments
gzvcf = args$gzvcf
foutput = paste0(args$out, '_FirstFilter')
keep_regexp = args$keep_regexp
exclude_bed = file.path(fd, args$exclude_bed)

# Second filter arguments
svcf = paste0(args$out, '_FirstFilter.recode.vcf')
ref_gff_file = file.path(fd, args$ref_gff)
soutput = paste0(args$out, '_SecondFilter')

# Load VCF to R and generate rGenome object

tvcf = paste0(args$out, '_SecondFilter.recode.vcf')
n = as.integer(args$n)
threshold = as.numeric(args$thres)

imagename = paste0(args$out, '_1st_2nd_filters.RData')



### Step 2: Low support (Read depth) and missing data----


# Check packages and functions----

source(file.path(fd,'load_libraries.R'))
source(file.path(fd,'functions.R'))
#sourceCpp(file.path(fd,'Rcpp_functions.cpp'))



#### 1.2 Load the rGenome objects:----

for(file in list.files('PerVen/Pviv/', pattern = '.+thres1.+.RData')){
  load(file.path('PerVen/Pviv/', file))
  rm(list = ls()[!grepl('PvP01_(\\d+|MIT|API)_v1_rGenome_object', ls())])
}

rGenome_objects = list(PvP01_01_v1_rGenome_object,
                       PvP01_02_v1_rGenome_object,
                       PvP01_03_v1_rGenome_object,
                       PvP01_04_v1_rGenome_object,
                       PvP01_05_v1_rGenome_object,
                       PvP01_06_v1_rGenome_object,
                       PvP01_07_v1_rGenome_object,
                       PvP01_08_v1_rGenome_object,
                       PvP01_09_v1_rGenome_object,
                       PvP01_10_v1_rGenome_object,
                       PvP01_11_v1_rGenome_object,
                       PvP01_12_v1_rGenome_object,
                       PvP01_13_v1_rGenome_object,
                       PvP01_14_v1_rGenome_object,
                       PvP01_MIT_v1_rGenome_object,
                       PvP01_API_v1_rGenome_object
)

rm(list = c('PvP01_01_v1_rGenome_object',
            'PvP01_02_v1_rGenome_object',
            'PvP01_03_v1_rGenome_object',
            'PvP01_04_v1_rGenome_object',
            'PvP01_05_v1_rGenome_object',
            'PvP01_06_v1_rGenome_object',
            'PvP01_07_v1_rGenome_object',
            'PvP01_08_v1_rGenome_object',
            'PvP01_09_v1_rGenome_object',
            'PvP01_10_v1_rGenome_object',
            'PvP01_11_v1_rGenome_object',
            'PvP01_12_v1_rGenome_object',
            'PvP01_13_v1_rGenome_object',
            'PvP01_14_v1_rGenome_object',
            'PvP01_MIT_v1_rGenome_object',
            'PvP01_API_v1_rGenome_object'
))

source('../functions_libraries/functions.R')
Pviv_rGenome_object = rGenome(gt = NULL,
                              loci_table = NULL,
                              metadata = NULL)

metadata = rGenome_objects[[1]]@metadata

for(obj in 1:length(rGenome_objects)){
  Pviv_rGenome_object@gt = rbind(Pviv_rGenome_object@gt, rGenome_objects[[obj]]@gt)
  Pviv_rGenome_object@loci_table = rbind(Pviv_rGenome_object@loci_table, rGenome_objects[[obj]]@loci_table)
  rGenome_objects[[obj]] = NULL
  
}

Pviv_rGenome_object@metadata = rbind(Pviv_rGenome_object@metadata, rGenome_objects[[obj]]@metadata)

rm(list = c('rGenome_objects', 'obj'))


## Metadata ----

# There are some spelling mistakes in the sample codes of the sequencing files, The table VENPER2019_2022.csv match the incorrect codes with the correct codes
external_metadata = read.csv('merged_metadata.csv')
PerVen_seq_codes = read.csv('VENPER2019_2022.csv')

# Remove Duplicated records
external_metadata = external_metadata[!duplicated(external_metadata$Sample_id),]
PerVen_seq_codes = PerVen_seq_codes[!duplicated(PerVen_seq_codes$PreferedSampleID),]

# merge incorrect and correct codes with metadata
external_metadata = merge(external_metadata, PerVen_seq_codes, by.x = 'Sample_id', by.y = 'AlternateSampleID', all.x = T)

external_metadata %<>% mutate(PreferedSampleID = case_when(
  is.na(PreferedSampleID) ~ Sample_id,
  !is.na(PreferedSampleID) ~ PreferedSampleID
))

# Add metadata to the rGenome object
Pviv_rGenome_object@metadata = merge(Pviv_rGenome_object@metadata, external_metadata[,c("PreferedSampleID", "Sample_id", "Study", "Country", 'Date_of_Collection', "Year_of_Collection",'Subnational_level0', 'Subnational_level1', "Subnational_level2")], by.x = 'Sample_id', by.y = 'PreferedSampleID', all.x = T)

# Sorting samples
rownames(Pviv_rGenome_object@metadata) = Pviv_rGenome_object@metadata$Sample_id
Pviv_rGenome_object@metadata = Pviv_rGenome_object@metadata[colnames(Pviv_rGenome_object@gt),]
Pviv_rGenome_object@metadata[is.na(Pviv_rGenome_object@metadata$Country),][['Country']] = 'Venezuela'








#### Setp 2.1: Site coverage (Read depth by site)----
# Calculate the total and mean read depth per each site per each country using the function Summarise_ReadDepth.

Read_Depth_Summ = NULL

for(w in 1:100){
  start = Sys.time()
  Read_Depth_Summ = rbind(Read_Depth_Summ, Summarise_ReadDepth(obj = Pviv_rGenome_object, by = 'Country', w = w, n = 100))
  end = Sys.time()
  print(paste0('Iteration ',w , ' took ', end - start, ' secs'))
}

# Transform the Read_Depth_Summ table from wide to long format

## Extract data from Peru
Read_Depth_Summ_Peru = data.frame(Chrom = 
                                    gsub('_\\d+$','',
                                         rownames(Read_Depth_Summ)),
                                  Pos = gsub('^.+_','',
                                             rownames(Read_Depth_Summ)),
                                  Country = 'Peru',
                                  Read_Depth_Summ[grepl('_Per',
                                                        colnames(Read_Depth_Summ))])

colnames(Read_Depth_Summ_Peru) = gsub('_Per.$','',
                                      colnames(Read_Depth_Summ_Peru))
## Extract data Venezuela
Read_Depth_Summ_Venezuela = data.frame(Chrom = 
                                         gsub('_\\d+$','',
                                              rownames(Read_Depth_Summ)),
                                       Pos = gsub('^.+_','',
                                                  rownames(Read_Depth_Summ)),
                                       Country = 'Venezuela',
                                       Read_Depth_Summ[grepl('_Ven',
                                                             colnames(Read_Depth_Summ))])

colnames(Read_Depth_Summ_Venezuela) = gsub('_Ven.+$','',
                                           colnames(Read_Depth_Summ_Venezuela))

## Extract data from the total Population
Read_Depth_Summ_Total = data.frame(Chrom = gsub('_\\d+$','',
                                                rownames(Read_Depth_Summ)),
                                   Pos = gsub('^.+_','',
                                              rownames(Read_Depth_Summ)),
                                   Country = 'Total',
                                   Read_Depth_Summ[grepl('_Total$',
                                                         colnames(Read_Depth_Summ))])

colnames(Read_Depth_Summ_Total) = gsub('_Total$','',
                                       colnames(Read_Depth_Summ_Total))


Read_Depth_Summ = rbind(Read_Depth_Summ_Peru,
                        Read_Depth_Summ_Venezuela,
                        Read_Depth_Summ_Total)

rm(list = c('Read_Depth_Summ_Peru',
            'Read_Depth_Summ_Venezuela',
            'Read_Depth_Summ_Total'))

Read_Depth_Summ$Pos = as.integer(Read_Depth_Summ$Pos)
Read_Depth_Summ$Country = factor(Read_Depth_Summ$Country,
                                 levels = c('Peru',
                                            'Venezuela',
                                            'Total'))



####4. Plot the total and mean read depth by site, by country, and by chromosome.----

Read_Depth_Summ %>% ggplot(aes(x = Pos, y = Total_ReadDepth, color = Country, group = Country)) +
  geom_line() +
  facet_grid(Chrom~Country, scales = 'free_y') +
  scale_color_manual(values = c('dodgerblue3', 'firebrick3', 'black')) +
  theme_bw() +
  theme(axis.text.x = element_blank())

Read_Depth_Summ %>% ggplot(aes(x = Pos, y = mean_ReadDepth, color = Country, group = Country)) +
  geom_line() +
  facet_grid(Chrom~Country, scales = 'free_y') +
  scale_color_manual(values = c('dodgerblue3', 'firebrick3', 'black')) +
  theme_bw() +
  theme(axis.text.x = element_blank())



#### Step 2.2: Remove samples with low amplification rate.----


# 6. Calculate the proportion of amplified loci (amplification rate), with a coverage equals or greater than 1, 2, 3, 4 and 5, of the samples.


# Sample amplification rate including all detected reads
SampAmpRate_thres0 = SampleAmplRate(Pviv_rGenome_object, update = FALSE, threshold = NULL)
SampAmpRate_thres0 = data.frame(Sample_id = names(SampAmpRate_thres0),
                                AmpRate = SampAmpRate_thres0,
                                Threshold = 1)
SampAmpRate_thres0 = merge(SampAmpRate_thres0,
                           Pviv_rGenome_object@metadata[,c('Sample_id', 'Country')],
                           by = 'Sample_id',
                           all.x = TRUE)

# Sample amplification rate including only sites supported by at least 2 reads
SampAmpRate_thres1 = SampleAmplRate(Pviv_rGenome_object, update = FALSE, threshold = 1, n = 100)
SampAmpRate_thres1 = data.frame(Sample_id = names(SampAmpRate_thres1),
                                AmpRate = SampAmpRate_thres1,
                                Threshold = 2)
SampAmpRate_thres1 = merge(SampAmpRate_thres1,
                           Pviv_rGenome_object@metadata[,c('Sample_id', 'Country')],
                           by = 'Sample_id',
                           all.x = TRUE)

# Sample amplification rate including only sites supported by at least 3 reads
SampAmpRate_thres2 = SampleAmplRate(Pviv_rGenome_object, update = FALSE, threshold = 2, n = 100)
SampAmpRate_thres2 = data.frame(Sample_id = names(SampAmpRate_thres2),
                                AmpRate = SampAmpRate_thres2,
                                Threshold = 3)
SampAmpRate_thres2 = merge(SampAmpRate_thres2,
                           Pviv_rGenome_object@metadata[,c('Sample_id', 'Country')],
                           by = 'Sample_id',
                           all.x = TRUE)

# Sample amplification rate including only sites supported by at least 4 reads
SampAmpRate_thres3 = SampleAmplRate(Pviv_rGenome_object, update = FALSE, threshold = 3, n = 100)
SampAmpRate_thres3 = data.frame(Sample_id = names(SampAmpRate_thres3),
                                AmpRate = SampAmpRate_thres3,
                                Threshold = 4)
SampAmpRate_thres3 = merge(SampAmpRate_thres3,
                           Pviv_rGenome_object@metadata[,c('Sample_id', 'Country')],
                           by = 'Sample_id',
                           all.x = TRUE)

# Sample amplification rate including only sites supported by at least 5 reads
SampAmpRate_thres4 = SampleAmplRate(Pviv_rGenome_object, update = FALSE, threshold = 4, n = 100)
SampAmpRate_thres4 = data.frame(Sample_id = names(SampAmpRate_thres4),
                                AmpRate = SampAmpRate_thres4,
                                Threshold = 5)
SampAmpRate_thres4 = merge(SampAmpRate_thres4,
                           Pviv_rGenome_object@metadata[,c('Sample_id', 'Country')],
                           by = 'Sample_id',
                           all.x = TRUE)

# Combining all vectors of SampAmpRate
SampAmpRate_summ = rbind(SampAmpRate_thres0,
                         SampAmpRate_thres1,
                         SampAmpRate_thres2,
                         SampAmpRate_thres3,
                         SampAmpRate_thres4)

rm(list = c('SampAmpRate_thres0',
            'SampAmpRate_thres1',
            'SampAmpRate_thres2',
            'SampAmpRate_thres3',
            'SampAmpRate_thres4'))

#7. Plot the proportion of samples that success to amplify with an specific amplification rate.

SampAmpRate_summ %>% group_by(Country, Threshold) %>%
  summarise(AmpRate5 = round(100*sum(AmpRate >= .05)/n(), 1),
            AmpRate10 = round(100*sum(AmpRate >= .10)/n(), 1),
            AmpRate15 = round(100*sum(AmpRate >= .15)/n(), 1),
            AmpRate20 = round(100*sum(AmpRate >= .20)/n(), 1),
            AmpRate25 = round(100*sum(AmpRate >= .25)/n(), 1),
            AmpRate30 = round(100*sum(AmpRate >= .30)/n(), 1),
            AmpRate35 = round(100*sum(AmpRate >= .35)/n(), 1),
            AmpRate40 = round(100*sum(AmpRate >= .40)/n(), 1),
            AmpRate45 = round(100*sum(AmpRate >= .45)/n(), 1),
            AmpRate50 = round(100*sum(AmpRate >= .50)/n(), 1),
            AmpRate55 = round(100*sum(AmpRate >= .55)/n(), 1),
            AmpRate60 = round(100*sum(AmpRate >= .60)/n(), 1),
            AmpRate65 = round(100*sum(AmpRate >= .65)/n(), 1),
            AmpRate70 = round(100*sum(AmpRate >= .70)/n(), 1),
            AmpRate75 = round(100*sum(AmpRate >= .75)/n(), 1),
            AmpRate80 = round(100*sum(AmpRate >= .80)/n(), 1),
            AmpRate85 = round(100*sum(AmpRate >= .85)/n(), 1),
            AmpRate90 = round(100*sum(AmpRate >= .90)/n(), 1),
            AmpRate95 = round(100*sum(AmpRate >= .95)/n(), 1),
            AmpRate100 = round(100*sum(AmpRate >= 1)/n(), 1)) %>%
  pivot_longer(cols = paste0('AmpRate', seq(5, 100, 5)),
               values_to = 'Percentage',
               names_to = 'AmpRate') %>%
  mutate(AmpRate = as.numeric(gsub('AmpRate','', AmpRate)))%>%
  ggplot(aes(x = AmpRate, y = Percentage, color = as.factor(Threshold), group = as.factor(Threshold))) +
  geom_line() +
  geom_vline(xintercept = 75, linetype = 2) +
  facet_grid(.~Country) +
  theme_bw() +
  labs(x = '% of amplified loci (amplification rate)', y = '% of Samples', color = 'Min Coverage')


9. Calculate the proportion of amplified samples per each loci using a coverage of 1 and 5 to call an allele.
```{r, eval = FALSE}
# Create an rGenome object per  each country

Peru_rGenome = filter_samples(obj = Pviv_rGenome_object,
                              v = Pviv_rGenome_object@metadata$Country == 'Perú')


Venezuela_rGenome = filter_samples(obj = Pviv_rGenome_object,
                                   v = Pviv_rGenome_object@metadata$Country == 'Venezuela')


Peru_LocusAmpRate_thres0 = LocusAmplRate(Peru_rGenome, update = FALSE, threshold = NULL)

Peru_LocusAmpRate_thres0 = 
  data.frame(Locus_id = names(Peru_LocusAmpRate_thres0),
             Chrom = gsub(
               '_\\d+$',
               '',
               names(Peru_LocusAmpRate_thres0)),
             Pos = as.integer(gsub(
               '^.+_',
               '',
               names(Peru_LocusAmpRate_thres0))),
             AmpRate = Peru_LocusAmpRate_thres0,
             Threshold = 1)

Venezuela_LocusAmpRate_thres0 = LocusAmplRate(Venezuela_rGenome, update = FALSE, threshold = NULL)

Venezuela_LocusAmpRate_thres0 = 
  data.frame(Locus_id = names(Venezuela_LocusAmpRate_thres0),
             Chrom = gsub(
               '_\\d+$',
               '',
               names(Venezuela_LocusAmpRate_thres0)),
             Pos = as.integer(gsub(
               '^.+_',
               '',
               names(Venezuela_LocusAmpRate_thres0))),
             AmpRate = Venezuela_LocusAmpRate_thres0,
             Threshold = 1)


Peru_LocusAmpRate_thres4 = LocusAmplRate(Peru_rGenome, update = FALSE, threshold = 4, n =100)

Peru_LocusAmpRate_thres4 = 
  data.frame(Locus_id = names(Peru_LocusAmpRate_thres4),
             Chrom = gsub(
               '_\\d+$',
               '',
               names(Peru_LocusAmpRate_thres4)),
             Pos = as.integer(gsub(
               '^.+_',
               '',
               names(Peru_LocusAmpRate_thres4))),
             AmpRate = Peru_LocusAmpRate_thres4,
             Threshold = 5)


Venezuela_LocusAmpRate_thres4 = LocusAmplRate(Venezuela_rGenome, update = FALSE, threshold = 4, n =100)

Venezuela_LocusAmpRate_thres4 = 
  data.frame(Locus_id = names(Venezuela_LocusAmpRate_thres4),
             Chrom = gsub(
               '_\\d+$',
               '',
               names(Venezuela_LocusAmpRate_thres4)),
             Pos = as.integer(gsub(
               '^.+_',
               '',
               names(Venezuela_LocusAmpRate_thres4))),
             AmpRate = Venezuela_LocusAmpRate_thres4,
             Threshold = 5)
```

#10. Manhattan plot the of the proportion of amplified samples by loci and by country.

ggdraw() + 
  draw_plot(
    Read_Depth_Summ %>%
      filter(grepl('_01_', Chrom), Country == 'Venezuela')%>% ggplot(aes(x = Pos, y = mean_ReadDepth)) + 
      geom_line() +
      facet_grid(.~Chrom)+
      theme_bw()+
      theme(axis.text.x = element_blank(),
            legend.position = 'none',
            axis.title.x = element_blank())+
      labs(y = 'Mean read depth)'),
    x = .5,
    y = .75,
    heigh = .25,
    width = .5
  )+
  draw_plot(
    Venezuela_LocusAmpRate_thres0 %>%
      filter(grepl('_01_', Chrom))%>% ggplot(aes(x = Pos, y = AmpRate)) + 
      geom_point(alpha = .3, size = .25) +
      theme_bw()+
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank())+
      labs(y = 'Amp. rate (thres = 1)'),
    x = .5,
    y = .5,
    heigh = .25,
    width = .5
  )+
  draw_plot(
    Venezuela_LocusAmpRate_thres4 %>%
      filter(grepl('_01_', Chrom))%>% ggplot(aes(x = Pos, y = AmpRate)) + 
      geom_point(alpha = .3, size = .25) +
      theme_bw()+
      theme(axis.text.x = element_blank())+
      labs(y = 'Amp. rate (thres = 5)',
           x = 'Venezuela'),
    x = .5,
    y = .25,
    heigh = .25,
    width = .5
  ) + 
  draw_plot(
    Read_Depth_Summ %>%
      filter(grepl('_01_', Chrom), Country == 'Peru')%>% ggplot(aes(x = Pos, y = mean_ReadDepth)) + 
      geom_line() +
      facet_grid(.~Chrom)+
      theme_bw()+
      theme(axis.text.x = element_blank(),
            legend.position = 'none',
            axis.title.x = element_blank())+
      labs(y = 'Mean read depth'),
    x = 0,
    y = .75,
    heigh = .25,
    width = .5
  )+
  draw_plot(
    Peru_LocusAmpRate_thres0 %>%
      filter(grepl('_01_', Chrom))%>% ggplot(aes(x = Pos, y = AmpRate)) + 
      geom_point(alpha = .3, size = .25) +
      theme_bw()+
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank())+
      labs(y = 'Amp. rate (thres = 1)'),
    x = 0,
    y = .5,
    heigh = .25,
    width = .5
  )+
  draw_plot(
    Peru_LocusAmpRate_thres4 %>%
      filter(grepl('_01_', Chrom))%>% ggplot(aes(x = Pos, y = AmpRate)) + 
      geom_point(alpha = .3, size = .25) +
      theme_bw()+
      theme(axis.text.x = element_blank())+
      labs(y = 'Amp. rate (thres = 5)',
           x = 'Peru'),
    x = 0,
    y = .25,
    heigh = .25,
    width = .5
  ) +
  draw_plot(data.frame(Country = 'Venezuela',
                       mean_ReadDepth = Read_Depth_Summ %>%
                         filter(grepl('_01_', Chrom), Country == 'Venezuela') %>%
                         select(mean_ReadDepth) %>% unlist(),
                       LocusAmpRate = Venezuela_LocusAmpRate_thres4 %>%
                         filter(grepl('_01_', Chrom)) %>%
                         select(AmpRate) %>% unlist()
  )%>%
    ggplot(aes(x = mean_ReadDepth, y = LocusAmpRate))+
    geom_point(alpha = .3)+
    theme_bw()+
    theme(legend.position = 'none') + 
    labs(y = 'Amp. rate (thres = 5)',
         x = 'Mean RD'),
  x = .5,
  y = 0,
  heigh = .25,
  width = .5) +
  draw_plot(data.frame(Country = 'Peru',
                       mean_ReadDepth = Read_Depth_Summ %>%
                         filter(grepl('_01_', Chrom), Country == 'Peru') %>%
                         select(mean_ReadDepth) %>% unlist(),
                       LocusAmpRate = Venezuela_LocusAmpRate_thres4 %>%
                         filter(grepl('_01_', Chrom)) %>%
                         select(AmpRate) %>% unlist()
  )%>%
    ggplot(aes(x = mean_ReadDepth, y = LocusAmpRate))+
    geom_point(alpha = .3)+
    theme_bw()+
    theme(legend.position = 'none') + 
    labs(y = 'Amp. rate (thres = 5)',
         x = 'Mean RD'),
  x = 0,
  y = 0,
  heigh = .25,
  width = .5)



#11. Mask alleles with less than 5 reads of coverage 


# rm(Peru_rGenome)
# rm(Venezuela_rGenome)
Pviv_rGenome_thres5 = Pviv_rGenome_object

Pviv_rGenome_thres5@gt = prune_alleles(Pviv_rGenome_thres5, threshold = 4, n = 100)





#12. Check the distribution of the amplification rate of the samples using a histogram


Pviv_rGenome_thres5@metadata$SampAmpRate = SampAmpRate_summ %>% filter(Threshold == 5) %>% select(AmpRate) %>% unlist

Pviv_rGenome_thres5@metadata %>% ggplot(aes(x = SampAmpRate, fill = Country))+
  geom_histogram(position = 'stack', binwidth = .05, color = 'gray40')+
  theme_bw()+
  labs(title = 'Sample Amplification Rate distribution',
       x = 'Amplification Rate')



#13 Remove samples with less than 50% of amplified loci

Pviv_rGenome_thres5 = filter_samples(obj = Pviv_rGenome_thres5, v =
                                       Pviv_rGenome_thres5@metadata$SampAmpRate >= .5)

#14. Remove monomorphic sites

#14.1 Update allele counts

allele_counts = NULL

test = get_AC(obj = Pviv_rGenome_thres5, w = w, n = 100)

for(w in 1:100){
  allele_counts = rbind(allele_counts, get_AC(obj = Pviv_rGenome_thres5, w = w, n = 100))
  print(w)
}

Pviv_rGenome_thres5@loci_table = cbind(Pviv_rGenome_thres5@loci_table,
                                       allele_counts)

#14.2 Filter monomorphic sites
Pviv_rGenome_thres5 = filter_loci(Pviv_rGenome_thres5,
                                  v = Pviv_rGenome_thres5@loci_table$Cardinality > 1)

#15. Remove loci with high missing data

#15.1 Calculate the amplification rate of each loci
Pviv_rGenome_thres5 = LocusAmplRate(Pviv_rGenome_thres5)

#15.2 Check the distribution of the amplification rate of the samples using a histogram

Pviv_rGenome_thres5@loci_table %>%
  mutate(CHROM2  = gsub('(^(\\d|P|v)+_|_v1)', '',CHROM))%>%
  ggplot(aes(x = POS, y = LocusAmplRate)) +
  geom_point(alpha = 0.7, size = .25) +
  facet_grid(.~ CHROM2)+
  theme_bw()+
  labs(y = 'Ampl. Rate', x = 'Chromosomal position')+
  theme(legend.position = 'right',
        axis.text.x = element_blank())

#15.3 Filter loci with amplification rate below .75
Pviv_rGenome_thres5 = filter_loci(Pviv_rGenome_thres5, v = Pviv_rGenome_thres5@loci_table$LocusAmplRate >= .75)

#16. Update alternative alleles

Pviv_rGenome_thres5@loci_table$ALT =
  ifelse(Pviv_rGenome_thres5@loci_table$REF !=
           gsub(',([ATCG]|\\*)+',
                '',
                gsub(':\\d+',
                     '',
                     Pviv_rGenome_thres5@loci_table$Alleles)),
         gsub(':\\d+',
              '',
              Pviv_rGenome_thres5@loci_table$Alleles),
         gsub('^([ATCG]|\\*)+,',
              '',
              gsub(':\\d+', '', Pviv_rGenome_thres5@loci_table$Alleles))
  )

### Step 3: Filter of potential PCR genotyping artifacts (Homopolymers & Short Tandem Repeats)

#### Step 3.1 Clasify SNPs, INDELS, Homopolymers and Short tandem repeats
#17. Differentiate between SNPs, INDELS, Homopolymers and Short tandem repeats

Pviv_rGenome_thres5@loci_table$TypeOf_Markers =
  TypeOf_Marker(Pviv_rGenome_thres5, w = 1, n = 1)

#### Step 3.2 Effect of PCR genotyping artifacts on heterozygosity

#18. Get Proportion of Heterozygous samples per site

Pviv_rGenome_thres5@loci_table$ObsHet = get_ObsHet(Pviv_rGenome_thres5, by = 'loci', w = 1, n = 1)

#19. Calculate fraction of Heterozygous samples per Alternative alleles per site

Pviv_rGenome_thres5@loci_table$frac_ofHet_pAlts =
  frac_ofHet_pAlt(Pviv_rGenome_thres5, w = 1, n = 1)

#20. Classify sites based on the proportion of alternative alleles in polyclonal samples per site

#20.1 Distribution of the fraction of alternative alleles in polyclonal samples

Pviv_rGenome_thres5@loci_table %>% ggplot(aes(x = frac_ofHet_pAlts))+
  geom_histogram(binwidth = 0.01)+
  labs(x = 'Fraction of heterozygous samples per alternative allele per site',
       y = 'Number of sites (Loci)')+
  theme_bw()

#20.2 Classify sites based on the proportion of alternative alleles in polyclonal samples per site
Pviv_rGenome_thres5@loci_table %<>% mutate(ALT_FILTER =case_when(
  frac_ofHet_pAlts == 1 ~ '100%',
  frac_ofHet_pAlts < 1 & frac_ofHet_pAlts > .5 ~ '50 - 99%',
  frac_ofHet_pAlts <= .5 ~ '<=50%',
))

#20.3 Distribution of observed heterozygosity per locus, by type of marker, by group (defined in previuos step 10.2)
Pviv_rGenome_thres5@loci_table %>%
  ggplot(aes(x = ObsHet,
             fill = factor(ALT_FILTER,
                           levels = c('<=50%', '50 - 99%', '100%'))))+
  geom_histogram(binwidth = .01)+
  scale_fill_manual(values = c('dodgerblue3', 'gold3', 'firebrick3'))+
  facet_wrap(.~factor(TypeOf_Markers, levels = c(
    'SNP',
    'INDEL',
    'INDEL:Homopolymer',
    'INDEL:Dinucleotide_STR',
    'INDEL:Trinucleotide_STR',
    'INDEL:Tetranucleotide_STR',
    'INDEL:Pentanucleotide_STR',
    'INDEL:Hexanucleotide_STR'
  )), scales = 'free_y', ncol = 4)+
  labs(y = 'Number of Loci',
       x = 'Observed Heterozygosity per locus',
       fill = 'Het/Alt')+
  theme_bw()

#21. Remove Homopolymers and DiSTRs

Pviv_rGenome_thres5 = filter_loci(Pviv_rGenome_thres5,
                                  v = !(Pviv_rGenome_thres5@loci_table$TypeOf_Markers %in%
                                          c('INDEL:Homopolymer', 'INDEL:Dinucleotide_STR')))

### Step 4: Filter of problematic genomic regions to map the reads

#### Step 4.1: Regions with excess of Heterozygosity
#22. Distribution of proportion of heterozygous per site by sites with \<= 50% Het/Alt, 50 - 99% Het/Alt and 100% of Het/Alt

# In that way we identify that maximum fraction of Heterozygous samples per loci is 0.4855491

Pviv_rGenome_thres5@loci_table %>%
  mutate(TypeOf_Markers2 = case_when(
    TypeOf_Markers == 'SNP' ~ 'SNP',
    TypeOf_Markers != 'SNP' ~ 'INDELs'),
    CHROM2  = gsub('(^(\\d|P|v)+_|_v1)', '',CHROM))%>%
  ggplot(aes(x = POS, y = ObsHet, color = factor(ALT_FILTER, levels = c('<=50%', '50 - 99%', '100%')))) +
  geom_point(alpha = 0.5, size = .25) +
  scale_color_manual(values = c('dodgerblue3', 'gold3', 'firebrick3'))+
  facet_grid(TypeOf_Markers2 ~ CHROM2)+
  theme_bw()+
  labs(y = 'Observed Heterozygosity', x = 'Chromosomal position', color = 'Het/Alt')+
  theme(legend.position = 'right',
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

#23. Calculate the average of heterozygous samples per site per gene

Pviv_rGenome_thres5@loci_table = mean_ObsHet(Pviv_rGenome_thres5, gff = '../reference/genes.gff')

#24. Define a threshold for trusted variants
ObsHet_Threshold = quantile(Pviv_rGenome_thres5@loci_table %>%
                              filter(ALT_FILTER == '<=50%', TypeOf_Markers == 'SNP') %>%
                              group_by(gene_id) %>%
                              dplyr::summarise(mean_ObsHet = max(mean_ObsHet)) %>%
                              select(mean_ObsHet) %>%
                              unlist, .95)

#25. Visualize the Average Observed Heteozygosity in SNPs per genomic position
Pviv_rGenome_thres5@loci_table %>%
  filter(ALT_FILTER == '<=50%', TypeOf_Markers == 'SNP') %>%
  group_by(gene_id) %>%
  dplyr::summarise(mean_ObsHet = max(mean_ObsHet)) %>%
  ggplot(aes(x = mean_ObsHet))+
  geom_histogram(binwidth = .005)+
  geom_vline(xintercept = ObsHet_Threshold)+
  labs(y = 'Number of genomic regions', x = 'Average Observed Heteozygosity in SNPs') +
  theme_bw()

#26. Differentiate genomic regions that PASS the threshold

Pviv_rGenome_thres5@loci_table %<>% mutate(ObsHet_Filter = case_when(
  mean_ObsHet > ObsHet_Threshold ~ 'OUT',
  mean_ObsHet <= ObsHet_Threshold ~ 'PASS'
))

#27. Visualize which genomic regions PASS the filter
Pviv_rGenome_thres5@loci_table %>%
  mutate(TypeOf_Markers2 = case_when(
    TypeOf_Markers == 'SNP' ~ 'SNP',
    TypeOf_Markers != 'SNP' ~ 'INDEL'),
    CHROM2  = gsub('(^(\\d|P|v)+_|_v1)', '',CHROM))%>%
  ggplot(aes(x = POS, y = ObsHet, color = ObsHet_Filter)) +
  geom_point(alpha = 0.5, size = .25) +
  geom_hline(yintercept = ObsHet_Threshold)+
  scale_color_manual(values = c('firebrick2', 'dodgerblue3'))+
  facet_grid(TypeOf_Markers2 ~ CHROM2)+
  theme_bw()+
  labs(y = 'Observed Heterozygosity',
       x = 'Chromosomal position',
       color = 'Mean ObsHet > Th')+
  theme(legend.position = 'right',
        axis.text.x = element_blank())

#### Step 4.2: Regions with excess of density of Polymorphism

#28. Calculate the density of SNPs per genomic region and define a threshold
Pviv_rGenome_thres5@loci_table = SNP_density(Pviv_rGenome_thres5, gff = '../reference/genes.gff')

SNP_density_threshold = quantile(Pviv_rGenome_thres5@loci_table %>%
                                   filter(TypeOf_Markers == 'SNP') %>%
                                   group_by(gene_id) %>%
                                   dplyr::summarise(SNP_density = max(SNP_density)) %>%
                                   select(SNP_density) %>%
                                   unlist, .925)

#29. Distribution of SNP density
Pviv_rGenome_thres5@loci_table %>%
  filter(TypeOf_Markers == 'SNP') %>%
  group_by(gene_id) %>%
  dplyr::summarise(SNP_density = max(SNP_density)) %>%
  ggplot(aes(x = SNP_density))+
  geom_histogram(binwidth = .001)+
  geom_vline(xintercept = SNP_density_threshold)+
  labs(y = 'Number of genomic regions', x = 'SNP density') +
  theme_bw()

#30. Differentiate/Filter genomic regions that PASS the threshold
Pviv_rGenome_thres5@loci_table %<>% 
  mutate(SNP_density_Filter = case_when(
    SNP_density > SNP_density_threshold ~ 'OUT',
    SNP_density <= SNP_density_threshold ~ 'PASS'
  ))

#31. Visualize which genomic regions PASS the filter
Pviv_rGenome_thres5@loci_table %>%
  mutate(TypeOf_Markers2 = case_when(
    TypeOf_Markers == 'SNP' ~ 'SNP',
    TypeOf_Markers != 'SNP' ~ 'INDEL'),
    CHROM2  = gsub('(^(\\d|P|v)+_|_v1)', '',CHROM))%>%
  ggplot(aes(x = POS, y = ObsHet, color = SNP_density_Filter)) +
  geom_point(alpha = 0.5, size = .25) +
  scale_color_manual(values = c('firebrick2', 'dodgerblue3'))+
  facet_grid(TypeOf_Markers2 ~ CHROM2)+
  theme_bw()+
  labs(y = 'Observed Heterozygosity', x = 'Chromosomal position', color = 'SNP density')+
  theme(legend.position = 'right',
        axis.text.x = element_blank())

#32. Differentiate of problematic genomic regions to map reads

Pviv_rGenome_thres5@loci_table %<>%
  mutate(Alignment_Filter = case_when(
    ObsHet_Filter == 'OUT' & SNP_density_Filter == 'OUT' ~ 'OUT',
    !(ObsHet_Filter == 'OUT' & SNP_density_Filter == 'OUT') ~ 'PASS'
  ))

#33. Visualize which genomic regions PASS the filter
Pviv_rGenome_thres5@loci_table %>%
  mutate(TypeOf_Markers2 = case_when(
    TypeOf_Markers == 'SNP' ~ 'SNP',
    TypeOf_Markers != 'SNP' ~ 'INDEL'),
    CHROM2  = gsub('(^(\\d|P|v)+_|_v1)', '',CHROM))%>%
  ggplot(aes(x = POS, y = ObsHet, color = Alignment_Filter)) +
  geom_point(alpha = 0.5, size = .25) +
  scale_color_manual(values = c('firebrick2', 'dodgerblue3'))+
  facet_grid(TypeOf_Markers2 ~ CHROM2)+
  theme_bw()+
  labs(y = 'Observed Heterozygosity', x = 'Chromosomal position', color = 'Alignment\nFilter')+
  theme(legend.position = 'right',
        axis.text.x = element_blank())

#34. Check which genomic regions or genes are removed
Pviv_rGenome_thres5@loci_table %>% 
  filter(Alignment_Filter == 'OUT') %>% 
  select(gene_id) %>% 
  unique()

#35. Remove problematic genomic regions to map reads
Pviv_rGenome_filtered = 
  filter_loci(Pviv_rGenome_thres5,
              v = Pviv_rGenome_thres5@loci_table$Alignment_Filter == 'PASS')

