# run_vcftools ----
run_vcftools = function(vcf = NULL,
                           gzvcf = NULL,
                           bcf = NULL,
                           out = NULL,
                           
                           keep_regexp = NULL, # regular expression pattern that identify samples
                           remove_regexp = NULL,
                           
                           keep = NULL,
                           remove = NULL,
                           
                           chr = NULL,
                           not_chr = NULL,
                           
                           bed = NULL,
                           exclude_bed = NULL,
                           
                           positions = NULL,
                           exclude_positions = NULL,
                        
                           keep_only_indels = FALSE,
                           remove_indels = FALSE,
                           
                           remove_filtered_all = FALSE,
                           
                           maf = NULL,
                           max_maf = NULL,
                           
                           non_ref_af = NULL,
                           max_non_ref_af = NULL,
                           non_ref_ac = NULL,
                           max_non_ref_ac = NULL,
                           
                           non_ref_af_any = NULL,
                           max_non_ref_af_any = NULL,
                           non_ref_ac_any = NULL,
                           max_non_ref_ac_any = NULL,
                           
                           mac = NULL,
                           max_mac = NULL,
                           
                           min_alleles = NULL,
                           max_alleles = NULL,
                           
                           # Output options
                           freq = FALSE,
                           counts = FALSE,
                           
                           depth = FALSE,
                           site_depth = FALSE,
                           site_mean_depth = FALSE,
                           geno_depth = FALSE,
                           
                           hap_r2 = FALSE,
                           geno_r2 = FALSE,
                           geno_chisq = FALSE,
                           hap_r2_positions = NULL,
                           geno_r2_positions = NULL,
                           ld_window = NULL,
                           ld_window_bp = NULL,
                           ld_window_min = NULL,
                           ld_window_bp_min = NULL,
                           min_r2 = NULL,
                           interchrom_hap_r2 = FALSE,
                           interchrom_geno_r2 = FALSE,
                           
                           TsTv = NULL,
                           TsTv_by_count = FALSE,
                           TsTv_by_qual = FALSE,
                           
                           site_pi = FALSE,
                           window_pi = NULL,
                           window_pi_step = NULL,
                           
                           weir_fst_pop = NULL,
                           fst_window_size = NULL,
                           fst_window_step = NULL,
                           
                           het = FALSE,
                           TajimaD = NULL,
                           relatedness = FALSE,
                           relatedness2 = FALSE,
                           
                           recode = FALSE,
                           recode_bcf = FALSE,
                           recode_INFO_all = FALSE
){
  
  vcf_run_file = c('#!/bin/bash',
                   'source /broad/software/scripts/useuse',
                   'use .vcftools-0.1.14',
                   'use Tabix')
  
  vcf_arguments = 'vcftools'
  
  if(!is.null(vcf)){vcf_arguments = paste0(vcf_arguments, ' --vcf ', vcf)}
  if(!is.null(gzvcf)){vcf_arguments = paste0(vcf_arguments, ' --gzvcf ', gzvcf)}
  if(!is.null(bcf)){vcf_arguments = paste0(vcf_arguments, ' --bcf ', bcf)}
  
  # Filters
  
  if(!is.null(keep_regexp)){
    system(paste0("zgrep '^#[A-Z]' ",  gzvcf, " > ", "temp_samples.indv"))
    temp_samples = as.character(read.csv("temp_samples.indv", header = FALSE, sep = '\t'))
    samples = temp_samples[grepl(keep_regexp,temp_samples)]
    write.table(samples, 'samples.indv', sep = '\t', quote = FALSE, row.names = FALSE, col.names = F)
    system(paste0('rm ', "temp_samples.indv"))
    vcf_arguments = paste0(vcf_arguments, ' --keep samples.indv')
  }
  
  if(!is.null(remove_regexp)){
    system(paste0("zgrep '^#[A-Z]' ",  gzvcf, " > ", "temp_samples.indv"))
    temp_samples = as.character(read.csv("temp_samples.indv", header = FALSE, sep = '\t'))
    samples = temp_samples[grepl(keep_regexp,temp_samples)]
    write.table(samples, 'rsamples.indv', sep = '\t', quote = FALSE, row.names = FALSE, col.names = F)
    system(paste0('rm ', "temp_samples.indv"))
    vcf_arguments = paste0(vcf_arguments, ' --remove rsamples.indv')
  }
  
  if(!is.null(keep)){
    vcf_arguments = paste0(vcf_arguments, ' --keep ', keep)
  }
  
  if(!is.null(remove)){
    vcf_arguments = paste0(vcf_arguments, ' --remove ', remove)
  }
  
  if(!is.null(chr)){
    for(chromosome in chr){
      vcf_arguments = paste0(vcf_arguments, ' --chr ', chromosome)
    }
  }
  
  if(!is.null(not_chr)){
    for(chromosome in not_chr){
      vcf_arguments = paste0(vcf_arguments, ' --not-chr ', chromosome)
    }
  }
  
  if(!is.null(bed)){vcf_arguments = paste0(vcf_arguments, ' --bed ', bed)}
  if(!is.null(exclude_bed)){vcf_arguments = paste0(vcf_arguments, ' --exclude-bed ', exclude_bed)}
  
  if(!is.null(positions)){vcf_arguments = paste0(vcf_arguments, ' --positions ', positions)}
  if(!is.null(exclude_positions)){vcf_arguments = paste0(vcf_arguments, ' --exclude-positions ', exclude_positions)}
  
  if(keep_only_indels){vcf_arguments = paste0(vcf_arguments, ' --keep-only-indels')}
  if(remove_indels){vcf_arguments = paste0(vcf_arguments, ' --remove-indels')}
  
  if(remove_filtered_all){vcf_arguments = paste0(vcf_arguments, ' --remove-filtered-all')}
  
  if(!is.null(maf)){vcf_arguments = paste0(vcf_arguments, ' --maf ', maf)}
  if(!is.null(max_maf)){vcf_arguments = paste0(vcf_arguments, ' --max-maf ', max_maf)}
  
  if(!is.null(non_ref_af)){vcf_arguments = paste0(vcf_arguments, ' --non-ref-af ', non_ref_af)}
  if(!is.null(max_non_ref_af)){vcf_arguments = paste0(vcf_arguments, ' --max-non-ref-af ', max_non_ref_af)}
  if(!is.null(non_ref_ac)){vcf_arguments = paste0(vcf_arguments, ' --non-ref-ac ', non_ref_ac)}
  if(!is.null(max_non_ref_ac)){vcf_arguments = paste0(vcf_arguments, ' --max-non-ref-ac ', max_non_ref_ac)}
  
  if(!is.null(non_ref_af_any)){vcf_arguments = paste0(vcf_arguments, ' --non-ref-af-any ', non_ref_af_any)}
  if(!is.null(max_non_ref_af_any)){vcf_arguments = paste0(vcf_arguments, ' --max-non-ref-af-any ', max_non_ref_af_any)}
  if(!is.null(non_ref_ac_any)){vcf_arguments = paste0(vcf_arguments, ' --non-ref-ac-any ', non_ref_ac_any)}
  if(!is.null(max_non_ref_ac_any)){vcf_arguments = paste0(vcf_arguments, ' --max-non-ref-ac-any ', max_non_ref_ac_any)}
  if(!is.null(mac)){vcf_arguments = paste0(vcf_arguments, ' --mac ', mac)}
  if(!is.null(max_mac)){vcf_arguments = paste0(vcf_arguments, ' --max-mac ', max_mac)}
  
  if(!is.null(min_alleles)){vcf_arguments = paste0(vcf_arguments, ' --min-alleles ', min_alleles)}
  if(!is.null(max_alleles)){vcf_arguments = paste0(vcf_arguments, ' --max-alleles ', max_alleles)}
  
  # Output options
  
  if(freq){vcf_arguments = paste0(vcf_arguments, ' --freq')}
  if(counts){vcf_arguments = paste0(vcf_arguments, ' --counts')}
  
  if(depth){vcf_arguments = paste0(vcf_arguments, ' --depth')}
  if(site_depth){vcf_arguments = paste0(vcf_arguments, ' --site-depth')}
  if(site_mean_depth){vcf_arguments = paste0(vcf_arguments, ' --site-mean-depth')}
  if(geno_depth){vcf_arguments = paste0(vcf_arguments, ' --geno_depth')}
  
  if(hap_r2){vcf_arguments = paste0(vcf_arguments, ' --hap-r2')}
  if(geno_r2){vcf_arguments = paste0(vcf_arguments, ' --geno-r2')}
  if(geno_chisq){vcf_arguments = paste0(vcf_arguments, ' --geno-chisq')}
  if(!is.null(ld_window)){vcf_arguments = paste0(vcf_arguments, ' --ld-window ', ld_window)}
  if(!is.null(ld_window_bp)){vcf_arguments = paste0(vcf_arguments, ' --ld-window-bp ', ld_window_bp)}
  if(!is.null(ld_window_min)){vcf_arguments = paste0(vcf_arguments, ' --ld-window-min ', ld_window_min)}
  if(!is.null(ld_window_bp_min)){vcf_arguments = paste0(vcf_arguments, ' --ld-window-bp-min ', ld_window_bp_min)}
  if(!is.null(min_r2)){vcf_arguments = paste0(vcf_arguments, ' --min-r2 ', min_r2)}
  if(interchrom_hap_r2){vcf_arguments = paste0(vcf_arguments, ' --interchrom-hap-r2')}
  if(interchrom_geno_r2){vcf_arguments = paste0(vcf_arguments, ' --interchrom-geno-r2')}
  
  if(!is.null(TsTv)){vcf_arguments = paste0(vcf_arguments, ' --TsTv ', TsTv)}
  if(TsTv_by_count){vcf_arguments = paste0(vcf_arguments, ' --TsTv-by-count')}
  if(TsTv_by_qual){vcf_arguments = paste0(vcf_arguments, ' --TsTv-by-qual')}
  
  if(site_pi){vcf_arguments = paste0(vcf_arguments, ' --site-pi')}
  if(!is.null(window_pi)){vcf_arguments = paste0(vcf_arguments, ' --window-pi ', window_pi)}
  if(!is.null(window_pi_step)){vcf_arguments = paste0(vcf_arguments, ' --window-pi-step ', window_pi_step)}
  
  if(!is.null(weir_fst_pop)){
    for(pop in weir_fst_pop){
      vcf_arguments = paste0(vcf_arguments, ' --weir-fst-pop ', pop)
    }
  }
  
  if(!is.null(fst_window_size)){vcf_arguments = paste0(vcf_arguments, ' --fst-window-size ', fst_window_size)}
  if(!is.null(fst_window_step)){vcf_arguments = paste0(vcf_arguments, ' --fst-window-step ', fst_window_step)}
  
  if(het){vcf_arguments = paste0(vcf_arguments, ' --het')}
  if(!is.null(TajimaD)){vcf_arguments = paste0(vcf_arguments, ' --TajimaD ', TajimaD)}
  if(relatedness){vcf_arguments = paste0(vcf_arguments, ' --relatedness')}
  if(relatedness2){vcf_arguments = paste0(vcf_arguments, ' --relatedness2')}
  
  if(recode){vcf_arguments = paste0(vcf_arguments, ' --recode')}
  if(recode_bcf){vcf_arguments = paste0(vcf_arguments, ' --recode-bcf')}
  if(recode_INFO_all){vcf_arguments = paste0(vcf_arguments, ' --recode-INFO-all')}
  
  if(!is.null(out)){vcf_arguments = paste0(vcf_arguments, ' --out ', out)}
  
  vcf_run_file = c(vcf_run_file, vcf_arguments)
  
  write.table(vcf_run_file, 'vcf_run_file.sh', row.names = FALSE, quote = FALSE, col.names = FALSE)
  
  system(paste0('chmod 777 ', 'vcf_run_file.sh'))
  system('./vcf_run_file.sh')
  
}



# load_vcf----
load_vcf = function(vcf = NULL,
                       gzvcf = NULL){
  
  if(!is.null(vcf)){
    system(paste0("grep -v '^##' ", vcf, " > temp_vfc.tsv"))
  }
  
  if(!is.null(gzvcf)){
    system(paste0("zgrep -v '^##' ", vcf, " > temp_vfc.tsv"))
  }
  
  system(paste0("grep '^#' temp_vfc.tsv > col_names.tsv"))
  
  col_names = strsplit(readLines('col_names.tsv'), '\t')[[1]]
  col_names[1] = 'CHROM'
  
  vcf = read.table('temp_vfc.tsv')
  names(vcf) = col_names
  
  system('rm temp_vfc.tsv')
  system('rm col_names.tsv')
  return(vcf)
}


# rGenome S4class and vcf2rGenome----

## rGenome S4 class

setClass('rGenome',
         representation = representation(data="environment"))

## rGenome constructor
rGenome = function(gt = NULL,
                   loci_table = NULL,
                   metadata = NULL){
  obj = new('rGenome')
  obj@data$gt = gt
  obj@data$loci_table = loci_table
  obj@data$metadata = metadata
  
  return(obj)
}

## set function vcf2rGenome

vcf2rGenome = function(vcf, n = 500, threshold = 5) {
            
            # Generate metadata
            metadata = data.frame(sample = names(vcf_object)[-1:-9])
            
            # generate loci_table
            loci_table = vcf[,c(1,2,4,5)]
            rownames(loci_table) = paste(loci_table$CHROM, loci_table$POS, sep = '_')
            
            # generate a genotype table (gt)
            
            gt = NULL
            for(w in 1:n){
              start = Sys.time()
              gt = rbind(gt, get_GTAD_matrix(vcf, w = w, n = n, threshold = threshold))
              end = Sys.time()
              print(w)
              print(end-start)
            }
            
            obj = rGenome(gt = gt, loci_table = loci_table, metadata = metadata)
          
            return(obj)
          }

# SampleAmplRate----

setGeneric("SampleAmplRate", function(obj) standardGeneric("SampleAmplRate"))

setMethod("SampleAmplRate", signature(obj = "rGenome"),
          
          function(obj) {
            
            obj@data$metadata[['SampleAmplRate']] =
              1 - colSums(is.na(obj@data$gt))/nrow(obj@data$gt)
            
            return(obj)
          }
)

# LocusAmplRate----

setGeneric("LocusAmplRate", function(obj) standardGeneric("LocusAmplRate"))

setMethod("LocusAmplRate", signature(obj = "rGenome"),
          
          function(obj) {
            
            obj@data$loci_table[['LocusAmplRate']] =
              1 - rowSums(is.na(obj@data$gt))/ncol(obj@data$gt)
            
            return(obj)
          }
)

# filter_samples----

setGeneric("filter_samples", function(obj, v = NULL) standardGeneric("filter_samples"))

setMethod("filter_samples", signature(obj = "rGenome"),
          
          function(obj, v = NULL) {
            
            obj@data$gt = obj@data$gt[,v]
            obj@data$metadata = obj@data$metadata[v,]
            
            
            return(obj)
          }
)

# filter_loci----

setGeneric("filter_loci", function(obj, v = NULL) standardGeneric("filter_loci"))

setMethod("filter_loci", signature(obj = "rGenome"),
          
          function(obj, v = NULL) {
            
            obj@data$gt = obj@data$gt[v,]
            obj@data$loci_table = obj@data$loci_table[v,]
            
            return(obj)
          }
)

# get_AC----

setGeneric("get_AC", function(obj = NULL, w = 1, n = 1,
                              update_alleles = TRUE,
                              monoclonals = NULL, polyclonals = NULL) standardGeneric("get_AC"))

setMethod("get_AC", signature(obj = "rGenome"),
          function(obj = NULL, w = 1, n = 1,
                  update_alleles = TRUE,
                  monoclonals = NULL, polyclonals = NULL){
  
            gt = obj@data$gt
            loci = obj@data$loci_table
            
            s = round(seq(1,nrow(gt)+1, length.out=n+1))
            low = s[w]
            high = s[w+1]-1
            
            gt = matrix(gsub(':\\d+', '', gt[low:high,]), nrow = high - low + 1, ncol = ncol(gt),
                        dimnames = list(rownames(gt)[low:high], colnames(gt)))
            
            if(is.null(monoclonals) & is.null(polyclonals)){
              gt1 = gsub('/\\d+', '', gt)
              gt2 = gsub('\\d+/', '', gt)
              gt2[gt2==gt1] = NA
              gt3 = cbind(gt1, gt2)
            }else{
              gt_mono = gt[,monoclonals]
              gt_mono = gsub('/\\d+', '', gt_mono)
              
              gt_poly = gt[,polyclonals]
              gt_poly1 = gsub('/\\d+', '', gt_poly)
              gt_poly2 = gsub('\\d+/', '', gt_poly)
              gt3 = cbind(gt_mono, gt_poly1, gt_poly2)
            }
            
            
            if(update_alleles){
              loci = loci[low:high,c('REF', 'ALT')]
            }
            
            alleles = t(sapply(1:nrow(gt3), function(locus){
              alleles = unique(gt3[locus,])
              alleles = sort(alleles[!is.na(alleles)])
              
              if(update_alleles){
                original_alleles = c(loci[locus,'REF'], strsplit(loci[locus,'ALT'], ',')[[1]])}
              nalleles = length(alleles)
              AC = sapply(alleles, function(allele){
                sum(gt3[locus,] == allele, na.rm = T)
              })
              
              AC = paste(paste(alleles, AC, sep = ':'), collapse = ',')
              
              if(update_alleles){
                alleles = paste(
                  paste(
                    original_alleles[as.integer(alleles)+1], alleles, sep = ':'), collapse = ',')}
              
              if(update_alleles){
                c(nalleles, alleles, AC)}else{
                  c(nalleles, AC)
                }
            }))
            
            alleles = as.data.frame(alleles)
            
            if(update_alleles){
              colnames(alleles) = c('Cardinality', 'Alleles', 'Allele_Counts')
            }else{
              colnames(alleles) = c('Cardinality', 'Allele_Counts')
            }
            rownames(alleles) = rownames(gt)
            alleles$Cardinality = as.integer(alleles$Cardinality)
            
            return(alleles)
            
          }
          )


# get_ExpHet----

setGeneric("get_ExpHet", function(obj = NULL, update_AC = FALSE, monoclonals = NULL, polyclonals = NULL) standardGeneric("get_ExpHet"))

setMethod("get_ExpHet", signature(obj = "rGenome"),
          function(obj = NULL, update_AC = FALSE, monoclonals = NULL, polyclonals = NULL){
  
  gt = obj@data$gt
  loci = obj@data$loci_table
  
  if(!update_AC){
    AC = sapply(1:nrow(loci), function(x){
      AC = strsplit(loci[x,'Allele_Counts'], ',')[[1]]
      gsub('^\\d+:', '', AC)
    })
  }else{
    AC = get_AC(obj = obj,w =1, n = 1, update_alleles = FALSE, monoclonals = monoclonals, polyclonals = polyclonals)
    AC = AC$Allele_Counts
    AC = sapply(AC, function(x){
      AC = strsplit(x, ',')[[1]]
      gsub('^\\d+:', '', AC)
    })
  }
  
  ExpHet = sapply(1:length(AC), function(pos){
    n = sum(as.numeric(AC[[pos]]))
    allele_counts = as.numeric(AC[[pos]])
    allele_freq = allele_counts/n
    sp2 = sum(allele_freq^2)
    ExpHet = n * (1 - sp2)/(n - 1)
  })
  
  return(ExpHet)
})

# get_EffCard----

setGeneric("get_EffCard", function(obj = NULL, update_AC = FALSE, monoclonals = NULL, polyclonals = NULL) standardGeneric("get_EffCard"))

setMethod("get_EffCard", signature(obj = "rGenome"),
          function(obj = NULL, update_AC = FALSE, monoclonals = NULL, polyclonals = NULL){
  
            gt = obj@data$gt
            loci = obj@data$loci_table
            
            if(!update_AC){
              AC = sapply(1:nrow(loci), function(x){
                AC = strsplit(loci[x,'Allele_Counts'], ',')[[1]]
                gsub('^\\d+:', '', AC)
              })
            }else{
              AC = get_AC(gt = gt,w =1, n = 1, update_alleles = FALSE, monoclonals = monoclonals, polyclonals = polyclonals)
              AC = AC$Allele_Counts
              AC = sapply(AC, function(x){
                AC = strsplit(x, ',')[[1]]
                gsub('^\\d+:', '', AC)
              })
            }
  
            EffCard = sapply(1:length(AC), function(pos){
              n = sum(as.numeric(AC[[pos]]))
              allele_counts = as.numeric(AC[[pos]])
              allele_freq = allele_counts/n
              sp2 = sum(allele_freq^2)
              ExpHet = n * (1 - sp2)/(n - 1)
              EffCard = 1/(1-ExpHet)
            })
            
            return(EffCard)
            }
          )

# get_ObsHet----

setGeneric("get_ObsHet", function(obj = NULL, by = 'loci', w = 1, n = 1) standardGeneric("get_ObsHet"))

setMethod("get_ObsHet", signature(obj = "rGenome"),
          function(obj = NULL, by = 'loci', w = 1, n = 1){
            
            gt = obj@data$gt
            s = round(seq(1,nrow(gt)+1, length.out=n+1))
            low = s[w]
            high = s[w+1]-1
            
            if(by == 'loci'){
              ObsHet = rowSums(matrix(grepl('/', gt), nrow = nrow(gt), ncol = ncol(gt)))/ncol(gt)
            }else if(by == 'sample'){
              ObsHet = colSums(matrix(grepl('/', gt), nrow = nrow(gt), ncol = ncol(gt)))/nrow(gt)
            }else{
              print('by argument should be loci or sample')
            }
            
            
            return(ObsHet)
          }
  
)

# Type of Marker----

setGeneric("TypeOf_Marker", function(obj, w = 1, n = 100) standardGeneric("TypeOf_Marker"))

setMethod("TypeOf_Marker", signature(obj = "rGenome"),
          function(obj, w = 1, n = 100){
            
            loci_table = obj@data$loci_table
            
            s = round(seq(1,nrow(loci_table)+1, length.out=n+1))
            low = s[w]
            high = s[w+1]-1
            
            loci_df = loci_table[low:high,c('REF', 'ALT')]
            
            alleles = apply(loci_df, 1, function(x){paste(x[1], x[2], sep = ',')})
            
            logical_vector = sapply(alleles,
                                    function(site){
                                      prod(nchar(strsplit(site, ',')[[1]]) == 1)
                                    })
            
            del_vector = sapply(alleles,
                                function(site){
                                  grepl('\\*', site)
                                })
            
            type_of_marker = ifelse(logical_vector == 1 & del_vector == FALSE, 'SNP', 'INDEL')
            
            indels = alleles[type_of_marker == 'INDEL']
            
            homopolymers = sapply(indels, function(site){
              ((length(strsplit(paste(gsub('^.','',strsplit(site, ',')[[1]]), collapse = ''), '')[[1]]) > 1)*
                 length(unique(strsplit(paste(gsub('^.','',strsplit(site, ',')[[1]]), collapse = ''), '')[[1]])) == 1|
                 (length(strsplit(paste(gsub('.$','',strsplit(site, ',')[[1]]), collapse = ''), '')[[1]]) > 1)*
                 length(unique(strsplit(paste(gsub('.$','',strsplit(site, ',')[[1]]), collapse = ''), '')[[1]])) == 1)
            })
            
            type_of_marker[type_of_marker == 'INDEL'] = ifelse(homopolymers, 'INDEL:Homopolymer', 'INDEL')
            
            # short tandem repeats
            
            indels = alleles[type_of_marker == 'INDEL']
            
            STRs = sapply(indels, function(site){
              str = strsplit(site, ',')[[1]]
              
              str = gsub('^.', '', str[which.max(sapply(str, function(allele) {nchar(allele)}))])
              str2 = gsub('.$', '', str[which.max(sapply(str, function(allele) {nchar(allele)}))])
              
              indel = NULL
              
              if(nchar(str) < 2){
                indel = 'INDEL'
              }
              
              if(nchar(str) > 2 & nchar(str)%%2 == 0 & is.null(indel)){
                dinucletide = substring(str, seq(1, nchar(str), 2), seq(1, nchar(str), 2) + 1)
                dinucletide2 = substring(str2, seq(1, nchar(str2), 2), seq(1, nchar(str2), 2) + 1)
                if(prod(grepl(dinucletide[1], dinucletide)) == 1 | prod(grepl(dinucletide2[1], dinucletide2)) == 1){
                  indel = 'INDEL:Dinucleotide_STR'
                }
              }
              
              if(nchar(str) > 3 & nchar(str)%%3 == 0 & is.null(indel)){
                trinucletide = substring(str, seq(1, nchar(str), 3), seq(1, nchar(str), 3) + 2)
                trinucletide2 = substring(str2, seq(1, nchar(str2), 3), seq(1, nchar(str2), 3) + 2)
                if(prod(grepl(trinucletide[1], trinucletide)) == 1| prod(grepl(trinucletide2[1], trinucletide2)) == 1)
                  indel = 'INDEL:Trinucleotide_STR'
              }
              
              if(nchar(str) > 4 & nchar(str)%%4 == 0 & is.null(indel)){
                tetranucletide = substring(str, seq(1, nchar(str), 4), seq(1, nchar(str), 4) + 3)
                
                tetranucletide2 = substring(str2, seq(1, nchar(str2), 4), seq(1, nchar(str2), 4) + 3)
                if(prod(grepl(tetranucletide[1], tetranucletide)) == 1| prod(grepl(tetranucletide2[1], tetranucletide2)) == 1){
                  indel = 'INDEL:Tetranucleotide_STR'
                }
              }
              
              if(nchar(str) > 5 & nchar(str)%%5 == 0 & is.null(indel)){
                pentanucletide = substring(str, seq(1, nchar(str), 5), seq(1, nchar(str), 5) + 4)
                pentanucletide2 = substring(str, seq(1, nchar(str2), 5), seq(1, nchar(str2), 5) + 4)
                if(prod(grepl(pentanucletide[1], pentanucletide)) == 1 | prod(grepl(pentanucletide2[1], pentanucletide2)) == 1){
                  indel = 'INDEL:Pentanucleotide_STR'
                }
              }
              
              if(nchar(str) > 6 & nchar(str)%%6 == 0 & is.null(indel)){
                hexanucletide = substring(str, seq(1, nchar(str), 6), seq(1, nchar(str), 6) + 5)
                hexanucletide2 = substring(str2, seq(1, nchar(str2), 6), seq(1, nchar(str2), 6) + 5)
                if(prod(grepl(hexanucletide[1], hexanucletide)) == 1 | prod(grepl(hexanucletide2[1], hexanucletide2)) == 1){
                  indel = 'INDEL:Hexanucleotide_STR'
                }
              }
              
              if(is.null(indel)){
                indel = 'INDEL'
              }
              
              indel
              
            })
            
            type_of_marker[type_of_marker == 'INDEL'] = STRs
            
            return(type_of_marker)
            
          }
)

# Fraction of heterozygous samples per alternative allele per site----

setGeneric("frac_ofHet_pAlt", function(obj = NULL, w = 1, n = 1) standardGeneric("frac_ofHet_pAlt"))

setMethod("frac_ofHet_pAlt", signature(obj = "rGenome"),
          
          function(obj = NULL, w = 1, n = 1){
            
            gt = obj@data$gt
            loci = obj@data$loci_table
            
            s = round(seq(1,nrow(gt)+1, length.out=n+1))
            low = s[w]
            high = s[w + 1] - 1
            
            alt = gsub('^\\d+,', '', gsub('(\\w+|\\*):', '', loci[low:high, 'Alleles']))
            gt = gsub(':\\d+', '',gt[low:high,])
            
            
            HetPos = matrix(grepl('/', gt), ncol = ncol(gt), nrow = nrow(gt))
            
            frac_ofHet_pAlt = sapply(1:nrow(gt), function(variant) {
              temp_gts = gt[variant,]
              alleles = strsplit(alt[variant], ',')[[1]]
              
              genotypes = sapply(alleles,
                                 function(allele){
                                   genotypes = grepl(allele, temp_gts)})
              
              het_genotypes = (genotypes == 1 & HetPos[variant,] == 1)    
              
              sum(het_genotypes, na.rm = T)/sum(genotypes, na.rm = T)
            })
            
            return(frac_ofHet_pAlt)
            
          }
          )

# Get genotype matrix----

get_GT_matrix = function(vcf, w = 1, n = 1, threshold = 5){
  
  s = round(seq(1,nrow(vcf)+1, length.out=n+1))
  low = s[w]
  high = s[w+1]-1
  
  w_gt_table = vcf[low:high,-1:-8]
  
  gt = t(sapply(1:nrow(w_gt_table), function(variant) {
    gt_pos = grep('GT',strsplit(w_gt_table[variant,1], ':')[[1]])
    ad_pos = grep('AD',strsplit(w_gt_table[variant,1], ':')[[1]])
    temp_gts = w_gt_table[variant,-1]
    
    sapply(1:length(temp_gts), function(sample){
      gt = strsplit(strsplit(as.character(temp_gts[sample]), ':')[[1]][gt_pos], '/')[[1]]
      ad = strsplit(strsplit(as.character(temp_gts[sample]), ':')[[1]][ad_pos], ',')[[1]]
      
      if(gt[1] != '.'){
        ad = ad[as.numeric(gt) + 1]
        gt_df = data.frame(gt = gt, ad = as.numeric(ad))
        gt_df = gt_df[order(gt_df$ad, decreasing = T),]
        gt_df = gt_df[gt_df$ad >= threshold,]
        gt = gt_df$gt
        gt = unique(gt)
        
        if(length(gt) > 0){
          gt = paste(gt, collapse = '/')
        }else{
          gt = NA
        }
        
      }else{
        gt = NA
      }
      
    })
  }))
  
  rownames(gt) = vcf[low:high,1:2] %>% mutate(Locus = paste(CHROM, POS, sep = '_')) %>% select(Locus) %>% unlist
  colnames(gt) = colnames(vcf)[-1:-9]
  
  return(gt)
  
}

get_GTAD_matrix = function(vcf, w = 1, n = 100, threshold = 5){
  start = Sys.time()
  s = round(seq(1,nrow(vcf)+1, length.out=n+1))
  low = s[w]
  high = s[w+1]-1
  
  w_gt_table = vcf[low:high,-1:-9]
  
  gt_temp = gsub(':.+', '', as.matrix(w_gt_table))
  ad_temp = gsub(':.+', '', gsub('^\\d+/\\d+:', '', as.matrix(w_gt_table)))
  
  # gt_allele1 = gsub('/(\\d|\\.)+', '', as.matrix(gt_temp))
  # gt_allele2 = gsub('^(\\d|\\.)/', '', as.matrix(gt_temp))
  # 
  
  # gt = matrix(NA, nrow = nrow(gt_temp), ncol = ncol(gt_temp))
  # 
  # gt[gt_allele1 == gt_allele2] = gt_allele1[gt_allele1 == gt_allele2]
  # 
  # ad = matrix(NA, nrow = nrow(gt_temp), ncol = ncol(gt_temp))
  

  gt = t(sapply(1:nrow(gt_temp), function(variant) {
    sapply(1:ncol(gt_temp), function(sample){
      gt = strsplit(as.character(gt_temp[variant, sample]), '/')[[1]]
      ad = strsplit(as.character(ad_temp[variant, sample]), ',')[[1]]
      
      if(gt[1] != '.'){
        ad = ad[as.numeric(gt) + 1]
        gt_df = data.frame(gt = gt, ad = as.numeric(ad))
        gt_df = gt_df[order(gt_df$ad, decreasing = T),]
        gt_df = gt_df[gt_df$ad >= threshold,]
        gt = paste(gt_df$gt, gt_df$ad, sep = ':')
        gt = unique(gt)
        
        if(length(gt) > 0){
          gt = paste(gt, collapse = '/')
        }else{
          gt = NA
        }
        
      }else{
        gt = NA
      }
      
    })
  }))
  end = Sys.time()
  end - start
  
  rownames(gt) = vcf[low:high,1:2] %>% mutate(Locus = paste(CHROM, POS, sep = '_')) %>% select(Locus) %>% unlist
  colnames(gt) = colnames(vcf)[-1:-9]
  
  return(gt)
  
}



# Calculate the average read depth by gene----

mean_ReadDepth = function(data, gff = 'genes.gff'){
  
  ref_gff = ape::read.gff('genes.gff')
  coding_regions = ref_gff[grepl('gene', ref_gff$type)&
                                   !grepl('^Transfer',ref_gff$seqid),c('seqid', 'start', 'end', 'attributes')]
  
  coding_regions$attributes = gsub('ID=','',str_extract(coding_regions$attributes, 'ID=PVP01_([0-9]+|MIT[0-9]+|API[0-9]+)'))
  
  coding_regions = coding_regions[order(coding_regions$start),]
  coding_regions = coding_regions[order(coding_regions$seqid),]
  rownames(coding_regions) = 1:nrow(coding_regions)
  
  coding_regions$mean_DP = NA
  
  for(gene in 1:nrow(coding_regions)){
    coding_regions[gene, ][['mean_DP']] = mean(data[data$CHROM == coding_regions[gene, ][['seqid']]&
            data$POS >= coding_regions[gene, ][['start']]&
            data$POS <= coding_regions[gene, ][['end']],][['DP']])
    
  }
  
  coding_regions = coding_regions[!is.na(coding_regions$mean_DP),]
  
  data$mean_DP = NA
  data$gene_id = NA
  
  for(gene in 1:nrow(coding_regions)){
    data[data$CHROM == coding_regions[gene, ][['seqid']]&
            data$POS >= coding_regions[gene, ][['start']]&
            data$POS <= coding_regions[gene, ][['end']],][['gene_id']] = coding_regions[gene, ][['attributes']]
    
    data[data$CHROM == coding_regions[gene, ][['seqid']]&
           data$POS >= coding_regions[gene, ][['start']]&
           data$POS <= coding_regions[gene, ][['end']],][['mean_DP']] = coding_regions[gene, ][['mean_DP']]
  }
  
  return(data)
}

# mean_ObsHet ----

setGeneric("mean_ObsHet", function(obj = NULL, gff = 'genes.gff') standardGeneric("mean_ObsHet"))

setMethod("mean_ObsHet", signature(obj = "rGenome"),
          
          function(obj = NULL, gff = 'genes.gff'){
  
            data = obj@data$loci_table
            
            ref_gff = ape::read.gff(gff)
            coding_regions = ref_gff[grepl('gene', ref_gff$type)&
                                       !grepl('^Transfer',ref_gff$seqid),
                                     c('seqid', 'start', 'end', 'attributes')]
            
            coding_regions$attributes = gsub('ID=','',str_extract(coding_regions$attributes, 'ID=PVP01_([0-9]+|MIT[0-9]+|API[0-9]+)'))
            
            coding_regions = coding_regions[order(coding_regions$start),]
            coding_regions = coding_regions[order(coding_regions$seqid),]
            rownames(coding_regions) = 1:nrow(coding_regions)
            
            coding_regions$mean_ObsHet = NA
            
            for(gene in 1:nrow(coding_regions)){
              coding_regions[gene, ][['mean_ObsHet']] = mean(data[data$CHROM == coding_regions[gene, ][['seqid']]&
                                                                    data$POS >= coding_regions[gene, ][['start']]&
                                                                    data$POS <= coding_regions[gene, ][['end']]&
                                                                    data$TypeOf_Marker == 'SNP',][['ObsHet']])
              
            }
            
            coding_regions2 = coding_regions[is.na(coding_regions$mean_ObsHet),]
            coding_regions = coding_regions[!is.na(coding_regions$mean_ObsHet),]
            
            data$mean_ObsHet = NA
            data$gene_id = NA
            
            for(gene in 1:nrow(coding_regions)){
              data[data$CHROM == coding_regions[gene, ][['seqid']]&
                     data$POS >= coding_regions[gene, ][['start']]&
                     data$POS <= coding_regions[gene, ][['end']],][['gene_id']] = coding_regions[gene, ][['attributes']]
              
              data[data$CHROM == coding_regions[gene, ][['seqid']]&
                     data$POS >= coding_regions[gene, ][['start']]&
                     data$POS <= coding_regions[gene, ][['end']],][['mean_ObsHet']] = coding_regions[gene, ][['mean_ObsHet']]
            }
            
            data[is.na(data$mean_ObsHet), ][['mean_ObsHet']] = 0
            
            for(pos in rownames(data[is.na(data$gene_id), ])){
              
              attrib = coding_regions2[coding_regions2[['seqid']] == data[pos,][['CHROM']] &
                                         coding_regions2[['start']] <= data[pos,][['POS']] &
                                         coding_regions2[['end']] >= data[pos,][['POS']],][['attributes']]
              
              data[pos,][['gene_id']] = ifelse(length(attrib) == 0, NA, attrib)
              
            }
            
            return(data)
          }
  
)

# Polymorphims density per target----

setGeneric("SNP_density", function(obj = NULL, gff = 'genes.gff') standardGeneric("SNP_density"))

setMethod("SNP_density", signature(obj = "rGenome"),
          
          function(obj = NULL, gff = 'genes.gff'){
            
            data = obj@data$loci_table
            
            ref_gff = ape::read.gff(gff)
            coding_regions = ref_gff[grepl('gene', ref_gff$type)&
                                       !grepl('^Transfer', ref_gff$seqid), c('seqid', 'start', 'end', 'attributes')]
            
            coding_regions$attributes = gsub('ID=', '', str_extract(coding_regions$attributes, 'ID=PVP01_([0-9]+|MIT[0-9]+|API[0-9]+)'))
            
            coding_regions = coding_regions[order(coding_regions$start),]
            coding_regions = coding_regions[order(coding_regions$seqid),]
            rownames(coding_regions) = 1:nrow(coding_regions)
            
            coding_regions$SNP_density = NA
            
            for(gene in 1:nrow(coding_regions)){
              coding_regions[gene, ][['SNP_density']] = nrow(data[data$CHROM == coding_regions[gene, ][['seqid']]&
                                                                    data$POS >= coding_regions[gene, ][['start']]&
                                                                    data$POS <= coding_regions[gene, ][['end']]&
                                                                    data$TypeOf_Marker == 'SNP',])/(coding_regions[gene, ][['end']] - coding_regions[gene, ][['start']] + 1)
              
            }
            
            data$SNP_density = NA
            
            for(gene in unique(data[!is.na(data$gene_id),][['gene_id']])){
              data[data$gene_id == gene &
                     !is.na(data$gene_id),][['SNP_density']] = coding_regions[coding_regions$attributes == gene, ][['SNP_density']]
            }
            
            data[is.na(data$SNP_density), ][['SNP_density']] = 0
            
            return(data)
          }
)

# get_allReaddepth----

get_allReaddepth = function(vcf, w = 1, n = 100){
  
  s = round(seq(1,nrow(vcf)+1, length.out=n+1))
  low = s[w]
  high = s[w+1]-1
  
  w_gt_table = vcf[low:high,-1:-8]
  
  allele_depth = NULL
  
  test = gsub('', '', w_gt_table)
  
  for(variant in 1:nrow(w_gt_table)) {
    gt_pos = grep('GT',strsplit(w_gt_table[variant,1], ':')[[1]])
    ad_pos = grep('AD',strsplit(w_gt_table[variant,1], ':')[[1]])
    temp_gts = w_gt_table[variant,-1]
    
    for(sample in 1:length(temp_gts)){
      gt = strsplit(strsplit(as.character(temp_gts[sample]),':')[[1]][gt_pos], '/')[[1]]
      ad = strsplit(strsplit(as.character(temp_gts[sample]),':')[[1]][ad_pos], ',')[[1]]
      
      gt = unique(gt)
      
      if(gt[1] == '.'){
        ad_df = data.frame(vcf_pos = variant + low - 1,
                           sample = names(temp_gts)[sample],
                           allele = NA,
                           allele_depth = NA,
                           ap_to_major = NA,
                           typeof_gt = NA,
                           typeof_all = NA)
      }else{
        
        ad_df = data.frame(vcf_pos = variant + low - 1,
                           sample = names(temp_gts)[sample],
                           allele = as.character(0:(length(ad)-1)),
                           allele_depth = as.numeric(ad))
        
        ad_df = ad_df[ad_df$allele_depth != 0,]
        
        if(nrow(ad_df) == 0){
          
          ad_df = data.frame(vcf_pos = variant + low - 1,
                             sample = names(temp_gts)[sample],
                             allele = 'uninformative alleles',
                             allele_depth = 0,
                             ap_to_major = NA,
                             typeof_gt = NA,
                             typeof_all = NA)
          
        }else{
          
          ad_df = ad_df[order(ad_df$allele_depth, decreasing = T),]
          
          ad_df$ap_to_major = ad_df$allele_depth/max(ad_df$allele_depth)
          
          ad_df$typeof_gt = ifelse(length(gt) == 1, 'Homozygous', 'Heterozygous')
          
          ad_df$typeof_all = NA
          
          if(length(gt) == 1){
            
            ad_df[1, ][['typeof_all']] = 'Homozygous'
            
            if(nrow(ad_df[-1,]) > 0 ){
              
              ad_df[-1,][['typeof_all']] = 'Excluded minor allele'
              
            }
            
          }else if(length(gt) == 2){
            
            ad_df[1, ][['typeof_all']] = 'Major allele'
            ad_df[2, ][['typeof_all']] = 'Minor allele'
            
            if(nrow(ad_df[-1:-2,]) > 0 ){
              
              ad_df[-1:-2,][['typeof_all']] = 'Excluded minor allele'
              
            }
          }
          
        }
        
      }
    
      allele_depth = rbind(allele_depth, ad_df)
      
      }
  }
  
  return(allele_depth)
  
}


# Fws----

setGeneric("get_Fws", function(obj = NULL, w = 1, n = 1) standardGeneric("get_Fws"))

setMethod("get_Fws", signature(obj = "rGenome"),
          function(obj = NULL, w = 1, n = 1){
  
            gt = obj@data$gt
            loci = obj@data$loci_table
            
            s = round(seq(1,nrow(gt)+1, length.out=n+1))
            low = s[w]
            high = s[w+1]-1
            
            gt = gt[low:high,]
            
            ExpHet = get_ExpHet(obj = obj, update_AC = TRUE)
            
            Hw = sapply(1:ncol(gt), function(sample){
              
              samp_alleles= gsub(':\\d+', '', gt[,sample])
              samp_allcounts = gsub('\\d+:', '', gt[,sample])
              
              samp_alleles1 = gsub('/\\d+$', '', samp_alleles)
              samp_alleles2 = gsub('^\\d+/', '', samp_alleles)
              
              samp_check = samp_alleles1 != samp_alleles2
              
              samp_alleles2[!samp_check] = NA
              
              samp_allcounts1 = gsub('/\\d+$', '', samp_allcounts)
              samp_allcounts2 = gsub('^\\d+/', '', samp_allcounts)
              
              samp_allcounts2[!samp_check] = NA
              
              samp_allcountsT = rowSums(cbind(as.integer(samp_allcounts1), as.integer(samp_allcounts2)), na.rm = T)
              
              samp_allfreq = cbind(as.integer(samp_allcounts1), as.integer(samp_allcounts2))/samp_allcountsT
              
              Hw = 1 - rowSums(samp_allfreq^2, na.rm = T)
              Hw[Hw==1] = NA
              Hw
            })
            
            Fws = 1 - (Hw/ExpHet)
            
            colMeans(Fws, na.rm = T)
            
          }
)

# get_genclon----

get_genclone = function(gt,
                    loci,
                    monoclonals,
                    polyclonals,
                    exclude_indels = T,
                    window = 150,
                    metadata){
  
  library(poppr)
  library(pegas)
  
  if(exclude_indels){
    gt = gt[loci$TypeOf_Markers == 'SNP',]
    loci %<>% filter(TypeOf_Markers == 'SNP')
    
  }
  
  gt = matrix(gsub(':\\d+', '', gt), nrow = nrow(gt), ncol = ncol(gt),
              dimnames = list(rownames(gt), colnames(gt)))
  
  if(is.null(monoclonals) & is.null(polyclonals)){
    gt1 = gsub('/\\d+', '', gt)
    gt2 = gsub('\\d+/', '', gt)
    gt2[gt2==gt1] = NA
    gt3 = cbind(gt1, gt2)
  }else{
    gt_mono = gt[,monoclonals]
    gt_mono = gsub('/\\d+', '', gt_mono)
    
    gt_poly = gt[,polyclonals]
    gt_poly1 = gsub('/\\d+', '', gt_poly)
    gt_poly2 = gsub('\\d+/', '', gt_poly)
    gt3 = cbind(gt_mono, gt_poly1, gt_poly2)
  }
  
  loci = loci[,c('CHROM', 'POS')]
  
  chrom_length = loci %>% group_by(CHROM) %>% summarise(length = max(POS))
  
  chrom_intervals = sapply(chrom_length$length, function(chrom){
    seq(1, chrom, window)
  })
  
  
  chrom_win = NULL
  
  for(chrom in 1:length(chrom_intervals)){
    chrom_win = rbind(chrom_win, data.frame(CHROM = chrom_length[chrom,][['CHROM']],
                                            start = chrom_intervals[[chrom]],
                                            end = chrom_intervals[[chrom]] - 1 + 150))
    
  }
  options(scipen=999)
  
  window_matrix = NULL
  
  for(i in 1:window){
    window_matrix = cbind(window_matrix,paste(chrom_win[['CHROM']], as.character(chrom_win[['start']] + i - 1), sep = '_'))
  }
  
  options(scipen=0)
  
  chrom_win$nVar = rowSums(matrix(window_matrix %in% rownames(loci),
                                  ncol = ncol(window_matrix),
                                  nrow = nrow(window_matrix)), na.rm = T)
  
  chrom_win %<>% filter(nVar != 0)
  
  chrom_win %<>% mutate(Filter = nVar >=3)
  
  Filtered_pos  = rep(chrom_win[['Filter']], chrom_win[['nVar']])
  
  chrom_win%<>%filter(nVar >=3)
  
  loci = loci[Filtered_pos,]
  
  gt3 = gt3[Filtered_pos,]
  
  collapsed_gt3 = NULL
  
  for(hap in 1:ncol(gt3)){
    start = Sys.time()
    collapsed_gt3 = rbind(collapsed_gt3,sapply(1:nrow(chrom_win), function(win){
      
      if(win > 1){
        bin = 1:chrom_win[win,][['nVar']] + sum(chrom_win[1:win - 1,][['nVar']])
      }else{
        bin = 1:chrom_win[win,][['nVar']]
      }
      paste(gt3[bin,hap], collapse = '')
    }))
    
    end = Sys.time()
    print(paste(hap, 'in', end - start))
  }
  
  collapsed_gt3[grepl('NA',collapsed_gt3)] = NA
  
  colnames(collapsed_gt3) = paste(chrom_win$CHROM, chrom_win$start, chrom_win$end, sep = '_')
  rownames(collapsed_gt3) = colnames(gt3)
  
  loci_format = pegas::as.loci(collapsed_gt3)
  
  genind_format = loci2genind(loci_format, ploidy = 1)
  
  genclone_format = as.genclone(genind_format)
  
  genclone_format@strata = data.frame(sample = rownames(collapsed_gt3))
  
  genclone_format@strata = left_join(genclone_format@strata,
                                     metadata,
                                     by = 'sample')
  
  return(genclone_format)
  
}

# filter_gt_matrix ----

filter_gt_matrix = function(gt, # genotype matrix
                            loci, # table with loci information
                            filter_table, # Table of filtered segments or positions
                            by = 'segments',
                            keep = TRUE,
                            exclude_indels = TRUE){
  
  if(exclude_indels){
    gt = gt[loci$TypeOf_Markers == 'SNP',]
    loci %<>% filter(TypeOf_Markers == 'SNP')
  }
  
  
  positions = NULL
  
  options(scipen=999)
  
  for(pos in 1:nrow(filter_table)){
    positions = c(positions, paste(filter_table[pos,][['CHROM']],
                                   filter_table[pos,][['start']]:filter_table[pos,][['end']], sep = '_'))
    
  }
  
  options(scipen=0)
  
  rownames(loci) %in% positions
 
  gt = gt[rownames(loci) %in% positions,]
  
  return(gt)
  
}

# handle_ploidy----

handle_ploidy = function(gt, monoclonals, polyclonals){
  
  gt = matrix(gsub(':\\d+', '', gt), nrow = nrow(gt), ncol = ncol(gt),
              dimnames = list(rownames(gt), colnames(gt)))
  
  if(is.null(monoclonals) & is.null(polyclonals)){
    gt1 = gsub('/\\d+', '', gt)
    gt2 = gsub('\\d+/', '', gt)
    gt2[gt2==gt1] = NA
    gt3 = cbind(gt1, gt2)
  }else{
    gt_mono = gt[,monoclonals]
    gt_mono = gsub('/\\d+', '', gt_mono)
    
    gt_poly = gt[,polyclonals]
    gt_poly1 = gsub('/\\d+', '', gt_poly)
    gt_poly2 = gsub('\\d+/', '', gt_poly)
    gt3 = cbind(gt_mono, gt_poly1, gt_poly2)
  }
  
  return(gt3)
  
}


# fastGRMcpp----
#' C++ implementation of a Genomic relationship matrix 'GRM'
#' 
#' @param X Matrix of the type 'MatrixXd' for which the GRM will be calculated.
#' 
#' @return Genomic relationship matrix (GRM).
#' 
#' @importFrom Rdpack reprompt
#' @references https://doi.org/10.1016/j.ajhg.2010.11.011
#' 
#' @examples
#' require(fastGRM)
#' Data = matrix(sample(0:1, 9000, TRUE, c(.9,.1)), 90)
#' X = grm(Data)
#' 
#' @export
#' 

grm = function(X){
  grmCpp(X)
}


#' C++ implementation of a fast singular value decomposition (SVD)
#' 
#' @param X Symmetric matrix of the type 'MatrixXd' for which the SVD will be calculated.
#' @param k Number of first k eigen vectors to return
#' @param q Auxiliary exponent
#' 
#' @return SVD matrix of size .
#' 
#' @importFrom Rdpack reprompt
#' @references https://doi.org/10.48550/arXiv.0909.4061
#' 
#' @examples
#' require(fastGRM)
#' Data = matrix(sample(0:1, 9000, TRUE, c(.9,.1)), 90)
#' X = grm(Data)
#' V = fastSVD(X, 2)
#' 
#' @export

fastSVD = function(X, k, q = 2){
  fastSVDCpp(X, k, q)
}

#' C++ implementation of a fast GRM function
#' 
#' @param X Matrix of the type 'MatrixXd' for which the fastGRM matrix will be calculated.
#' @param k Number of first k eigen vectors to return
#' @param q Auxiliary exponent
#' 
#' @return SVD matrix of size .
#' 
#' @importFrom Rdpack reprompt
#' @references https://doi.org/10.1016/j.ajhg.2010.11.011
#' @references https://doi.org/10.48550/arXiv.0909.4061
#' 
#' @examples
#' require(fastGRM)
#' Data = matrix(sample(0:1, 9000, TRUE, c(.9,.1)), 90)
#' X = fastGRM(Data, 2)
#' 
#' @export
#' 

fastGRM = function(X, k, q = 2){
  fastGRMCpp(X, k, q)
}




# merge_rGenome----

# Locus_info_temp1 = data.frame(locus_id = rownames(Locus_info), Locus_info)
# names(Locus_info_temp1) = c('locus_id',paste0(names(Locus_info_temp1)[-1], '.temp1'))
# 
# Locus_info_temp2 = data.frame(locis_id = rownames(Locus_info_Pv4), Locus_info_Pv4)
# names(Locus_info_temp2) = c('locus_id',paste0(names(Locus_info_temp2)[-1], '.temp2'))
# 
# gt_temp1 = genotypes_AD
# gt_temp2 = genotypes_AD_Pv4
# 
# Locus_info_merged = merge(Locus_info_temp1, Locus_info_temp2, by = 'locus_id', all = T)
# 
# 
# Locus_info_merged$ALT = sapply(1:nrow(Locus_info_merged), function(pos){
#   ALT = unique(c(unlist(str_split(gsub(':\\d+','',Locus_info_merged[pos,][['Alleles.temp1']]), ',', simplify = T)),
#            unlist(str_split(gsub(':\\d+','',Locus_info_merged[pos,][['Alleles.temp2']]), ',', simplify = T))))
#   
#   ALT = ALT[!is.na(ALT) & ALT != Locus_info_merged[pos,][['REF.temp1']]]
# 
#   
#   paste(ALT, collapse = ',')
#   
# }, simplify = T)
# 
# Locus_info_merged$Alleles = sapply(1:nrow(Locus_info_merged), function(pos){
#   paste(paste(c(Locus_info_merged[pos,][['REF.temp1']],
#                 str_split(Locus_info_merged[pos,][['ALT']], ',', simplify = T)) ,
#               0:(length(c(Locus_info_merged[pos,][['REF.temp1']],
#                          str_split(Locus_info_merged[pos,][['ALT']], ',', simplify = T))) - 1), sep = ':'), collapse = ',')
# })
# 
# View(Locus_info_merged[, c("Alleles.temp1", "Alleles.temp2", "Alleles")])
# 
# for(locus_id in Locus_info_merged$locus_id){
#   
#   if(!is.na(Locus_info_merged[Locus_info_merged$locus_id == locus_id,][['Alleles.temp2']])){
#    
#     if(Locus_info_merged[Locus_info_merged$locus_id == locus_id,][['Alleles.temp1']] != Locus_info_merged[Locus_info_merged$locus_id == locus_id,][['Alleles.temp2']]){
#       
#       Alleles = Locus_info_merged[Locus_info_merged$locus_id == locus_id,][['Alleles']]
#       Alleles = data.frame(Alleles = t(gsub(':\\d+','',str_split(Alleles, ',', simplify = T))),
#                            Codes = t(gsub('([A-Z]+|\\*):','',str_split(Alleles, ',', simplify = T))))
#       
#       if(locus_id %in% rownames(gt_temp1)){
#         Alleles.temp1 = Locus_info_merged[Locus_info_merged$locus_id == locus_id,][['Alleles.temp1']]
#         Alleles.temp1 = data.frame(Alleles = t(gsub(':\\d+','',str_split(Alleles.temp1, ',', simplify = T))),
#                                    Codes = t(gsub('([A-Z]+|\\*):','',str_split(Alleles.temp1, ',', simplify = T))))
#         
#         Alleles.temp1 = Alleles.temp1[order(Alleles.temp1$Codes, decreasing = T),]
#         rownames(Alleles.temp1) = 1:nrow(Alleles.temp1)
#         
#         for(allele in Alleles.temp1$Alleles){
#           if(Alleles.temp1[Alleles.temp1$Alleles == allele, ][['Codes']] != Alleles[Alleles$Alleles == allele, ][['Codes']]){
#             
#             gt_temp1[locus_id,] = gsub(paste0('^',Alleles.temp1[Alleles.temp1$Alleles == allele, ][['Codes']],':'),
#                                        paste0(Alleles[Alleles$Alleles == allele, ][['Codes']],':')
#                                        , gt_temp1[locus_id,])
#             
#             gt_temp1[locus_id,] = gsub(paste0('/',Alleles.temp1[Alleles.temp1$Alleles == allele, ][['Codes']],':'),
#                                        paste0('/',Alleles[Alleles$Alleles == allele, ][['Codes']],':')
#                                        , gt_temp1[locus_id,])
#           }
#         }
#         
#       }
#       
#       if(locus_id %in% rownames(gt_temp2)){
#         Alleles.temp2 = Locus_info_merged[Locus_info_merged$locus_id == locus_id,][['Alleles.temp2']]
#         Alleles.temp2 = data.frame(Alleles = t(gsub(':\\d+','',str_split(Alleles.temp2, ',', simplify = T))),
#                                    Codes = t(gsub('([A-Z]+|\\*):','',str_split(Alleles.temp2, ',', simplify = T))))
#         
#         Alleles.temp2 = Alleles.temp2[order(Alleles.temp2$Codes, decreasing = T),]
#         rownames(Alleles.temp2) = 1:nrow(Alleles.temp2)
#         
#         for(allele in Alleles.temp2$Alleles){
#           if(Alleles.temp2[Alleles.temp2$Alleles == allele, ][['Codes']] != Alleles[Alleles$Alleles == allele, ][['Codes']]){
#             
#             gt_temp2[locus_id,] = gsub(paste0('^',Alleles.temp2[Alleles.temp2$Alleles == allele, ][['Codes']],':'),
#                                        paste0(Alleles[Alleles$Alleles == allele, ][['Codes']],':')
#                                        , gt_temp2[locus_id,])
#             
#             gt_temp2[locus_id,] = gsub(paste0('/',Alleles.temp2[Alleles.temp2$Alleles == allele, ][['Codes']],':'),
#                                        paste0('/',Alleles[Alleles$Alleles == allele, ][['Codes']],':')
#                                        , gt_temp2[locus_id,])
#             
#           }
#         }
#       }
#     }
#   }
# }
# 
# 
# vcf[vcf$CHROM == Locus_info_merged[Locus_info_merged$locus_id == locus_id,][['CHROM.temp1']]&
#       vcf$POS == Locus_info_merged[Locus_info_merged$locus_id == locus_id,][['POS.temp1']], 1:8]
# 
# vcf_Pv4[vcf_Pv4$CHROM == Locus_info_merged[Locus_info_merged$locus_id == locus_id,][['CHROM.temp1']]&
#       vcf_Pv4$POS == Locus_info_merged[Locus_info_merged$locus_id == locus_id,][['POS.temp1']], 1:8]
# 


# filter_loci_table----

filter_loci_table = function(loci_table = NULL, # table with loci information
                             vcf_object = NULL, # table with loci information in vcf format
                            filter_table, # Table of filtered segments or positions
                            by = 'segments'){
  
  
  if(!is.null(loci_table)){
    
    positions = NULL
    
    options(scipen=999)
    
    for(pos in 1:nrow(filter_table)){
      positions = c(positions, paste(filter_table[pos,][['CHROM']],
                                     filter_table[pos,][['start']]:filter_table[pos,][['end']], sep = '_'))
      
    }
    
    options(scipen=0)
    
    loci_table = loci_table[rownames(loci_table) %in% positions,]
  
  }else if(!is.null(vcf_table)){
    
    rownames(vcf_object) = paste(vcf_object[['CHROM']], vcf_object[['POS']], sep = '_')
    
    positions = NULL
    
    options(scipen=999)
    
    for(pos in 1:nrow(filter_table)){
      positions = c(positions, paste(filter_table[pos,][['CHROM']],
                                     filter_table[pos,][['start']]:filter_table[pos,][['end']], sep = '_'))
      
    }
    
    options(scipen=0)
    
    vcf_object = vcf_object[rownames(vcf_object) %in% positions,]
    
    loci_table = vcf_object[,c(1,2,4,5)]
    
    gt = NULL
    
    for(w in 1:10){
      gt = rbind(gt, get_GTAD_matrix(vcf = vcf_object, w = w, n = 10, threshold = 5))
    }
    
    loci_table = cbind(loci_table, get_AC(gt, loci_table, w = 1, n = 1))
    
    loci_table$TypeOf_Markers = TypeOf_Marker(loci_table, w = 1, n = 1)
    
  }
  
  return(loci_table)
  
}


# get_MHAP_fragment ----

get_MHAP_fragments = function(MHAPs = NULL,
                             upstream_buffer = 50,
                             downstream_buffer = 50,
                             path_to_reference_gff = 'genes.gff',
                             path_to_reference_genome = 'PvP01.v1.fasta',
                             loci_table = Locus_info
                             ){
  library(Biostrings)
  
  reference_gff = ape::read.gff(path_to_reference_gff)
  reference_genome = Biostrings::readDNAStringSet(path_to_reference_genome, format = 'fasta')
  
  filtered_loci_table = NULL
  target_sequences = NULL
  
  for(MHAP in 1:nrow(MHAPs)){
    
    temp_filtered_loci_table = filter_loci_table(loci_table = loci_table, filter_table = MHAPs[MHAP,])
    
    
    upstream_filter_table = data.frame(CHROM = MHAPs[MHAP,][['CHROM']],
                                       start = MHAPs[MHAP,][['start']] - upstream_buffer,
                                       end = min(temp_filtered_loci_table$POS - 1))
    
    temp_filtered_upstream_loci_table = filter_loci_table(loci_table = loci_table, filter_table = upstream_filter_table)
    
    
    downstream_filter_table = data.frame(CHROM = MHAPs[MHAP,][['CHROM']],
                                         start = max(temp_filtered_loci_table$POS) + 1,
                                         end = MHAPs[MHAP,][['end']] + downstream_buffer)
    
    temp_filtered_downstream_loci_table = filter_loci_table(loci_table = loci_table, filter_table = downstream_filter_table)
    
    
    temp_filtered_loci_table$polymorphism_location = 'Target sequence'
    if(nrow(temp_filtered_upstream_loci_table) > 0){
      temp_filtered_upstream_loci_table$polymorphism_location = 'Upstream region'
    }
    if(nrow(temp_filtered_downstream_loci_table)>0){
      temp_filtered_downstream_loci_table$polymorphism_location = 'Downstream region'
    }
    
    temp_filtered_loci_table = rbind(
      temp_filtered_upstream_loci_table,
      temp_filtered_loci_table,
      temp_filtered_downstream_loci_table
    )
    
    temp_filtered_loci_table = temp_filtered_loci_table[,c('CHROM', 'POS', 'Alleles', 'Allele_Counts','TypeOf_Markers', 'polymorphism_location')]
    
    target_sequences[[MHAPs[MHAP,][['Locus']]]] = as.character(subseq(reference_genome[grep(MHAPs[MHAP,][['CHROM']], names(reference_genome))],
                                                                      start = MHAPs[MHAP,][['start']] - upstream_buffer,
                                                                      end = MHAPs[MHAP,][['end']] + downstream_buffer))
    
    target_sequences[[MHAPs[MHAP,][['Locus']]]] = unlist(str_split(target_sequences[[MHAPs[MHAP,][['Locus']]]], ''))
    
    for(polymorphism in 1:nrow(temp_filtered_loci_table)){
      
      if(temp_filtered_loci_table[polymorphism,][['TypeOf_Markers']] == 'SNP'){
        
        target_sequences[[MHAPs[MHAP,][['Locus']]]][temp_filtered_loci_table[polymorphism,][['POS']] - (MHAPs[MHAP,][['start']] - upstream_buffer) + 1] = mergeIUPACLetters(
          paste(unlist(str_split(gsub(':\\d+','',temp_filtered_loci_table[polymorphism,][['Alleles']]),',')), collapse = ''))
        
        
      }else{
        target_sequences[[MHAPs[MHAP,][['Locus']]]][temp_filtered_loci_table[polymorphism,][['POS']] - (MHAPs[MHAP,][['start']] - upstream_buffer) + 1] = paste(
          rep('-',max(nchar(unlist(str_split(gsub(':\\d+','',temp_filtered_loci_table[polymorphism,][['Alleles']]),','))))), collapse = '')
        
      }
      
    }
    
    
    upstream_region = target_sequences[[MHAPs[MHAP,][['Locus']]]][1:(min(temp_filtered_loci_table[temp_filtered_loci_table$polymorphism_location == 'Target sequence',][['POS']]) -
                                                                       (MHAPs[MHAP,][['start']] - upstream_buffer))]
    
    targeted_region = target_sequences[[MHAPs[MHAP,][['Locus']]]][(min(temp_filtered_loci_table[temp_filtered_loci_table$polymorphism_location == 'Target sequence',][['POS']]) -
                                                                     (MHAPs[MHAP,][['start']] - upstream_buffer) + 1):(max(temp_filtered_loci_table[temp_filtered_loci_table$polymorphism_location == 'Target sequence',][['POS']]) -
                                                                                                                         (MHAPs[MHAP,][['start']] - upstream_buffer) + 1)]
    
    downstream_region = target_sequences[[MHAPs[MHAP,][['Locus']]]][(max(temp_filtered_loci_table[temp_filtered_loci_table$polymorphism_location == 'Target sequence',][['POS']]) -
                                                                       (MHAPs[MHAP,][['start']] - upstream_buffer) + 2):length(target_sequences[[MHAPs[MHAP,][['Locus']]]])]
    
    
    target_sequences[[MHAPs[MHAP,][['Locus']]]] = paste(paste(upstream_region, collapse = ''),
                                                        '[',
                                                        paste(targeted_region, collapse = ''),
                                                        ']',
                                                        paste(downstream_region, collapse = ''), sep = '')
    
    
    temp_filtered_loci_table = data.frame(MHAPs[MHAP,c('Locus', 'start', 'end')], temp_filtered_loci_table)
    
    temp_filtered_loci_table$start = temp_filtered_loci_table$start - upstream_buffer
    temp_filtered_loci_table$end = temp_filtered_loci_table$end + downstream_buffer
    
    filtered_loci_table = rbind(filtered_loci_table,
                                temp_filtered_loci_table)
    
  }
  
  return(target_sequences)
  
}


# get_MHAPs_for_geneTarget ----

get_MHAPs_for_geneTarget = function(gene_id,
                                    upstream_buffer = 50,
                                    downstream_buffer = 50,
                                    path_to_reference_gff = 'genes.gff',
                                    path_to_reference_genome = 'PvP01.v1.fasta',
                                    loci_table = Locus_info,
                                    max_length){
  
  library(Biostrings)
  
  reference_gff = ape::read.gff(path_to_reference_gff)
  reference_genome = Biostrings::readDNAStringSet(path_to_reference_genome, format = 'fasta')
  
  selected_CDS = reference_gff[grepl(gene_id,reference_gff$attributes)&reference_gff$type == 'CDS',]
  
  selected_CDS = selected_CDS[,c('seqid', 'start', 'end')]
  names(selected_CDS) = c('CHROM', 'start', 'end')
  
  selelcted_CDS_loci_table = filter_loci_table(loci_table = loci_table,
                                               filter_table = selected_CDS,
                                               by = 'segments')
  
  selelcted_CDS_loci_table$dist = c(selelcted_CDS_loci_table$POS[-1] -
                                      selelcted_CDS_loci_table$POS[-length(selelcted_CDS_loci_table$POS)],
                                    Inf)
  
  MHAPs = NULL
  pos = 1
  start_pos = 1
  end_pos = 1
  
  temp_MHAP = data.frame(CHROM = selelcted_CDS_loci_table$CHROM[pos],
                         start = selelcted_CDS_loci_table$POS[pos],
                         end = selelcted_CDS_loci_table$POS[pos])
  
  while(pos <= length(selelcted_CDS_loci_table$POS)){
    
    if(temp_MHAP$end - temp_MHAP$start + selelcted_CDS_loci_table$dist[pos] <= max_length){
      pos = pos + 1
      end_pos = pos
      
      temp_MHAP = data.frame(CHROM = selelcted_CDS_loci_table$CHROM[pos],
                             start = selelcted_CDS_loci_table$POS[start_pos],
                             end = selelcted_CDS_loci_table$POS[end_pos])
      
    }else{
      pos = pos + 1
      start_pos = pos
      end_pos = pos
      
      MHAPs = rbind(MHAPs, temp_MHAP)
      
      temp_MHAP = data.frame(CHROM = selelcted_CDS_loci_table$CHROM[pos],
                             start = selelcted_CDS_loci_table$POS[start_pos],
                             end = selelcted_CDS_loci_table$POS[end_pos])
    }
    
  }
  
  MHAPs = data.frame(Locus = paste(MHAPs$CHROM, MHAPs$start, MHAPs$end, sep = '_'), MHAPs)
  
  target_sequences = get_MHAP_fragments(MHAPs = MHAPs)
  
  return(target_sequences)
  
}


# find_sequence----

find_DNAsequence = function(sequences = NULL,
                         reference_genome = NULL 
                         ){
  
  sequence_location = NULL
  
  for(DNAsequence in sequences){
    
    sequence_length = nchar(DNAsequence)
    
    print(paste0('searching DNA sequence ', DNAsequence))
    
    grep(DNAsequence, as.character(reference_genome)) 
    
    chrom = grep(DNAsequence, as.character(reference_genome)) 
    
    if(length(chrom) == 0){
      DNAsequence = as.character(reverseComplement(DNAString(DNAsequence)))
      chrom = grep(DNAsequence, as.character(reference_genome)) 
    }
    
    print(paste0('searching DNA sequence in chrom ', chrom))
    
    temp_chrom = unlist(str_split(as.character(reference_genome[[chrom]]),''))
    
    w = 100000
    
    low = 1
    high = nchar(reference_genome[[chrom]])
    
    while(w > sequence_length){
      
      s = c(seq(low, high, w), high + 1)
      
      temp_search = data.frame(CHROM =names(reference_genome)[chrom], start = s[-length(s)],
                               end = s[-1] -1)
      
      seqs = apply(temp_search, 1, function(x){paste(temp_chrom[seq(x['start'], x['end'],1)], collapse = '')})
      
      i = grep(DNAsequence, seqs)
      
      while(length(i) == 0){
        s = s + 1
        
        temp_search = data.frame(CHROM =names(reference_genome)[chrom], start = s[-length(s)],
                                 end = s[-1] -1)
        
        seqs = apply(temp_search, 1, function(x){paste(temp_chrom[seq(x['start'], x['end'],1)], collapse = '')})
        
        i = grep(DNAsequence, seqs)
        
        
      }
      
      low = s[i]
      high = s[i+1]
      
      w = w/10
      
    }
    
    temp_search = data.frame(CHROM =names(reference_genome)[chrom], start = low:(high - (sequence_length - 1)),
                             end = (low + sequence_length -1):high)
    
    seqs = apply(temp_search, 1, function(x){paste(temp_chrom[seq(x['start'], x['end'],1)], collapse = '')})
    
    sequence_location = rbind(sequence_location, data.frame(DNAsequence = DNAsequence, temp_search[grep(DNAsequence, seqs),]))
    
  }
  
  rownames(sequence_location) = 1:nrow(sequence_location)
  
  return(sequence_location)
  
}



# load Rcpp functions----

if(!require(Rcpp)){
  install.packages('Rcpp')
  library(Rcpp)
}

#sourceCpp('Rcpp_functions.cpp')
