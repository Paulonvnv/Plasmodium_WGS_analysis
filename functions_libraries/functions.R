# run_vcftools ----
run_vcftools = function(vcf = NULL,
                        gzvcf = NULL,
                        bcf = NULL,
                        
                        bash_file = NULL, 
                        
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
  
  print(vcf_arguments)
  
  write.table(vcf_run_file, bash_file, row.names = FALSE, quote = FALSE, col.names = FALSE)
  
  system(paste0('chmod 777 ', bash_file))
  system(paste0('./', bash_file))
  
}

# load_vcf----
load_vcf = function(vcf = NULL,
                    gzvcf = NULL,
                    na.rm = TRUE
                    ){
  
  if(!is.null(vcf)){
    
    temp_tsv_name = gsub('vcf', 'tsv', vcf)
    system(paste0("grep -v '^##' ", vcf, " > ", temp_tsv_name))
  }
  
  if(!is.null(gzvcf)){
    system(paste0("zgrep -v '^##' ", vcf, " > ", temp_tsv_name))
  }
  
  col_names_file = paste0('col_names_', temp_tsv_name)
  
  system(paste0("grep '^#' ", temp_tsv_name, " > ", col_names_file))
  
  col_names = strsplit(readLines(col_names_file), '\t')[[1]]
  col_names[1] = 'CHROM'
  
  vcf = read.table(temp_tsv_name)
  names(vcf) = col_names
  
  if(na.rm){
    vcf %<>% filter(ALT != '.')
  }
  
  system(paste0('rm ', temp_tsv_name))
  system(paste0('rm ', col_names_file))
  return(vcf)
}


# rGenome S4class and vcf2rGenome----

## rGenome S4 class

setClass('rGenome', slots = c(
  gt = "ANY",
  loci_table = "ANY",
  metadata = "ANY"
))

## rGenome constructor
rGenome = function(gt = NULL,
                   loci_table = NULL,
                   metadata = NULL){
  obj = new('rGenome')
  obj@gt = gt
  obj@loci_table = loci_table
  obj@metadata = metadata
  
  return(obj)
}

## set function vcf2rGenome

vcf2rGenome = function(vcf, n = 500, threshold = 5) {
            
            # Generate metadata
            metadata = data.frame(Sample_id = names(vcf)[-1:-9])
            rownames(metadata) = metadata[['Sample_id']]
            
            # generate loci_table
            loci_table = vcf[,c(1,2,4,5)]
            rownames(loci_table) = paste(loci_table$CHROM, loci_table$POS, sep = '_')
            
            # generate a haplotype table (gt)
            
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

setGeneric("SampleAmplRate", function(obj, update = TRUE, threshold = NULL, n = 100) standardGeneric("SampleAmplRate"))

setMethod("SampleAmplRate", signature(obj = "rGenome"),
          
          function(obj, update = TRUE, threshold = NULL, n = 100) {
            
            obj2 = obj
            
            if(!is.null(threshold)){
              obj2@gt = prune_alleles(obj = obj2, threshold = threshold, n = n)
            }
            
            if(update){
              obj2@metadata[['SampleAmplRate']] =
                1 - colSums(is.na(obj2@gt))/nrow(obj2@gt)
              return(obj2)
            }else{
              result = 1 - colSums(is.na(obj2@gt))/nrow(obj2@gt)
              return(result)
            }
            
          }
)

# LocusAmplRate----

setGeneric("LocusAmplRate", function(obj, update = TRUE, threshold = NULL, n = 100) standardGeneric("LocusAmplRate"))

setMethod("LocusAmplRate", signature(obj = "rGenome"),
          
          function(obj, update = TRUE, threshold = NULL, n = 100) {
            
            obj2 = obj
            
            if(!is.null(threshold)){
              obj2@gt = prune_alleles(obj = obj2, threshold = threshold, n = n)
            }
            
            if(update){
              obj2@loci_table[['LocusAmplRate']] =
                1 - rowSums(is.na(obj2@gt))/ncol(obj2@gt)
              return(obj2)
            }else{
              result = 1 - rowSums(is.na(obj2@gt))/ncol(obj2@gt)
              return(result)
            }
            
          }
)

# filter_samples----


setGeneric("filter_samples", function(obj, v = NULL) standardGeneric("filter_samples"))

setMethod("filter_samples", signature(obj = "rGenome"),
          
          function(obj, v = NULL) {
            
            obj2 = rGenome(gt = obj@gt[,v],
                           loci_table = obj@loci_table,
                           metadata = obj@metadata[v,])
            
            return(obj2)
          }
)

# filter_loci----

setGeneric("filter_loci", function(obj, v = NULL) standardGeneric("filter_loci"))

setMethod("filter_loci", signature(obj = "rGenome"),
          
          function(obj, v = NULL) {
            
            obj2 = obj
            obj2@gt = obj2@gt[v,]
            obj2@loci_table = obj2@loci_table[v,]
            
            if(is.null(nrow(obj2@gt))){
              
              obj2@gt = matrix(obj2@gt, nrow = 1, ncol = length(obj2@gt),
                               dimnames = list(
                                 rownames(obj2@loci_table),
                                 names(obj2@gt)
                               ))
              
            }
            
            return(obj2)
          }
)

# handle_ploidy----

handle_ploidy = function(gt, monoclonals, polyclonals, w = 1, n = 1){
  
  s = round(seq(1,nrow(gt)+1, length.out=n+1))
  low = s[w]
  high = s[w+1]-1
  
  if(sum(grepl(':', gt)) > 0){
      gt = matrix(gsub(':\\d+', '', gt[low:high,]), nrow = high - low + 1, ncol = ncol(gt),
                  dimnames = list(rownames(gt)[low:high], colnames(gt)))
  }
  
  # if(sum(grepl('/', gt)) == 0){
  #   
  #   gt3 = gt
  #   
  # }else
  
  if(is.null(monoclonals) & is.null(polyclonals)){
    gt1 = gsub('/\\d+', '', gt)
    gt2 = gsub('\\d+/', '', gt)
    gt2[gt2==gt1] = NA
    gt3 = cbind(gt1, gt2)
  }else{
    gt_mono = matrix(gt[,monoclonals],
                     nrow = nrow(gt),
                     ncol = length(monoclonals),
                     dimnames = list(rownames(gt), monoclonals))
    
    gt_mono = gsub('/\\d+', '', gt_mono)
    
    
    gt_poly = matrix(gt[,polyclonals],
                     nrow = nrow(gt),
                     ncol = length(polyclonals),
                     dimnames = list(rownames(gt), polyclonals)
                     )
    gt_poly1 = gsub('/\\d+', '', gt_poly)
    
    gt_poly2 = gsub('\\d+/', '', gt_poly)
    
    if(!is.null(polyclonals)){
      colnames(gt_poly1) = paste(colnames(gt_poly1), 'C1', sep = '_')
      colnames(gt_poly2) = paste(colnames(gt_poly2), 'C2', sep = '_')
    }
    
    gt3 = cbind(gt_mono, gt_poly1, gt_poly2)
  }
  
  return(gt3)
  
}

# get_AC----

setGeneric("get_AC", function(obj = NULL, w = 1, n = 1,
                              update_alleles = TRUE,
                              monoclonals = NULL, polyclonals = NULL) standardGeneric("get_AC"))

setMethod("get_AC", signature(obj = "rGenome"),
          function(obj = NULL, w = 1, n = 1,
                  update_alleles = TRUE,
                  monoclonals = NULL, polyclonals = NULL){
  
            gt = obj@gt
            loci = obj@loci_table
            
            s = round(seq(1,nrow(gt)+1, length.out=n+1))
            low = s[w]
            high = s[w+1]-1
            
            gt3 = handle_ploidy(gt = gt, w = w, n = n, monoclonals = monoclonals, polyclonals = polyclonals)
            
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
            rownames(alleles) = rownames(gt3)
            alleles$Cardinality = as.integer(alleles$Cardinality)
            
            return(alleles)
            
          }
          )


# get_ExpHet----

setGeneric("get_ExpHet", function(obj = NULL, update_AC = FALSE, monoclonals = NULL, polyclonals = NULL, by = NULL) standardGeneric("get_ExpHet"))

setMethod("get_ExpHet", signature(obj = "rGenome"),
          function(obj = NULL, update_AC = FALSE, monoclonals = NULL, polyclonals = NULL, by = NULL){
  
  gt = obj@gt
  loci = obj@loci_table
  metadata = obj@metadata
  
  if(!is.null(by)){
    
    populations = t(table(metadata[[by]]))
    populations = data.frame(population = colnames(populations), nsamples = populations[1,])
    
    ExpHet = NULL
    
    for(pop in populations$population){
      
      if(populations[pop,][['nsamples']] >= 2){
        
        samples = metadata[metadata[[by]] == pop,][['Sample_id']]
        temp_pop = filter_samples(obj = obj, v = samples)
  
        temp_monoclonals = monoclonals[monoclonals %in% samples]
        if(length(temp_monoclonals) == 0){
          temp_monoclonals = NULL
        }
        temp_polyclonals = polyclonals[polyclonals %in% samples]
        if(length(temp_polyclonals) == 0){
          temp_polyclonals = NULL
        }
        
        temp_AC = get_AC(obj = temp_pop, w =1, n = 1, update_alleles = FALSE, monoclonals = temp_monoclonals, polyclonals = temp_polyclonals)
        temp_AC = temp_AC$Allele_Counts
        temp_AC = sapply(temp_AC, function(x){
          temp_AC = strsplit(x, ',')[[1]]
          gsub('^\\d+:', '', temp_AC)
        })
        
        temp_ExpHet = sapply(1:length(temp_AC), function(pos){
          n = sum(as.numeric(temp_AC[[pos]]))
          allele_counts = as.numeric(temp_AC[[pos]])
          allele_freq = allele_counts/n
          sp2 = sum(allele_freq^2)
          ExpHet = n*(1 - sp2)/(n-1)
        })
        
      }else{
        
        # fill with NA's if there is only one sample in the population
        temp_ExpHet = rep(NA, nrow(temp_pop@gt))
        
      }
      
      ExpHet = cbind(ExpHet, temp_ExpHet)
      
    }
    
    # Calculate Heterozygosity for the whole population
    
    AC = get_AC(obj = obj, w =1, n = 1, update_alleles = FALSE, monoclonals = monoclonals, polyclonals = polyclonals)
    AC = AC$Allele_Counts
    AC = sapply(AC, function(x){
      AC = strsplit(x, ',')[[1]]
      gsub('^\\d+:', '', AC)
    })
    
    Total = unlist(sapply(1:length(AC), function(pos){
      n = sum(as.numeric(AC[[pos]]))
      allele_counts = as.numeric(AC[[pos]])
      allele_freq = allele_counts/n
      sp2 = sum(allele_freq^2)
      ExpHet = n*(1 - sp2)/(n-1)
    }))
    
    ExpHet = cbind(ExpHet, Total)
    
    colnames(ExpHet) = c(populations$population, 'Total')
    rownames(ExpHet) = rownames(gt)
    
    ExpHet = cbind(loci[,c('CHROM','POS','Cardinality','TypeOf_Markers')], ExpHet)
    
  }else{
    
    if(!update_AC){
      AC = sapply(1:nrow(loci), function(x){
        AC = strsplit(loci[x,'Allele_Counts'], ',')[[1]]
        gsub('^\\d+:', '', AC)
      })
    }else{
      
      AC = get_AC(obj = obj, w =1, n = 1, update_alleles = FALSE, monoclonals = monoclonals, polyclonals = polyclonals)
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
      ExpHet = n*(1 - sp2)/(n-1)
    })
    
  }
  
  return(ExpHet)
})


# get_nuc_div----

setGeneric("get_nuc_div", function(obj = NULL, monoclonals = NULL, polyclonals = NULL, gff = NULL, type_of_region = NULL, window = NULL, by = NULL) standardGeneric("get_nuc_div"))

setMethod("get_nuc_div", signature(obj = "rGenome"),
          function(obj = NULL, monoclonals = NULL, polyclonals = NULL, gff = NULL, type_of_region = NULL,  window = NULL, by = NULL){
            
            loci = obj@loci_table
            metadata = obj@metadata
            
            if(is.object(gff) & is.null(window)){
              
              dna_regions = gff
              
            }else if(file.exists(gff) & is.null(window)){
              
              ref_gff = ape::read.gff(gff)
              dna_regions = ref_gff[grepl(type_of_region, ref_gff$type)&
                                      !grepl('^Transfer',ref_gff$seqid),
                                    c('seqid', 'start', 'end', 'attributes')]
              
              dna_regions = dna_regions[order(dna_regions$start),]
              dna_regions = dna_regions[order(dna_regions$seqid),]
              rownames(dna_regions) = 1:nrow(dna_regions)
              
              dna_regions = cbind(dna_regions, as.data.frame(t(sapply(1:nrow(dna_regions), function(gene){
                attributes = strsplit(dna_regions[gene,][['attributes']], ';')[[1]]
                c(gene_id = gsub('^ID=','',attributes[grep('^ID=', attributes)]),
                  gene_description = gsub('^description=','',attributes[grep('^description=', attributes)]))
              }))))
              
              dna_regions = dna_regions[,c('seqid', 'start', 'end', 'gene_id', 'gene_description')]
              
            }else if(!is.null(window)){
              
              ref_gff = ape::read.gff(gff)
              ref_gff = ref_gff[grepl(type_of_region, ref_gff$type)&
                                  !grepl('^Transfer',ref_gff$seqid),
                                c('seqid', 'start', 'end', 'attributes')]
              
              ref_gff = ref_gff[order(ref_gff$start),]
              ref_gff = ref_gff[order(ref_gff$seqid),]
              rownames(ref_gff) = 1:nrow(ref_gff)
              
              ref_gff = cbind(ref_gff, as.data.frame(t(sapply(1:nrow(ref_gff), function(gene){
                attributes = strsplit(ref_gff[gene,][['attributes']], ';')[[1]]
                c(gene_id = gsub('^ID=','',attributes[grep('^ID=', attributes)]),
                  gene_description = gsub('^description=','',attributes[grep('^description=', attributes)]))
              }))))
              
              chrom_length = loci %>% group_by(CHROM) %>% summarise(length = max(POS))
              
              chrom_intervals = sapply(chrom_length$length, function(chrom){
                seq(1, chrom, window)
              })
              
              
              dna_regions = NULL
              
              for(chrom in 1:length(chrom_intervals)){
                dna_regions = rbind(dna_regions, data.frame(seqid = chrom_length[chrom,][['CHROM']],
                                                            start = chrom_intervals[[chrom]],
                                                            end = chrom_intervals[[chrom]] - 1 + window))
                
              }
              
              
              dna_regions$gene_ids = sapply(1:nrow(dna_regions),function(bin){
                paste(ref_gff[ref_gff[['seqid']] == dna_regions[bin,][['seqid']] &
                                ((ref_gff[['start']] > dna_regions[bin,][['start']] &
                                    ref_gff[['start']] < dna_regions[bin,][['end']])|
                                   
                                   (ref_gff[['end']] > dna_regions[bin,][['start']] &
                                      ref_gff[['end']] < dna_regions[bin,][['end']])|
                                   
                                   (dna_regions[bin,][['start']] > ref_gff[['start']] &
                                      dna_regions[bin,][['start']] < ref_gff[['end']])|
                                   
                                   (dna_regions[bin,][['end']] > ref_gff[['start']] &
                                      dna_regions[bin,][['end']] < ref_gff[['end']])
                                ),][['gene_id']], collapse = ',')}, simplify = T)
              
              ###### Add gene_description
              dna_regions$genes_description = sapply(1:nrow(dna_regions),function(bin){
                paste(ref_gff[ref_gff[['seqid']] == dna_regions[bin,][['seqid']] &
                                ((ref_gff[['start']] > dna_regions[bin,][['start']] &
                                    ref_gff[['start']] < dna_regions[bin,][['end']])|
                                   
                                   (ref_gff[['end']] > dna_regions[bin,][['start']] &
                                      ref_gff[['end']] < dna_regions[bin,][['end']])|
                                   
                                   (dna_regions[bin,][['start']] > ref_gff[['start']] &
                                      dna_regions[bin,][['start']] < ref_gff[['end']])|
                                   
                                   (dna_regions[bin,][['end']] > ref_gff[['start']] &
                                      dna_regions[bin,][['end']] < ref_gff[['end']])
                                ),][['gene_description']], collapse = ',')}, simplify = T)
              
            }else if(is.null(window) & is.null(window)){
              
              print('You must provide at least a gff file')
              
            }
            
            
            if(!is.null(by)){
              
              populations = t(table(metadata[[by]]))
              populations = data.frame(population = colnames(populations), nsamples = populations[1,])
              
              
              for(pop in populations$population){
                
                if(populations[pop,][['nsamples']] >= 2){
                  
                  samples = metadata[metadata[[by]] == pop,][['Sample_id']]
                  temp_pop = filter_samples(obj = obj, v = samples)
                  
                  temp_monoclonals = monoclonals[monoclonals %in% samples]
                  if(length(temp_monoclonals) == 0){
                    temp_monoclonals = NULL
                  }
                  temp_polyclonals = polyclonals[polyclonals %in% samples]
                  if(length(temp_polyclonals) == 0){
                    temp_polyclonals = NULL
                  }
                  
                  gt = temp_pop@gt
                  
                  gt3 = handle_ploidy(gt, monoclonals = temp_monoclonals, polyclonals = temp_polyclonals)
                  gt3 = as.data.frame(gt3)
                  
                  pi = NULL
                  var = NULL
                  
                  for(region in 1:nrow(dna_regions)){
                    positions = paste(dna_regions[region, ][['seqid']],dna_regions[region, ][['start']]:dna_regions[region, ][['end']], sep = '_')
                    
                    region_length = dna_regions[region, ][['end']] - dna_regions[region, ][['start']] + 1
                    
                    temp_gt = gt3[rownames(gt3) %in% positions,]
                    
                    n = ncol(temp_gt)
                    
                    if(nrow(temp_gt)>0){
                      
                      temp_loci = loci[rownames(loci) %in% positions,]
                      
                      temp_indels = temp_loci[temp_loci$TypeOf_Markers != 'SNP',]
                      
                      temp_gt = temp_gt[temp_loci$TypeOf_Markers == 'SNP',] # Keep only SNPs
                      
                      if(nrow(temp_gt)>0){
                        
                        if(length(temp_indels$ALT) > 0){
                          
                          ALTs = temp_indels$ALT
                          
                          ALTs = strsplit(ALTs, ',')
                          
                          gaps =  nchar(temp_indels$REF) - sapply(ALTs, function(alt){
                            max(nchar(alt))
                          })
                          
                          gaps[gaps < 0] = 0
                          
                          region_length = region_length - sum(gaps) # remove deletions from the total size of the region
                          
                        }
                        
                        haplotypes_counts =  summary(as.factor(sapply(1:ncol(temp_gt), function(x){paste(temp_gt[,x], collapse = '')})))
                        
                        if(length(haplotypes_counts) > 1){
                          
                          names(haplotypes_counts) = gsub('NA', "_", names(haplotypes_counts))
                          
                          haplotypes_freqs = haplotypes_counts/sum(haplotypes_counts)
                          haplotypes = names(haplotypes_counts)
                          
                          combinations = combn(1:length(haplotypes),2)
                          
                          temp_pi = NULL
                          
                          for(comb in 1:ncol(combinations)){
                            
                            x_i = haplotypes_freqs[combinations[1,comb]]
                            x_j = haplotypes_freqs[combinations[2,comb]]
                            
                            seq_i = unlist(strsplit(haplotypes[combinations[1,comb]], ''))
                            seq_i[seq_i == '_'] = NA
                            seq_j = unlist(strsplit(haplotypes[combinations[2,comb]], ''))
                            seq_j[seq_j == '_'] = NA
                            
                            pi_ij = sum(seq_i != seq_j, na.rm = T)/(region_length - sum(is.na(seq_i != seq_j)))
                            temp_pi = c(temp_pi, 2*x_i*x_j*pi_ij)
                          }
                          
                          
                          pi = c(pi, (n/(n-1))*(sum(temp_pi)))
                          var = c(var, (n + 1)*sum(temp_pi)/(3*(n - 1)*region_length) + 2*(n^2 + n + 3)*sum(temp_pi)^2/(9*n*(n - 1)))
                          
                        }else{
                          
                          pi = c(pi, 0)
                          var = c(var, 0)
                        }
                        
                      }else{
                        
                        pi = c(pi, NA)
                        var = c(var, NA)
                      }
                      
                    }else{
                      
                      pi = c(pi, NA)
                      var = c(var, NA)
                      
                    }
                    
                  }
                  
                }else{
                  
                  print(paste0(pop, ' will be excluded because sample size is les than 2 individuals')) 
                  
                }
                
                temp_dna_regions_pi = data.frame(pi, var)
                
                names(temp_dna_regions_pi) = c(paste0(pop, '_pi'),
                                               paste0(pop, '_pi_var'))
                
                dna_regions = cbind(dna_regions, temp_dna_regions_pi)
              }
              
              gt = obj@gt
              
              gt3 = handle_ploidy(gt, monoclonals = monoclonals, polyclonals = polyclonals)
              gt3 = as.data.frame(gt3)
              
              pi = NULL
              var = NULL
              
              for(region in 1:nrow(dna_regions)){
                positions = paste(dna_regions[region, ][['seqid']],dna_regions[region, ][['start']]:dna_regions[region, ][['end']], sep = '_')
                
                region_length = dna_regions[region, ][['end']] - dna_regions[region, ][['start']] + 1
                
                temp_gt = gt3[rownames(gt3) %in% positions,]
                
                n = ncol(temp_gt)
                
                if(nrow(temp_gt)>0){
                  
                  temp_loci = loci[rownames(loci) %in% positions,]
                  
                  temp_indels = temp_loci[temp_loci$TypeOf_Markers != 'SNP',]
                  
                  temp_gt = temp_gt[temp_loci$TypeOf_Markers == 'SNP',] # Keep only SNPs
                  
                  if(nrow(temp_gt)>0){
                    
                    if(length(temp_indels$ALT) > 0){
                      
                      ALTs = temp_indels$ALT
                      
                      ALTs = strsplit(ALTs, ',')
                      
                      gaps =  nchar(temp_indels$REF) - sapply(ALTs, function(alt){
                        max(nchar(alt))
                      })
                      
                      gaps[gaps < 0] = 0
                      
                      region_length = region_length - sum(gaps) # remove deletions from the total size of the region
                      
                    }
                    
                    haplotypes_counts =  summary(as.factor(sapply(1:ncol(temp_gt), function(x){paste(temp_gt[,x], collapse = '')})))
                    
                    if(length(haplotypes_counts) > 1){
                      
                      names(haplotypes_counts) = gsub('NA', "_", names(haplotypes_counts))
                      
                      haplotypes_freqs = haplotypes_counts/sum(haplotypes_counts)
                      haplotypes = names(haplotypes_counts)
                      
                      combinations = combn(1:length(haplotypes),2)
                      
                      temp_pi = NULL
                      
                      for(comb in 1:ncol(combinations)){
                        
                        x_i = haplotypes_freqs[combinations[1,comb]]
                        x_j = haplotypes_freqs[combinations[2,comb]]
                        
                        seq_i = unlist(strsplit(haplotypes[combinations[1,comb]], ''))
                        seq_i[seq_i == '_'] = NA
                        seq_j = unlist(strsplit(haplotypes[combinations[2,comb]], ''))
                        seq_j[seq_j == '_'] = NA
                        
                        pi_ij = sum(seq_i != seq_j, na.rm = T)/(region_length - sum(is.na(seq_i != seq_j)))
                        temp_pi = c(temp_pi, 2*x_i*x_j*pi_ij)
                      }
                      
                      
                      pi = c(pi, (n/(n-1))*(sum(temp_pi)))
                      var = c(var, (n + 1)*sum(temp_pi)/(3*(n - 1)*region_length) + 2*(n^2 + n + 3)*sum(temp_pi)^2/(9*n*(n - 1)))
                      
                    }else{
                      
                      pi = c(pi, 0)
                      var = c(var, 0)
                    }
                    
                  }else{
                    
                    pi = c(pi, NA)
                    var = c(var, NA)
                  }
                  
                }else{
                  
                  pi = c(pi, NA)
                  var = c(var, NA)
                  
                }
                
              }
              
              dna_regions = cbind(dna_regions, data.frame(Total_pi = pi, Total_pi_var = var))
              
              
            }else{
              
              gt = obj@gt
              
              gt3 = handle_ploidy(gt, monoclonals = monoclonals, polyclonals = polyclonals)
              gt3 = as.data.frame(gt3)
              
              pi = NULL
              var = NULL
              
              for(region in 1:nrow(dna_regions)){
                positions = paste(dna_regions[region, ][['seqid']],dna_regions[region, ][['start']]:dna_regions[region, ][['end']], sep = '_')
                
                region_length = dna_regions[region, ][['end']] - dna_regions[region, ][['start']] + 1
                
                temp_gt = gt3[rownames(gt3) %in% positions,]
                
                n = ncol(temp_gt)
                
                if(nrow(temp_gt)>0){
                  
                  temp_loci = loci[rownames(loci) %in% positions,]
                  
                  temp_indels = temp_loci[temp_loci$TypeOf_Markers != 'SNP',]
                  
                  temp_gt = temp_gt[temp_loci$TypeOf_Markers == 'SNP',] # Keep only SNPs
                  
                  if(nrow(temp_gt)>0){
                    
                    if(length(temp_indels$ALT) > 0){
                      
                      ALTs = temp_indels$ALT
                      
                      ALTs = strsplit(ALTs, ',')
                      
                      gaps =  nchar(temp_indels$REF) - sapply(ALTs, function(alt){
                        max(nchar(alt))
                      })
                      
                      gaps[gaps < 0] = 0
                      
                      region_length = region_length - sum(gaps) # remove deletions from the total size of the region
                      
                    }
                    
                    haplotypes_counts =  summary(as.factor(sapply(1:ncol(temp_gt), function(x){paste(temp_gt[,x], collapse = '')})))
                    
                    if(length(haplotypes_counts) > 1){
                      
                      names(haplotypes_counts) = gsub('NA', "_", names(haplotypes_counts))
                      
                      haplotypes_freqs = haplotypes_counts/sum(haplotypes_counts)
                      haplotypes = names(haplotypes_counts)
                      
                      combinations = combn(1:length(haplotypes),2)
                      
                      temp_pi = NULL
                      
                      for(comb in 1:ncol(combinations)){
                        
                        x_i = haplotypes_freqs[combinations[1,comb]]
                        x_j = haplotypes_freqs[combinations[2,comb]]
                        
                        seq_i = unlist(strsplit(haplotypes[combinations[1,comb]], ''))
                        seq_i[seq_i == '_'] = NA
                        seq_j = unlist(strsplit(haplotypes[combinations[2,comb]], ''))
                        seq_j[seq_j == '_'] = NA
                        
                        pi_ij = sum(seq_i != seq_j, na.rm = T)/(region_length - sum(is.na(seq_i != seq_j)))
                        temp_pi = c(temp_pi, 2*x_i*x_j*pi_ij)
                      }
                      
                      
                      pi = c(pi, (n/(n-1))*(sum(temp_pi)))
                      var = c(var, (n + 1)*sum(temp_pi)/(3*(n - 1)*region_length) + 2*(n^2 + n + 3)*sum(temp_pi)^2/(9*n*(n - 1)))
                      
                    }else{
                      
                      pi = c(pi, 0)
                      var = c(var, 0)
                    }
                    
                  }else{
                    
                    pi = c(pi, NA)
                    var = c(var, NA)
                  }
                  
                }else{
                  
                  pi = c(pi, NA)
                  var = c(var, NA)
                  
                }
                
              }
              
              dna_regions$pi = pi
              dna_regions$pi_var = var
              
            }
            
            return(dna_regions)
          })


# get_EffCard----

setGeneric("get_EffCard", function(obj = NULL, update_AC = FALSE, monoclonals = NULL, polyclonals = NULL) standardGeneric("get_EffCard"))

setMethod("get_EffCard", signature(obj = "rGenome"),
          function(obj = NULL, update_AC = FALSE, monoclonals = NULL, polyclonals = NULL){
  
            gt = obj@gt
            loci = obj@loci_table
            
            if(!update_AC){
              AC = sapply(1:nrow(loci), function(x){
                AC = strsplit(loci[x,'Allele_Counts'], ',')[[1]]
                gsub('^\\d+:', '', AC)
              })
            }else{
              AC = get_AC(gt = obj,w =1, n = 1, update_alleles = FALSE, monoclonals = monoclonals, polyclonals = polyclonals)
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
            
            gt = obj@gt
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
            
            loci_table = obj@loci_table
            
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
            
            gt = obj@gt
            loci = obj@loci_table
            
            s = round(seq(1,nrow(gt)+1, length.out=n+1))
            low = s[w]
            high = s[w + 1] - 1
            
            alt = gsub('^\\d+,', '', gsub('(\\w+|\\*):', '', loci[low:high, 'Alleles']))
            gt = gsub(':\\d+', '',gt[low:high,])
            
            
            HetPos = matrix(grepl('/', gt), ncol = ncol(gt), nrow = nrow(gt))
            
            frac_ofHet_pAlt = sapply(1:nrow(gt), function(variant) {
              temp_gts = gt[variant,]
              alleles = strsplit(alt[variant], ',')[[1]]
              
              haplotypes = sapply(alleles,
                                 function(allele){
                                   haplotypes = grepl(allele, temp_gts)})
              
              het_haplotypes = (haplotypes == 1 & HetPos[variant,] == 1)    
              
              sum(het_haplotypes, na.rm = T)/sum(haplotypes, na.rm = T)
            })
            
            return(frac_ofHet_pAlt)
            
          }
          )

# Get haplotype matrix----

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
  
            data = obj@loci_table
            
            ref_gff = ape::read.gff(gff)
            coding_regions = ref_gff[grepl('gene', ref_gff$type)&
                                       !grepl('^Transfer',ref_gff$seqid),
                                     c('seqid', 'start', 'end', 'attributes')]
            
            coding_regions$attributes = gsub('ID=','',str_extract(coding_regions$attributes, 'ID=(PVP01|PF3D7)_([0-9]+|MIT[0-9]+|API[0-9]+)'))
            
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



# get_gene_description ----

get_gene_description = function(obj = NULL, gff = 'genes.gff'){
            
            data = obj@loci_table
            data$POS_end = data$POS + nchar(data$REF) - 1
            
            ref_gff = ape::read.gff(gff)
            coding_regions = ref_gff[grepl('gene', ref_gff$type)&
                                       !grepl('^Transfer',ref_gff$seqid),
                                     c('seqid', 'start', 'end', 'attributes')]
            
            coding_regions$gene_id = gsub('ID=','',str_extract(coding_regions$attributes, 'ID=(PVP01|PF3D7)_([0-9]+|MIT[0-9]+|API[0-9]+)'))
            
            coding_regions$gene_description = gsub('description=','',str_extract(coding_regions$attributes, '(description=.+$|description=.+;)'))
            
            coding_regions = coding_regions[order(coding_regions$start),]
            coding_regions = coding_regions[order(coding_regions$seqid),]
            rownames(coding_regions) = 1:nrow(coding_regions)
            
            
            data$gene_id = NA
            data$gene_description = NA
            
            
            for(gene in 1:nrow(coding_regions)){
              
              if(nrow(data[data$CHROM == coding_regions[gene, ][['seqid']]&
                           data$POS >= coding_regions[gene, ][['start']]&
                           data$POS <= coding_regions[gene, ][['end']],]) != 0){
                
                data[data$CHROM == coding_regions[gene, ][['seqid']]&
                       data$POS >= coding_regions[gene, ][['start']]&
                       data$POS <= coding_regions[gene, ][['end']],][['gene_id']] = coding_regions[gene, ][['gene_id']]
                
                data[data$CHROM == coding_regions[gene, ][['seqid']]&
                       data$POS >= coding_regions[gene, ][['start']]&
                       data$POS <= coding_regions[gene, ][['end']],][['gene_description']] = coding_regions[gene, ][['gene_description']]
              }
              
            }
            
            # indels that start out of the gene region but end 
            
            for(pos in rownames(data[is.na(data$gene_id), ])){
              
              if(nrow(coding_regions[coding_regions[['seqid']] == data[pos,][['CHROM']] &
                                     coding_regions[['start']] <= data[pos,][['POS_end']] &
                                     coding_regions[['end']] >= data[pos,][['POS_end']],]) != 0){
                data[pos,][['gene_id']] = coding_regions[coding_regions[['seqid']] == data[pos,][['CHROM']] &
                                                           coding_regions[['start']] <= data[pos,][['POS_end']] &
                                                           coding_regions[['end']] >= data[pos,][['POS_end']],][['gene_id']]
                
                data[pos,][['gene_description']] = coding_regions[coding_regions[['seqid']] == data[pos,][['CHROM']] &
                                                                    coding_regions[['start']] <= data[pos,][['POS_end']] &
                                                                    coding_regions[['end']] >= data[pos,][['POS_end']],][['gene_description']]
              }
              
            }
            
            return(data[, c('gene_id', 'gene_description')])
          }
          

# Polymorphims density per target----

setGeneric("SNP_density", function(obj = NULL, gff = 'genes.gff') standardGeneric("SNP_density"))

setMethod("SNP_density", signature(obj = "rGenome"),
          
          function(obj = NULL, gff = 'genes.gff'){
            
            data = obj@loci_table
            
            ref_gff = ape::read.gff(gff)
            coding_regions = ref_gff[grepl('gene', ref_gff$type)&
                                       !grepl('^Transfer', ref_gff$seqid), c('seqid', 'start', 'end', 'attributes')]
            
            coding_regions$attributes = gsub('ID=', '', str_extract(coding_regions$attributes, 'ID=(PVP01|PF3D7)_([0-9]+|MIT[0-9]+|API[0-9]+)'))
            
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

# Summarise_ReadDepth----

Summarise_ReadDepth = function(obj, by = NULL, w = 1, n = 100){
  
  s = round(seq(1,nrow(obj@gt)+1, length.out=n+1))
  low = s[w]
  high = s[w+1]-1
  
  metadata = obj@metadata
  
  if(!is.null(by)){
    
    populations = t(table(metadata[[by]]))
    populations = data.frame(population = colnames(populations), nsamples = populations[1,])
    
    ReadDepth_Summ = NULL
    
    # For each subpopulation
    
    for(pop in populations$population){
      
      samples = metadata[metadata[[by]] == pop & !is.na(metadata[[by]]),][['Sample_id']]
      temp_pop = filter_samples(obj = obj, v = samples)
      
      ad = gsub('\\d+:', '', as.matrix(temp_pop@gt[low:high,]))
      
      ad1 = matrix(as.integer(gsub('/\\d+','',ad)), nrow = nrow(ad), ncol = ncol(ad))
      ad2 = matrix(as.integer(gsub('^\\d+|^\\d+/','',ad)), nrow = nrow(ad), ncol = ncol(ad))
      
      ad3 = matrix(rowSums(cbind(c(ad1),c(ad2)), na.rm = T), nrow = nrow(ad), ncol = ncol(ad))
      
      temp_ReadDepth_Summ = data.frame(Total_ReadDepth = rowSums(ad3, na.rm = T),
                                       mean_ReadDepth = apply(ad3, 1, function(x){mean(x, na.rm = T)})#,
                                       #sd_ReadDepth = apply(ad3, 1, function(x){sd(x, na.rm = T)}),
                                       #median_ReadDepth = apply(ad3, 1, function(x){median(x, na.rm = T)}),
                                       #quantile25_ReadDepth = apply(ad3, 1, function(x){quantile(x, probs = .25, na.rm = T)}),
                                       #quantile75_ReadDepth = apply(ad3, 1, function(x){quantile(x, probs = .75, na.rm = T)}),
                                       #iqr_ReadDepth = apply(ad3, 1, function(x){IQR(x, na.rm = T)})
                                       )
      
      names(temp_ReadDepth_Summ) = paste(names(temp_ReadDepth_Summ), pop, sep = "_")
      
      if(is.null(ReadDepth_Summ)){
        ReadDepth_Summ = temp_ReadDepth_Summ
      }else{
        ReadDepth_Summ = cbind(ReadDepth_Summ, temp_ReadDepth_Summ)
      }
      
    }
    
    # For the total population

    ad = gsub('\\d+:', '', as.matrix(obj@gt[low:high,]))
    
    ad1 = matrix(as.integer(gsub('/\\d+','',ad)), nrow = nrow(ad), ncol = ncol(ad))
    ad2 = matrix(as.integer(gsub('^\\d+|^\\d+/','',ad)), nrow = nrow(ad), ncol = ncol(ad))
    
    ad3 = matrix(rowSums(cbind(c(ad1),c(ad2)), na.rm = T), nrow = nrow(ad), ncol = ncol(ad))
    
    temp_ReadDepth_Summ = data.frame(Total_ReadDepth = rowSums(ad3, na.rm = T),
                                     mean_ReadDepth = apply(ad3, 1, function(x){mean(x, na.rm = T)})#,
                                     #sd_ReadDepth = apply(ad3, 1, function(x){sd(x, na.rm = T)}),
                                     #median_ReadDepth = apply(ad3, 1, function(x){median(x, na.rm = T)}),
                                     #quantile25_ReadDepth = apply(ad3, 1, function(x){quantile(x, probs = .25, na.rm = T)}),
                                     #quantile75_ReadDepth = apply(ad3, 1, function(x){quantile(x, probs = .75, na.rm = T)}),
                                     #iqr_ReadDepth = apply(ad3, 1, function(x){IQR(x, na.rm = T)})
                                     )
    
    names(temp_ReadDepth_Summ) = paste(names(temp_ReadDepth_Summ), 'Total', sep = "_")
    
    ReadDepth_Summ = cbind(ReadDepth_Summ, temp_ReadDepth_Summ)
    rownames(ReadDepth_Summ) = rownames(ad)
    
  }else{
    
    ad = gsub('\\d+:', '', as.matrix(obj@gt[low:high,]))
    
    ad1 = matrix(as.integer(gsub('/\\d+','',ad)), nrow = nrow(ad), ncol = ncol(ad))
    ad2 = matrix(as.integer(gsub('^\\d+|^\\d+/','',ad)), nrow = nrow(ad), ncol = ncol(ad))
    
    ad3 = matrix(rowSums(cbind(c(ad1),c(ad2)), na.rm = T), nrow = nrow(ad), ncol = ncol(ad))
    
    ReadDepth_Summ = data.frame(Total_ReadDepth = rowSums(ad3, na.rm = T),
                                mean_ReadDepth = apply(ad3, 1, function(x){mean(x, na.rm = T)})#,
                                #sd_ReadDepth = apply(ad3, 1, function(x){sd(x, na.rm = T)}),
                                #median_ReadDepth = apply(ad3, 1, function(x){median(x, na.rm = T)}),
                                #quantile25_ReadDepth = apply(ad3, 1, function(x){quantile(x, probs = .25, na.rm = T)}),
                                #quantile75_ReadDepth = apply(ad3, 1, function(x){quantile(x, probs = .75, na.rm = T)}),
                                #iqr_ReadDepth = apply(ad3, 1, function(x){IQR(x, na.rm = T)})
                                )
    
    rownames(ReadDepth_Summ) = rownames(ad)

  }
  
  return(ReadDepth_Summ)
  
}

# prune_alleles----

prune_alleles = function(obj, threshold = 4, n = 100){
  library(svMisc)
  s = round(seq(1,nrow(obj@gt)+1, length.out=n+1))
  gt6 = NULL
  
  for(w in 1:n){
    low = s[w]
    high = s[w+1]-1
  
    gt = gsub(':\\d+', '', as.matrix(obj@gt[low:high,]))
    gt1 = matrix(as.integer(gsub('/\\d+','',gt)), nrow = nrow(gt), ncol = ncol(gt))
    gt2 = matrix(as.integer(gsub('^\\d+|^\\d+/','',gt)), nrow = nrow(gt), ncol = ncol(gt))
  
    ad = gsub('\\d+:', '', as.matrix(obj@gt[low:high,]))

    ad1 = matrix(as.integer(gsub('/\\d+','',ad)), nrow = nrow(ad), ncol = ncol(ad))
    ad2 = matrix(as.integer(gsub('^\\d+|^\\d+/','',ad)), nrow = nrow(ad), ncol = ncol(ad))
  
    # prune alleles
    gt1[ad1 <= threshold] = NA
    ad1[ad1 <= threshold] = NA
  
    gt2[ad2 <= threshold] = NA
    ad2[ad2 <= threshold] = NA
  
    gt3 = matrix(paste(gt1, ad1, sep = ':'), nrow = nrow(gt), ncol = ncol(gt))
    gt3[gt3 == 'NA:NA'] = NA
  
    gt4 = matrix(paste(gt2, ad2, sep = ':'), nrow = nrow(gt), ncol = ncol(gt))
    gt4[gt4 == 'NA:NA'] = NA
  
    gt5 = matrix(paste(gt3, gt4, sep = '/'), nrow = nrow(gt), ncol = ncol(gt))
    gt5[gt5 == 'NA/NA'] = NA
    gt5 = gsub('/NA', '', gt5)
  
    colnames(gt5) = colnames(gt)
    rownames(gt5) = rownames(gt)
    
    gt6 = rbind(gt6, gt5)
    
    progress(round(100*w/n))
    
  }
  
  return(gt6)
  
  }

# Fws----

setGeneric("get_Fws", function(obj = NULL, w = 1, n = 1) standardGeneric("get_Fws"))

setMethod("get_Fws", signature(obj = "rGenome"),
          function(obj = NULL, w = 1, n = 1){
  
            gt = obj@gt
            loci = obj@loci_table
            
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
                                            end = chrom_intervals[[chrom]] - 1 + window))
    
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

filter_gt_matrix = function(gt, # haplotype matrix
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

fastGRM = function(obj, k, monoclonals = NULL, polyclonals = NULL, Pop = NULL, q = 2){
  
  gt = obj@gt
  loci = obj@loci_table
  metadata = obj@metadata
  
  
  X = handle_ploidy(gt = gt, w = 1, n = 1, monoclonals = monoclonals, polyclonals = polyclonals)
  
  
  X = matrix(as.numeric(X), ncol = ncol(X),
                                        nrow = nrow(X), 
                                        dimnames = list(
                                          rownames(X),
                                          colnames(X)
                                        ))
  
  X[is.na(X)] = 0
  
  evector = fastGRMCpp(X, k, q)
  
  #### Add metadata to the PCA----
  Pop_col = merge(data.frame(Sample_id = gsub('_C[1,2]$','',colnames(X)),
                             order = 1:ncol(X)), metadata[,c('Sample_id', Pop)], by = 'Sample_id', all.x = T)
  
  Pop_col = Pop_col[order(Pop_col$order),]
  
  evector = data.frame(Pop_col, evector)
  names(evector) = c(colnames(Pop_col), paste0(rep('PC', k), 1:k))
  
  return(evector)
  
}

# merge_rGenome----

merge_rGenome = function(obj1 = NULL, obj2 = NULL){
  
  # Check if the column Alleles is not present in the first data set
  if(sum(grepl('Alleles', names(obj1@loci_table))) == 0){ # if not generate Alleles and Allele_Counts
    # count Alleles
    Locus_info_temp1 = get_AC(obj1, update_alleles = T)
    Locus_info_temp1 = cbind(obj1@loci_table, Locus_info_temp1)
  }else{
    Locus_info_temp1 = obj1@loci_table[, c('CHROM',
                                           'POS',
                                           'REF', 
                                           'ALT',
                                           'Alleles',
                                           'Allele_Counts')]
  }
  
  # Check if the column Alleles is not present in the second data set
  if(sum(grepl('Alleles', names(obj2@loci_table))) == 0){ # if not generate Alleles and Allele_Counts
    # count Alleles
    Locus_info_temp2 = get_AC(obj2, update_alleles = T)
    Locus_info_temp2 = cbind(obj2@loci_table, Locus_info_temp2)
  }else{
    Locus_info_temp2 = obj2@loci_table[, c('CHROM',
                                           'POS',
                                           'REF', 
                                           'ALT',
                                           'Alleles',
                                           'Allele_Counts')]
  }
  
  # Relabel columns of the first data set
  
  Locus_info_temp1 = data.frame(locus_id = rownames(Locus_info_temp1), Locus_info_temp1)
  names(Locus_info_temp1) = c('locus_id',paste0(names(Locus_info_temp1)[-1], '.temp1'))
  
  # Relabel columns of the second data set
  Locus_info_temp2 = data.frame(locis_id = rownames(Locus_info_temp2), Locus_info_temp2)
  names(Locus_info_temp2) = c('locus_id',paste0(names(Locus_info_temp2)[-1], '.temp2'))
  
  
  # Store the genotype tables of both data sets
  gt_temp1 = obj1@gt
  gt_temp2 = obj2@gt
  
  # Merge the loci_table of both data sets
  Locus_info_merged = merge(Locus_info_temp1, Locus_info_temp2, by = 'locus_id', all = T)
  
  # Check if the reference allele is the same for both data sets in each position
  Locus_info_merged %<>% mutate(
    REF = case_when(
      REF.temp1 == REF.temp2 ~ REF.temp1,
      !is.na(REF.temp1) & is.na(REF.temp2) ~ REF.temp1,
      is.na(REF.temp1) & !is.na(REF.temp2) ~ REF.temp2,
      REF.temp1 != REF.temp2 ~ 'ERROR'
    )
  )
  
  # Check if the alternative alleles  and their labels are the same for both data sets in each position
  Locus_info_merged$Alleles = apply(Locus_info_merged, 1, function(pos){
    if(!is.na(pos['Alleles.temp1'])&is.na(pos['Alleles.temp2'])){ # If the second data set is all missing data
      paste0(pos['Alleles.temp1'], ';Only First data set') 
    }else if(is.na(pos['Alleles.temp1'])&!is.na(pos['Alleles.temp2'])){ # If the first data set is all missing data
      paste0(pos['Alleles.temp2'], ';Only Second data set') 
    }else if(is.na(pos['Alleles.temp1'])&is.na(pos['Alleles.temp2'])){ # If the first and second data sets are all missing data
      NA
    }else if(pos['Alleles.temp1'] == pos['Alleles.temp2']){ # If the alternative alleles  and their labels are the same for both data sets
      pos['Alleles.temp1']
    }else if(sum(!(unlist(str_split(pos['Alleles.temp2'], ',', simplify = T)) %in% unlist(str_split(pos['Alleles.temp1'], ',', simplify = T)))) == 0){ # If the alternative alleles  and their labels of second data set are contained in the first data set
      pos['Alleles.temp1']
    }else if(sum(!(unlist(str_split(pos['Alleles.temp1'], ',', simplify = T)) %in% unlist(str_split(pos['Alleles.temp2'], ',', simplify = T)))) == 0){ # If the alternative alleles  and their labels of first data set are contained in the second data set
      pos['Alleles.temp2']
    }else{'ERROR'} # If at least one alternative allele  and its label does not coincide between both data sets
  })
  
  # Split the loci table between:
  
  ## loci which reference alleles do not coincide between data sets
  Locus_info_merged_w_REFdiscrepancies = Locus_info_merged %>% filter(REF == 'ERROR')
  
  ## loci which reference alleles coincide but the aleternative alleles do not coincide between data sets
  Locus_info_merged_w_ALTdiscrepancies = Locus_info_merged %>% filter(Alleles == 'ERROR', REF != 'ERROR')
  
  ## loci which reference and alternative alleles coincide between data sets and do not have missing data
  Locus_info_merged_good = Locus_info_merged %>% filter(Alleles != 'ERROR' & REF != 'ERROR' & !grepl('Only',Alleles))
  
  ## loci with no amplification in one data set
  Locus_info_merged_w_missing = Locus_info_merged %>% filter(grepl('Only',Alleles) & REF != 'ERROR')
  
  # Add row names for each splited set of loci
  rownames(Locus_info_merged_w_REFdiscrepancies) = Locus_info_merged_w_REFdiscrepancies$locus_id
  rownames(Locus_info_merged_w_ALTdiscrepancies) = Locus_info_merged_w_ALTdiscrepancies$locus_id
  rownames(Locus_info_merged_good) = Locus_info_merged_good$locus_id
  rownames(Locus_info_merged_w_missing) = Locus_info_merged_w_missing$locus_id
  
  # Relabel each position in the set of loci that reference alleles coincide but the aleternative alleles do not coincide between data sets
  for(pos in rownames(Locus_info_merged_w_ALTdiscrepancies)){
    
    # Get Alleles, allele labels, and allele counts for the first data set
    Alleles.temp1 = as.character(unlist(str_split(gsub(':\\d+','',Locus_info_merged_w_ALTdiscrepancies[pos,][['Alleles.temp1']]), ',', simplify = T)))
    Allele_labels.temp1 = unlist(str_split(gsub('([ATGC]+|\\*):','',Locus_info_merged_w_ALTdiscrepancies[pos,][['Alleles.temp1']]), ',', simplify = T))
    Allele_counts.temp1 = as.integer(unlist(str_split(gsub('\\d+:','',Locus_info_merged_w_ALTdiscrepancies[pos,][['Allele_Counts.temp1']]), ',', simplify = T)))
    names(Alleles.temp1) = Allele_labels.temp1
    names(Allele_counts.temp1) = Allele_labels.temp1
    
    # Get Alleles, allele labels, and allele counts for the first data set
    Alleles.temp2 = as.character(unlist(str_split(gsub(':\\d+','',Locus_info_merged_w_ALTdiscrepancies[pos,][['Alleles.temp2']]), ',', simplify = T)))
    Allele_labels.temp2 = unlist(str_split(gsub('([ATGC]+|\\*):','',Locus_info_merged_w_ALTdiscrepancies[pos,][['Alleles.temp2']]), ',', simplify = T))
    Allele_counts.temp2 = as.integer(unlist(str_split(gsub('\\d+:','',Locus_info_merged_w_ALTdiscrepancies[pos,][['Allele_Counts.temp2']]), ',', simplify = T)))
    names(Alleles.temp2) = Allele_labels.temp2
    names(Allele_counts.temp2) = Allele_labels.temp2
    
    # Get the reference and unique alternative alleles observed in both data sets
    REF.Allele = unique(c(Alleles.temp1[names(Alleles.temp1)[names(Alleles.temp1) == '0']], Alleles.temp2[names(Alleles.temp2)[names(Alleles.temp2) == '0']]))
    ALT.Alleles = unique(c(Alleles.temp1[names(Alleles.temp1)[names(Alleles.temp1) != '0']], Alleles.temp2[names(Alleles.temp2)[names(Alleles.temp2) != '0']]))
    
    # If reference allele is not present in either set
    if(length(REF.Allele) == 0){
      REF.Allele = Locus_info_merged_w_ALTdiscrepancies[pos,][['REF']]
    }
    
    # Calculate the total number of samples that have each ealternative allele in both data sets
    ALT.Allele_counts = NULL
    
    for(allele in ALT.Alleles){ # For each alternative allele
      
      ALT.Allele_counts = c(ALT.Allele_counts,
                            
                            if(allele %in% Alleles.temp1 & allele %in% Alleles.temp2){ # if the allele is present in both data sets
                              Allele_counts.temp1[names(Alleles.temp1[Alleles.temp1 == allele])] +  
                                Allele_counts.temp2[names(Alleles.temp2[Alleles.temp2 == allele])]
                            }else if(allele %in% Alleles.temp1 & !(allele %in% Alleles.temp2)){ # if allele is present only in the first data set
                              Allele_counts.temp1[names(Alleles.temp1[Alleles.temp1 == allele])]
                            }else if(!(allele %in% Alleles.temp1) & allele %in% Alleles.temp2){ # if allele is present only in the second data set
                              Allele_counts.temp2[names(Alleles.temp2[Alleles.temp2 == allele])]
                            }
      )
    }
    
    # name (index) the allele count with the allele
    names(ALT.Allele_counts) = ALT.Alleles
    
    # sort the alternative alleles based on their allele count
    ALT.Allele_counts = sort(ALT.Allele_counts, decreasing = T)
    
    # Calculate the total number of samples that have the reference allele in both data sets
    REF.Allele_count = if(REF.Allele %in% Alleles.temp1 & REF.Allele %in% Alleles.temp2){ # if the allele is present in both data sets
      Allele_counts.temp1[names(Alleles.temp1[Alleles.temp1 == REF.Allele])] +  
        Allele_counts.temp2[names(Alleles.temp2[Alleles.temp2 == REF.Allele])]
    }else if(REF.Allele %in% Alleles.temp1 & !(REF.Allele %in% Alleles.temp2)){ # if allele is present only in the first data set
      Allele_counts.temp1[names(Alleles.temp1[Alleles.temp1 == REF.Allele])]
    }else if(!(REF.Allele %in% Alleles.temp1) & REF.Allele %in% Alleles.temp2){ # if allele is present only in the second data set
      Allele_counts.temp2[names(Alleles.temp2[Alleles.temp2 == REF.Allele])]
    }
    
    if(!is.null(REF.Allele_count)){
      
      # name (index) the allele count with the allele
      names(REF.Allele_count) = REF.Allele
      # Combine the reference and the alternative alleles in one object
      Allele_counts = c(REF.Allele_count, ALT.Allele_counts)
      Alleles = names(Allele_counts)
      
    }else{
      
      Allele_counts = ALT.Allele_counts
      Alleles = names(Allele_counts)
      
    }
    
    # RELABEL the alleles for the combined data set
    
    if(!is.null(REF.Allele_count)){
      Allele_labels = as.character(0:(length(Alleles) - 1))
      names(Alleles) = Allele_labels
    }else{
      
      Allele_labels = as.character(1:length(Alleles))
      names(Alleles) = Allele_labels
    }
    
    # If the alternative alleles  and their labels of First data set are NOT contained in the merged data set
    if(sum(!(paste(Alleles.temp1, Allele_labels.temp1, sep = ':') %in% paste(Alleles, Allele_labels, sep = ':'))) != 0){
      
      # Identify and select the alleles which lables have changed
      Changed.Alleles.temp1 = Alleles.temp1[!(paste(Alleles.temp1, Allele_labels.temp1, sep = ':') %in% paste(Alleles, Allele_labels, sep = ':'))]
      
      # Identify the new labels of the alleles
      Relabeled.Alleles.temp1 = Alleles[Alleles %in% Changed.Alleles.temp1]
      
      # For each changed label, modify the label in the genotype table
      for(clabel in names(Changed.Alleles.temp1)){
        
        rlabel = names(Relabeled.Alleles.temp1[Relabeled.Alleles.temp1 == Changed.Alleles.temp1[clabel]])
        
        gt_temp1[pos,] = gsub(paste0(clabel,':'), paste0(rlabel,'R:'), gt_temp1[pos,]) # The R: denotes that that label is been modify and avoids overwriting
        
      }
      
      # Once screening of all alleles is been completed, delete the R: identifier
      gt_temp1[pos,] = gsub('R:', ':', gt_temp1[pos,])
      
    }
    
    # If the alternative alleles  and their labels of Second data set are NOT contained in the merged data set
    if(sum(!(paste(Alleles.temp2, Allele_labels.temp2, sep = ':') %in% paste(Alleles, Allele_labels, sep = ':'))) != 0){
      
      # Identify and select the alleles which lables have changed
      Changed.Alleles.temp2 = Alleles.temp2[!(paste(Alleles.temp2, Allele_labels.temp2, sep = ':') %in% paste(Alleles, Allele_labels, sep = ':'))]
      
      # Identify the new labels of the alleles
      Relabeled.Alleles.temp2 = Alleles[Alleles %in% Changed.Alleles.temp2]
      
      # For each changed label, modify the label in the genotype table
      for(clabel in names(Changed.Alleles.temp2)){
        
        rlabel = names(Relabeled.Alleles.temp2[Relabeled.Alleles.temp2 == Changed.Alleles.temp2[clabel]])
        
        gt_temp2[pos,] = gsub(paste0(clabel,':'), paste0(rlabel,'R:'), gt_temp2[pos,]) # The R: denotes that that label is been modify and avoids overwriting
        
      }
      
      # Once screening of all alleles is been completed, delete the R: identifier
      gt_temp2[pos,] = gsub('R:', ':', gt_temp2[pos,])
      
    }
    
    # Update Alleles in the Locus_info_merged_w_ALTdiscrepancies data.frame
    Locus_info_merged_w_ALTdiscrepancies[pos,][['Alleles']] = paste(paste(Alleles, Allele_labels,sep = ':'), collapse = ',')
    
  }
  
  # Filter positions that were able to merge
  Locus_info_merged_final = rbind(Locus_info_merged_good, Locus_info_merged_w_ALTdiscrepancies)
  
  # Update Alternative Alleles (ALT)
  Locus_info_merged_final$ALT  = apply(Locus_info_merged_final, 1, function(ALT){
    gsub(':\\d+', '', gsub(paste0('^',paste(ALT['REF'], '0', sep = ':'), ','), '', ALT['Alleles']))
  })
  
  # Sort positions by position and chromosome
  Locus_info_merged_final = Locus_info_merged_final[order(Locus_info_merged_final$POS.temp1),]
  Locus_info_merged_final = Locus_info_merged_final[order(Locus_info_merged_final$CHROM.temp1),]
  
  # Select columns that will be returned
  Locus_info_merged_final = Locus_info_merged_final[,c('CHROM.temp1',
                                                       'POS.temp1',
                                                       'REF',
                                                       'ALT',
                                                       'Alleles')]
  
  colnames(Locus_info_merged_final) = c('CHROM',
                                        'POS',
                                        'REF',
                                        'ALT',
                                        'Alleles')
  
  # Filter and merge the genotype tables
  gt_final = cbind(gt_temp1[rownames(Locus_info_merged_final),], gt_temp2[rownames(Locus_info_merged_final),])
  
  # Merge the metadata
  merged_metadata = merge(obj1@metadata, obj2@metadata, by = 'Sample_id', all = T)
  
  rownames(merged_metadata) = merged_metadata$Sample_id
  
  merged_metadata = merged_metadata[colnames(gt_final),]
  
  # Create the merged rGenome object
  merged_rGenome = rGenome(gt = gt_final,
                           loci_table = Locus_info_merged_final,
                           metadata = merged_metadata)
  
  return(merged_rGenome)
  
}

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

# get_haplotypes_respect_to_reference----

get_haplotypes_respect_to_reference = function(obj,
                                               gene_ids = NULL,
                                               gene_labels = NULL,
                                               gff_file = NULL,
                                               fasta_file = NULL,
                                               monoclonals = NULL,
                                               polyclonals = NULL,
                                               variables = NULL){
  library(ape)
  library(Biostrings)
  library(msa)
  
  # Call reference genome and its corresponding anotation in the gff file
  reference_gff = read.gff(gff_file)
  reference_genome = readDNAStringSet(fasta_file)
  
  # Get gene ids, names, and descriptions 
  gene_names = data.frame(gene_id = gsub(';.+$','',gsub('ID=','',reference_gff %>% filter(grepl(paste(gene_ids, collapse = '|'), attributes),  grepl('gene', type)) %>% select(attributes) %>% unlist())),
                          gene_name = gsub(';.+$','',gsub('^.+Name=','',reference_gff %>% filter(grepl(paste(gene_ids, collapse = '|'), attributes),  grepl('gene', type)) %>% select(attributes) %>% unlist())),
                          gene_description = gsub(';.+$','',gsub('^.+description=','',reference_gff %>% filter(grepl(paste(gene_ids, collapse = '|'), attributes),  grepl('gene', type)) %>% select(attributes) %>% unlist())))
  
  # Split gene coordinates by coding sequences (CDS)
  genes_gff = reference_gff %>% filter(grepl(paste(gene_ids, collapse = '|'), attributes) & type == 'CDS')
  
  # Get gene ids
  genes_gff %<>% mutate(gene_id = gsub(';.+$|\\..+$','',gsub('^.+gene_id=|^ID=','',attributes)))
  
  # Merge CDS coordinates with gene names and descriptions by gene_id
  genes_gff = merge(genes_gff, gene_names, by = 'gene_id', all.x = T)
  
  # For each gene:
  
  Sample_id = obj@metadata$Sample_id
  
  monoclonals_ids = Sample_id[Sample_id %in% monoclonals]
  polyclonals_ids = Sample_id[Sample_id %in% polyclonals]
  
  haplotypes_ids = c(monoclonals_ids,
                    paste0(polyclonals_ids, '_C1'),
                    paste0(polyclonals_ids, '_C2'))
  
  haplotypes = matrix(NA, ncol = length(gene_ids), nrow = length(haplotypes_ids),
                     dimnames = list(haplotypes_ids, gene_ids))
  
  for(gene in unique(genes_gff$gene_id)){
    
    # Filter polymorphic positions located in the gene of interest (GOI)
    ## Get the cds coordinates of the goi
    cds_gff = genes_gff[genes_gff$gene_id == gene,]
    ## Get the chromosome sequence where the goi is located
    ref_seq = reference_genome[[grep(unique(cds_gff$seqid), names(reference_genome))]]
    
    # Generate a vector with all DNA coordinates in the chromosome
    dna_regions = NULL
    gene_seq = NULL
    
    for(cds in 1:nrow(cds_gff)){
      
      # Get nucleotide coordinates
      dna_regions = c(dna_regions,
                      paste(cds_gff[cds,][['seqid']],
                            cds_gff[cds,][['start']]:cds_gff[cds,][['end']],
                            sep = '_'))
      
      # Get DNA sequence
      gene_seq = paste0(gene_seq, as.character(subseq(ref_seq, start = cds_gff[cds,][['start']], end = cds_gff[cds,][['end']])))
      
    }
    
    names(gene_seq) = 'reference'
    
    if(sum(rownames(obj@loci_table) %in% dna_regions) > 0){
      
      # Filter the rGenome object based on the gene coordinates
      
      gene_obj = filter_loci(obj, v = rownames(obj@loci_table) %in% dna_regions)
      
      ## Remove read abundance
      
      gt = handle_ploidy(gt = gene_obj@gt, monoclonals = monoclonals, polyclonals = polyclonals)
      #gt = handle_ploidy(gene_obj@gt, monoclonals = NULL, polyclonals = c(monoclonals, polyclonals))
      
      # Tarnsform haplotype codes to nucleotides
      nuc_gt = NULL
      for(pos in 1:nrow(gt)){
        temp_nuc_gt = gt[pos,]
        allele_codes = unique(unlist(strsplit(gt[pos,], '/')))
        allele_codes = min(as.integer(allele_codes), na.rm = T):max(as.integer(allele_codes), na.rm = T)
        alleles = unlist(strsplit(gsub(':\\d+','',gene_obj@loci_table[pos,"Alleles"]), ','))
        for(allele in allele_codes + 1){
          temp_nuc_gt = gsub(allele_codes[allele],alleles[allele],temp_nuc_gt)
        }
        
        temp_nuc_gt = gsub('\\*', '',temp_nuc_gt)
        nuc_gt = rbind(nuc_gt, temp_nuc_gt)
      }
      rownames(nuc_gt) = rownames(gt)
      
      ref_polymorphims_length = nchar(gene_obj@loci_table$REF)
      
      # Define the coordinates of the polymorphims in the gene sequence
      polymorphic_positions = which(dna_regions %in% rownames(nuc_gt))
      
      # For each sample for each codon define the aminoacid changes in the gene
      
      sample_seqs = NULL
      
      for(sample in 1:ncol(nuc_gt)){
        
        sample_seq = gene_seq
        
        for(pos in 1:length(polymorphic_positions)){
          
          substr(sample_seq, polymorphic_positions[pos], polymorphic_positions[pos] + ref_polymorphims_length[pos] - 1) = nuc_gt[pos,sample]
          
        }
        
        sample_seqs = c(sample_seqs, sample_seq)
        
      }
      
      names(sample_seqs) = colnames(nuc_gt)
      
      sample_seqs = sample_seqs[!is.na(sample_seqs)]
      
      if(cds_gff$strand[1] == '+'){
        
        aa_seqs = AAStringSet(c(as.character(translate(DNAString(gene_seq))),
                                as.character(translate(DNAStringSet(sample_seqs)))))
        
      }else if(cds_gff$strand[1] == '-'){
        
        aa_seqs = AAStringSet(c(as.character(translate(reverseComplement(DNAString(gene_seq)))),
                                as.character(translate(reverseComplement(DNAStringSet(sample_seqs))))))
        
      }
      
      names(aa_seqs) = c('reference', names(sample_seqs))
      
      aa_alignment = msa(aa_seqs)
      
      aa_alignment_matrix = matrix(unlist(strsplit(as.character(aa_alignment@unmasked), '')), ncol = length(aa_alignment@unmasked),
                                   nrow = nchar(as.character(aa_alignment@unmasked[['reference']])),
                                   dimnames = list(
                                     1:nchar(as.character(aa_alignment@unmasked[['reference']])),
                                     names(aa_alignment@unmasked)
                                   )
      )
      
      aa_polymorphic_positions = which(rowSums(aa_alignment_matrix[,'reference'] != aa_alignment_matrix) > 0)
      
      if(length(aa_polymorphic_positions) > 0){
        
        aa_haplotypes = NULL
        
        for(sample in names(sample_seqs)){
          
          ref_aa_aligned = aa_alignment_matrix[,'reference']
          
          samp_aa_aligned = aa_alignment_matrix[,sample]
          
          haplotype = paste(paste0(ref_aa_aligned[aa_polymorphic_positions], aa_polymorphic_positions, samp_aa_aligned[aa_polymorphic_positions]), collapse = ' ')
          
          names(haplotype) = sample
          
          aa_haplotypes = c(aa_haplotypes, haplotype)
          
        }
        
        haplotypes[rownames(haplotypes) %in% names(aa_haplotypes), gene] = aa_haplotypes
        
        
      }else{
        
        haplotypes[rownames(haplotypes) %in% names(sample_seqs), gene] = '.'
        
      }
      
    }else{
      
      haplotypes[, gene] = '.'
      
    }
    
  }
  
  haplo_freqs = data.frame(Sample_id = rownames(haplotypes), haplotypes)
  
  metadata = obj@metadata[,c('Sample_id', variables)]
  colnames(metadata) = c('Sample_id', 'Var1', 'Var2')
  
  metadata_monoclonals = metadata[monoclonals,]
  metadata_polyclonals = metadata[polyclonals,]
  metadata_polyclonals1 = metadata_polyclonals
  metadata_polyclonals1$Sample_id = paste0(metadata_polyclonals1$Sample_id, '_C1')
  metadata_polyclonals2 = metadata_polyclonals
  metadata_polyclonals2$Sample_id = paste0(metadata_polyclonals2$Sample_id, '_C2')
  
  metadata = rbind(metadata_monoclonals, metadata_polyclonals1, metadata_polyclonals2)
  
  haplo_freqs = merge(haplo_freqs, metadata, by = 'Sample_id')
  
  haplo_freqs %<>% pivot_longer(cols = colnames(haplotypes), names_to = 'Gene_id', values_to = 'Haplotype')
  
  gene_labels = data.frame(gene_ids, gene_labels)
  
  haplo_freqs$gene_label = NA
  
  for(gene in gene_labels$gene_ids){
    
    haplo_freqs[haplo_freqs$Gene_id == gene,][['gene_label']] = gene_labels[gene_labels$gene_ids == gene,][['gene_labels']]
    
  }
  
  haplo_freqs %<>% 
    mutate(Haplotype = paste(gene_label, Haplotype, sep = ':'))%>%
    group_by(Var1, Var2, gene_label) %>%
    summarise(Sample_id = Sample_id,
              Haplotype = Haplotype,
              Sample_size = n()
              ) %>% group_by(Var1, Var2, gene_label, Haplotype) %>%
    summarise(Count = n(),
              Freq = n()/unique(Sample_size))
  
  haplo_freq_plot = haplo_freqs%>%
    ggplot(aes(y = Freq, x = Var2, fill = Haplotype))+
    geom_col(position = 'stack')+
    facet_grid(Var1~gene_label)+
    labs(x = variables[2])+
    theme(legend.position = 'bottom')
  
  
  return(list(Haplotypes = haplotypes,
              haplo_freqs = haplo_freqs,
              haplo_freq_plot = haplo_freq_plot))
  
}

