## Functions for code to replicate results in 
# "Quantifying the global film festival circuit: Networks, diversity, and public value creation"
# by Vejune Zemaityte, Andres Karjus, Ulrike Rohn, Maximilian Schich, Indrek Ibrus
# Codebase (c) Andres Karjus 2023, andreskarjus.github.io

# This file is meant to be sourced from the main script file.

## Note: internal diversity is the same label as in the paper
# but external is called "contributing diversity" in the paper now,
# and system is referred to as "(external) circuit diversity" in the paper.



#### Load packages; will install packages if needed ####
# Make sure reticulate works (requires a Python distro; will attempt to install miniconda), 
# otherwise can't do language diversity

p = c("colorspace","scales" ,"umap","text2vec" ,"shadowtext","geosphere" , "stringr", "Hmisc" ,"igraph", "ggraph", "ggbeeswarm" ,"maps","rworldmap" ,"boot","ISOcodes" ,"svs","tidyverse" , "ggrepel","patchwork" ,"reticulate")
px = setdiff(p, installed.packages()[,"Package"])
if(length(px)>0){ install.packages(px) } 
invisible(lapply(p, require, character.only = TRUE)); try(rm(p, px))


# Try to install miniconda if not present - debug if needed
if(is.null(reticulate::conda_version())){install_miniconda()}
# Try to install required lang2vec package if not present
if(!reticulate::py_module_available("lang2vec")){
  conda_install(envname = 'r-reticulate', c('lang2vec'), pip = TRUE, python_version="3.6")
}
tryCatch({use_condaenv('r-reticulate', required = T ); py_run_string("import lang2vec.lang2vec")}, error=function(e){print("Failed to install and/or load Conda and/or lang2vec; so the language section below won't run. Please try debugging in the scripts file, or comment out this section in the scripts file. Other parts don't require Python.")})


#### Package versions ####

# We used the following versions:
# boot          1.3.28
# colorspace     2.0.3
# geosphere     1.5.14
# ggbeeswarm     0.6.0
# ggraph         2.1.0
# ggrepel        0.9.1
# Hmisc          5.0.1
# igraph         1.3.4
# ISOcodes   2022.9.29
# maps           3.4.0
# patchwork      1.1.1
# reticulate      1.25
# rworldmap      1.3.6
# scales         1.2.0
# shadowtext     0.1.2
# stringr        1.5.0
# svs            3.0.0
# text2vec       0.6.1
# tidygraph      1.2.2
# tidyverse      2.0.0
# umap         0.2.8.0



# --------------------------------------------
# Functions

#### Data wrangling, latent space and diversity calculation functions ####

disttomean=function(x, m, w, maxd=NULL, geo=F, divmult=2){
  # with bootstrapping
  if(!geo){
    d = dist2(x %>% select(starts_with("V")) %>% as.matrix,
              m,   #  festival or grand mean
              method = "euclidean", norm = "none") %>% 
      {.[,1]/maxd*divmult} # scale by max latent dist, indiv multiply*2 so scale is [0,1]
  } else { # geo
    d = distGeo(x %>% select(long, lat), m)/1000 # distance from geo-mean centre
  }
  div = weighted.mean(d,w, na.rm=T)
  ## bootstrapped confidence intervals
  if(length(unique(d))==1){
    bi=c(div,div) # if no variance/constant then bootstrap would fail anyway.
  } else {
    bi=
      boot(d, function(dx,y){weighted.mean(dx[y], w=w[y])}, R=1000) %>% # multicore won't work on Win
      boot.ci(type = "norm") %>%
      {.$normal[2:3]}
  }
  return(c(div, bi))
}

# for change calc
meanvecs=function(x, inds, ems, loopindex, maxdist){
  xt = x[inds,] %>% filter(EventID_VZ==ems$EventID_VZ[loopindex]) %>% 
    select(starts_with("V")) %>%
    colMeans() %>% matrix(nrow=1)
  xp = x[inds,] %>% filter(EventID_VZ==ems$EventID_VZ[loopindex-1]) %>%  # prev year
    select(starts_with("V")) %>% 
    colMeans() %>% matrix(nrow=1)
  return(dist(rbind(xt, xp))[1]/maxdist) 
  # not multiplied by 2, as it's comparison to an external mean, like the external diversity calculation
}


# round_any = function(x, accuracy, f=round){f(x/ accuracy) * accuracy}



dolatentspaces = function(metatype, filmgen=NULL){
  
  if(metatype=="genre"){
    require(svs)
    # the only inferred latent space so far
    genmat = filmgen %>% filter(txtKind!="unknown", !is.na(txtKind)) %>% 
      select(refFilm, txtKind) %>% mutate(tmp=1) %>% 
      pivot_wider(values_from=tmp, values_fill = 0, id_cols = "refFilm", names_from = "txtKind") %>% 
      column_to_rownames("refFilm") %>% as.matrix() %>% 
      as.dfm() %>% fcm() %>% 
      #as( "generalMatrix") %>%   #
      as.matrix() %>%  # deprecated
      {.[lower.tri(., diag = F)]=t(.)[lower.tri(., diag = F)];.}
    diag(genmat) = (filmgen %>% group_by(refFilm) %>% filter(n()==1, NoOrder=="1") %>% 
                      ungroup() %>% count(txtKind) %>% 
                      column_to_rownames("txtKind"))[rownames(genmat),"n"] #  counts of single-tags i.e. occurs with itself, single-genre films. helps give meaning to e.g. drama.
    pgenmat=pmi(genmat) %>% {.[.<0]=0;.} # Positive PMI
    # fcm method gives only upper triangle, which is fine below for fest-fest (no double links needed), but for vectors need full matrix
    m = LatentSemanticAnalysis$new(10)
    lat = m$fit_transform(pgenmat) 
  }
  
  # Geo.
  # all distances need to be haversine/great circle or ellipsoid analogue.
  # festival is film-countries weighted geo-average location geomean()
  # internal diversity is weighted mean of distance matrix of geo distances (in meters). weighs are from films (a country can occur multiple times in the dist matrix), for the w.mean the weights are multiplied (so 0.5 and 0.5 become a 0.25 for a given comparison).
  # weights %o% outer product
  # distm with default distGeo is spherical distance on ellipsoid (more accurate than haversine which assumes perfect sphere).
  if(metatype=="geo"){
    require(maps) # for cities
    lat = world.cities %>% # from package
      mutate(txtKind=tolower(country.etc), name=tolower(name)) %>% 
      filter(capital==1 |  name %in% c("montenegro", "los angeles", "gazzah", "pristina", "podgorica")) %>%
      arrange(-pop) %>%   # LA pop > Washington, so keeps that
      filter(!duplicated(txtKind)) %>%   
      select(txtKind, long, lat)  %>% 
      column_to_rownames( var = "txtKind")
  }
  
  if(metatype=="lang"){
    require(reticulate)
    # conda_install(envname = 'r-reticulate', c('lang2vec'), pip = TRUE, python_version="3.6")
    py_run_string("
import lang2vec.lang2vec as l2v
import numpy as np
np.load.__defaults__=(None, True, True, 'ASCII') # fix numpy bug
"
    )
    py_run_string("langs = list(l2v.URIEL_LANGUAGES)")
    py_run_string(
      "
feats = l2v.get_features(list(l2v.URIEL_LANGUAGES), l2v.fs_concatenation(['syntax_knn', 'phonology_knn','inventory_knn','fam', 'geo']), header=False, minimal=False)
# featsl = l2v.get_features(r.isos, 'learned', header=False, minimal=False)
"
    )
    feats = do.call(rbind,  py$feats)
    rn=rownames(feats)
    feats=apply(feats #[intersect(filmlang2$iso, rn),]
                ,2, as.numeric) 
    rownames(feats)=rn # intersect(filmlang2$iso, rn)
    feats2 = feats %>% .[,apply(.,2,function(x) !all(x==0) & !all(x==1) )]; dim(feats2)  # 91 x 793
    
    # feats2 %>% reshape2::melt() %>% ggplot(aes(Var2, Var1, fill=value))+geom_tile()
    # d=dist(feats2); d[is.na(d)]=max(d, na.rm=T)
    # plot(hclust(d))
    # sort(as.matrix(d)["fra",]) %>% head
    # works!
    
    # reduce very long mostly binary feature vecs to a latent space
    m = LatentSemanticAnalysis$new(20)
    lat = m$fit_transform(feats2) 
    
    # sim2( langlat, langlat["est",,drop=F])[,1] %>% sort(T) %>% head(10)
    # dist2( langlat, langlat["est",,drop=F], method="euclidean", norm="none")[,1] %>% {names(.)=rownames(langlat);.} %>% sort %>% head(10)
    # dist2( genlat, genlat["Comedy",,drop=F], method="euclidean", "none") %>% {rownames(.)=rownames(genlat);.} %>% .[,1] %>% sort
    
  }
  
  if(metatype %in% c("lang", "genre")){
    latmap = umap(lat)$layout %>% as.data.frame() %>%
      mutate(txtKind=rownames(lat))
    if(metatype=="genre"){
      latmap=latmap %>% left_join(count(filmgen, txtKind), by="txtKind")
    }
    attr(lat, "latmap")=latmap
  }
  return(lat)
  
}



countryfixes = 
  tribble(
    ~txtKind, ~txtKind2,
    "united kingdom"  , "uk",
    "korea (south)", "korea south",
    "korea (north)", "korea north",
    "birmanie", "myanmar",
    "hong kong (china)", "china",
    "trinidad  & tobago", "trinidad and tobago",
    "central african rep.", "central african republic",
    "acores", "portugal",
    "azores" , "portugal",
    "macao", "china",
    "guernsey", "uk",
    "canary islands" , "spain",
    "kurdistan" , "iraq" , # their capital is on iraqi territory; so for matching purposes relabeling
    "curacao" , "venezuela" # dutch territory island but next to VE
  ) 
langfixes=
  tribble(
    ~txtKind, ~txtKind2,
    "belarussian"  , "belarusian",
    "ukranian", "ukrainian",
    "albanian" , "tosk albanian", # need one
    "flemish","dutch",
    "chinese", "mandarin chinese",
    "mandarin", "mandarin chinese",
    "mandarin (traditional)","mandarin chinese",
    "taiwanese", "mandarin chinese",
    "mexican","spanish",
    "farsi" , "western farsi",
    "slovakian" , "czech",  # close enough
    "greek" , "modern greek (1453-)",
    "cantonese" , "yue chinese",
    "persian","iranian persian",
    "hazaragi","iranian persian",
    "swahili","swahili (individual language)",
    "australian","english",
    "creole","haitian",
    "gaelic","scottish gaelic",
    "scottish" , "scots",
    # "inuktitut" , "inuinnaqtun", # no vector
    "malay","standard malay",
    "malaysian","standard malay",
    "sinhale","sinhala",
    "sesotho", "southern sotho",
    "berber","tarifit",  # it's a family, just picking one
    "tunisian","tunisian arabic",
    "azerbaijani","south azerbaijani",
    "mongole" , "mongolia buriat",
    "algerian", "algerian arabic",
    "moroccan","moroccan arabic",
    "moor" , "egyptian arabic",
    "moore" , "egyptian arabic",
    "arabic" , "egyptian arabic",  # need to pick one
    "arabic-egypt", "egyptian arabic",
    "arabic-lebanese", "north levantine arabic",
    "uzbek","northern uzbek",
    "serbo-croatian" , "serbian",
    "filipino" , "tagalog",
    "guarani" , "paraguayan guaraní",
    "fanagalo" , "zulu",  # a zulu-based pidgin
    "pashtu", "central pashto",
    "hassanya", "hassaniyya"
  ) 
# any(duplicated(langfixes[,1])) # ok

fixdat = function(dat, metatype, filmkinds, fixes=NULL){
  # might need to rewrite to be able to include all 2 meta kinda in a single object; for now it's ok.
  
  tmp = dat %>%  # full fest * film * genre vecs for exdiversity calc
    ungroup() %>% 
    left_join(filmkinds, by="refFilm", multiple = "all") %>%    # ! intentionally multiplies rows for each genre per film
    filter(!is.na(txtKind), !(txtKind %in% c("unknown", "Unknown"))) %>%    # remove empty and unk if any
    mutate(originalKind = txtKind)  # keep original
  
  if(metatype %in% c("geo", "lang")){ # need fixing and adding of external labels
    
    tmp2 = tmp %>% 
      mutate(txtKind=tolower(txtKind)) %>% 
      left_join(fixes, by="txtKind") %>%   # just joins fixes, no duplicates here, 2 is new
      # adds already ok tags after fixing the ones that needed fixing:
      mutate(txtKind2=case_when(is.na(txtKind2)~txtKind, T~txtKind2 )) # prepare for matching
    
    # event country
    if(metatype=="geo"){
      tmp3 = tmp2 %>% 
        mutate(txtKind3 = txtKind2) %>%  # already matches lat object
        # but add event country too:
        mutate(txtCountryEvent_VZ=tolower(txtCountryEvent_VZ)) %>% 
        left_join(fixes %>% rename(txtCountryEvent_VZ=txtKind, 
                                   eventcountry2=txtKind2), 
                  by="txtCountryEvent_VZ") %>%    # just joins fixes, no duplicates here
        
        mutate(eventcountry3=case_when(is.na(eventcountry2)~txtCountryEvent_VZ, 
                                       T~eventcountry2 )) 
    }
    
    # add labels from external iso object to be able to match languages
    if(metatype=="lang"){
      require(ISOcodes)
      tmp3 = tmp2 %>% 
        mutate(txtKind3 = ISO_639_3$Id[pmin(match(txtKind2, tolower(ISO_639_3$eng)), 
                                            match(txtKind2, tolower(ISO_639_3$Name)), na.rm = T)]) %>% 
        mutate(isoname = ISO_639_3$Name[pmin(match(txtKind2, tolower(ISO_639_3$eng)), 
                                             match(txtKind2, tolower(ISO_639_3$Name)), na.rm = T)]) %>% 
        filter(!is.na(txtKind3), !(txtKind3 %in% c("iku", "osi")))  # won't match
    }
  } 
  if(metatype %in% c("genre", "role") ){
    tmp3=tmp %>% mutate(txtKind2=txtKind, txtKind3=txtKind) # already matches, but for compatibility
  }
  
  return(tmp3)
}

dofestvecs=function(fesfilm2, metatype, sortvar, lat, filmkinds, fixes=NULL, 
                    minfilms=15, customlabs=NULL){
  # make sure no variable starts with "V", as this is reserved to filter latent vectors
  dat = fesfilm2 %>% 
    rename(nsort:={{sortvar}}) %>%  # variable passing now works
    mutate(tmp=paste0(TitleVA, "_", YearOfProduction )) %>% 
    group_by(EventID_VZ) %>%   # for duplicate film removal
    arrange(desc(nsort), is.na(YearOfProduction), .by_group = T) %>%  # prefer more complete entries
    filter(!duplicated(tmp)) %>%  # remove duplicate films within event (works, w ungroup fewer rows)
    group_by(EventID_VZ) %>% filter(n()>=minfilms) %>% ungroup()  # keep events with 15 or more films; but metadata specific analyses will lower that a bit, as not all films have metadata.
  
  # role/gender works differently, no fixing and space needed, already a diversity measure.
  # does need weighting, doing it here (mostly 1-4 but some have 20-30 dir/prod credits)
  if(metatype=="role"){
    d = dat %>% group_by(EventID_VZ) %>% filter(n()>=minfilms) # just to make faster
    # instead of joining everything, do calc first
    gratios = filmkinds %>% 
      filter(refFilm %in% d$refFilm) %>% 
      group_by(refFilm) %>% 
      summarise(
        ratiof=sum(txtKind=="F")/sum(txtKind %in% c("F", "M")),
        ratiofmax = sum(txtKind %in% c("F", "other")) / n(), # max if all unk would be F
        ratiofmin = sum(txtKind=="F") / n() ,   # min if all unk would be M
        nother=round(sum(txtKind=="other")/n(), 5)  # ratio of unknowns per film; round just in case to avoid any chance of floating point error, for filter below
      )
    
    dratios = filmkinds %>% 
      filter(refFilm %in% d$refFilm) %>% 
      filter(txtType == "Director") %>% 
      group_by(refFilm) %>% 
      summarise(
        ratiodf=sum(txtKind=="F")/sum(txtKind %in% c("F", "M")),
        ratiodfmax = sum(txtKind %in% c("F", "other")) / n(), # max if all unk would be F
        ratiodfmin = sum(txtKind=="F") / n()    # min if all unk would be M
        #dnother=round(sum(txtKind=="other")/n(), 5)  # ratio of unknowns per film; round just in case to avoid any chance of floating point error, for filter below
      )
    
    
    pratios = filmkinds %>% 
      filter(refFilm %in% d$refFilm) %>% 
      filter(txtType == "Producer") %>% 
      group_by(refFilm) %>% 
      summarise(
        ratiopf=sum(txtKind=="F")/sum(txtKind %in% c("F", "M")),
        ratiopfmax = sum(txtKind %in% c("F", "other")) / n(), # max if all unk would be F
        ratiopfmin = sum(txtKind=="F") / n()    # min if all unk would be M
        #pnother=round(sum(txtKind=="other")/n(), 5)  # ratio of unknowns per film; round just in case to avoid any chance of floating point error, for filter below
      )
    
    filmfestvecs = d %>% 
      left_join(gratios, by = "refFilm") %>% # does not multiply rows
      left_join(dratios, by = "refFilm") %>% 
      left_join(pratios, by = "refFilm") %>% 
      filter(!is.na(ratiof), nother<1) %>%   # some films don't have crew data, some only unk genders
      group_by(EventID_VZ) %>% 
      mutate(nf=n_distinct(refFilm)) %>% 
      filter(nf>=minfilms) # final filter here
    
    festvecs=filmfestvecs %>% # filter(EventID_VZ=="1") %>% # test  
      # grouped by event
      summarise(
        rfmaxboot = boot(ratiofmax, function(dx,y){mean(dx[y])}, R=1000) %>% boot.ci(type = "norm") %>% {.$normal[3]},  # max upper 
        rfminboot = boot(ratiofmin, function(dx,y){mean(dx[y])}, R=1000) %>% boot.ci(type = "norm") %>% {.$normal[2]},  # min lower
        
        across(starts_with("ratio"), mean,na.rm=T),
        across(c( "refEvent", "lstTypeEvent", "LibelleEvent_NEW", "LibelleEvent_NEW_VZ", "YearEvent", "libelleFestival_NEW_VZ", "ab", "PFIAF category", "nf"), first)
        
      ) %>% 
      mutate(tool=paste0(
        LibelleEvent_NEW_VZ, 
        ", FMmean=",  ratiof  %>% round(2),
        ", producersmean=", ratiopf %>% round(2),
        ", directorsmean=",  ratiodf %>% round(2),
        ", nfilms=", nf
      )) 
    
    
    # system averages
    sysdivs0 = d %>% 
      left_join(gratios, by = "refFilm") %>% 
      filter(!is.na(ratiof), nother<1) %>%   
      group_by(EventID_VZ) %>% 
      mutate(nf=n_distinct(refFilm)) %>% 
      filter(nf>=minfilms, YearEvent>=2012) %>%  # final filter + fixed year filter
      group_by(EventID_VZ) %>% 
      mutate(ws = 1/n_distinct(refFilm)) %>% 
      ungroup()
    
    sysdivs = rbind(sysdivs0, sysdivs0 %>% mutate(ab="Both")) %>% 
      group_by(ab, YearEvent) %>% 
      summarize(
        rfmaxboot = boot(ratiofmax, function(dx,y){weighted.mean(dx[y], w=ws[y])}, R=1000) %>% boot.ci(type = "norm") %>% {.$normal[3]},  # max upper 
        rfminboot = boot(ratiofmin, function(dx,y){weighted.mean(dx[y], w=ws[y])}, R=1000) %>% boot.ci(type = "norm") %>% {.$normal[2]},  # min lower
        across(starts_with("ratio"), weighted.mean, w=ws, na.rm=T), 
        .groups = "drop"
      ) %>% 
      rename(sysdiv=ratiof, supper=rfmaxboot, slower=rfminboot)
    
    info = paste(length(unique(filmfestvecs$refFilm)), "unique films," ,
                 length(unique(filmfestvecs$EventID_VZ)), "unique fests,", 
                 nrow(filmfestvecs), "film*fest pairs,",
                 "after final filter")
    info %>% print
    
    attr(festvecs, "filmfestvecs")=filmfestvecs
    attr(festvecs, "info") = info
    attr(festvecs, "sysdivs") = sysdivs
    
    return(festvecs)
  }
  
  # integrated fixes and join with metadata labels here, in prep for matching to spaces
  dat = fixdat(dat, metatype, filmkinds, fixes=fixes)
  
  paste(nrow(dat), "fest*films pre meta filter,", paste(unique(dat$EventID_VZ) %>% length, "fests pre meta filter")) %>% print # 
  
  filmfestvecs0 = 
    dat %>%   
    # final filter, as some have missing meta
    group_by(EventID_VZ) %>% filter(n_distinct(refFilm)>=minfilms) %>%  
    # weights
    group_by(refFilm) %>% mutate(weight = 1/n_distinct(txtKind))   # film meta weights
  if(metatype=="genre"){
    filmfestvecs=filmfestvecs0 %>% mutate(tmpweight=1)  # but not for genre
  } else {
    filmfestvecs=filmfestvecs0 %>% mutate(tmpweight=weight)
  }
  filmfestvecs = filmfestvecs %>%  
    group_by(EventID_VZ) %>% 
    mutate(exweight = tmpweight*(1/n_distinct(refFilm))) %>% 
    select(-tmpweight) %>% 
    
    # add festival main type, using original label
    #group_by(EventID_VZ, txtKind) %>% 
    # add_count(name = "nc") %>% # won't work with weighted
    # -festival main now added below
    # mutate(nc=sum(weight)) %>% 
    # group_by(EventID_VZ) %>% 
    # mutate(firstkind = txtKind[nc == max(nc)][1]) %>% 
    # select(-nc) %>% 
    ungroup() %>% 
    
    # film main type
    mutate(ord=suppressWarnings(as.numeric(NoOrder))) %>%    # will throw warning, is ok
    mutate(ord= case_when(is.na(ord)~100,T~ord)) %>% 
    group_by(refFilm) %>% mutate(firstfilmkind = txtKind[ord == min(ord)][1]) %>%  
    #select(-ord) %>% 
    ungroup() %>% 
    
    # add latent vectors
    cbind(as.data.frame(lat)[.$txtKind3,]) 
  
  # rank top discrete labels, assign main label to fest
  tops0 = filmfestvecs %>% group_by(EventID_VZ, txtKind) %>% 
    summarise(s=sum(weight), no=sum(1/ord), .groups = "drop_last" ) %>% 
    mutate(s=(s/sum(s)*100)) %>% 
    arrange(-s, -no, .by_group = T)
  tops=tops0 %>% 
    slice(1:5) %>% 
    filter(s>0) %>% 
    mutate(s=round(s)) %>% # round after exact ranking
    mutate(l = paste0(txtKind,":",s,"%"), toptest=txtKind[s==max(s)][1]) %>%  
    summarise(firstkind=toptest[1], tops=paste(l,collapse=" "), .groups = "drop")
  filmfestvecs = left_join(filmfestvecs, tops, by="EventID_VZ")
  # add event lat vector for geo
  if(metatype=="geo"){
    filmfestvecs = filmfestvecs %>% 
      cbind(as.data.frame(lat)[.$eventcountry3,] %>% 
              rename(eventlat=lat, eventlong=long)
      )
  }
  rownames(filmfestvecs)=NULL 
  
  
  # for some reason this is slower in new dplyr 1.1.0? rewriting as 2 cbinded operations, much faster (but prone to errors if not matching; checked now, works)
  if(metatype=="genre"){ # latent genres are averaged, but countries and languages not.
    filmfestvecs = 
      cbind(
        filmfestvecs %>% group_by(EventID_VZ, refFilm) %>% select(-starts_with("V", ignore.case=F)) %>% slice(1) %>% ungroup() ,
        filmfestvecs %>% group_by(EventID_VZ, refFilm) %>% summarise(across(starts_with("V", ignore.case=F),mean), .groups = "drop") %>% select(-c(EventID_VZ, refFilm))
      ) %>% 
      # filmfestvecs %>% group_by(EventID_VZ, refFilm) %>% 
      # summarise(across(!starts_with("V"), first),
      #           across(starts_with("V"),mean)) %>% 
      ungroup() %>% 
      mutate(weight=1)  # important, set weights to 1 now, since 1 vec per film!
    
  }
  info = paste(length(unique(filmfestvecs$refFilm)), "unique films," ,
               length(unique(filmfestvecs$EventID_VZ)), "unique fests,", 
               nrow(filmfestvecs), "film*fest pairs,",
               "after final filter")
  info %>% print
  #filmfestvecs<<-filmfestvecs  # save tmp object in global for aother calc
  
  exfilms = filmfestvecs %>% 
    filter(nchar(gsub("[^ ]", "", tmp))<=3, nchar(tmp)<25 ) %>% 
    group_by(tmp) %>% mutate(n=n()) %>% arrange(-n, .by_group = T) %>% 
    group_by(firstfilmkind) %>% slice_max(n,with_ties = F) %>% ungroup %>%  
    select(tmp, n, firstfilmkind) %>%  
    rename(firstkind=firstfilmkind) %>% 
    filter(!duplicated(tmp)) %>% arrange(-n) 
  
  
  # scaling factor and latent centre calculations
  maxdist=NULL
  if(metatype %in% c("lang", "genre")){
    # scaling value (not in geo)
    maxdist = dist(lat) %>% max
    
    ## Festival ecosystem grand mean vector, 
    # weighted by meta and festival size to balance things out
    fmv = filmfestvecs %>% 
      summarise(across(starts_with("V"), weighted.mean, w=exweight, na.rm=T)) %>% 
      as.matrix()
    maxexdist = dist2(lat, fmv, "euclidean", "none") %>% max # max to mean vector here!
    
    # fmv1 = tmpfestvecs %>% summarise(across(starts_with("V"), mean)) %>% select(starts_with("V")) %>% colMeans() %>% matrix(nrow=1) # initial simpler mean of means variant
  } else {
    if(metatype=="geo"){
      # no scaling factor required for geo
      # geo grand mean vec:
      fmv=geomean(filmfestvecs %>% select(long, lat),w = filmfestvecs$exweight)
    }
  }
  
  # System internal diversity
  # uses external diversity weights, because need to account for both n tags per films AND n films per fests, as the MAD is calculated across all films.
  sysdivs = filmfestvecs %>% 
    rbind(filmfestvecs %>% mutate(ab="Both")) %>% 
    filter(YearEvent>=2012) %>% # hard filter here, early years have too little data
    group_by(ab, YearEvent) %>% 
    group_split() %>% 
    lapply(function(x){
      if(metatype=="geo"){
        fmean=geomean(x %>% filter(!is.na(lat)) %>% select(long, lat), x$exweight) 
        colnames(fmean)=c("long", "lat")
      } else {
        fmean = summarise(x, across( starts_with("V"), weighted.mean, w=x$exweight, na.rm=T)) %>% as.matrix # year mean
      }
      #  system internal diversity + bootstrapped confidence
      di = disttomean(x, fmean, w=x$exweight, maxdist, geo=metatype=="geo", divmult=2)
      return(tibble(sysdiv=di[1], slower=di[2], supper=di[3], 
                    YearEvent=x$YearEvent[1], ab=x$ab[1], nf=nrow(x)  ))
    }
    ) %>% bind_rows()
  
  
  # do averaged festival vectors with bootstrapped confidence intervals
  festvecs = filmfestvecs %>% 
    group_by(EventID_VZ) %>% 
    group_split() %>% 
    lapply(function(x){
      x=ungroup(x)
      # inter-distance of all films as metric
      #indiv_interdist = dist2(x %>% select(starts_with("V")) %>% as.matrix, method = "euclidean", norm = "none") %>% .[lower.tri(.)] %>% mean
      
      # distance from mean but excluding self; 
      # xx=x %>% select(starts_with("V")) %>% as.matrix
      # d=rep(NA, nrow(xx))
      # for(i in 1:nrow(xx)){ 
      #   d[i]=dist2( matrix(colMeans(xx[-i,]),nrow=1), # festival mean
      #               xx[i,,drop=F],
      #         method = "euclidean", norm = "none")[1,1]
      # }
      # indiv = mean(d, na.rm=T)
      
      
      # internal diversity, distance from mean incl self, essentially the MAD
      # including self in mean might bias, but it's a simpler calc
      # and MAD seems pretty widely used, also includes self
      if(metatype=="geo"){
        fmean=geomean(x %>% filter(!is.na(lat)) %>% select(long, lat), x$weight) # is matrix
        colnames(fmean)=c("long", "lat")
      } else {
        fmean = summarise(x, across( starts_with("V"), weighted.mean, w=x$weight, na.rm=T)) %>% 
          as.matrix # festival mean
        #  genres: equal, averaged above already (constant weight); 
        #  for country and lang need weighted
      }
      
      
      # internal and external diversity + bootstrapped confidence
      # mean, lower, upper
      di=disttomean(x, fmean, w=x$weight, maxdist, geo=metatype=="geo", divmult=2)
      de=disttomean(x, fmv, w=x$weight, maxexdist, geo=metatype=="geo", divmult=1)
      
      res=cbind(
        summarise(.data=x,
                  across(c(
                    EventID_VZ,
                    LibelleEvent_NEW_VZ,
                    libelleFestival_NEW_VZ,
                    idFestival_NEW,
                    YearEvent,
                    txtCountryEvent_VZ, 
                    ab,
                    "PFIAF category",
                    originalKind,
                    # txtKind,
                    # txtKind2,
                    # txtKind3,
                    firstkind,
                    tops
                  ), first),
                  nf=n()
        ),
        #indiv_sdold= x %>% ungroup() %>% select(starts_with("V")) %>% apply(2,sd,na.rm=T) %>% mean(na.rm=T),
        #summarise(x, across(starts_with("V"), sd, na.rm=T, .names="sd_{.col}")) %>% 
        # indiv_interdist,
        indiv=di[1], ilower=di[2], iupper=di[3],       # this is now scaled (above)
        exdiv=de[1], elower=de[2], eupper=de[3], 
        fmean  # averaged festival latent vector
      ) 
      
      
      
      if(metatype=="geo"){
        # additional comparison, geo distance between an event coords and all film coords
        dev = disttomean(x,  # film longs, lats
                         cbind(x$eventlong[1], x$eventlat[1]), # compared to event longlat
                         w=x$weight, maxd=NA, geo=T)
        gx = cbind(
          eventlong=x$eventlong[1], eventlat=x$eventlat[1], 
          eventdist=dev[1], edlower=dev[2], edupper=dev[3] )
        res=cbind(res, gx) %>% mutate(
          tool=paste0(
            LibelleEvent_NEW_VZ, 
            ", int.div=", round(indiv), 
            ", contr.div=", round(exdiv),
            ", eventdist=", round(eventdist),
            ", nfilms=", nf,
            "\ntop: ",tops
          )
        )
      } else {
        res= res %>% mutate(
          tool=paste0(
            LibelleEvent_NEW_VZ, 
            ", int.div=", round(indiv,2), 
            ", contr.div=", round(exdiv,2),
            ", nfilms=", nf,
            "\ntop: ", tops
          ))
      }
      
      return(res)
      
    }) %>% 
    bind_rows()
  
  if(is.null(customlabs)){
    festvecs = festvecs %>% 
      mutate(colorgenre =   # top ones across all fests; for plotting later
               case_when(firstkind %in% head(names(sort(table(.$firstkind),decreasing = T)),6)~firstkind, T~"other" ) %>% as.factor()) %>% 
      mutate(colorgenre = fct_relevel(colorgenre, 
                                      c(names(sort(table(colorgenre),decreasing = T)) %>% .[.!="other"], "other") ))
  } else {
    festvecs = festvecs %>% 
      mutate(colorgenre =   # top ones across all fests; for plotting later
               case_when(firstkind %in% customlabs ~firstkind, T~"other" ) %>% 
               as.factor()) %>% 
      mutate(colorgenre = fct_relevel(colorgenre, customlabs))
  }
  
  
  paste(nrow(festvecs), "fest vectors") %>% print
  # genre: 576 out of 622 left.
  attr(festvecs, "exfilms")=exfilms
  attr(festvecs, "latentcentre")=fmv
  attr(festvecs, "maxdist")=maxdist
  attr(festvecs, "filmfestvecs")=filmfestvecs
  attr(festvecs, "info")=info
  attr(festvecs, "tops")=tops0
  attr(festvecs, "sysdivs")=sysdivs
  
  
  return(festvecs)
}





# change over time calc
dogenrechangevals = function(tmpfestvecs, festvecs, ex){
  maxdist=attr(festvecs, "maxdist")
  genrechangevals = tmpfestvecs %>% 
    filter(libelleFestival_NEW_VZ %in% ex) %>% 
    group_by(libelleFestival_NEW_VZ) %>% 
    #filter(n_distinct(EventID_VZ)>=minevents) %>%  # replace with target filter
    # ungroup %>% filter(libelleFestival_NEW_VZ==libelleFestival_NEW_VZ[1])#debug
    group_split() %>% 
    lapply(function(x){
      x=ungroup(x)
      emeans = x %>% group_by(EventID_VZ) %>% 
        summarise(across(starts_with("V"), mean),
                  YearEvent=mean(YearEvent)) %>% arrange(YearEvent)
      # fmean = emeans %>% ungroup() %>% select(starts_with("V")) %>% colMeans() %>% matrix(nrow=1) # comparison to mean is not really informative, correlates strongly with diversity, which makes sense.
      tmp = tibble(
        # change_interdist= dist2(xt, xr, method = "euclidean", norm = "none") %>% 
        #   .[lower.tri(.,diag = F)] %>% mean,
        EventID_VZ=emeans$EventID_VZ,
        dif = NA,
        diflower=NA,
        difupper=NA
      )
      for(i in 2:nrow(emeans)){
        b=boot(x %>% filter(EventID_VZ %in% emeans$EventID_VZ[c(i, i-1)]),
               meanvecs, R=200,
               ems=emeans, loopindex=i, 
               maxdist=maxdist)
        bci = boot.ci(b, type = "norm") %>%
          {.$normal[2:3]}
        tmp$dif[i] = b$t0
        tmp$diflower[i] = bci[1]
        tmp$difupper[i] = bci[2]
        
      }
      return(tmp)
    }) %>% do.call(rbind, .) %>% 
    left_join(., festvecs,  by="EventID_VZ")
  return(genrechangevals)
}



dogenrecommonspace = function(festvecs, filmfestvecs){ # hacky helper function
  allgenvecs = bind_rows(
    # genres
    genlat %>% as.data.frame() %>% rownames_to_column("firstkind") %>%  
      mutate(YearEvent=0, type="genre",event=firstkind %>% toupper()) %>% 
      filter(firstkind %in% names(tail(sort(table(filmgen$txtKind)),22))
      ) %>% 
      left_join(count(festvecs, firstkind) %>% rename( sizegen=n), by="firstkind") %>%
      mutate(sizegen=case_when(is.na(sizegen)~2L, sizegen<5L~3L, T~sizegen))
    ,
    # films
    filmfestvecs %>% filter(tmp %in% attr(festvecs, "exfilms")$tmp) %>% 
      filter(!duplicated(tmp)) %>% select(starts_with("V"), firstfilmkind, YearEvent, tmp) %>% 
      rename(firstkind=firstfilmkind) %>%  
      mutate(type="film",
             event= tolower(tmp) %>%  #gsub("[_ ]","\n",event) %>% gsub("^|$", '"',.), 
               gsub("_([0-9]{4,4})", '"\n(\\1)',.) %>% gsub("^", '"',.),
             sizegen=1,
             tool=tmp
      ) ,
    # fests
    festvecs %>% select(starts_with("V"), YearEvent,LibelleEvent_NEW_VZ,firstkind, ab, tool) %>% rename(event=LibelleEvent_NEW_VZ) %>% mutate(type="festival")
  ) %>% 
    #, tool=paste0(event, "\nmain: ", firstkind)) )  %>% 
    mutate(colorgenre=case_when(firstkind %in% levels(festvecs$colorgenre)~firstkind, T~"other")
    ) %>%   
    mutate(colorgenre = fct_relevel(as.factor(colorgenre), levels(festvecs$colorgenre))) 
  return(allgenvecs)
}


doUMAP = function(festvecs){
  u = umap(festvecs %>% select(starts_with("V")), 
           method = "umap-learn", min_dist=0.6, n_neighbors=20)$layout %>% 
    as.data.frame() %>% 
    cbind(festvecs %>% select(-starts_with("V"))) 
  # mutate(colorgenre = 
  #          case_when(firstkind %in% head(names(sort(table(festvecs$firstkind),decreasing = T)),6)~firstkind, T~"other" ) %>% as.factor()) %>% 
  # mutate(colorgenre = fct_relevel(colorgenre, 
  #                                 c(names(sort(table(colorgenre),decreasing = T)) %>% .[.!="other"], "other") )) %>%
  if(!("sizegen" %in% colnames(u))){
    u=u %>%  group_by(colorgenre) %>% 
      mutate(sizegen = n()) %>% ungroup()
  }
  
  return(u)
}



do_networkobjects=function(fesfilm3){
  # this creates the various objects needed for network plots and analyses
  # currently saves some intermediate objects directly into global environment; should refactor at some point
  ffmat = fesfilm3 %>% select(EventID_VZ, filmyear) %>% mutate(tmp=1) %>% 
    pivot_wider(values_from=tmp, values_fill = 0, id_cols = "filmyear", names_from = "EventID_VZ") %>% 
    column_to_rownames("filmyear") %>% as.matrix() %>% 
    as.dfm() %>% fcm()
  all((ffmat %>% as.matrix() %>% diag)==0)
  dim(ffmat) # 643
  edgelist = reshape2::melt(ffmat %>% as.matrix()) %>% 
    .[.[,3]>0,] # remove 0-weight links (ie no link, byproduct of the fcm)
  colnames(edgelist)=c("from", "to", "weight")
  edgelist = edgelist %>%   mutate(to=as.character(to), from=as.character(from)) 
  dim(edgelist)  #  14420 
  rownames(edgelist)=NULL
  
  e2 = edgelist %>%  
    left_join(fesfilm3 %>% filter(!duplicated(EventID_VZ)) %>% 
                mutate(tolibelle=paste0(YearEvent," ", libelleFestival_NEW_VZ)) %>% 
                rename(to=EventID_VZ, toyear=YearEvent) %>%
                select(to, toyear, tolibelle), by="to" ) %>% 
    left_join(fesfilm3 %>% filter(!duplicated(EventID_VZ)) %>% 
                mutate(fromlibelle=paste0(YearEvent," ", libelleFestival_NEW_VZ)) %>% 
                rename(from=EventID_VZ, fromyear=YearEvent) %>% 
                select(from, fromyear, fromlibelle), by="from" ) %>% 
    mutate(dist=abs(toyear-fromyear), coloryear = pmax(toyear, fromyear)) 
  
  heatmaplist = ffmat %>% as.matrix() %>% 
    {.[lower.tri(., diag = F)]=t(.)[lower.tri(., diag = F)];.} %>% 
    {diag(.)=NA;.} %>% 
    reshape2::melt() %>% {colnames(.)=c("from", "to", "weight");.} %>% 
    mutate(to=as.character(to), from=as.character(from), weight=na_if(weight, 0)) %>% 
    left_join(fesfilm3 %>% count(EventID_VZ,name="ton") %>% rename(to=EventID_VZ), by="to") %>% 
    left_join(fesfilm3 %>% count(EventID_VZ,name="fon") %>% rename(from=EventID_VZ), by="from") %>% 
    mutate(overlap=weight/fon*100) %>%  # how much of programme is covered by another
    left_join(fesfilm3 %>% filter(!duplicated(EventID_VZ)) %>% 
                mutate(tolibelle=paste0(YearEvent," ", libelleFestival_NEW_VZ)) %>% 
                rename(to=EventID_VZ, toyear=YearEvent) %>%
                select(to, toyear, tolibelle), by="to" ) %>% 
    left_join(fesfilm3 %>% filter(!duplicated(EventID_VZ)) %>% 
                mutate(fromlibelle=paste0(YearEvent," ", libelleFestival_NEW_VZ)) %>% 
                rename(from=EventID_VZ, fromyear=YearEvent) %>% 
                select(from, fromyear, fromlibelle), by="from" ) %>% 
    mutate(dist=abs(toyear-fromyear), reldist=(toyear-fromyear), coloryear = pmax(toyear, fromyear)) 
  
  nodelist = fesfilm3 %>% 
    mutate(EventID_VZ==fct_inorder(EventID_VZ)) %>%  
    group_by(EventID_VZ) %>% 
    summarise(nf=n(), YearEvent=mean(YearEvent, na.rm=T), 
              LibelleEvent_NEW_VZ=LibelleEvent_NEW_VZ[1],
              libelleFestival_NEW_VZ = libelleFestival_NEW_VZ[1],
              ab=ab[1])  %>% 
    # left_join(u, by="EventID_VZ") %>% 
    mutate(x=YearEvent+runif(nrow(.), -0.2, 0.2)) 
  d=degree(graph_from_data_frame(e2, directed = F, vertices = nodelist)) # all(nodelist$EventID_VZ==names(d))
  # b=betweenness(graph_from_data_frame(e2, directed = F, vertices = nodelist))
  s=strength(graph_from_data_frame(e2, directed = F, vertices = nodelist) )
  nodelist=nodelist %>% 
    mutate(degree=d[.$EventID_VZ]) %>% 
    #mutate(ndegree=degree/nf) %>% 
    mutate(lndegree=degree/log10(nf)) %>% 
    mutate(#betweenness=b, nbetweenness = b[.$EventID_VZ]/nf, 
      strength=s[.$EventID_VZ], nstrength=s[.$EventID_VZ]/log10(nf))  %>% 
    mutate(tool=paste0(LibelleEvent_NEW_VZ,"\n nfilms=", nf, " degree=",degree, " lognormdegree=",round(lndegree,2), " strength=",round(strength,2) )) #, " lognormstrength=",round(nstrength,2) ))
  
  attr(nodelist, "info") =
    paste(fesfilm3$EventID_VZ %>% unique %>% length, "festivals,", 
          fesfilm3 %>% group_by(libelleFestival_NEW_VZ) %>% filter(n()>1) %>% pull(libelleFestival_NEW_VZ) %>% unique %>% length, "series with 2+ events",
          fesfilm3$refFilm %>% unique %>% length , "unique films,", nrow(fesfilm3), "festival-film pairs")  # 
  
  # saves into env
  ffmat<<-ffmat
  heatmaplist<<-heatmaplist
  e2<<- e2
  return(nodelist)
}




#### Plotting functions ####

set_plottingvars = function(){
  convertplots <<- F
  
  abcols <<- divergingx_hcl(6, palette="Zissou 1")[c(1,4)] %>% c("gray2") %>% desaturate(0.6) #%>% show_col
  abshapes <<- c(16,4)
  
  sixgcols <<- c(divergingx_hcl(6, palette="Zissou 1") %>% {.[6]=darken(.[6],0.3);.},  "gray80") %>%
    {.[1]=lighten(.[1],0.2);.} %>% {.[3] = "#03408f";.} %>% 
    {.[c(2)]=darken(.[c(2)],0.1);.} %>% 
    {.[c(1,4,2,5,3,6,7)]} # %>%  show_col
  sixcols <<- c(divergingx_hcl(6, palette="Zissou 1"), "gray10",  "gray80") %>% {.[c(1)]=lighten(.[c(1)],0.5);.} %>% {.[c(2)]=darken(.[c(1)],0.5);.} %>% {.[c(4)]=lighten(.[c(4)],0.2);.} %>% {.[7]=darken(.[1],0.75);.} # %>% show_col
  langcols <<- c(divergingx_hcl(6, palette="Zissou 1") %>% 
                {.[c(1)]=darken(.[c(1)],0);.} %>% 
                {.[7]=darken(.[1],0.6);.} ,  "gray80")[c(1,3,5,6,7,8)] #%>% show_col
  rolecols <<- divergingx_hcl(31, palette="Zissou 1")  %>% {c(.[1:2], "gray75", .[30:31])}
  matcols <<- (divergingx_hcl(8, palette="Zissou 1")[c(3,5,7)] %>% {.[1]=lighten(.[1],0.4);.} %>%  {.[3]=darken(.[3],0.95);.}) # %>% show_col
  netcols <<- (divergingx_hcl(5, palette="Zissou 1")[c(1,1,2,3)] %>% {.[1]=darken(.[1],0.5);.} %>%  {.[4]=lighten(.[4],0.45);.} ) %>% {c("#01033a", .)} #  %>% show_col
  
  
  exgenfests <<- c("Festival de Cannes", "Sundance", "Hot Docs Festival", "Sitges FF")
  exgeofests<<-c("Festival de Cannes", "Sundance", "Tokyo IFF", "Zurich FF")
  exlangfests<<-c("Festival de Cannes", "Sundance", "Busan IFF", "BAFICI")
  exrolefests <<-c("Festival de Cannes", "Sundance","SXSW",  "Berlinale - Berlin IFF")
}


tags=function(g, x="", m=margin(0,-9,-10,0)){
  g+labs(tag=x)+
    theme(plot.tag=element_text(face = "bold",hjust=0,vjust=1, margin=m))
}


plotbees = function(vecs, tag, tlab="", col=NULL, cex=3, vw=1.7, comp=F, colvals=abcols){
  g=ggplot(vecs 
           , 
           aes(x=as.factor(ab) %>% as.numeric , group=as.factor(ab) %>% as.numeric, color=ab,y=target))+
    #geom_boxplot(fill=NA, color="gray30", width=0.2)+
    geom_violin(fill=NA, color="gray50", size=0.1,scale="count", width=vw)
  
  if(is.null(col)){
    if(!comp) {g=g+geom_beeswarm(aes(color=ab), cex = cex,priority = "random",size=0.5)}
    if(comp) {g=g+geom_beeswarm(aes(color=comp), cex = cex,priority = "random",size=0.5)}
  } else {
    g=g+geom_beeswarm(color=col, cex = cex,priority = "random",size=0.5)
  }
  
  g=g+stat_summary(na.rm=T, shape="-", color="black", fun = median, size=2)+
    #scale_shape_manual(values=abshapes)+
    scale_color_manual(values=colvals)+
    scale_size_manual(values=c(0.6,0.3))+
    coord_cartesian(xlim=c(0.6,2.6), expand=T)+
    labs(title = expr(paste(bold(!!tag), !!tlab)))+
    scale_x_continuous(breaks=1:2, labels=c("A", "B"))+
    scale_y_continuous(expand=expansion(0.01,0))+
    theme_bw()+
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.title = element_blank(),
          plot.margin = margin(1,1,1,1),
          plot.title = element_text(size=8, margin=margin(0,0,0,0)),
          legend.position = "none",
          axis.text.x = element_text(size=6, margin=margin(-1,0,0,0), vjust=1),
          axis.text.y = element_text(size=6, margin=margin(0,-1,0,0), hjust=1),
          axis.ticks = element_blank()
    )
}





dopies = function(vecs, tops=NULL){
  mains=levels(vecs$colorgenre) 
  if( !attr(vecs, "tops") %>% is.null){
    tops = attr(vecs, "tops") 
  }
  pies = tops %>% 
    left_join(vecs %>% filter(!duplicated(EventID_VZ)), by="EventID_VZ") %>% 
    select(txtKind, s, EventID_VZ, libelleFestival_NEW_VZ, YearEvent, LibelleEvent_NEW_VZ)   %>% 
    mutate(txtKind=case_when(txtKind%in%mains~txtKind, T~"other")) %>% 
    group_by(EventID_VZ, txtKind) %>% 
    summarise(s=sum(s), 
              across(c(LibelleEvent_NEW_VZ, libelleFestival_NEW_VZ, YearEvent), first), )    %>% 
    ungroup() %>% 
    mutate(colorgenre2=fct_relevel(txtKind, mains)) %>% select(-txtKind) %>% 
    arrange(colorgenre2)
  
  pies2 = pivot_wider(pies , names_from = colorgenre2, values_from=s, values_fill=0 ) %>%
    mutate(event=LibelleEvent_NEW_VZ) %>% 
    group_by(libelleFestival_NEW_VZ) %>% 
    arrange(YearEvent, .by_group = T) %>% 
    ungroup
  attr(pies2, "mains")=mains
  return(pies2)
}





div2plot = function(vecs, ex=NULL, convert=T,  colorlab="Main", xlab="Contributing diversity", ylab="Internal diversity", genrecols=sixcols, pies){
  # if(isgenre){
  #     vecs = vecs %>% arrange(YearEvent) %>% 
  #   mutate(colorgenre = 
  #          case_when(firstkind %in% head(names(sort(table(.$firstkind),decreasing = T)),6)~firstkind, T~"other" ) %>% as.factor()) %>% 
  #   mutate(colorgenre = fct_relevel(colorgenre, 
  #                                 c(names(sort(table(colorgenre),decreasing = T)) %>% .[.!="other"], "other") ))
  #     colorlab="Main\ngenre"
  # } else {
  #   vecs = vecs %>% arrange(YearEvent) %>% 
  #   mutate(colorgenre = 
  #          case_when(firstcountry %in% head(names(sort(table(.$firstcountry),decreasing = T)),6)~firstcountry, T~"other" ) %>% as.factor()) %>% 
  #   mutate(colorgenre = fct_relevel(colorgenre, 
  #                                 c(names(sort(table(colorgenre),decreasing = T)) %>% .[.!="other"], "other") ))
  
  # colorlab="Main\nprod.\ncountry"
  #}
  
  
  g={ggplot(vecs , 
            aes(exdiv, indiv, shape=ab, color=colorgenre, size=ab, text=tool))+
      geom_vline(xintercept = median(vecs$exdiv, na.rm=T) # mean(c(min(vecs$exdiv,na.rm=T), max(vecs$exdiv,na.rm=T)))
                 , color="gray25")+
      geom_hline(yintercept = median(vecs$indiv, na.rm=T)
                 #   mean(c(min(vecs$indiv,na.rm=T), max(vecs$indiv,na.rm=T)))
                 , color="gray25")+
      geom_point(#data=vecs %>% filter(!(libelleFestival_NEW_VZ%in%ex)), 
        alpha=0.7) +
      #geom_point(data=vecs %>% filter((libelleFestival_NEW_VZ%in%ex)), alpha=1)+
      #geom_line(aes(YearEvent, indiv), data=festvecs %>% group_by(YearEvent) %>% summarise(indiv=mean(indiv,na.rm=T)), inherit.aes = F)+
      # geom_path(aes(group=libelleFestival_NEW_VZ), 
      #           vecs %>% filter(libelleFestival_NEW_VZ%in%ex) ,
      #           color="black", 
      #           size=0.2)+
      #scale_color_manual(values=abcols)+
      scale_color_manual(values=genrecols, name=colorlab)+
      scale_shape_manual(values=abshapes, name="Festival")+
      scale_size_manual(values=c(2, 1.2), guide = "none")+
      labs(x=xlab, y=ylab)+
      theme_bw()+
      theme(#legend.position = 'left', 
        #panel.border = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.background = element_rect(color="white", fill="white")
      )} 
  if(convert){
    ggplotly(g,tooltip="text")
  } else {
    g
  }
}

divyearplot = function(vecs, ex, shape=16, yl="", splitab=T, convert=T, networkstats=T, ylims=NULL, genrecols=sixgcols, pies=NULL, nrowmod=6, contrdiv=T){
  if(splitab){
    vecs = vecs %>% mutate(YearEvent=case_when(ab=="A"~YearEvent-0.15,
                                               ab=="B"~YearEvent+0.15 )) %>% 
      arrange(libelleFestival_NEW_VZ) %>% 
      mutate(series = substr(libelleFestival_NEW_VZ, 1,15))
    if(length(unique(vecs$libelleFestival_NEW_VZ)) != 
       length(unique(vecs$series))){ stop("shortening mismatch")}
    vex=vecs
    g=ggplot(vecs , 
             aes(YearEvent, target, color=ab, shape=ab,  text=tool))
  } else {
    vex = vecs %>% 
      filter(libelleFestival_NEW_VZ %in% ex, !is.na(target)) %>% 
      mutate(libelleFestival_NEW_VZ=fct_relevel(libelleFestival_NEW_VZ, ex))
    g=ggplot(vex , 
             aes(YearEvent, target, color=ab, shape=ab,  text=tool))
  }
  
  
  if(splitab){
    
    col2 = vecs  %>% 
      filter(!duplicated(libelleFestival_NEW_VZ)) %>% 
      mutate(col=case_when(ab=="A"~ abcols[1], T~abcols[2])) %>% 
      pull(col)
    
    if(networkstats){
      fromy=2014.5; toy=2019.5
    } else{
      fromy=2013.5; toy=2021.5
    }
    g= g+
      stat_summary(aes(y = target, group=1), 
                   data=vecs %>% filter(YearEvent>fromy, YearEvent<toy, ab=="A"),
                   fun=mean, na.rm=T,  geom="line", size=2, color=abcols[1])+
      stat_summary(aes(y = target, group=1), 
                   data=vecs %>% filter(YearEvent>fromy, YearEvent<toy, ab=="B",),
                   fun=mean, na.rm=T,  geom="line", size=1.5, color=abcols[2])+
      scale_x_continuous(breaks=2009:2021)
    
    if(convertplots){
      g=g+geom_beeswarm(aes(color=series), # creates legend in ggplotly
                        alpha=1, cex=0.8,size=1.5, priority = "random") +
        scale_color_manual(values=col2)
    } else {
      g=g+geom_beeswarm(aes(color=ab),
                        alpha=1, cex=0.8,size=1.5, priority = "random") +
        scale_color_manual(values=abcols)
    }
    
  } else {
    g=g+scale_color_manual(values=abcols)+
      theme(legend.position = 'none')
  }
  
  g=g+  
    scale_shape_manual(values=abshapes)+
    #scale_size( range=c(0.5*2.2, 1.5*3.2), guide = "none")+
    theme_bw()+
    theme(
      #legend.position = 'none', 
      #panel.border = element_blank(),
      
      panel.grid.major.x  = element_blank(),
      panel.grid.minor.x  = element_line(color="gray80"),
      plot.background = element_rect(color="white", fill="white"),
      axis.title.x = element_blank(),
      axis.ticks.length = unit(0, "inch")
    )+
    labs(y=yl)
  # if(nrow( vecs %>% filter(libelleFestival_NEW_VZ%in%ex))>0){
  #   g=g+geom_line(aes(group=libelleFestival_NEW_VZ), 
  #             vecs %>% filter(libelleFestival_NEW_VZ%in%ex) ,
  #             color="gray40", 
  #             size=0.2)
  # }
  
  if(!splitab & !is.null(ex)){
    if(networkstats){
      fromy=2014.5; toy=2019.5
    } else{
      fromy=2008; toy=2022
      g=g+geom_errorbar(aes(ymin=lower, ymax=upper), color="gray30", width=0.1, size=0.45)
    }
    g=
      g+facet_wrap(~libelleFestival_NEW_VZ, nrow = ceiling(length(ex)/nrowmod))+
      # geom_point( aes(YearEvent, target),data=vecs %>% select(-libelleFestival_NEW_VZ), inherit.aes = F,
      #             color="gray90", size=0.9)+
      geom_line(data=vex %>% filter(YearEvent>fromy, YearEvent<toy) %>% 
                  group_by(YearEvent, libelleFestival_NEW_VZ) %>% 
                  summarise(target=mean(target, na.rm=T), tool=NA, ab=ab[1]), 
                color="gray70", size=0.3
      )+
      geom_point(data=vex, alpha=1, size=1.5, shape=16)+
      theme(
        panel.grid.major.x  = element_line(color="gray95"),
        panel.grid.minor.x  = element_blank(),
        panel.grid.minor.y  = element_blank(),
        legend.position = 'none',
        strip.text = element_text(hjust=0, margin=margin(0.5,0,0.5,1))
      )
    
    if(!is.null(pies)){
      p = dopies(vex, pies)
      m = attr(p, "mains")
      vp = vex %>%  select(EventID_VZ, target,YearEvent,libelleFestival_NEW_VZ) %>% 
        left_join(p %>% select(EventID_VZ, {{m}}),by="EventID_VZ")
      g = g+ geom_pie_glyph(aes(YearEvent, target), data = vp, slices=m, radius=0.2, inherit.aes = F)+
        scale_fill_manual(values=genrecols, breaks=m)
      
    }
    if(contrdiv){
      g = g+geom_errorbar(aes(ymin=elower, ymax=eupper), color="gray60", width=0, size=0.45,
                          position=position_nudge(0.45))+
        geom_point(aes(y=exdiv), color="gray60", size=0.5, shape=15, 
                   position=position_nudge(0.45))+
        theme(axis.title.y.right = element_text(color="gray50"), 
              axis.ticks.length.y.right = unit(0,"in"),
              axis.text.y.right = element_blank())
    }
  }
  
  minx = min(vex$YearEvent, na.rm=T)
  if(!is.null(ylims)){
    g=g+
      coord_cartesian(ylim = ylims, expand=F)+
      scale_x_continuous(breaks=seq(minx+1,2020,3), limits=c(minx-0.4,2021+0.4 ) )
  } else {
    g=g+ scale_x_continuous(breaks=seq(minx+1,2020,3), expand = expansion(0.1,0) )+
      scale_y_continuous(expand=c(0.05,0)
                         #limits=c( min(c(vecs$target, vecs$lower), na.rm=T), max(c(vecs$target, vecs$upper), na.rm=T )*1.001 )
      )
    if(contrdiv){
      g=g+scale_y_continuous(expand=c(0.05,0),sec.axis = sec_axis(~., name="Contributing diversity"))
    }
  }
  
  
  if(convert){
    
    return( ggplotly(g,tooltip="text") %>% plotly::config(doubleClickDelay=500) 
    )
  } else {
    return(g+theme(legend.position = "none"))
  }
  
}


sysplots = function(sysdiv, yl){
  
  g=  ggplot(sysdiv,aes(YearEvent, sysdiv, color=ab))+
    #geom_pointline(position=position_dodge(width = 0.25), distance=unit(10, "pt") )+
    geom_line(data=sysdiv %>% filter(ab=="Both"), position=position_nudge(x = 0.25/3))+
    geom_point(aes(shape=ab, size=ab), position=position_dodge(width = 0.25))+
    geom_errorbar(aes(ymin=slower, ymax=supper),width=0.1, linewidth=0.45,position=position_dodge(width = 0.25) )+
    scale_x_continuous(breaks=2012:2021)+
    scale_color_manual(values=abcols)+
    scale_shape_manual(values=abshapes %>% c(1))+
    scale_size_manual(values=c(0.6, 0.3, 0.7))+
    theme_bw()+
    theme(legend.position = 'none')+
    theme(
      panel.grid.minor.x  = element_blank(),
      panel.grid.minor.y  = element_blank(),
      plot.background = element_rect(color="white", fill="white"),
      axis.title.x = element_blank()
    )+
    labs(y=yl)
  return(g)
}


gg=function(g, do=T){
  if(do){
    ggplotly(g, tooltip="text") %>%  plotly::config(doubleClickDelay=500)
  } else{
    g
  }
}

pspacer=ggplot()+theme_void()+theme(plot.margin = margin(0,0,0,0))


