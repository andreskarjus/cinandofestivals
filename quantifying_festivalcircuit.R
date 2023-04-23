#### Introduction ##########
#
# This is an R script to replicate results in 
# "Quantifying the global film festival circuit: Networks, diversity, and public value creation"
# by Vejune Zemaityte, Andres Karjus, Ulrike Rohn, Maximilian Schich, Indrek Ibrus
# Codebase (c) Andres Karjus 2023, andreskarjus.github.io
# This software is distributed as is. No warranty expressed or implied, 
# nor guarantee for accuracy or applicability.
# This code may be freely used for non-commercial purposes, with attribution (to the above paper).
# This code has been developed in R 4.2.1, using tidyverse 2.0.0 (incl dplyr 1.1.0).

# Follow the instructions in this file and run everything line by line to reproduce
# the data and results as presented in the paper. If using RStudio, enable the 
# Outline (button top right of this script pane) for a handy overview.



#### 1. Make files accessible ####

# Do one of the following:
# - set working directory, i.e. where the downloaded data files and script file are found:
setwd()
# - or specify full paths below for the files.
# - or just place all the files in the current working directory reported by this function:
getwd()


#### 2. Load packages, scripts, data ####

## Load scripts (also loads packages)
source("quantifying_festivalcircuit_scripts.R")
#
# Effort has been made to automate package installation, including that of
# reticulate and minoconda (necessary for the Python part of this codebase,
# which is limited to the language part).
# If this fails, go to this functions file and attempt to debug. 


## Load data
# First unzip the  cinandofestivals_datasample.zip  data package into your folder of choice.
#
# This is a cleaned sample of the Cinando database used in the paper.
# All Cinando internal reference and technical ID values have been anonymized.
# All personal and identifying data has been removed/anonymized.
# Load all of these csv files:
fesfilm2 = read_csv("festivals_films.csv")  # Main table, festival programmes
filmrole = read_csv("filmrole.csv")         # Film crews (anonymized) metadata
filmlang = read_csv("filmlang.csv")         # Film language metadata
filmcountry = read_csv("filmcountry.csv")   # Film production country metadata
filmgen = read_csv("filmgen.csv")           # Film thematic/content/genre metadata

# If you get a "file not found" error, redo step 1
# If you get a "function not found", make sure package installation and loading worked.

#### ........................... ####




#### 3. Run analyses ####

# Run all of this to create the data objects.


##### Gender diversity ####

festrolevecs = dofestvecs(fesfilm2, metatype="role", sortvar="nroles", 
                          lat=NULL, filmkinds=filmrole, minfilms=15) 



##### Film and festival thematic vecs ####
# thematic/content type, subsumes genre and other aspects

genlat = dolatentspaces(metatype="genre", filmgen=filmgen)
festgenrevecs = dofestvecs(fesfilm2, metatype="genre", sortvar="ngens", 
                           lat=genlat, filmkinds=filmgen, minfilms=15) 
# 25833 unique films, 576 unique fests, 35118 film*fest pairs, after final filter.

# umap common space for genre example
ugenre = dogenrecommonspace(festgenrevecs, attr(festgenrevecs, "filmfestvecs")) %>% doUMAP()

genrechangevals = dogenrechangevals(attr(festgenrevecs, "filmfestvecs"), festgenrevecs, 
                                    ex=c("Berlinale - Berlin IFF", "San Sebastian FF", "Sundance",
                                         "Venice - Biennale", "Busan IFF", "Toronto - TIFF", 
                                         "Tallinn Black Nights IFF"))



##### Latent countries ####

geolat = dolatentspaces(metatype="geo", filmgen=NULL)
festcountryvecs = dofestvecs(fesfilm2, metatype="geo", sortvar="ncountries", 
                             lat=geolat, filmkinds=filmcountry, minfilms=15, 
                             fixes = countryfixes, 
                             customlabs = c("usa", "united kingdom",   
                                            "france",  
                                            "argentina",  "spain", 
                                            "germany", 
                                            "japan",
                                            "other" ) ) 
# get map for plots
mp = joinCountryData2Map(countryExData, joinCode = "ISO3", nameJoinColumn = "ISO3V10", mapResolution = "low") %>% fortify()



##### Latent language space ####

langlat = dolatentspaces(metatype="lang", filmgen=NULL)
festlangvecs = dofestvecs(fesfilm2, metatype="lang", sortvar="nlangs", 
                          lat=langlat, filmkinds=filmlang, minfilms=15, 
                          fixes = langfixes,
                          customlabs = 
                            c("english",  "french",  "spanish", "german", "japanese", "other"))
ulang = doUMAP(festlangvecs)


##### Networks ####


fesfilm3 = fesfilm2 %>% 
  filter(!is.na(TitleVA)) %>% 
  mutate(filmyear=paste0(TitleVA, "_", YearOfProduction)) %>%
  group_by(EventID_VZ) %>% filter(!duplicated(filmyear)) %>% 
  group_by(EventID_VZ) %>% filter(n()>=15) %>% ungroup()

nodelist = do_networkobjects(fesfilm3) # creates various necessary network objects

#### ........................... ####





#### 4. Recreate plots and analyses ####

# Running everything here will reproduce the plots in the paper and save them in the working directory:
getwd()

# Start here
set_plottingvars() # creates color and shape presets

### Networks
# Heatmap
hm=heatmaplist %>% 
  #group_by(fromyear) %>%  arrange(-overlap, .by_group = T) %>%  mutate(y = as.factor((100*fromyear)+(1:n())) ) %>% 
  mutate(fillweight=log10(weight)) %>% 
  mutate(overlap2 = case_when(overlap==0 ~ NA, T~overlap)) %>% 
  mutate(tool=paste0(fromlibelle,"\n",  
                     tolibelle, 
                     "\nn shared=",weight,
                     "\n% overlap=",overlap, "\ntime dist=",dist ))

hmy = hm %>% # year labels
  arrange(tolibelle) %>% 
  filter(toyear>2011, fromyear>2011) %>% 
  group_by(toyear ) %>% 
  filter(1:n()==round(n()/2)) %>% 
  mutate(fromlibelle=tolibelle)

fig_heatmap =  ggplot(hm %>% mutate(w=rescale(overlap2,to=c(1,1.3))) , 
                      aes(tolibelle, fromlibelle, 
                          fill=overlap2, text=tool))+
  coord_cartesian(expand=F)+
  geom_tile(aes(height=w, width=w), size=0)+
  geom_text(aes(label=toyear, x=0), data=hmy, hjust=0,vjust=1, angle=90,size=2.5)+ 
  #color=viridis_pal(end = 0.9)(13))+
  geom_text(aes(label=toyear,x=fromlibelle, y=0), data=hmy,vjust=0, size=2.5)+
  
  #scale_fill_viridis_c(end=0.9, na.value = "white")+
  #scale_fill_gradientn(colors = cls, na.value = "white", trans="log10", name="log\nstrength")+
  # scale_fill_viridis_c(direction=-1,na.value = "white", option="magma", trans="log10", name="log\nstrength")+
  scale_fill_gradientn(colors = matcols, na.value = "transparent",  name="Overlap %")+
  annotate("text", -Inf, Inf, hjust=-2, vjust=1, label="C", fontface='bold')+
  theme_bw()+
  theme(legend.position = c(0.94,0.22), 
        panel.border = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_rect(color="white", fill="white"),
        axis.title = element_blank(),
        axis.text = element_blank(), # element_text(size=3),
        axis.ticks = element_blank()
  )



# big force directed network
gr = as_tbl_graph(graph_from_data_frame(e2, directed = F, vertices = nodelist %>% filter(EventID_VZ %in% c(e2$from,e2$to)) %>% mutate(nsize=case_when(ab=="A"~lndegree*1.3, T~lndegree)) )) # adjust for glyph difference
fig_network = ggraph(gr  ) +
  geom_node_point(aes(color=YearEvent, size=nsize, shape=ab)) +
  #geom_node_text(aes(label=name)) +
  geom_edge_arc(aes(color=coloryear, alpha=dist),  width=0.1, strength=0.04) +
  #scale_color_viridis(end=0.9)+
  #scale_edge_color_viridis(end=0.9)+
  scale_color_gradientn(colors=netcols, name="Year", 
                        breaks=seq(2009,2021,4), 
                        labels=seq(2009,2021,4) %>% substr(3,4) %>% paste0("'",.)
  )+
  scale_edge_color_gradientn(colors=netcols, guide="none")+
  scale_edge_alpha(range=c(0.15, 0.5), guide="none")+
  scale_shape_manual(values=abshapes, guide="none")+
  scale_size(range=c(0.2,1.2),  guide="none")+
  scale_x_continuous(expand=expansion(0.02,0))+
  scale_y_continuous(expand=expansion(0.01,0))+
  #coord_flip()+
  #annotate("text", -Inf, Inf, hjust=0, vjust=1, label="A", fontface='bold')+
  theme_void()+
  labs(tag="A")+
  theme(plot.tag=element_text(face = "bold",hjust=0,vjust=1, margin=margin(0,-9,-10,0)))+
  theme(legend.position = c(0.94,0.22), 
        plot.margin = margin(0,0,0,0),
        plot.background = element_rect(color="white", fill="white")
  )

#ggsave("fig_heatmap.pdf", fig_network+fig_heatmap,width=10,height=5, scale=1)


fig_netbee = ggplot(nodelist %>% filter(YearEvent%in%(2015:2019)), 
                    aes(x=as.factor(ab) %>% as.numeric , group=as.factor(ab) %>% as.numeric, color=ab,shape=ab,size=ab, y=lndegree))+
  #geom_boxplot(fill=NA, color="gray30", width=0.2)+
  geom_violin(fill=NA, color="gray50", size=0.1,scale="count", width=1.7, bw=0.9)+
  geom_beeswarm(cex = 3,priority = "random", color="gray50")+
  stat_summary(na.rm=T, shape="-", color="black", fun = median, size=2)+
  scale_shape_manual(values=abshapes)+
  scale_color_manual(values=abcols)+
  scale_size_manual(values=c(0.6,0.3))+
  coord_cartesian(xlim=c(0.6,2.4), expand=T)+
  scale_x_continuous(breaks=1:2, labels=c("A", "B"))+
  labs(title = expression(paste(bold("B"), " Norm. degree")))+
  scale_y_continuous(expand=expansion(0.01,0))+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(1,0,1,1),
        plot.title = element_text(size=8, margin=margin(0,0,0,0)),
        legend.position = "none",
        axis.text.x = element_text(size=6, margin=margin(-1,0,0,0), vjust=1),
        axis.text.y = element_text(size=6, margin=margin(0,-1,0,0), hjust=1),
        axis.ticks = element_blank()
  )

ggsave("fig_heatmap.pdf",
       (fig_network+inset_element(fig_netbee, 0.79,0.69,1,1))+(pspacer)+fig_heatmap+plot_layout(nrow=1, widths=c(1,0.015,1))+plot_annotation( theme = theme(plot.margin = margin(0,1,0,0))),
       width=10,height=4.2, scale=1)


# ggsave("yearnetwork2.png", g, width=5000, height=5000, scale=0.6, units="px")


### Latent spaces

# Big UMAP for thematic space + change plots, and the second, language & geo plot
pies2 = dopies(festgenrevecs) %>% left_join(ugenre %>% filter(type=="festival") %>% select(V1,V2,event,ab), by="event")

exchange =  c("Berlinale - Berlin IFF", "Sundance", "Venice - Biennale", "Busan IFF", "Toronto - TIFF", "Tallinn Black Nights IFF")

plot.new()# just makes the spline work
fig_thematic = ggplot(pies2, aes(V1, V2))+
  geom_pie_glyph(radius=0.16, slices=attr(pies2, "mains"))+
  geom_point(aes(text=tool), data=ugenre %>% filter(type=="festival", ab=="A"), 
             size=0.5, fill=NA, alpha=0.9, color="black", shape=16)+
  geom_path(data=pies2 %>% ungroup %>%  filter(libelleFestival_NEW_VZ==exchange[1]) %>% select(V1, V2) %>% xspline(shape=-0.2, lwd=2, draw=F) %>% data.frame() %>% rename(V1=x,V2=y) , linewidth=0.5, arrow=arrow(angle=25, length=unit(0.06, "inch"), type="closed"), color="gray10", alpha=0.5)+
  #geom_point(aes(color=colorgenre), data=ugenre %>% filter(type=="film"), size=0.2, fill=NA,show.legend = FALSE)+
  geom_text_repel(aes(label=event,color=colorgenre, size=sizegen %>% log),
                  data=ugenre %>%
                    filter(type %in% c("genre", "film")) %>%
                    group_by(type) %>%
                    filter(type=="genre" | (type=="film" & ( (1:n()) %in% c(15,7, 8, 10, 23, 5)) ) )
                  #mutate(event=case_when(type=="film"~as.character(1:n()), T~type))
                  ,
                  lineheight=0.55,
                  alpha=0.9, hjust=1,
                  bg.color=alpha("white", 0.7),
                  show.legend = FALSE, max.overlaps=1000,
                  box.padding = 0.2, point.padding = 10, min.segment.length = 999,
                  force_pull=0.1,force=2,  max.time = 2)+
  scale_size(range=c(2.8,3.6), guide="none")+
  #scale_radius_manual(values=c(0.2, 0.14), guide="none")+
  #scale_shape_manual(values=abshapes, name="Festival\ntype")+
  
  #annotate("text", -Inf, Inf, hjust=0, vjust=1, label="A", fontface='bold')+
  
  scale_fill_manual(values=sixgcols, name="", breaks=attr(pies2, "mains") )+
  scale_color_manual(values=sixgcols, guide="none", breaks=attr(pies2, "mains") )+
  scale_x_continuous(expand=c(0.01,0))+
  scale_y_continuous(expand=c(0.01,0))+
  theme_void()+
  labs(tag="A")+
  theme(plot.tag=element_text(face = "bold",hjust=0,vjust=1, margin=margin(0,-9,-10,0)))+
  theme(legend.position = c(0.75,0.2),
        plot.background = element_rect( fill="white")
  )


fig_change = divyearplot(genrechangevals %>% 
                           mutate(target=dif, lower=diflower, upper=difupper), 
                         ex=exchange, 
                         yl="Yearly change", ylims = c(-0.005,0.09),
                         convert=convertplots, splitab=F, networkstats = F, 
                         genrecols = sixgcols, 
                         pies=attr(festgenrevecs, "tops"),
                         nrowmod=1, contrdiv = F
)+
  scale_y_continuous(breaks=c(0,0.04,0.08))+
  geom_path(aes(YearEvent, target), data=genrechangevals %>%
              filter(libelleFestival_NEW_VZ==exchange[1]) %>%  
              mutate(target=dif), inherit.aes=F)+
  labs(tag="B")+
  theme(plot.tag=element_text(face = "bold",hjust=1,vjust=1, margin=margin(-5,-9,-11,0)), 
        plot.margin = margin(1,1,0,0))+theme(axis.title.y=element_blank())


ggsave("fig_thematic_map.pdf", fig_thematic+pspacer+fig_change+plot_layout(widths  = c(0.8,0.01,0.2)), width=10,height=5)



# other 2 maps

piesg = dopies(festcountryvecs) %>% left_join(festcountryvecs %>%  select(long,lat,EventID_VZ,ab), by="EventID_VZ")

fig_world = ggplot()+
  #geom_text(aes(V1, V2, color=colorgenre,label=libelleFestival_NEW_VZ), data=ulang %>% sample_n(30), size=3, show.legend = F, inherit.aes = F, alpha=0.9)+
  geom_polygon(mapping= aes(long, lat, group = group), data=mp, inherit.aes = F, color="gray80", fill=NA, size=0.1)+
  geom_pie_glyph(aes(long, lat), radius=0.13, data=piesg, slices=attr(piesg, "mains"))+
  geom_point(aes(long, lat), data=festcountryvecs %>% filter(ab=="A")
             , size=0.4, color="black")+
  #scale_size(range=c(0.3,0.1), guide="none")+
  #scale_radius_manual(values=c(0.16, 0.12), guide="none")+
  #scale_shape_manual(values=c(16,3), name="Festival\ntype")+
  scale_fill_manual(values=sixcols, 
                    name="", breaks=levels(festcountryvecs$colorgenre)
  )+ # add breaks if multiple geoms
  coord_cartesian(xlim=range(festcountryvecs$long)*1.01, ylim=range(festcountryvecs$lat)*1.03,expand=F)+
  theme_void()+
  labs(tag="A")+
  theme(plot.tag=element_text(face = "bold",hjust=0,vjust=1, margin=margin(0,-9,-10,0)))+
  theme(legend.position = c(0.08,0.27),
        legend.background = element_rect(color="transparent", fill="transparent"),
        plot.background = element_rect( fill="white"),
        axis.text = element_text(size=3.8)
  )+
  NULL
fig_eu = fig_world+coord_cartesian(xlim=c(-15, 25), ylim=c(30,51.5))+
  theme(#axis.text = element_blank(), 
    legend.position = "none")+labs(tag=NULL)

fig_eventmap = ggplot(festcountryvecs %>%  # small event location inset
                        group_by(paste(eventlong, eventlat)) %>% 
                        mutate(eventlong=eventlong+runif(n(), -log(n()), log(n()) )) %>% 
                        mutate(eventlat=eventlat+runif(n(), -log(n()), log(n()) ))
)+
  coord_cartesian(xlim=range(festcountryvecs$eventlong)*1.05, ylim=range(festcountryvecs$eventlat)*1.07,expand=F)+
  geom_polygon(mapping= aes(long, lat, group = group), data=mp, inherit.aes = F, color="gray80", fill=NA, size=0.2)+
  geom_point(aes(eventlong, eventlat), festcountryvecs,alpha=0.2, position=position_jitter(3,3), color="black", size=1)+
  theme_void()+
  labs(tag="B")+
  theme(plot.tag=element_text(face = "bold",hjust=0,vjust=1, margin=margin(0,-9,-10,0)))+
  theme(legend.position = 'none',
        plot.background = element_rect(fill="white")
  )


fig_geomap = fig_world+inset_element(fig_eu, 0.55,0,1,0.5, align_to = "full")+
  inset_element(fig_eventmap, 0.8,0.8,1,1, align_to = "full")




piesl = dopies(festlangvecs) %>% left_join(ulang %>%  select(V1,V2,EventID_VZ,ab), by="EventID_VZ")

fig_langmap = ggplot(piesl, mapping=aes(V1, V2))+
  #geom_text(aes(V1, V2, color=colorgenre,label=libelleFestival_NEW_VZ), data=ulang %>% sample_n(30), size=3, show.legend = F, inherit.aes = F, alpha=0.9)+
  geom_pie_glyph(radius=0.16,data=piesl, slices=attr(piesl, "mains"))+
  
  geom_point(aes( text=tool), data=ulang %>% filter(ab=="A") , size=0.6, color="black")+
  #scale_size(range=c(0.3,0.1), guide="none")+
  #scale_radius_manual(values=c(0.16, 0.12), guide="none")+
  #scale_shape_manual(values=c(16,3), name="Festival\ntype")+
  scale_fill_manual(values=langcols, 
                    name="", breaks=levels(ulang$colorgenre)
  )+ # add breaks if multiple geoms
  theme_void()+
  labs(tag="C")+
  theme(plot.tag=element_text(face = "bold",hjust=0,vjust=1, margin=margin(0,-9,-10,0)))+
  scale_x_continuous(expand=c(0.005,0))+
  theme(legend.position = c(0.07,0.28),
        legend.background = element_rect(color="transparent", fill="transparent"),
        plot.background = element_rect( fill="white")
  )+
  NULL

ggsave("fig_geolang_map.pdf", 
       fig_geomap+
         (pspacer)+
         fig_langmap+plot_layout(heights = c(0.55,0.005, 0.45)), width=10,height=7)


## Diversity

fig_g1=ggplot(festrolevecs, mapping=aes(ratiopf, ratiodf, color=ratiof, shape=ab, size=ab, text=tool))+
  geom_vline(xintercept = 0.5)+
  geom_hline(yintercept = 0.5)+
  geom_point()+
  scale_size_manual(values =c(1.4,0.8)*1.2, guide="none")+
  scale_shape_manual(values=abshapes, name="Festival\ntype")+
  scale_color_gradientn(colours =rolecols, name="Fraction\nwomen", limits=c(0,1))+
  theme_bw()+
  theme(#legend.position = 'none',
    plot.background = element_rect( fill="white"),
    plot.margin = margin(1,2,-5,0)
  )+
  labs(x="Producers", y="Directors")

fig_gendersys = sysplots(attr(festrolevecs,"sysdivs"), yl="Circuit-wide weighted fraction of women")+lims(y=c(0.2,0.45))+ theme(legend.position = c(0.75,0.1),legend.direction = "horizontal", legend.title = element_blank())

fig_gex = divyearplot(festrolevecs  %>% mutate(target=ratiof, lower=rfminboot, upper=rfmaxboot),
                      ex=exrolefests, yl="% of women", convert=convertplots, splitab=F, networkstats = F, ylims=c(0,0.59), nrowmod = 2,  contrdiv = F)+theme(axis.title = element_blank())

fig_gender = tags(fig_g1, "A", margin(0,-11,-10,0)) + 
  (tags(fig_gendersys,"B")+inset_element(
    plotbees(festrolevecs %>% mutate(target=ratiof), "C", "") , 
    0.02,0.6,0.4,0.98 ))+
  tags(fig_gex, "D")+plot_annotation(theme=theme(plot.margin = margin(0,0,-5,0)))+
  plot_layout(widths = c(0.8,1.1, 0.9))

ggsave("fig_gender.pdf", fig_gender, width=10,height=4)




tbees = 
  plotbees(festgenrevecs %>% mutate(target=indiv), "C", tlab=" Internal", col="gray50", cex=2.8, vw=1.5) +
  plotbees(festgenrevecs %>% mutate(target=exdiv), "D", tlab=" Contributing", col="gray50", cex=2)+plot_layout(nrow=2)

tdiv2 = (div2plot(festgenrevecs, ex=NULL, convert=F, genrecols = sixgcols %>% darken(0.1))+
           theme(legend.position = "right")+
           guides(color = guide_legend(keyheight = unit(0.6, "lines"), override.aes = list(size = 4)),
                  shape = guide_legend(keyheight = unit(0.6, "lines"))
           )
) %>% tags("E")

fig_thematic_div = 
  (tags(sysplots(attr(festgenrevecs,"sysdivs"), yl="Circuit thematic diversity")+
          theme(legend.position = c(0.85,0.22), legend.title = element_blank()),"A") +
     tags(divyearplot(festgenrevecs  %>% mutate(target=indiv, lower=ilower, upper=iupper), 
                      ex=exgenfests, yl="Internal thematic diversity", convert=F, splitab=F, networkstats = F, nrowmod = 2, pies = pies2, genrecols = sixgcols), "B") + plot_layout(widths = c(0.4,0.65))
  ) /
  ((tbees | tdiv2 | pspacer)+plot_layout(widths = c(0.3,0.7,0), guides = "keep"))


ggsave("fig_thematic_div.pdf", fig_thematic_div, width=10,height=6)


# geo

tbees = 
  plotbees(festcountryvecs %>% mutate(target=indiv), "C", tlab=" Internal", col="gray50", cex=3, vw = 1.4) +
  plotbees(festcountryvecs %>% mutate(target=exdiv), "D", tlab=" Contributing", col="gray50", cex=2.9, vw=1.35)+plot_layout(nrow=2)

tdiv2 = (div2plot(festcountryvecs, ex=NULL, convert=F, genrecols = sixcols %>% darken(0.05))+
           theme(legend.position = "right")+
           guides(color = guide_legend(keyheight = unit(0.6, "lines"), override.aes = list(size = 4)),
                  shape = guide_legend(keyheight = unit(0.6, "lines"))
           )
) %>% tags("E")

fig_geo_div = 
  (tags(sysplots(attr(festcountryvecs,"sysdivs"), yl="Circuit geographic diversity")+
          theme(legend.position = c(0.6,0.22), legend.title = element_blank()),"A") +
     tags(divyearplot(festcountryvecs  %>% mutate(target=indiv, lower=ilower, upper=iupper), 
                      ex=exgeofests, yl="Internal geographic diversity", convert=F, splitab=F, networkstats = F, nrowmod = 2, pies=piesg, genrecols = sixcols), "B") + plot_layout(widths = c(0.35,0.65))
  ) /
  ((tbees | tdiv2 | pspacer)+plot_layout(widths = c(0.3,0.7,0), guides = "keep"))+
  plot_annotation(theme=theme(plot.margin=margin(0,0,0,0)))


ggsave("fig_geo_div.pdf", fig_geo_div, width=10,height=6)


# language

tbees = 
  plotbees(festlangvecs %>% mutate(target=indiv), "C", tlab=" Internal", col="gray50", cex=3, vw = 1.7) +
  plotbees(festlangvecs %>% mutate(target=exdiv), "D", tlab=" Contributing", col="gray50", cex=2.3, vw=1.5)+plot_layout(nrow=2)

tdiv2 = (div2plot(festlangvecs, ex=NULL, convert=F, genrecols = langcols)+
           theme(legend.position = "right")+
           guides(color = guide_legend(keyheight = unit(0.6, "lines"), override.aes = list(size = 4)),
                  shape = guide_legend(keyheight = unit(0.6, "lines"))
           )
) %>% tags("E")

fig_lang_div = 
  (tags(sysplots(attr(festlangvecs,"sysdivs"), yl="Circuit language diversity")+
          theme(legend.position = c(0.6,0.1),legend.direction = "horizontal", legend.title = element_blank()),"A") +
     tags(divyearplot(festlangvecs  %>% mutate(target=indiv, lower=ilower, upper=iupper), 
                      ex=exlangfests, yl="Internal language diversity", convert=F, splitab=F, networkstats = F, nrowmod = 2, pies=piesl, genrecols = langcols), "B")+ plot_layout(widths = c(0.4,0.65))
  ) /
  ((tbees | tdiv2 | pspacer)+plot_layout(widths = c(0.3,0.7,0), guides = "keep"))


ggsave("fig_lang_div.pdf", fig_lang_div, width=10,height=6)





##### Run regression models ####

# Consult the paper for details.

# intercept only prod-dir
festrolevecs %>% mutate(dif=ratiopf-ratiodf) %>% lm(dif~1, data=., weights=nf) %>% summary

attr(festrolevecs, "sysdivs") %>% filter(YearEvent %in% c(2012,2021))

lm(ratiof~ab, festrolevecs, weights=nf) %>% summary() # gender

festrolevecs %>% filter(grepl("Sundance", libelleFestival_NEW_VZ) ) %>% as.data.frame() %>% arrange(YearEvent) %>% select(rfminboot, rfmaxboot)
festrolevecs %>% filter(grepl("SXSW", libelleFestival_NEW_VZ) ) %>% as.data.frame() %>% arrange(YearEvent) %>% select(rfminboot, rfmaxboot)

# language, AB
lm(indiv~ab, festlangvecs, weights=nf) %>% summary()











