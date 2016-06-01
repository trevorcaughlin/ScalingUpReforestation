# ScalingUpReforestation
Models to explore when seed limitation matters for scaling up reforestation from patches to landscapes

Our modeling framework applies to deforested patches that are initially covered by a dense layer of shade-intolerant vegetation ("grass layer") with little or no tree cover. Our primary objective is to predict time-to-canopy-closure – the moment at which tree crowns completely cover the patch and the grass layer is no longer a competitive threat to tree establishment. We model the crown area of tree cohorts over time, dividing cohorts into those beneath and those above the top of the grass layer. Tree height (at the top of the crown) determines whether a cohort is in the grass or canopy layer and is related to crown area by an allometric equation. We assume that the spatial arrangement of individual tree crowns is “perfectly plastic,” filling any available space in the horizontal before overlapping with other tree crowns (Strigul et al. 2008). 
The R code here presents a simple version of the model where tree recruitment in an early-successional patch depends solely on a fixed rate of seed rain from sources external to the patch. The simplicity of this model enables an exact solution for time-to-canopy-closure that users with a wide range of quantitative expertise can explore interactively via a web app (https://t-trevorcaughlin.shinyapps.io/ReforestationDynamics). The code for this webapp is available in the "ui.R" and "server.R" files. R code for a multispecies version of this simple model is available in the file labelled "Appendix S3 multispecies model.R"
