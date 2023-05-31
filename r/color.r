library(colorRamps)
library(RColorBrewer)
# standardise colors

info("Loading default colors")
default.cols = function(n){
    #info(sprintf("Getting %s default colors", n))
    if(n<=20){
        #print(n)
        #info("Using 'Kelly' cols")
        kelly.cols(n)
    }else{
        warn("More than 20 requested, using 'Distinct' cols")
        distinct.cols(n)
    }
} 

#https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
brewer_all = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

image_cols = colorRampPalette(c("black", viridis::plasma(6)))(20)

wyrb.heat = colorRampPalette(c("white", "yellow3", "red2", "black"))(20)

# Color maps
brewer16 = c(brewer.pal(9, "Set1"), brewer.pal(7, "Set2"))
brewer16[6] = "khaki2"
brewer16[8] = "lightskyblue2"
isc.subset.cols = rev(brewer.pal(3, "Set1")) #colorRampPalette(material.700[9:12])(3)

brewer20 = c(brewer16, brewer.pal(12, "Set3")[c(3, 9, 8)], "violetred4")
brewer20[3] = brewer.pal(9, "Greens")[7]
brewer20[14] = brewer.pal(9, "Greens")[4]


perceptual_rainbow_12 = c('#873B61', '#8F489D',
                          '#7966CF', '#677BDC', '#5492DF', '#45AAD7', '#3BC0C5', '#3CD2AC',
                          '#47DF91', '#5DE578', '#A1E35F', '#E9D575')

perceptual_rainbow_16 = c('#873B61', '#8F407F', '#8F489D', '#8755B9', 
    '#7966CF', '#677BDC', '#5492DF', '#45AAD7', '#3BC0C5', '#3CD2AC', 
    '#47DF91', '#5DE578', '#7CE767', '#A1E35F', '#C6DC64', '#E9D575')

tol14rainbow=c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C")
tol15rainbow=c("#114477", "#4477AA", "#77AADD", "#117755", "#44AA88", "#99CCBB", "#777711", "#AAAA44", "#DDDD77", "#771111", "#AA4444", "#DD7777", "#771144", "#AA4477", "#DD77AA")
tol18rainbow=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
# ...and finally, the Paul Tol 21-color salute
tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")


## Use these colours for all plots for the Atlas paper.
celltype.cols = list(
    "Goblet"=brewer20[1], 
    "Paneth"=brewer20[18], 
    "Stem"=brewer20[3], 
    "Endocrine"=brewer20[20],
    "Tuft"=brewer20[5],
    "Tuft-1"=brewer20[5],
    "TA.Early"=brewer20[15], 
    "TA.G1"=brewer20[7],
    "TA.G2"=brewer20[8],
    "TA"=brewer20[14],
    "Enterocyte.Progenitor.Early"=brewer20[9], 
    "Enterocyte.Progenitor.Late"= brewer20[10], 
    "Enterocyte.Immature.Proximal"=brewer20[13], 
    "Enterocyte.Progenitor"=brewer20[11], 
    "Enterocyte.Immature.Distal"=brewer20[12], 
    "Enterocyte.Mature.Proximal"=brewer20[4], 
    "Enterocyte"=brewer20[2], 
    "Enterocyte.Mature.Distal"=brewer20[2], 
    "Enterocyte.Immature"=brewer20[6],
    "Enterocyte.Mature"=brewer20[16],
    "Tuft-2"=brewer20[19],
    "M.cell"=perceptual_rainbow_16[9],
    "M-cell"=perceptual_rainbow_16[9],
    "hc-ISC"=isc.subset.cols[3],
    "lc-ISC"=isc.subset.cols[1],
    "mc-ISC"=isc.subset.cols[2],
    "rmv"="black"
)

weather_heat <- function(n)
{
    sp = brewer.pal(11, "Spectral")
    wh = c(sp[2:5], sp[7:10], "lightblue", "grey95")
    rev(colorRampPalette(wh)(n))
}

# Qualitative color schemes by Paul Tol
 tol1qualitative=c("#4477AA")
 tol2qualitative=c("#4477AA", "#CC6677")
 tol3qualitative=c("#4477AA", "#DDCC77", "#CC6677")
 tol4qualitative=c("#4477AA", "#117733", "#DDCC77", "#CC6677")
 tol5qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677")
 tol6qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677","#AA4499")
 tol7qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77", "#CC6677","#AA4499")
 tol8qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677","#AA4499")
 tol9qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499")
 tol10qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")
 tol11qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")
 tol12qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", "#882255", "#AA4499")

# MONOCHROME PALETTES
# sort(brewer.pal(8,"Greens"))
redmono = c("#99000D", "#CB181D", "#EF3B2C", "#FB6A4A", "#FC9272", "#FCBBA1", "#FEE0D2", "#FFF5F0")
greenmono = c("#005A32", "#238B45", "#41AB5D", "#74C476", "#A1D99B", "#C7E9C0", "#E5F5E0", "#F7FCF5")
bluemono = c("#084594", "#2171B5", "#4292C6", "#6BAED6", "#9ECAE1", "#C6DBEF", "#DEEBF7", "#F7FBFF")
grey8mono = c("#000000","#252525", "#525252", "#737373", "#969696", "#BDBDBD", "#D9D9D9", "#F0F0F0")
grey6mono = c("#242424", "#494949", "#6D6D6D", "#929292", "#B6B6B6", "#DBDBDB")

# EQUAL WEIGHT
# Generated with rainbow(12, s = 0.6, v = 0.75)
rainbow12equal = c("#BF4D4D", "#BF864D", "#BFBF4D", "#86BF4D", "#4DBF4D", "#4DBF86", "#4DBFBF", "#4D86BF", "#4D4DBF", "#864DBF", "#BF4DBF", "#BF4D86")
rainbow10equal = c("#BF4D4D", "#BF914D", "#A8BF4D", "#63BF4D", "#4DBF7A", "#4DBFBF", "#4D7ABF", "#634DBF", "#A84DBF", "#BF4D91")
rainbow8equal = c("#BF4D4D", "#BFA34D", "#86BF4D", "#4DBF69", "#4DBFBF", "#4D69BF", "#864DBF", "#BF4DA3")
rainbow6equal = c("#BF4D4D", "#BFBF4D", "#4DBF4D", "#4DBFBF", "#4D4DBF", "#BF4DBF")
 
# Generated with package "gplots" function rich.colors(12)
rich12equal = c("#000040", "#000093", "#0020E9", "#0076FF", "#00B8C2", "#04E466", "#49FB25", "#E7FD09", "#FEEA02", "#FFC200", "#FF8500", "#FF3300")
rich10equal = c("#000041", "#0000A9", "#0049FF", "#00A4DE", "#03E070", "#5DFC21", "#F6F905", "#FFD701", "#FF9500", "#FF3300")
rich8equal = c("#000041", "#0000CB", "#0081FF", "#02DA81", "#80FE1A", "#FDEE02", "#FFAB00", "#FF3300")
rich6equal = c("#000043", "#0033FF", "#01CCA4", "#BAFF12", "#FFCC00", "#FF3300")
 
# Generated with package "fields" function tim.colors(12), which is said to emulate the default matlab colorset
tim12equal = c("#00008F", "#0000EA", "#0047FF", "#00A2FF", "#00FEFF", "#5AFFA5", "#B5FF4A", "#FFED00", "#FF9200", "#FF3700", "#DB0000", "#800000")
tim10equal = c("#00008F", "#0000FF", "#0070FF", "#00DFFF", "#50FFAF", "#BFFF40", "#FFCF00", "#FF6000", "#EF0000", "#800000")
tim8equal = c("#00008F", "#0020FF", "#00AFFF", "#40FFBF", "#CFFF30", "#FF9F00", "#FF1000", "#800000")
tim6equal = c("#00008F", "#005AFF", "#23FFDC", "#ECFF13", "#FF4A00", "#800000")
 
# Generated with sort(brewer.pal(8,"Dark2")) #Dark2, Set2
dark8equal = c("#1B9E77", "#666666", "#66A61E", "#7570B3", "#A6761D", "#D95F02", "#E6AB02", "#E7298A")
dark6equal = c("#1B9E77", "#66A61E", "#7570B3", "#D95F02", "#E6AB02", "#E7298A")
set8equal = c("#66C2A5", "#8DA0CB", "#A6D854", "#B3B3B3", "#E5C494", "#E78AC3", "#FC8D62", "#FFD92F")
set6equal = c("#66C2A5", "#8DA0CB", "#A6D854", "#E78AC3", "#FC8D62", "#FFD92F")

## from Sam Rs compoHeatMap.R
## Get good colors for use in the heatmap. A wrapper function for
## color.palette(). Allows steps to be rescaled so that middle color
## corresponds to a given value in the range.
### ARGS:
## steps: vector of colors.
### n.steps.between: integer vector of #steps between each color given
## range.val: numeric vector of length 2 giving lower and upper limits
## of range of values that will be plotted;
### mid.val: numeric value in range.val that should be represented by
### the middle index (rounding up) in the steps vector; ignored if
### range.val is NULL.
## n.steps.final: the number of colors desired in the output vector.
### ...: Unspecified arguments are sent to color.palette().
## RETURNS:
### a vector of length n.steps.final.
get.hmap.col <- function(steps=c("blue", "cyan", "yellow", "red"), n.steps.between=c(9,1,10),
                         range.val=NULL, mid.val=NULL, n.steps.final=30,...) {
    if (!is.null(range.val)) {
        if (is.null(mid.val)) {
            mid.val=(range.val[2]-range.val[1])/2.0
        }
        mid.index=ceiling(length(n.steps.between)/2)
        ## fractional steps in the low vs. high ranges
        frac.steps.low=n.steps.between[1:mid.index]/sum(n.steps.between[1:mid.index])
        frac.steps.high=n.steps.between[(mid.index+1):length(n.steps.between)]/sum(n.steps.between[(mid.index+1):length(n.steps.between)])
        ## fraction of actual values in the low vs. high ranges
        frac.low=(mid.val-range.val[1])/(range.val[2]-range.val[1])
        frac.high=(range.val[2]-mid.val)/(range.val[2]-range.val[1])
        ## Get the right resolution and scale:
        ## n.steps is the total number of steps that will be used
        n.steps=max(255, ceiling(10^abs(log10(min(frac.low*frac.steps.low)))), ceiling(10^abs(log10(min(frac.high*frac.steps.high)))))
        ## sum(frac.high*frac.steps.high)+sum(frac.low*frac.steps.low) == 1
        n.steps.low=round(frac.low*frac.steps.low*n.steps)
        n.steps.high=round(frac.high*frac.steps.high*n.steps)
        n.steps.between=c(n.steps.low, n.steps.high)
    }
    hmcol=color.palette(steps=steps, n.steps.between=n.steps.between, ...)(n.steps.final)
    return(hmcol)
}

## Wrapper function for colorRampPalette based on
## http://stackoverflow.com/questions/13327326/r-image-function-in-r
## It allows for the definition of the number of intermediate colors
## between the main colors.  Using this option, one can stretch out
## colors that should predominate the palette spectrum. Additional
## arguments of colorRampPalette can also be added regarding the type
## and bias of the subsequent interpolation.
### ARGS:
## steps: integer.
### n.steps.between: NULL or integer.
## ...: Unspecified arguments sent to colorRampPalette().
### RETURNS:
## a color palette function, as returned by colorRampPalette.
### Usage:
## Compare pal.1 <- colorRampPalette(c("blue", "cyan", "yellow",
### "red"), bias=1)
## with
### pal.2 <- color.palette(c("blue", "cyan", "yellow", "red"),
### n.steps.between=c(10,1,10))
color.palette <- function(steps, n.steps.between=NULL, ...){
    if(is.null(n.steps.between)) n.steps.between <- rep(0, (length(steps)-1))
    if(length(n.steps.between) != length(steps)-1) stop("Must have one less n.steps.between value than steps")
    fill.steps <- cumsum(rep(1, length(steps))+c(0,n.steps.between))
    RGB <- matrix(NA, nrow=3, ncol=fill.steps[length(fill.steps)])
    RGB[,fill.steps] <- col2rgb(steps)
    for(i in which(n.steps.between>0)){
        col.start=RGB[,fill.steps[i]]
        col.end=RGB[,fill.steps[i+1]]
        for(j in seq(3)){
            vals <- seq(col.start[j], col.end[j], length.out=n.steps.between[i]+2)[2:(2+n.steps.between[i]-1)]
            RGB[j,(fill.steps[i]+1):(fill.steps[i+1]-1)] <- vals
        }
    }
    new.steps <- rgb(RGB[1,], RGB[2,], RGB[3,], maxColorValue = 255)
    pal <- colorRampPalette(new.steps, ...)
    return(pal)
}

# https://www.materialui.co/colors
material.cols <- c("#f44336", #red
    "#E91E63", #pink
    "#9C27B0", #purple
    "#673AB7", #deep purple 
    "#3F51B5",  # indigo
    "#2196F3", # blue
    "#03A9F4", # light blue
    "#00BCD4", #cyan
    "#009688", # teal
    "#4CAF50", #green
    "#8BC34A", #light green
    "#CDDC39", # lime
    "#FFEB3B", #yellow
    "#FFC107", # amber
    "#FF9800", # organe
    "#FF5722", # deep orange
    "#795548", #brown
    "#9E9E9E", # grey
    "#607D8B" #blue grey
    )

isc.subset.cols = rev(brewer.pal(3, "Set1")) #colorRampPalette(material.700[9:12])(3)

# reverse engineered from the darksky weather app
darksky <- function(n)
{
    colorRampPalette(c(rgb(110/255, 41/255, 132/255, 1), # magenta
    rgb(33/255, 46/255, 115/255), # navy
    rgb(25/255, 96/255, 155/255), # blue
    #rgb(50/255, 150/255, 86/255), # blue 2
    rgb(80/255, 170/255, 183/255), # seafoam
    rgb(127/255, 200/255, 178/255), # teal
    rgb(235/255, 238/255, 207/255), # goldenrod
    rgb(246/255, 226/255, 155/255), # light orange
    rgb(247/255, 170/255, 86/255), # orange
    rgb(240/255, 92/255, 38/255), #scarlet
    rgb(142/255, 40/255, 11/255), #maroon
    rgb(99/255, 27/255, 7/255)) # deep red
    )(n)
}

material.heat <- function(n)
{

    mh = c(
        #"#607D8B", #blue grey
        "#283593", #indigo 800
        "#3F51B5",  #indigo
        "#2196F3", # blue
        #"#03A9F4", # light blue
        "#00BCD4", #cyan
        #"#009688", # teal
        "#4CAF50", #green
        "#8BC34A", #light green
        "#CDDC39", # lime
        "#FFEB3B", #yellow
        "#FFC107", # amber
        "#FF9800", # organe
        "#FF5722", # deep orange)
        "#f44336")
    colorRampPalette(mh)(n)
}

material.heat.new <- function(n)
{

    mh = c(
        "black",
        #"grey10",
        "grey20",
        "#2D3B79", #blue grey
        "#283593", #indigo 800
        "#3F51B5",  #indigo
        "#2196F3", # blue
        "#00BCD4", #cyan
        "#009688", # teal
        "#4CAF50", #green
        "#8BC34A", #light green
        "#CDDC39", # lime
        "#FFEB3B", #yellow
        "#FFC107", # amber
        "#FF9800", # organe
        "#FF5722", # deep orange)
        "#f44336")
    colorRampPalette(mh)(n)
}

material.700 = c("#d32f2f",
                    "#C2185B",
                    "#7B1FA2",
                    "#512DA8",
                    "#303F9F",
                    "#1976D2",
                    "#0288D1",
                    "#0097A7",
                    "#00796B",
                    "#388E3C",
                    "#689F38",
                    "#AFB42B",
                    "#FBC02D",
                    "#FFA000",
                    "#F57C00",
                    "#E64A19",
                    "#5D4037",
                    "#616161",
                    "#455A64")

material.800.heat <- function(n)
{
    m8h = c("#37474F",
    "#283593",
    "#1565C0",
    "#0277BD",
    "#00838F",
    "#00695C",
    "#2E7D32",
    "#558B2F",
    "#9E9D24",
    "#F9A825",
    "#FF8F00",
    "#EF6C00",
    "#D84315")
    colorRampPalette(m8h)(n)
}


flat.cols <- function(n)
{
    fc = c("#34495e", #wet asphalt
    "#9b59b6",  #amythest 
    "#3498db",  #peter river
    "#2ecc71",  # emerald
    #"#1abc9c",  #turquiose
    "#f1c40f",  # sunflower 
    "#e67e22",   # carrot
    "#e74c3c")   # alizarin
    #"#ecf0f1",   #clouds
    #"#95a5a6")   # concrete
    return(colorRampPalette(fc)(n))

}

heat.cols.adam <- function(n)
{
    rev(colorRampPalette(c("honeydew2", "lightgoldenrod1", "darkgoldenrod1", "firebrick1"))(n))
}

hmap.cols <- function(min.val, mid.val, max.val)
{
    get.hmap.col(mid.val=mid.val, range.val=c(min.val, max.val), 
                                steps=c("midnightblue", "blue", "cyan", "yellow", "red"), 
                                n.steps.between=c(30, 60, 60, 60))
}

## a combination of set1 and set2 with similar pinks and ugly yellow replaced 
brewer16 = c(brewer.pal(9, "Set1"), brewer.pal(7, "Set2"))
brewer16[6] = "khaki2"
brewer16[8] = "lightskyblue2"

brewer20 = c(brewer16, brewer.pal(12, "Set3")[c(3, 9, 8)], "violetred4")
brewer20[3] = brewer.pal(9, "Greens")[7]
brewer20[14] = brewer.pal(9, "Greens")[4]


distinct.cols <- function(n)
{
    return (intense.100 [1:n])
}

kelly.cols <- function(n)
{
    if(n <= 20)
    {
        return(kelly[1:n])
    }else
    {
        warn("Only 20 kelly colours available")
        return(kelly)
    }
}

# cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# color.blind.friendly.cols <- function(n)
# {
#      if (n < length(cbbPalette))
#     {
#         return(cbbPalette[1:n])
#     }else
#     {
#         warn(sprintf("Dont have that many colourblind friendly colours. Returning %s maximally contrasting colors.", n))
#         return (c(cbbPalette[1:length(cbbPalette)], cols[1:n-length(cbbPalette)])
#     }
# }

intense.cols <-function(n)
{
    if (n < length(intense))
    {
        return(intense[1:n])
    }else if( n< length(intense.100))
    {
        return(intense.100[1:n])
    }else
    {
        warn(sprintf("Dont have that many intense colours. Returning %s maximally contrasting colors.", n))
        return (distinct.cols(n))
    }
    
}

intense.100 = c( '#6CC839','#0DA1E8','#F243F5','#3D3D12','#512257','#2BAF98','#E17516',
    '#F6A1D7','#4751C5','#D48B79','#025466','#A1AF65','#70091B','#A82F9A','#2F6401',
    '#FA3707','#A18D06','#95B0C8','#155D44','#E3A359','#F282F0','#5FC981','#C3051E',
    '#5C2D37','#993700','#E11F90','#B3C424','#94A0F5','#07326B','#0499B3','#A977F9',
    '#F6774F','#FD7C7E','#A3691D','#B33FCE','#2F3442','#C5AADA','#5485F7','#7B1E52',
    '#DC9EB2','#819722','#2176C8','#0E9164','#338D2E','#09602A','#EC61AB','#0B3D37',
    '#3B3280','#725515','#7FBAA8','#5C9A18','#611D6E','#F7A77C','#8B0938','#FB1A2A',
    '#F59A9E','#95150E','#EF5447','#35B4E5','#FC2C4F','#1D4260','#EB93E2','#6BC968',
    '#66CCA2','#1F7684','#C82EB2','#F1905E','#ED127F','#2B4A2F','#7EBCC0','#644037',
    '#ACB4DE','#81C1E0','#CCA660','#5B1E3C','#FA3774','#81A200','#7C79EB','#87C875',
    '#8F1121','#473142','#DA5BEE','#751D61','#EC63D9','#437C04','#9059DF','#A3B235',
    '#56450F','#13A367','#E47EFC','#E1AD4F','#236E42','#3875B8','#1B5118','#3F490F',
    '#F1153C','#5A2F29','#0F6050','#1E6976')

# http://tools.medialab.sciences-po.fr/iwanthue/
intense = c("#F41AA7",
    "#12802C",
    "#F56700",
    "#4371B9",
    "#5C270E",
    "#B9BD09",
    "#D99B68",
    "#157D71",
    "#40334F",
    "#D80F35",
    "#CFA6DA",
    "#B8085C",
    "#B633C3",
    "#3D4326",
    "#6447AD",
    "#76B57E",
    "#F5A83E",
    "#F570BE",
    "#851122",
    "#31360C")

cols = c(
        "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
        "#5A0007", "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
        "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
        "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
        "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
        "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
        "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
        
        "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
        "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
        "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
        "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
        "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
        "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800",
        "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
        "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58",
        
        "#7A7BFF", "#D68E01", "#353339", "#78AFA1", "#FEB2C6", "#75797C", "#837393", "#943A4D",
        "#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A", "#001325", "#02525F", "#0AA3F7", "#E98176",
        "#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", "#02684E", "#962B75", "#8D8546", "#9695C5",
        "#E773CE", "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B", "#E704C4",
        "#00005F", "#A97399", "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", "#643127", "#513A01",
        "#6B94AA", "#51A058", "#A45B02", "#1D1702", "#E20027", "#E7AB63", "#4C6001", "#9C6966",
        "#64547B", "#97979E", "#006A66", "#391406", "#F4D749", "#0045D2", "#006C31", "#DDB6D0",
        "#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9", "#FFFFFE", "#C6DC99", "#203B3C",

        "#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", "#CCAA35", "#374527", "#8BB400", "#797868",
        "#C6005A", "#3B000A", "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C", "#B88183",
        "#AA5199", "#B5D6C3", "#A38469", "#9F94F0", "#A74571", "#B894A6", "#71BB8C", "#00B433",
        "#789EC9", "#6D80BA", "#953F00", "#5EFF03", "#E4FFFC", "#1BE177", "#BCB1E5", "#76912F",
        "#003109", "#0060CD", "#D20096", "#895563", "#29201D", "#5B3213", "#A76F42", "#89412E",
        "#1A3A2A", "#494B5A", "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#958A9F",
        "#BDC9D2", "#9FA064", "#BE4700", "#658188", "#83A485", "#453C23", "#47675D", "#3A3F00",
        "#061203", "#DFFB71", "#868E7E", "#98D058", "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66",
        
        "#2F5D9B", "#6C5E46", "#D25B88", "#5B656C", "#00B57F", "#545C46", "#866097", "#365D25",
        "#252F99", "#00CCFF", "#674E60", "#FC009C", "#92896B")

kelly = c(
    "#00538A", # Strong Blue
    "#C10020", # Vivid Red
    "#007D34", # Vivid Green
    "#FFB300", # Vivid Yellow
    "#803E75", # Strong Purple
    "#FF6800", # Vivid Orange
    "#A6BDD7", # Very Light Blue
    "#CEA262", # Grayish Yellow
    "#817066", # Medium Gray
    "#F6768E", # Strong Purplish Pink
    "#FF7A5C", # Strong Yellowish Pink
    "#53377A", # Strong Violet
    "#FF8E00", # Vivid Orange Yellow
    "#B32851", # Strong Purplish Red
    "#F4C800", # Vivid Greenish Yellow
    "#7F180D", # Strong Reddish Brown
    "#93AA00", # Vivid Yellowish Green
    "#593315", # Deep Yellowish Brown
    "#F13A13", # Vivid Reddish Orange
    "#232C16" # Dark Olive Green
    )

color.kelly <- function(n)
{
   if(n < length(kelly))
   {
	return (sample(kelly, n))
   }else{
	return (colorRampPalette(kelly)(n))
   }
}

# plot some colors
 pal <- function(col, border = "light gray", ...)
 {
     n <- length(col)
     plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
     axes = FALSE, xlab = "", ylab = "", ...)
     rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
 }

 rainbow.hcl <- function(n=10)
 {
    library(colorspace)
    rainbow_hcl(n, c = 50, l = 70, start = 0, end = 360*(n-1)/n,
        gamma = NULL, fixup = TRUE, alpha = 1)
 }

 heat.hcl <- function(n=10)
 {
    library(colorspace)
    rev(heat_hcl(n, c = c(80, 30), l = c(30, 90), power = c(1/5, 2)))
 }

intense.16 = c("#FFBC6B", "#727EFF", "#8FE300", "#FF075F", "#01EB50", "#9A5D9E", "#F1DE3B", "#85B6FF", "#9A6A0F", "#25D0FF", "#BB4F67", "#00E59F", "#FFC6C7", "#00873A", "#00A1B2",
"#2D8163")
 

intense.16.2 = toupper(c("#bfca41", "#ab48f9", "#019f25", "#b5318e", "#bf9600", "#0190f3", "#ff8e39", "#eb98ff", "#007a54", "#ff649e", "#94ccc6", "#ca232d", "#bcc0da", "#a9466f", "#bfc68f", "#81624d"))

ad.cubehelix.old <- function(n)
{
    library(rje)
    cubeHelix(n, gamma = 0.4, hue=1.6, r=-0.9, start = -3)[ceiling(n/10):n]
}

ad.cubehelix <- function(n)
{
    library(rje)
    cubeHelix(n, gamma = 0.4, hue=1.6, r=-0.9, start = -3)[ceiling(n/10):floor(n*9/10)]
}

distinct.phrogz <- function(n)
{
    x = c("#BF3030","#FFA280","#A68A53","#EEF2B6","#00BF80","#206C80","#262A33","#7340FF","#D9A3CE","#331A1A","#B2622D","#FFF240","#829926","#2D594A","#0088CC",
    "#535EA6","#655673","#E63995","#A67C7C","#593E2D","#333226","#41F200","#7CA698","#002999","#BFC8FF","#B630BF","#592D3E","#F24100","#F29D3D","#555916",
    "#24661A","#3DE6F2","#101D40","#110040","#6D1D73","#66001B")
    return(x[1:n])
}

# from python package 'palettable'
cubehelix.classic16 = c('#000000', '#160A22', '#182044', '#103E53', 
    '#0E5E4A', '#237433', '#507D23', '#8A7A2D', '#BE7555', '#DA7991', 
    '#DB8ACB', '#CCA7F0', '#BFC9FB', '#C3E5F4', '#DCF6EF', '#FFFFFF')

cubehelix1.16 = c('#000000', '#1B0F00', '#411704', '#681B20', 
    '#85214B', '#932D7E', '#9042AF', '#8160D2', '#6F83E3', 
    '#63A6E2', '#65C5D3', '#78DBC2', '#99E9B9', '#C1F0BF', '#E6F5D8', '#FFFFFF')

cubehelix2.16 = c('#000000', '#001C0E', '#00332F', '#07415B', 
    '#234787', '#4E48A8', '#8148B8', '#B14DB5', '#D65AA5', 
    '#EB718F', '#EE8E80', '#E6AF7F', '#DBCE90', '#D8E7B2', '#E2F7DB', '#FFFFFF')

jim_special_16 = c('#000000', '#160A22', '#251944', '#2F2B63', 
    '#34417D', '#375892', '#3B70A0', '#4089A9', '#4AA0AD', '#59B5AF',
     '#6DC7B1', '#86D6B4', '#A3E3BD', '#C3EDCB', '#E2F6E1', '#FFFFFF')



red_16 = c('#000000', '#130C23', '#2C1641', '#49205A', 
    '#68296B', '#863576', '#A2437C', '#B9537E', '#CC6680', 
    '#D87C82', '#E19488', '#E5AC93', '#E8C4A4', '#ECDBBD', '#F2EEDB', '#FFFFFF')



#http://godsnotwheregodsnot.blogspot.com/2012/09/color-distribution-methodology.html
distinct64 =c("#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
"#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
"#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
"#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
"#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
"#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
"#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
"#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
"#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
"#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
"#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
"#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
"#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")

color.brewer <- function(n)
{
    library(RColorBrewer)
    if(n < 9)
    {
        return (brewer.pal(n, "Set1"))
    }else if(n < 18)
    {
        return (c(brewer.pal(9, "Set1"), brewer.pal(19-n, "Set2")))
    }else if(n < 27)
    {
        return (c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), c(brewer.pal(26-n, "Dark2"))))
    }
}
