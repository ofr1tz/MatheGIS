# Vorbedingungen

## Installiert Tidyverse falls nötig
if (!require(tidyverse)) {
      install.packages("tidyverse")
      library(tidyverse)
}

# Hilfsfunktionen

## Konvertiert Grad oder Neugrad in Radiant
in_radiant <- function(
      angle=0,
      angle_unit) {
      
      if(angle_unit=="grad" | angle_unit=="gon") rad <- angle*pi/200
      if(angle_unit=="deg") rad <- angle*pi/180
      if(angle_unit=="rad") rad <- angle
      rad
}

## Bestimmt den orientierten Normalenvektor (onv) aus einem directionalVector
## Der onv steht im Drehsin des Koordinatensystems orthogonal zum directionalVector
## Drehsinn geodaetisch -> mit dem Uhrzeigersinn, kartesisch -> gegen die Uhr
onv <- function(directionalVector=c(0,0)) {
      c(-directionalVector[2], directionalVector[1])
}

## Berechnet die Laenge der Strecke zwischen zwei gegebenen Punkten P1 und P2
strecke <- function(
      P1=c(0,0),
      P2=c(0,0)) {
      
      sqrt(streckenquadrat(P1, P2))
}

## Berechnet das Quadrat der Strecke zwischen zwei gegebenen Punkten P1 und P2
streckenquadrat <- function(
      P1=c(0,0),
      P2=c(0,0)) {
      
      (P2[1]-P1[1])^2+(P2[2]-P1[2])^2
}

## Berechnet den Streckenparameter t (p/s3) als Funktion von s1, s2, und s3^2
param_ps3 <- function(s1, s2, s3square) {
    0.5*(1+(s2^2/s3square)-(s1^2/s3square))
}

## Berechnet den Streckenparameter u (h/s3) als Funktion von s1, s2, s3^2 und p/s3
param_hs3 <- function(s1, s2, s3square, ps3) {
    sqrt((s2^2/s3square)-ps3^2)
}

## Berechnet den Lotfusspunkt von P3 als Funktion von P1, P2, s1 und s2
lotfuss <- function(P1, P2, s1, s2) {
    P1+param_ps3(s1, s2, streckenquadrat(P1, P2))*(P2-P1)
}


# Schnittverhalten (3.3)


## Lösung geodätischer Probleme mithilfe orientierter Wegbeschreibung (3.5)

### Polaraufnahme
### Die Punkte P1 und P2 bekannt, die Strecke s2 und der Winkel alpha1 sind gemessen.
### Gesucht ist der Punkt P3.
polaraufnahme <- function(
      P1=c(0,0),
      P2=c(0,0),
      alpha1=0,
      s2=0,
      angle_unit="grad",
      orientation="geodetic",
      full_report=FALSE) {
      
      output <- list()
      
      if(full_report==TRUE) {
            if(angle_unit=="grad" | angle_unit=="gon") max_angle <- 400
            else if(angle_unit=="deg") max_angle <- 360
            else if(angle_unit=="rad") max_angle <- 2*pi
            else stop("Invalid angle unit provided. Unit must be deg for degrees, grad for gradiant or rad for radiant.")
            
            
            if(alpha1<0 | alpha1>max_angle) stop("Invalid angle value provided.")
            
            # Andersrum, wenn kartesisch?
            if(alpha1<max_angle*.25 | alpha1>max_angle*.75) output$P3_direction <- "in front of P1 relative to P1P2"
            if(alpha1<max_angle*.25 | alpha1>max_angle*.75) output$P3_direction <- "behind P1 relative to P1P2"
            if(alpha1<max_angle*.5) output$P3_position <- "right"
            if(alpha1>max_angle*.5) output$P3_position <- "left"
      }
      
      output$s3 <- strecke(P1, P2)
      output$cos_alpha1 <- cos(in_radiant(alpha1, angle_unit))
      output$sin_alpha1 <- sin(in_radiant(alpha1, angle_unit))
      output$directionalVector_P2P1 <- P2-P1
      output$orientedNormalVector_P2P1 <- onv(output$directionalVector_P2P1)
      output$P3 <- P1+
            (s2/output$s3)*
            (output$cos_alpha1*output$directionalVector_P2P1+
             output$sin_alpha1*output$orientedNormalVector_P2P1)
      if(full_report==TRUE) output else output$P3
}


### Vorwaertsschnitt
### Die Punkte P1 und P2 bekannt, die Winkel alpha1 und alpha2 sind gemessen.
### Gesucht ist der Punkt P3.


### Bogenschnitt
### Die Punkte P1 und P2 bekannt, die Strecken s1 und s2 sind gemessen.
### Gesucht ist der Punkt P3.





# Berechnet den Punkt P3 aus P1, P2, s1 und s2
bogenschnitt_P3 <- function(P1, P2, s1, s2, orientation, position) {
    if (position=="right" & orientation=="geodetic" |
        position=="left" & orientation=="cartesian") sign <- 1
    if (position=="left" & orientation=="geodetic" |
        position=="right" & orientation=="cartesian") sign <- -1
        
    s3square <- streckenquadrat(P1, P2)
    t <- param_ps3(s1, s2, s3square)
    u <- param_hs3(s1, s2, s3square, t)
    
    P1+t*(P2-P1)+sign*u*onv(P2-P1)
}


# Returns a complete analysis of the bogenschnitt (-> RMarkdown document)
report_bogenschnitt <- function(
      P1=c(0,0),
      P2=c(0,0),
      s1=0,
      s2=0,
      orientation="geodetic",
      P3_position,
      full_report=FALSE,
      plot=FALSE) {

      output <- list()
      output$s3square <- streckenquadrat(P1, P2)
      output$p_divided_by_s3 <- param_ps3(s1, s2, output$s3square)
      output$h_divided_by_s3 <- param_hs3(s1, s2, output$s3square, output$p_divided_by_s3)
      output$directionalVector_P2P1 <- P2-P1
      output$orientedNormalVector_P2P1 <- onv(output$directionalVector_P2P1)
      output$P3_right <- bogenschnitt_P3(P1, P2, s1, s2, orientation=orientation, position="right")
      output$P3_left <- bogenschnitt_P3(P1, P2, s1, s2, orientation=orientation, position="left")
      
      if(plot==TRUE) plot_bogenschnitt(
            P1=P1, 
            P2=P2, 
            P3_right=output$P3_right, 
            P3_left=output$P3_left, 
            s1=s1, 
            s2=s2,
            orientation=orientation)
            
      
      if(full_report==TRUE) output
      else if(P3_position=="rechts" | P3_position=="right") output$P3_right
      else if(P3_position=="links" | P3_position=="left") output$P3_left
      else list(output$P3_left, output$P3_right)
      
}

### Visualisiert den Bogenschnitt
plot_bogenschnitt <- function(
    P1=c(0,0),
    P2=c(0,0),
    P3_left=NULL,
    P3_right=NULL,
    s1=0,
    s2=0,
    orientation="geodetic") {
    
    if(is.null(P3_left)) P3_left <- bogenschnitt_P3(P1, P2, s1, s2, orientation=orientation, position="left")
    if(is.null(P3_right)) P3_right <- bogenschnitt_P3(P1, P2, s1, s2, orientation=orientation, position="right")
    
    dat <- tribble(
        ~point, ~x, ~y, 
        "P1", P1[1], P1[2], 
        "P2", P2[1], P2[2], 
        "P3 links", P3_left[1], P3_left[2], 
        "P3 rechts", P3_right[1], P3_right[2])
    nudge_x <- (range(dat["x"])[2]-range(dat["x"])[1])/20
    PL <- lotfuss(P1, P2, s1, s2)
    
    
    g <- ggplot()+
        # circle around P1 with radius s2
        annotate(
            "path",
            x=P1[1]+s2*cos(seq(0,2*pi,length.out=100)),
            y=P1[2]+s2*sin(seq(0,2*pi,length.out=100)),
            col="grey")+
        # circle around P2 with radius s1
        annotate(
            "path",
            x=P2[1]+s1*cos(seq(0,2*pi,length.out=100)),
            y=P2[2]+s1*sin(seq(0,2*pi,length.out=100)),
            col="grey")+
        # P1 to P2 (s3)
        geom_segment(aes(
            x=dat[1,"x"],
            y=dat[1,"y"],
            xend=dat[2,"x"],
            yend=dat[2,"y"]),
            arrow=arrow(
                length=unit(12, "points"),
                angle=20,
                type="closed"),
            col="darkgrey")+
        # first path segment from P1 to P3 (p)
        geom_segment(aes(
            x=dat[1,"x"],
            y=dat[1,"y"],
            xend=PL[1],
            yend=PL[2]),
            arrow=arrow(
                length=unit(12, "points"),
                angle=20,
                type="closed"),
            linetype="dashed",
            col="red")+
        # label p
        geom_text(aes(x=mean(c(P1[1], PL[1])),
                      y=mean(c(P1[2], PL[2]))),
                  label="p",
                  nudge_x=nudge_x,
                  hjust=0,
                  col="red",
                  alpha=.5)+
        # second path segment from P1 to P3_left (h)
        geom_segment(aes(
            x=PL[1],
            y=PL[2],
            xend=dat[3, "x"],
            yend=dat[3, "y"]),
            arrow=arrow(
                length=unit(12, "points"),
                angle=20,
                type="closed"),
            linetype="dashed",
            col="red")+
        # second path segment from P1 to P3_right (h)
        geom_segment(aes(
            x=PL[1],
            y=PL[2],
            xend=dat[4, "x"],
            yend=dat[4, "y"]),
            arrow=arrow(
                length=unit(12, "points"),
                angle=20,
                type="closed"),
            linetype="dashed",
            col="red")+
        # labels h
        geom_text(aes(x=mean(c(PL[1], P3_left[1])),
                      y=mean(c(PL[2], P3_left[2]))),
                  label="h",
                  nudge_x=nudge_x,
                  hjust=0,
                  col="red",
                  alpha=.5)+
        geom_text(aes(x=mean(c(PL[1], P3_right[1])),
                      y=mean(c(PL[2], P3_right[2]))),
                  label="h",
                  nudge_x=nudge_x,
                  hjust=0,
                  col="red",
                  alpha=.5)+
        # s2
        annotate(
            "path",
            x=dat[c(3,1,4), "x"],
            y=dat[c(3,1,4), "y"], 
            col="grey")+
        # labels s2
        geom_text(aes(x=mean(c(P1[1], P3_left[1])),
                      y=mean(c(P1[2], P3_left[2]))),
                  label="s2",
                  nudge_x=nudge_x,
                  hjust=0,
                  col="grey",
                  alpha=.67)+
        geom_text(aes(x=mean(c(P1[1], P3_right[1])),
                      y=mean(c(P1[2], P3_right[2]))),
                  label="s2",
                  nudge_x=nudge_x,
                  hjust=0,
                  col="grey",
                  alpha=.67)+
        # s1
        annotate(
            "path",
            x=dat[c(3,2,4), "x"],
            y=dat[c(3,2,4), "y"], 
            col="grey")+
        geom_point(
            data=dat, 
            mapping=aes(x=x, y=y))+
        # labels s1
        geom_text(aes(x=mean(c(P2[1], P3_left[1])),
                      y=mean(c(P2[2], P3_left[2]))),
                  label="s1",
                  nudge_x=nudge_x,
                  hjust=0,
                  col="grey")+
        geom_text(aes(x=mean(c(P2[1], P3_right[1])),
                      y=mean(c(P2[2], P3_right[2]))),
                  label="s1",
                  nudge_x=nudge_x,
                  hjust=0,
                  col="grey")+
        # axis labels
        labs(x="x", y="y")+
        # point labels
        geom_text(
            data=dat, 
            aes(x=x, y=y, label=point),
            nudge_x=nudge_x,
            hjust=0,
            alpha=.5)+
        geom_text(
            aes(x=PL[1], y=PL[2], label="PL"),
            nudge_x=nudge_x,
            hjust=0,
            col="red",
            alpha=.5
        )+
        # fix x and y axis and equal distance between units
        coord_equal()+
        theme_minimal()
    
    # flip coordinate system if orientation is geodetic instead of cartesian
    if(orientation=="geodetic") {
        xlim=ggplot_build(g)$layout$panel_ranges[[1]]$x.range
        ylim=ggplot_build(g)$layout$panel_ranges[[1]]$y.range
        g <- g+coord_flip(xlim=xlim, ylim=ylim)+theme(aspect.ratio=1)
    }
    else if(orientation!="cartesian") warning("Invalid orientation paramenter. Coordinate system defaulted to cartesian")
    
    print(g)
}