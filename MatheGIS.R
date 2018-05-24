# Vorbedingungen

## Installiert Tidyverse falls nötig
if (!require(tidyverse)) {
      install.packages("tidyverse")
      library(tidyverse)
}

# Hilfsfunktionen

## Rechnet Winkelmass um
umrechnenWinkel <- function(
    alpha=0,
    quelle=c("deg", "grad", "gon", "rad"),
    ziel=c("deg", "grad", "gon", "rad")) {
    
    if(quelle==ziel) return(alpha)
    
    halbkreis <- list(deg=180, grad=200, gon=200, rad=pi)
    alpha*halbkreis[[ziel]]/halbkreis[[quelle]]
}

## Berechnet den Kotangens
cot <- function(alpha) 1/tan(alpha)

## Berechnet Winkel aus drei Seiten mithilfe des Kosinussatzes
## und gibt ihn in gewünschter Einheit aus
winkel_kosinus_seiten <- function(s1, s2, s3, angle_unit=c("deg", "grad", "gon", "rad")) {
    umrechnenWinkel(acos((s1^2-s2^2-s3^2)/(2*s2*s3)), 
                    "rad", 
                    angle_unit)
}

## Berechnet Winkel an Punkt 1 aus drei Punkten mithilfe des Kosinussatzes
## und gibt ihn in gewünschter Einheit aus
winkel_kosinus_punkte <- function(P1, P2, P3, angle_unit=c("deg", "grad", "gon", "rad")) {
    
    umrechnenWinkel(acos(-1*(skalarprodukt(P1-P3, P2-P1)/(strecke(P3, P1)*strecke(P1, P2)))),
                    "rad",
                    angle_unit)
}

## Berechnet Winkel an Punkt 1 aus drei Punkten mithilfe zweier Skalarprodukte und Kotangens
## und gibt ihn in gewünschter Einheit aus
winkel_kotangens_punkte <- function(P1, P2, P3, angle_unit=c("deg", "grad", "gon", "rad")) {
    
    umrechnenWinkel(atan(1/(c(1,-1)*(skalarprodukt(P3-P1, P2-P1)/skalarprodukt(P3-P1, onv(P2-P1))))),
                    "rad",
                    angle_unit)
}

## Berechnet das Skalarprodukt zweier Vektoren aus Koordinaten
skalarprodukt <- function(V1,V2) {
    V1[1]*V2[1]+V1[2]*V2[2]
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
param_hs3 <- function(s2, s3square, ps3) {
    sqrt((s2^2/s3square)-ps3^2)
}

## Berechnet den Lotfusspunkt von P3 als Funktion von P1, P2 und dem Streckenparameter t
lotfuss <- function(P1, P2, t) {
    P1+t*(P2-P1)
}


# Schnittverhalten (3.3)


## Lösung geodätischer Probleme mithilfe orientierter Wegbeschreibung (3.5)

#### Bestimmt die Lage von P3 relativ zu P1 und der Strecke P1P2 aus alpha1 (Polaraufnahme und Vorwaertsschnitt)
quadrant_P3 <- function(alpha1, angle_unit) {
    if(angle_unit!="rad") alpha1 <- umrechnenWinkel(alpha1, angle_unit, "rad")
    ceiling(alpha1/(pi/2)) 
    
    # quadrant 2 und 3: hinter P1 relativ zu P1P2
    # quadrant 1 und 4: vor P1 relativ zu P1P2
    # quadrant 1 und 2: kartesisch links von P1P2, geodätisch rechts
    # quadrant 3 und 4: kartesisch rechts von P1P2, geodätisch links
    
}

### Polaraufnahme
### Die Punkte P1 und P2 bekannt, die Strecke s2 und der Winkel alpha1 sind gemessen.
### Gesucht ist der Punkt P3.

#### Berechnet den Punkt P3 aus P1, P2, alpha1 und s2
polaraufnahme_P3 <- function(
    P1=c(0,0),
    P2=c(0,0),
    alpha1=0,
    s2=0,
    angle_unit=c("deg","grad","gon","rad")) {
    
    alpha1 <- umrechnenWinkel(alpha1, angle_unit, "rad")

    P1+(s2/strecke(P1,P2))*(cos(alpha1)*(P2-P1)+sin(alpha1)*onv(P2-P1))
}

#### Visualisiert die Polaraufnahme
## Hier noch hinzufügen: Visualisierung von alpha1!
plot_polaraufnahme <- function(
    P1=c(0,0),
    P2=c(0,0),
    P3=NULL,
    alpha1=0,
    s2=0,
    angle_unit=c("deg","grad","gon","rad"),
    orientation=c("geodetic", "cartesian")) {
    
    if(is.null(P3)) P3 <- polaraufnahme_P3(P1, P2, alpha1, s2, angle_unit)
    alpha1 <- umrechnenWinkel(alpha1, angle_unit, "rad")
    
    dat <- tribble(
        ~point, ~x, ~y, 
        "P1", P1[1], P1[2], 
        "P2", P2[1], P2[2], 
        "P3", P3[1], P3[2])
    nudge_x <- (range(dat["x"])[2]-range(dat["x"])[1])/20
    PL <- lotfuss(P1, P2, (s2*cos(alpha1))/strecke(P1,P2))
    Ponv <- PL+onv(P2-P1)
    
    g <- ggplot()+
        # circle around P1 with radius s2
        annotate(
            "path",
            x=P1[1]+s2*cos(seq(0,2*pi,length.out=100)),
            y=P1[2]+s2*sin(seq(0,2*pi,length.out=100)),
            col="grey")+
        # PL + orientierter Normalenvektor (P2-P1)°
        geom_segment(aes(
            x=PL[1],
            y=PL[2],
            xend=Ponv[1],
            yend=Ponv[2]),
            arrow=arrow(
                length=unit(12, "points"),
                angle=20,
                type="closed"),
            col="darkgrey")+
        # label (P2-P1)°
        geom_text(aes(x=mean(c(PL[1], Ponv[1])),
                      y=mean(c(PL[2], Ponv[2]))),
                  label="(P2-P1)^o",
                  nudge_x=nudge_x,
                  hjust=0,
                  col="grey",
                  alpha=.67)+
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
        # label s3
        geom_text(aes(x=mean(c(P1[1], P2[1])),
                      y=mean(c(P1[2], P2[2]))),
                  label="s3",
                  nudge_x=nudge_x,
                  hjust=0,
                  col="grey",
                  alpha=.67)+
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
        # second path segment from P1 to P3 (h)
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
        # label h
        geom_text(aes(x=mean(c(PL[1], P3[1])),
                      y=mean(c(PL[2], P3[2]))),
                  label="h",
                  nudge_x=nudge_x,
                  hjust=0,
                  col="red",
                  alpha=.5)+
        # s2
        annotate(
            "path",
            x=dat[c(3,1), "x"],
            y=dat[c(3,1), "y"], 
            col="grey")+
        # label s2
        geom_text(aes(x=mean(c(P1[1], P3[1])),
                      y=mean(c(P1[2], P3[2]))),
                  label="s2",
                  nudge_x=nudge_x,
                  hjust=0,
                  col="grey",
                  alpha=.67)+
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

### Vorwaertsschnitt
### Die Punkte P1 und P2 bekannt, die Winkel alpha1 und alpha2 sind gemessen.
### Gesucht ist der Punkt P3.

#### Berechnet den Punkt P3 aus P1, P2, alpha1 und s2
vorwaertsschnitt_P3 <- function(
    P1=c(0,0),
    P2=c(0,0),
    alpha1=0,
    alpha2=0,
    angle_unit=c("deg","grad","gon","rad")) {
    
    alpha1 <- umrechnenWinkel(alpha1, angle_unit, "rad")
    alpha2 <- umrechnenWinkel(alpha2, angle_unit, "rad")
    
    P1+1/(cot(alpha1)+cot(alpha2))*(cot(alpha1)*(P2-P1)+onv(P2-P1))
}

#### Visualisiert des Vorwaertsschnitts
## Hier noch hinzufügen: Visualisierung von alpha1 und alpha2!
plot_vorwaertsschnitt <- function(
    P1=c(0,0),
    P2=c(0,0),
    P3=NULL,
    alpha1=0,
    alpha2=0,
    angle_unit=c("deg","grad","gon","rad"),
    orientation=c("geodetic", "cartesian")) {
    
    if(is.null(P3)) P3 <- vorwaertsschnitt_P3(P1, P2, alpha1, alpha2, angle_unit)
    alpha1 <- umrechnenWinkel(alpha1, angle_unit, "rad")
    alpha2 <- umrechnenWinkel(alpha2, angle_unit, "rad")
    
    dat <- tribble(
        ~point, ~x, ~y, 
        "P1", P1[1], P1[2], 
        "P2", P2[1], P2[2], 
        "P3", P3[1], P3[2])
    nudge_x <- (range(dat["x"])[2]-range(dat["x"])[1])/20
    PL <- lotfuss(P1, P2, (1/tan(alpha1))/(1/tan(alpha1)+1/tan(alpha2)))
    Ponv <- PL+onv(P2-P1)

        g <- ggplot()+
        # PL + orientierter Normalenvektor (P2-P1)°
        geom_segment(aes(
            x=PL[1],
            y=PL[2],
            xend=Ponv[1],
            yend=Ponv[2]),
            arrow=arrow(
                length=unit(12, "points"),
                angle=20,
                type="closed"),
            col="darkgrey")+
        # label (P2-P1)°
        geom_text(aes(x=mean(c(PL[1], Ponv[1])),
                      y=mean(c(PL[2], Ponv[2]))),
                  label="(P2-P1)^o",
                  nudge_x=nudge_x,
                  hjust=0,
                  col="grey",
                  alpha=.67)+
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
        # label s3
        geom_text(aes(x=mean(c(P1[1], P2[1])),
                      y=mean(c(P1[2], P2[2]))),
                  label="s3",
                  nudge_x=nudge_x,
                  hjust=0,
                  col="grey",
                  alpha=.67)+
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
        # second path segment from P1 to P3 (h)
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
        # label h
        geom_text(aes(x=mean(c(PL[1], P3[1])),
                      y=mean(c(PL[2], P3[2]))),
                  label="h",
                  nudge_x=nudge_x,
                  hjust=0,
                  col="red",
                  alpha=.5)+
        # s2
        annotate(
            "path",
            x=dat[c(3,1), "x"],
            y=dat[c(3,1), "y"], 
            col="grey")+
        # label s2
        geom_text(aes(x=mean(c(P1[1], P3[1])),
                      y=mean(c(P1[2], P3[2]))),
                  label="s2",
                  nudge_x=nudge_x,
                  hjust=0,
                  col="grey",
                  alpha=.67)+
        # s1
        annotate(
            "path",
            x=dat[c(3,2), "x"],
            y=dat[c(3,2), "y"], 
            col="grey")+
        # label s1
        geom_text(aes(x=mean(c(P2[1], P3[1])),
                      y=mean(c(P2[2], P3[2]))),
                  label="s1",
                  nudge_x=nudge_x,
                  hjust=0,
                  col="grey",
                  alpha=.67)+
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
        g <- g+coord_flip(xlim=xlim, ylim=ylim)+
        theme(aspect.ratio=1)
    }
    else if(orientation!="cartesian") warning("Invalid orientation paramenter. Coordinate system defaulted to cartesian")
    
    print(g)
}

### Bogenschnitt
### Die Punkte P1 und P2 bekannt, die Strecken s1 und s2 sind gemessen.
### Gesucht ist der Punkt P3.

#### Berechnet den Punkt P3 aus P1, P2, s1 und s2
bogenschnitt_P3 <- function(P1=c(0,0), 
                            P2=c(0,0), 
                            s1=0, 
                            s2=0, 
                            orientation=c("geodetic", "cartesian"), 
                            position=c("left", "right")) {
    if (position=="right" & orientation=="geodetic" | position=="left" & orientation=="cartesian") {
        sign <- 1
    } else if (position=="left" & orientation=="geodetic" | position=="right" & orientation=="cartesian") {
        sign <- -1
    } else stop("Invalid position or orientation.")

    s3square <- streckenquadrat(P1, P2)
    t <- param_ps3(s1, s2, s3square)
    u <- param_hs3(s2, s3square, t)
    
    P1+t*(P2-P1)+sign*u*onv(P2-P1)
}

#### Visualisiert den Bogenschnitt
plot_bogenschnitt <- function(
    P1=c(0,0),
    P2=c(0,0),
    P3_left=NULL,
    P3_right=NULL,
    s1=0,
    s2=0,
    orientation=c("geodetic", "cartesian")) {
    
    if(is.null(P3_left)) P3_left <- bogenschnitt_P3(P1, P2, s1, s2, orientation=orientation, position="left")
    if(is.null(P3_right)) P3_right <- bogenschnitt_P3(P1, P2, s1, s2, orientation=orientation, position="right")
    
    dat <- tribble(
        ~point, ~x, ~y, 
        "P1", P1[1], P1[2], 
        "P2", P2[1], P2[2], 
        "P3 links", P3_left[1], P3_left[2], 
        "P3 rechts", P3_right[1], P3_right[2])
    nudge_x <- (range(dat["x"])[2]-range(dat["x"])[1])/20
    PL <- lotfuss(P1, P2, param_ps3(s1, s2, streckenquadrat(P1, P2)))
    Ponv <- PL+onv(P2-P1)
    
    
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
        # PL + orientierter Normalenvektor (P2-P1)°
        geom_segment(aes(
            x=PL[1],
            y=PL[2],
            xend=Ponv[1],
            yend=Ponv[2]),
            arrow=arrow(
                length=unit(12, "points"),
                angle=20,
                type="closed"),
            col="darkgrey")+
        # label (P2-P1)°
        geom_text(aes(x=mean(c(PL[1], Ponv[1])),
                      y=mean(c(PL[2], Ponv[2]))),
                  label="(P2-P1)^o",
                  nudge_x=nudge_x,
                  hjust=0,
                  col="grey",
                  alpha=.67)+
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