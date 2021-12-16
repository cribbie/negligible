#' Proportional Distance Function (post hoc function - not to be used independently)
#'
#' @param effect observed effect
#' @param PD proportional distance for effect
#' @param EIsign equivalence interval value of the same sign as the effect
#' @param PDcil lower bound of the CI for the proportional distance
#' @param PDciu upper bound of the CI for the proportional distance
#' @param cil lower bound of the CI for the effect
#' @param ciu upper bound of the CI for the effect
#' @param Elevel 1-2alpha CI for the effect
#' @param Plevel 1-alpha CI for the PD
#' @param save Whether to save the plot or not
#'
#' @return nothing is returned
#' @export
#'
#' @examples
#' \dontrun{
#' 1+1
#' }

ProportionalDistance<-function(effect, PD, EIsign,
        PDcil, PDciu,  cil, ciu, Elevel, Plevel, save) {

  EIp<-c(EIsign,0) #actual value
  PDp<-c(effect,1) #not actual value
  effectp<-c(effect,0) #actual value
  start1<-c(0,1) #orgin for PD scale
  start2<-c(0,0) #orgin for effect scale

  #calculate 1/4 of the total 100*(1-x$alpha) distance to create the length of the scales (note that this is an arbritrary number)
  onefourth<-(abs((cil/2)+ (PDcil*effect/PD)) + abs((ciu/2)+ (PDciu*effect/PD)))/4
  extraL<-c(-3*onefourth,0)
  extraU<-c(3*onefourth,0)
  PDextraL<-c(-3*onefourth,1)
  PDextraU<-c(3*onefourth,1)

  #individual scale values after the (0,0) co-ordinate
  first<-onefourth
  second <- first + onefourth
  third<-second + onefourth
  first<-c(first,0)
  second<-c(second,0)
  third<-c(third,0)
  PDfirst<-onefourth
  PDsecond <-PDfirst+onefourth
  PDthird<-PDsecond+onefourth
  PDfirst<-c(PDfirst,1)
  PDsecond<-c(PDsecond,1)
  PDthird<-c(PDthird,1)

  #individual scale values before the (0,0) co-ordinate
  Bfirst<-onefourth
  Bsecond <- Bfirst + onefourth
  Bthird<-Bsecond + onefourth
  Bfirst<-c(-Bfirst,0)
  Bsecond<-c(-Bsecond,0)
  Bthird<-c(-Bthird,0)

  BPDfirst<-onefourth
  BPDsecond <-BPDfirst+onefourth
  BPDthird<-BPDsecond+onefourth
  BPDfirst<-c(-BPDfirst,1)
  BPDsecond<-c(-BPDsecond,1)
  BPDthird<-c(-BPDthird,1)

  #create a dataframe containing all the co-ordinates of interest

  dat<-rbind(BPDthird, BPDsecond, BPDfirst, start1, PDp, PDextraL, PDextraU, PDfirst, PDsecond, PDthird, Bthird, Bsecond, Bfirst, start2, effectp, EIp, extraL, extraU, first, second, third)
  dat<- data.frame(dat)
  print(
   ggplot2::ggplot(data=dat, ggplot2::aes(x=X1,y=X2)) +
     ggplot2::geom_point(shape=c(3,3,3,3,23,3,3,3,3,3,3,3,3,3,21,8,16,16,3,3,3), fill=c("black" ,"black", "black","black", "darkred", "black", "black" ,"black", "black" ,"black", "black" ,"black", "black" ,"black" ,"purple", "blue", "black", "black","black", "black","black"), color=c("black" ,"black", "black","black", "darkred", "black", "black" ,"black", "black" ,"black", "black" ,"black", "black" ,"black" ,"purple", "blue", "black", "black","black", "black","black"), size=c(3,3,3,3,2,3,3,3,3,3,3,3,3,3,2,2,0.1,0.1,3,3,3), stroke = 0.7) +
     ggplot2::ylim(-0.5, 3) + ggplot2::annotate("rect", xmin = cil, xmax = ciu, ymin = -0.05, ymax = 0.05, alpha = .2,fill = "blue") +
     ggplot2::annotate("rect", xmin = PDcil*effect/PD, xmax = PDciu*effect/PD, ymin = 0.95, ymax = 1.05, alpha = .2,fill = "red") +
     ggplot2::geom_step(data=dat[6:7,], color = "black") +
     ggplot2::geom_step(data=dat[4:5,], color = "darkred") +
     ggplot2::geom_step(data=dat[17:18,], color = "black") +
     ggplot2::geom_step(data=dat[14:15,], color = "blue") +
     ggplot2::annotate("text",x=effect, y=1.12, label=format(round(PD, 2), nsmall = 2), size=3, colour = "darkred") +
     ggplot2::annotate("text",x=PDextraU[1], y=0.88, label=format(round((PDextraU[1]*PD)/effect, 2), nsmall = 2), size=3, colour = "black") +
     ggplot2::annotate("text",x=PDextraL[1], y=0.88, label=format(round((PDextraL[1]*PD)/effect, 2), nsmall = 2), size=3, colour = "black") +
     ggplot2::annotate("text",x=effect, y=0.12, label=format(round(effect, 2), nsmall = 2), size=3, colour = "purple") +
     ggplot2::annotate("text",x=cil, y=0.12, label=format(round(cil, 2), nsmall = 2), size=3, colour = "purple") +
     ggplot2::annotate("text",x=ciu, y=0.12, label=format(round(ciu, 2), nsmall = 2), size=3, colour = "purple") +
     ggplot2::annotate("text",x=PDcil*effect/PD, y=1.12, label=format(round(PDcil, 2), nsmall = 2), size=3, colour = "darkred") +
     ggplot2::annotate("text",x=PDciu*effect/PD, y=1.12, label=format(round(PDciu, 2), nsmall = 2), size=3, colour = "darkred") +
     ggplot2::annotate("text",x=extraU[1], y=-0.12, label=format(round(extraU[1], 2), nsmall = 2), size=3, colour = "black") +
     ggplot2::annotate("text",x=extraL[1], y=-0.12, label=format(round(extraL[1], 2), nsmall = 2), size=3, colour = "black") +
     ggplot2::annotate("text",x=onefourth, y=-0.12, label=format(round(onefourth, 2), nsmall = 2), size=3, colour = "black") +
     ggplot2::annotate("text",x=2*onefourth, y=-0.12, label=format(round(2*onefourth, 2), nsmall = 2), size=3, colour = "black") +
     ggplot2::annotate("text",x=3*onefourth, y=-0.12, label=format(round(3*onefourth, 2), nsmall = 2), size=3, colour = "black") +
     ggplot2::annotate("text",x=-onefourth, y=-0.12, label=format(round(-onefourth, 2), nsmall = 2), size=3, colour = "black") +
     ggplot2::annotate("text",x=-2*onefourth, y=-0.12, label=format(round(-2*onefourth, 2), nsmall = 2), size=3, colour = "black") +
    ggplot2::annotate("text",x=-3*onefourth, y=-0.12, label=format(round(-3*onefourth, 2), nsmall = 2), size=3, colour = "black") +
    ggplot2::annotate("text",x=0, y=-0.12, label=0, size=3, colour = "black") +
    ggplot2::annotate("text",x=0, y=0.88, label=0, size=3, colour = "black") +
    ggplot2::annotate("text",x=-onefourth, y=0.88, label=format(round((-onefourth*PD)/effect, 2), nsmall = 2), size=3, colour = "black") +
    ggplot2::annotate("text",x=-2*onefourth, y=0.88, label=format(round((-2*onefourth*PD)/effect, 2), nsmall = 2), size=3, colour = "black") +
    ggplot2::annotate("text",x=onefourth, y=0.88, label=format(round((onefourth*PD)/effect, 2), nsmall = 2), size=3, colour = "black") +
    ggplot2::annotate("text",x=2*onefourth, y=0.88, label=format(round((2*onefourth*PD)/effect, 2), nsmall = 2), size=3, colour = "black") +
    ggplot2::annotate("text",x=-3*onefourth, y=0.4, label="Effect (Original Units)", size=3, hjust = 0) +
    ggplot2::annotate("text",x=-3*onefourth, y=1.4, label="The Proportional Distance (Effect / |EI|)", size=3, hjust = 0) +
    ggplot2::annotate("text",x=-3*onefourth, y=2.2, label=(paste0("    Observed effect and its ",Elevel,"% CI")), size=3.5, colour = "purple", hjust=0) +
    ggplot2::geom_point(ggplot2::aes(x=-3*onefourth, y=2.2), shape = 22, colour="blue", size = 4, fill = "grey") +
    ggplot2::geom_point(ggplot2::aes(x=-3*onefourth, y=2.2), shape = 16, colour="purple", size = 3) +
    ggplot2::annotate("text",x=-3*onefourth, y=2.6, label="    Equivalence Interval (EI) Boundary", size=3.5, colour = "blue4", hjust = 0) +
    ggplot2::geom_point(ggplot2::aes(x=-3*onefourth, y=2.6), shape = 8, colour="blue", size = 3) +
    ggplot2::annotate("text",x=-3*onefourth, y=3, label=(paste0("    Proportional Distance and its ",Plevel, "% CI")), size=3.5, colour = "darkred", hjust = 0) +
    ggplot2::geom_point(ggplot2::aes(x=-3*onefourth, y=3), shape = 22, colour="red", size = 4, fill = "grey") +
    ggplot2::geom_point(ggplot2::aes(x=-3*onefourth, y=3), shape = 18, colour="darkred", size = 3) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), axis.title = ggplot2::element_blank(),
           panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_blank(), axis.ticks= ggplot2::element_blank(), axis.text = ggplot2::element_blank())
   )

  if (!(save == FALSE)) {
    if (save == "png") {
      ggplot2::ggsave(filename = "Proportional_Distance_Plot.png" ,width = 9, height = 4, device='png', dpi=700)
    }
    if (save == "jpeg") {
      ggplot2::ggsave(filename = "Proportional_Distance_Plot.jpeg" ,width = 9, height = 4, device='png', dpi=700)
    }
  }
}






