#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(tidyr)
library(tidygraph) 
library(particles) 
#library(magicaxis)

# Reproducing Holmberg in 2021 <!--{{{-->
#vtot<-function(gr) {
#  v<-velocity(gr)
#  v<-sqrt(v[,1]^2+v[,2]^2)
#  v<-v-min(v)
#  v<-v/max(v)
#  v<-v*2/3
#  return=v
#}
map<-function(C,minmax=NULL) {
  if (!is.null(minmax)) { 
    #print(paste("Minmax:",paste(collapse=' ',minmax)))
    #print(range(C),minmax)
    C[which(C<minmax[1])]<-minmax[1]
    C[which(C>minmax[2])]<-minmax[2]
    C<-C-minmax[1]
    C<-C/(minmax[2]-minmax[1])
  } else { 
    #print(paste("Range:",paste(collapse=' ',range(C))))
    C<-C-min(C)
    C<-C/max(C)
  }
  C<-(1-C)*2/3
  # print(range(C))
  return=C
}
#map<-function(X) magicaxis::magmap(X,flip=TRUE)$map
pos.plot<-function(gr) 
  plot(pch=20,cex=cex,col=c(rep(cols[1],N.part),rep(cols[2],N.part)),
          position(gr),xlim=c(-box,box),ylim=c(-box,box))
velpos<-function(gr,npart,ngal) { 
  vx<-velocity(gr)[,1]
  vy<-velocity(gr)[,2]
  pos<-position(gr)
  vx.per<-vy.per<-NULL
  for (i in 1:ngal) { 
    gal.index<-(i-1)*npart + 1:npart
    vx.per<-c(vx.per,mean(vx[gal.index]))
    vy.per<-c(vy.per,mean(vy[gal.index]))
    vx[gal.index]<-vx[gal.index]-vx.per[i]
    vy[gal.index]<-vy[gal.index]-vy.per[i]
  }
  return=list(xpos=pos[,1],ypos=pos[,2],xvel=vx,yvel=vy,xvel.per=vx.per,yvel.per=vy.per)
}
arr.plot<-function(gr,npart,ngal,box,...) {  
  parm<-par(mar=c(2,2,2,1))
  vp<-velpos(gr,npart=npart,ngal=ngal)
  #print(range(sqrt(vx^2+vy^2)))
  #print(range(map(sqrt(vx^2+vy^2))))
  with(vp,{
    plot(xpos,ypos,col=hsv(map(sqrt(xvel^2+yvel^2))),pch=20,cex=cex,
            xlim=c(-box,box),ylim=c(-box,box),asp=1,
            #mtline = 1,
            xlab=expression(Delta*x),ylab=expression(Delta*y),
            ...)
            #labels=c(F,F,F,F),side=1:4,...)
            #labels=c(F,F,F,F),side=1:4,...)
    
    
    suppressWarnings({
      arrows(xpos,ypos,
             xpos+xvel*fact,
             ypos+yvel*fact,
             length=0.05,col=hsv(map(sqrt(xvel^2+yvel^2)),a=0.7)
      )
      for (i in 1:ngal) { 
        gal.index<-(i-1)*npart + 1:npart
        mean.x<-mean(xpos[gal.index])
        min.y<-min(ypos[gal.index])
        width.y<-abs(diff(range(ypos[gal.index])))
        arrows(mean.x,min.y-width.y/2,
               mean.x+xvel.per[i]*fact,
               min.y-width.y/2+yvel.per[i]*fact,
               length=0.05,lwd=2
        )
      }
    })
  })
  par(parm)
  return=range(sqrt(vp$xvel^2+vp$yvel^2))
}

cex<-1.0
box<-20
N.step<-10
N.part<-37
N.gal<-2
radius<-1
f.grav<-0.06
v.scale<-1
x.offset<-  -(5)
y.offset<-   (5)
vel.pec<-  1.5
vx.offset<- 0.3
v.mult<-sqrt(f.grav*v.scale)
fact<-2e-1/v.mult




cex<-1.0
box<-20
N.step<-1
N.part<-37
N.gal<-2
radius<-1
f.grav<-0.06
v.scale<-1
x.offset<-  -(5)
y.offset<-   (5)
vel.pec<-  1.5
vx.offset<- 0.3
v.mult<-sqrt(f.grav*v.scale)
fact<-2e-1/v.mult

tmp<-create_empty(N.part) %>%
  simulate(alpha=1,alpha_target=1,alpha_decay=0,velocity_decay = 0,
           setup = phyllotactic_genesis(radius=radius))

##Single "galaxy"
#pos.x<-position(tmp)[,1]
#pos.y<-position(tmp)[,2]
##phi = pi-atan(pos.y/pos.x)
#phi = atan2(y=pos.y,x=pos.x)
#vel.tot = sqrt(pos.x^2+pos.y^2)*v.mult
#vel.x = -vel.tot*sin(phi)*ifelse(pos.x<0,-1,1)
#vel.y = -vel.tot*cos(phi)*ifelse(pos.x<0,-1,1)
#
##Two "Galaxies" 
#if (N.gal>1) { 
#  vy.offset<-vel.pec*v.scale
#  pos.x<-c(pos.x+x.offset,
#           pos.x-x.offset)
#  pos.y<-c(pos.y+y.offset,
#           pos.y-y.offset)
#  vel.x<-c(vel.x-vx.offset,
#           (vel.x+vx.offset))
#  vel.y<-c(vel.y-vy.offset,
#           (vel.y+vy.offset))
#} else { 
#  vx.offset<-vy.offset<-0
#}
#
#holmberg.anticlockwise <- create_empty(N.part*N.gal) %>%
#  simulate(alpha=1,alpha_target=1,alpha_decay=0,velocity_decay = 0,
#           setup = predefined_genesis(pos.x,pos.y,vel.x,vel.y)) %>%
#  wield(manybody_force,f.grav,theta=0.01,min.dist=0.01,max.dist=0.1)

#minmax.clock<-arr.plot(holmberg.clockwise)
#minmax<-minmax.clock

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Holmberg 1943 in your browser!"),
    p(a("Read the lecture: Simulation",
        href = "https://anguswright.github.io/LectureNotes/IntroductionToStatistics_Section2d_html.html")),
    p(a("Explore uncertain galaxy encounters in the statistical companion",
        href = "../ForwardModelling/")),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          actionButton("goButton", "Go/Stop"),
          actionButton("reset", "Reset"),
            sliderInput("N.gal",
                      "Number of Galaxies:",
                      min = 1,
                      max = 15,
                      step=1,
                      value = 2),
          selectInput("rotation", "Rotation Direction:",
                      c("Clockwise" = "clockwise",
                        "Counter-Clockwise" = "counter-clockwise")),
            sliderInput("N.part",
                        "Number of particles per galaxy:",
                        min = 1,
                        max = 500,
                        value = 50),
            sliderInput("N.step",
                        "Number of simulation steps:",
                        min = 1,
                        max = 50,
                        value = 15),
            sliderInput("radius",
                        withMathJax("Galaxy radius scaling factor $$(R=10^x \\times R_0)$$"),
                        min = -1, 
                        max = 1,
                        step = 0.01, 
                        value = 0), 
            sliderInput("grav.scale",
                        withMathJax("Gravity scaling factor $$(g=10^x g_0)$$"),
                        min = -1, 
                        max = 1,
                        step = 0.01, 
                        value = 0), 
            sliderInput("v.scale",
                        expression("Velocity scaling factor $$(v=10^x v_0)$$"),
                        min = -1,
                        max = 1,
                        step = 0.01, 
                        value = 0),
            sliderInput("xoff.scale",
                        expression("x-position offset scaling factor $$(x_{\\rm off}=10^x x_{\\rm off, 0})$$"),
                        min = -1,
                        max = 1,
                        step = 0.01, 
                        value = 0),
            sliderInput("yoff.scale",
                        expression("y-position offset scaling factor $$(y_{\\rm off}=10^x\\times y_{\\rm off,0})$$"),
                        min = -1,
                        max = 1,
                        step = 0.01, 
                        value = 0),
            sliderInput("vpec.scale",
                        expression("peculiar velocity scaling factor $$(v_{\\rm pec} = 10^x \\times v_{\\rm pec,0})$$"),
                        min = -1,
                        max = 1,
                        step = 0.01, 
                        value = 0),
            sliderInput("vpec.phi",
                        expression("peculiar velocity direction:"),
                        min = -180,
                        max = 180,
                        step = 1, 
                        value = 90)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot",width="100%")
        )
    )
)

counter<-1 

# Define server logic required to draw a histogram
server <- function(input, output,session) {
  
    ts <- reactiveValues(N.part=isolate(input$N.part),
                         N.gal=isolate(input$N.gal),
                         radius=isolate(input$radius),
                         xoff.scale=isolate(input$xoff.scale),
                         yoff.scale=isolate(input$yoff.scale),
                         vpec.scale=isolate(input$vpec.scale),
                         vpec.phi=isolate(input$vpec.phi),
                         counter=0)
    
    #radius<-1
    #f.grav<-0.06
    #v.scale<-1
    #x.offset<-  -(5)
    #y.offset<-   (5)
    #vel.pec<-  1.5
    #vx.offset<- 0.3
    output$distPlot <- renderPlot({
      #if (length(ts$holmberg.clockwise)!=0) { 
      if (nrow(position(ts$holmberg))!=input$N.part*input$N.gal | 
          input$radius != ts$radius | 
          input$vpec.phi != ts$vpec.phi |
          input$xoff.scale != ts$xoff.scale |
          input$yoff.scale != ts$yoff.scale |
          input$vpec.scale != ts$vpec.scale | 
          input$rotation != ts$rotation) { 
        tmp<-create_empty(input$N.part) %>%
          simulate(alpha=1,alpha_target=1,alpha_decay=0,velocity_decay = 0,
                   setup = phyllotactic_genesis(radius=10^input$radius))
        
        #Single "galaxy"
        pos.x<-position(tmp)[,1]
        pos.y<-position(tmp)[,2]
        #phi = pi-atan(pos.y/pos.x)
        phi = pi-atan2(x=pos.y,y=pos.x)
        v.mult<-sqrt(f.grav*v.scale*10^input$v.scale)
        vel.tot = sqrt(pos.x^2+pos.y^2)*v.mult
        vel.x = vel.tot*cos(phi)#*ifelse(pos.x<0,-1,1)
        vel.y = vel.tot*sin(phi)#*ifelse(pos.y<0,-1,1)
        if (input$rotation=='clockwise') { 
          vel.x <- -vel.x
          vel.y <- -vel.y
        }
        
        if (input$N.gal==2) { 
          #Two "Galaxies" 
          vpec<-10^input$vpec.scale*vel.pec*v.scale
          vx.offset = vpec*cos(input$vpec.phi*2*pi/360)#*ifelse(pos.x<0,-1,1)
          vy.offset = vpec*sin(input$vpec.phi*2*pi/360)#*ifelse(pos.y<0,-1,1)
          pos.x<-c(pos.x+10^input$xoff.scale*x.offset,
                   pos.x-10^input$xoff.scale*x.offset)
          pos.y<-c(pos.y+10^input$yoff.scale*y.offset,
                   pos.y-10^input$yoff.scale*y.offset)
          vel.x<-c(vel.x+vx.offset,
                   vel.x-vx.offset)
          vel.y<-c(vel.y-vy.offset,
                   vel.y+vy.offset)
        } else if (input$N.gal>1) { 
          #Two "Galaxies" 
          vpec<-10^input$vpec.scale*vel.pec*v.scale
          vx.offset = vpec*cos(input$vpec.phi*2*pi/360)#*ifelse(pos.x<0,-1,1)
          vy.offset = vpec*sin(input$vpec.phi*2*pi/360)#*ifelse(pos.y<0,-1,1)
          single<-list(pos.x=pos.x,pos.y=pos.y,vel.x=vel.x,vel.y=vel.y)
          pos.x<-pos.y<-vel.x<-vel.y<-NULL
          if (input$N.gal<5) { 
            loc<-expand.grid(c(-1,1),c(-1,1)) 
          } else { 
            loc<-matrix(runif(input$N.gal*2,min=-1.5,1.5),ncol=2) 
          } 
          for (i in 1:input$N.gal) { 
            pos.x<-c(pos.x,single$pos.x+loc[i,1]*10^input$xoff.scale*x.offset)
            pos.y<-c(pos.y,single$pos.y+loc[i,2]*10^input$yoff.scale*y.offset)
            vel.x<-c(vel.x,single$vel.x+loc[i,1]/abs(loc[i,1])*vx.offset)
            vel.y<-c(vel.y,single$vel.y+loc[i,2]/abs(loc[i,2])*vy.offset)
          }
        } else if (input$N.gal==1) { 
          vy.offset<-10^input$vpec.scale*vel.pec*v.scale
          pos.x<-c(pos.x+10^input$xoff.scale*x.offset)
          pos.y<-c(pos.y+10^input$yoff.scale*y.offset)
          vel.x<-c(vel.x-10^input$vpec.scale*vx.offset)
          vel.y<-c(vel.y-vy.offset)
        }
        
        ts$holmberg <- create_empty(input$N.part*input$N.gal) %>%
          simulate(alpha=1,alpha_target=1,alpha_decay=0,velocity_decay = 0,
                   setup = predefined_genesis(pos.x,pos.y,vel.x,vel.y)) %>%
          wield(manybody_force,10^(input$grav.scale)*f.grav,theta=0.01,min.dist=0.01,max.dist=0.1)
        ts$radius<-input$radius
        ts$xoff.scale<-input$xoff.scale
        ts$yoff.scale<-input$yoff.scale
        ts$vpec.scale<-input$vpec.scale
        ts$vpec.phi<-input$vpec.phi
        ts$box<-max(abs(position(ts$holmberg)))*1.5
        ts$rotation<-input$rotation
        arr.plot(isolate(ts$holmberg),
                 npart=input$N.part,
                 ngal=input$N.gal,box=ts$box,
                 main=paste('INITIALISATION')) 
      } else { 
        arr.plot(isolate(ts$holmberg),
                 npart=input$N.part,
                 ngal=input$N.gal,box=ts$box,
                 main=paste('Iteration',ts$counter)) 
      }
      #} else { 
      #  plot(1,1,type='n',xlab='',ylab='',axes=F)
      #  text(1,1,pos=NULL,labels="Press 'Go!' to start",cex=2)
      #  text(1,0.8,pos=NULL,labels=paste("(counter=",ts$counter,")"))
      #}
    },height=1000,width=1000)
      
    
    observe({
      isolate({
        if (ts$counter==0 | length(ts$holmberg)==0){
          tmp<-create_empty(input$N.part) %>%
            simulate(alpha=1,alpha_target=1,alpha_decay=0,velocity_decay = 0,
                     setup = phyllotactic_genesis(radius=10^input$radius))
      
          #Single "galaxy"
          pos.x<-position(tmp)[,1]
          pos.y<-position(tmp)[,2]
          #phi = pi-atan(pos.y/pos.x)
          phi = pi-atan2(x=pos.y,y=pos.x)
          v.mult<-sqrt(f.grav*v.scale*10^input$v.scale)
          vel.tot = sqrt(pos.x^2+pos.y^2)*v.mult
          vel.x = vel.tot*cos(phi)#*ifelse(pos.x<0,-1,1)
          vel.y = vel.tot*sin(phi)#*ifelse(pos.y<0,-1,1)
          if (input$rotation=='clockwise') { 
            vel.x <- -vel.x
            vel.y <- -vel.y
          }
          
          if (input$N.gal==2) { 
            #Two "Galaxies" 
            vy.offset<-10^input$vpec.scale*vel.pec*v.scale
            pos.x<-c(pos.x+10^input$xoff.scale*x.offset,
                     pos.x-10^input$xoff.scale*x.offset)
            pos.y<-c(pos.y+10^input$yoff.scale*y.offset,
                     pos.y-10^input$yoff.scale*y.offset)
            vel.x<-c(vel.x-10^input$vpec.scale*vx.offset,
                     vel.x+10^input$vpec.scale*vx.offset)
            vel.y<-c(vel.y-vy.offset,
                     vel.y+vy.offset)
          } else if (input$N.gal>1) { 
            #Two "Galaxies" 
            vpec<-10^input$vpec.scale*vel.pec*v.scale
            vx.offset = vpec*cos(input$vpec.phi*2*pi/360)#*ifelse(pos.x<0,-1,1)
            vy.offset = vpec*sin(input$vpec.phi*2*pi/360)#*ifelse(pos.y<0,-1,1)
            single<-list(pos.x=pos.x,pos.y=pos.y,vel.x=vel.x,vel.y=vel.y)
            pos.x<-pos.y<-vel.x<-vel.y<-NULL
            if (input$N.gal<5) { 
              loc<-expand.grid(c(-1,1),c(-1,1)) 
            } else { 
              loc<-matrix(runif(input$N.gal*2,min=-1.5,1.5),ncol=2) 
            } 
            for (i in 1:input$N.gal) { 
              pos.x<-c(pos.x,single$pos.x+loc[i,1]*10^input$xoff.scale*x.offset)
              pos.y<-c(pos.y,single$pos.y+loc[i,2]*10^input$yoff.scale*y.offset)
              vel.x<-c(vel.x,single$vel.x+loc[i,1]/abs(loc[i,1])*vx.offset)
              vel.y<-c(vel.y,single$vel.y+loc[i,2]/abs(loc[i,2])*vy.offset)
            }
          } else if (input$N.gal==1) { 
            vy.offset<-10^input$vpec.scale*vel.pec*v.scale
            pos.x<-c(pos.x+10^input$xoff.scale*x.offset)
            pos.y<-c(pos.y+10^input$yoff.scale*y.offset)
            vel.x<-c(vel.x-10^input$vpec.scale*vx.offset)
            vel.y<-c(vel.y-vy.offset)
          }
          ts$holmberg <- create_empty(input$N.part*input$N.gal) %>%
            simulate(alpha=1,alpha_target=1,alpha_decay=0,velocity_decay = 0,
                     setup = predefined_genesis(pos.x,pos.y,vel.x,vel.y)) %>%
            wield(manybody_force,10^(input$grav.scale)*f.grav,
                  theta=10^(input$grav.scale)*f.grav/10, #0.01,
                  min.dist=10^(input$grav.scale)*f.grav/10, #0.01,
                  max.dist=10^(input$grav.scale)*f.grav*10) #0.1)
          ts$radius<-input$radius
          ts$xoff.scale<-input$xoff.scale
          ts$yoff.scale<-input$yoff.scale
          ts$vpec.scale<-input$vpec.scale
          ts$vpec.phi<-input$vpec.phi
          ts$box<-max(abs(position(ts$holmberg)))*1.5
          ts$rotation<-input$rotation
        }
        
        if (ts$counter > 1){
          ts$holmberg<-ts$holmberg %>% evolve(1)
        }
        
        ts$counter=ts$counter+1    
      })
      
      if (!is.null(input$goButton)) { 
        if (((isolate(ts$counter) < isolate(input$N.step))) & (input$goButton %% 2 > 0)){
          invalidateLater(200, session)
        }
      }
      
    })
    
    observe({
      if (input$reset > 0){
        #invalidateLater(200, session)
        ts$counter <<- 0
      }
    })
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
